#!/usr/bin/env python3
"""Independent nodal P1 weak form, Gauss assembly, BE or own BDF2/Picard.

Does not import the primary mesh, RHS, face fluxes, boundary recovery,
reconstruction, or time integrator. Robin conditions act directly on endpoint
nodes. Consistent mass matrices are retained. Appendix constants and raw
observations are the shared physical inputs, not a shared discretization.
"""
import argparse
import json
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
from numpy.polynomial.legendre import leggauss
from openpyxl import load_workbook
from scipy.linalg import solve_banded

from src.data_io import ROOT, digest, load_config, save_json, verify_sources


class IndependentFEM:
    def __init__(self, elements, config):
        self.config = config
        self.p = config["property_sets"]["q23"]
        self.R = config["geometry"]["reference_radius_m"]
        s = np.linspace(0., 1., elements+1)
        self.nodes = self.R*(2*s-s*s)
        self.dx = np.diff(self.nodes)
        z, w = leggauss(3)
        self.N = np.column_stack(((1-z)/2, (1+z)/2))
        self.rq = self.nodes[:-1, None]*self.N[:, 0]+self.nodes[1:, None]*self.N[:, 1]
        self.wq = self.dx[:, None]/2*w*self.rq
        self.mass = self.mass_matrix(np.ones_like(self.rq))
        self.h = config["boundary"]["heat_transfer_W_m2_K"]
        self.hm = config["boundary"]["normalized_mass_transfer_m_s"]
        wb = load_workbook(Path(config["source_directory"])/"附件1.xlsx", read_only=True, data_only=True)
        self.air = np.asarray(list(wb.active.values)[1:], float)
        wb.close()

    def gauss_state(self, nodal):
        return nodal[:-1, None]*self.N[:, 0]+nodal[1:, None]*self.N[:, 1]

    def coefficients(self, C, T):
        p = self.p
        rho = p["rho_kg_m3"]["offset"]+p["rho_kg_m3"]["slope_C"]*C
        cp = p["cp_J_kg_K"]["offset"]+p["cp_J_kg_K"]["fraction_C_over_1_plus_C"]*C/(1+C)
        k = p["k_W_m_K"]["offset"]+p["k_W_m_K"]["fraction_C_over_1_plus_C"]*C/(1+C)
        d = p["D_m2_s"]
        D = d["prefactor"]*np.exp(-d["moisture_exponent"]/C-d["thermal_exponent_K"]/(T+273.15))
        return rho*cp, k, D

    def mass_matrix(self, coefficient):
        w = self.wq*coefficient
        left = np.sum(w*self.N[:, 0]**2, axis=1)
        right = np.sum(w*self.N[:, 1]**2, axis=1)
        off = np.sum(w*self.N[:, 0]*self.N[:, 1], axis=1)
        diagonal = np.r_[left, 0.]+np.r_[0., right]
        return diagonal, off

    def stiffness(self, coefficient):
        value = np.sum(self.wq*coefficient, axis=1)/self.dx**2
        return np.r_[value, 0.]+np.r_[0., value], -value

    @staticmethod
    def multiply(matrix, vector):
        diagonal, off = matrix
        out = diagonal*vector
        out[:-1] += off*vector[1:]
        out[1:] += off*vector[:-1]
        return out

    def linear_step(self, mass, stiffness, old, dt, robin, ambient):
        diagonal = mass[0]/dt+stiffness[0]
        off = mass[1]/dt+stiffness[1]
        diagonal = diagonal.copy()
        diagonal[-1] += self.R*robin
        # Solve for excess above the prescribed ambient. The stiffness applied
        # to a constant is exactly zero analytically. Removing that constant
        # before the linear solve avoids cancellation on very thin elements;
        # it changes neither the weak form nor the boundary condition.
        rhs = self.multiply(mass, old-ambient)/dt
        band = np.zeros((3, len(diagonal)))
        band[0, 1:], band[1], band[2, :-1] = off, diagonal, off
        return ambient+solve_banded((1, 1), band, rhs, check_finite=False)

    def step(self, t, dt, oldT, oldC, older=None, previous_dt=None):
        ta = np.interp(t, self.air[:, 0], self.air[:, 1])
        ce = self.config["boundary"]["beta"]*np.interp(t, self.air[:, 0], self.air[:, 2])
        T, C = oldT.copy(), oldC.copy()
        effectiveT, effectiveC, effective_dt = oldT, oldC, dt
        if older is not None:
            ratio = dt/previous_dt
            a0 = (1+2*ratio)/(1+ratio)
            effectiveT = ((1+ratio)*oldT-ratio**2/(1+ratio)*older[0])/a0
            effectiveC = ((1+ratio)*oldC-ratio**2/(1+ratio)*older[1])/a0
            effective_dt = dt/a0
        for iteration in range(50):
            Cg, Tg = self.gauss_state(C), self.gauss_state(T)
            capacity, k, _ = self.coefficients(Cg, Tg)
            newT = self.linear_step(self.mass_matrix(capacity), self.stiffness(k), effectiveT, effective_dt, self.h, ta)
            _, _, D = self.coefficients(Cg, self.gauss_state(newT))
            newC = self.linear_step(self.mass, self.stiffness(D), effectiveC, effective_dt, self.hm, ce)
            if np.min(newC) <= 0 or not np.isfinite(newT).all():
                raise RuntimeError("Nonphysical independent FEM iterate")
            errT, errC = np.max(np.abs(newT-T)), np.max(np.abs(newC-C))
            T, C = newT, newC
            if errT < 2e-10 and errC < 2e-12:
                return T, C, iteration+1
        raise RuntimeError(f"Independent Picard failed at t={t}, dt={dt}: dT={errT}, dC={errC}")


def run(elements, max_step, label=None, method="backward_euler"):
    started = time.perf_counter()
    verify_sources()
    cfg = load_config()
    model = IndependentFEM(elements, cfg)
    if method not in ("backward_euler", "bdf2"):
        raise ValueError("Unknown independent time method")
    paper_times = np.arange(1800., 10801., 1800.)
    times = np.r_[0., 1., 10., 100., 300., 600., 900., 1200., 1500., paper_times]
    cuts = np.unique(np.r_[times, model.air[model.air[:, 0] <= 10800, 0]])
    T = np.full(elements+1, cfg["initial"]["temperature_C"])
    C = np.full(elements+1, cfg["initial"]["moisture_dry_basis"])
    Trows, Crows = [T.copy()], [C.copy()]
    collected = [0.]
    t, steps, iterations = 0., 0, 0
    older, previous_dt = None, None
    for stop in cuts[1:]:
        while t < stop-1e-10:
            # Refine the incompatible initial water boundary at every level;
            # refinements also halve this early-time cap, not only the late cap.
            dt = min(max_step, .04*(max_step/30)*max(t, .01), stop-t)
            if method == "bdf2" and previous_dt is not None:
                dt = min(dt, 1.8*previous_dt)
            newT, newC, count = model.step(t+dt, dt, T, C,
                older=older if method == "bdf2" else None, previous_dt=previous_dt)
            older, previous_dt = (T, C), dt
            T, C = newT, newC
            t += dt
            steps += 1
            iterations += count
        t = float(stop)
        if np.any(np.isclose(times, stop, rtol=0, atol=1e-8)):
            collected.append(t)
            Trows.append(T.copy())
            Crows.append(C.copy())
    collected = np.asarray(collected)
    assert np.array_equal(collected, times)
    radii = np.arange(21)*.001
    cache = {"time_s": times, "radius_m": radii, "nodes_m": model.nodes,
             "temperature_nodes": np.asarray(Trows), "moisture_nodes": np.asarray(Crows),
             "temperature": np.array([np.interp(radii, model.nodes, row) for row in Trows]),
             "moisture": np.array([np.interp(radii, model.nodes, row) for row in Crows])}
    label = label or f"fem_{method}_n{elements}_dt{max_step:g}"
    folder = ROOT/"results/cache/q2_independent"
    folder.mkdir(parents=True, exist_ok=True)
    path = folder/(label+".npz")
    np.savez_compressed(path, **cache)
    meta = {"label": label, "method": "P1_consistent_mass_"+method+"_Picard",
            "elements": elements, "max_step_s": max_step, "steps": steps, "picard_iterations": iterations,
            "early_step_rule": "min(max_step,0.04*(max_step/30)*max(t,0.01),next_knot-t)",
            "bdf2_startup": "one_BE_step; subsequent_variable_step_BDF2; growth_limited_to_1.8",
            "initialization": "fresh", "property_set": "q23", "end_s": 10800.,
            "script_sha256": digest(__file__), "config_sha256": digest(ROOT/"config/model_config.json"),
            "cache_sha256": digest(path), "elapsed_s": time.perf_counter()-started}
    save_json(folder/(label+".json"), meta)
    print(json.dumps(meta), flush=True)
    return cache, meta


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--elements", type=int, default=160)
    parser.add_argument("--max-step", type=float, default=30.)
    parser.add_argument("--method", choices=["backward_euler", "bdf2"], default="backward_euler")
    args = parser.parse_args()
    run(args.elements, args.max_step, method=args.method)
