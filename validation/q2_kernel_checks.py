#!/usr/bin/env python3
"""Coupled Jacobian, conservation identities, boundary and Q1 degeneration."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np

from scripts.export_question1 import read_case as read_q1
from scripts.question2 import fingerprint
from src.boundary import Environment
from src.coupled.model import CoupledFVM
from src.coupled.integrate import integrate
from src.data_io import ROOT, digest, load_config, save_json, verify_sources
from src.fvm import RadialMesh
from src.properties import PropertyLaw


def main():
    verify_sources()
    config = load_config()
    env = Environment(config)
    report = {"source_fingerprint": fingerprint(), "script_sha256": digest(__file__), "jacobian": []}
    for group, options in (("q1", {}), ("q23", {}),
                           ("q23", {"freeze_D_temperature": True}),
                           ("q23", {"freeze_thermal_moisture": True})):
        for grading in (1., 2.):
            model = CoupledFVM(RadialMesh(9, .02, grading), config, env, group, options)
            state = np.empty(18)
            state[::2], state[1::2] = np.linspace(32., 37., 9), np.linspace(2.5, 1.7, 9)
            numerical = np.column_stack([(model.rhs(600., state+1e-5*v)-model.rhs(600., state-1e-5*v))/2e-5 for v in np.eye(18)])
            analytic = model.jacobian(600., state).toarray()
            error = float(np.max(np.abs(analytic-numerical))/np.max(np.abs(analytic)))
            assert error < 1e-7
            report["jacobian"].append({"group": group, "grading": grading, "options": options, "relative_error": error})
            b = model.thermal(state[1::2])[0]
            _, _, _, _, qt, qc = model.fluxes(600., state)
            rhs = model.rhs(600., state)
            assert abs(np.sum(model.mesh.weights*b*rhs[::2])+2*qt/.02) < 1e-7
            assert abs(np.sum(model.mesh.weights*rhs[1::2])+2*qc/.02) < 1e-14
    model = CoupledFVM(RadialMesh(320, .02, 2.), config, env)
    C, T = np.array([2.55, 1.3, .2]), np.array([28., 42., 50.])
    b, k, D = PropertyLaw(config, "q23").evaluate(C, T)
    assert np.allclose(model.thermal(C)[0], b, rtol=1e-14)
    assert np.allclose(model.thermal(C)[1], k, rtol=1e-14)
    assert np.allclose(model.diffusion(C, T)[0], D, rtol=1e-14)
    test_times = np.array([1., 1000., 10000.])
    ts, cs = model.surface(test_times, T, C)
    ta, ce = model.ambient(test_times)
    km = model.thermal((C+cs)/2)[1]
    am = model.A*np.exp(-model.B/((T+ts)/2+273.15))
    heat_residual = km*(T-ts)/model.mesh.half_width-model.h*(ts-ta)
    water_residual = am*model.primitive.integral(cs, C)/model.mesh.half_width-model.hm*(cs-ce)
    report["boundary_heat_residual_W_m2"] = float(np.max(np.abs(heat_residual)))
    report["boundary_water_residual_m_s"] = float(np.max(np.abs(water_residual)))
    assert np.max(np.abs(heat_residual)) < 2e-7
    assert np.max(np.abs(water_residual)) < 1e-14
    model.h = model.hm = 0.
    state = np.tile([31., 1.2], model.mesh.n)
    assert np.max(np.abs(model.rhs(10., state))) < 1e-10
    rejected = False
    state[1] = -.1
    try:
        model.rhs(10., state)
    except ValueError:
        rejected = True
    assert rejected
    report["constant_insulated_state_and_nonphysical_rejection"] = True
    old, _ = read_q1("g2_n320")
    q1_model = CoupledFVM(RadialMesh(320, .02, 2.), config, env, "q1")
    cache, stats = integrate(q1_model, end=1800., rtol=1e-10, atol_temperature=1e-12,
                             atol_moisture=1e-12, max_step=60.)
    report["q1_degenerate_comparison"] = {
        f: float(np.max(np.abs(cache[f]-old[f]))) for f in ("temperature", "moisture")}
    assert all(value < 1e-7 for value in report["q1_degenerate_comparison"].values())
    report["q1_degenerate_statistics"] = stats
    report["status"] = "passed_coupled_kernel_checks_not_full_Q23_validation"
    save_json(ROOT/"reports/q2_kernel_checks.json", report)
    print(report)


if __name__ == "__main__":
    main()
