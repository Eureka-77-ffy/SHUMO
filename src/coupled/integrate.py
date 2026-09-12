"""Segmented coupled integration, streaming reconstruction and balance checks."""
import time

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.integrate import solve_ivp
from scipy.interpolate import PchipInterpolator

from src.solver import InitializedBDF


def reconstruct(model, times, state, radii):
    """PCHIP has no overshoot inside each interval; axis uses an even quadratic."""
    mesh = model.mesh
    T, C = state[0::2].T, state[1::2].T
    surfaces = model.surface(times, T[:, -1], C[:, -1])
    x0, x1 = mesh.centres[:2]
    points = np.r_[0., mesh.centres, mesh.radius]
    fields, extended = [], []
    for cells, surface, initial in zip((T, C), surfaces, (model.T0, model.C0)):
        axis = (x1*x1*cells[:, 0]-x0*x0*cells[:, 1])/(x1*x1-x0*x0)
        nodes = np.column_stack((axis, cells, surface))
        sampled = PchipInterpolator(points, nodes, axis=1)(radii)
        near = radii < x0
        sampled[:, near] = axis[:, None]+((cells[:, 1]-cells[:, 0])/(x1*x1-x0*x0))[:, None]*radii[near]**2
        sampled[np.asarray(times) == 0] = initial
        nodes[np.asarray(times) == 0] = initial
        fields.append(sampled)
        extended.append(nodes)
    # The central even quadratic and intervalwise PCHIP are bounded by endpoints.
    maximum_index = np.argmax(extended[1], axis=1)
    diagnostics = {
        "temperature_mean": np.sum(T*mesh.weights, axis=1),
        "moisture_mean": np.sum(C*mesh.weights, axis=1),
        "moisture_global_max": np.max(extended[1], axis=1),
        "moisture_argmax_radius_m": points[maximum_index],
        "moisture_radial_increase_max": np.maximum(0., np.max(np.diff(extended[1], axis=1), axis=1)),
    }
    return fields, diagnostics


def integrate(model, end=10800., rtol=1e-9, atol_temperature=1e-10,
              atol_moisture=1e-12, max_step=10., method="BDF", output_step=1.):
    started = time.perf_counter()
    n = model.mesh.n
    initial = np.empty(2*n)
    initial[0::2], initial[1::2] = model.T0, model.C0
    state = initial.copy()
    times = np.unique(np.r_[np.arange(0., end+1e-9, output_step), end])
    radii = np.arange(21)*.001
    cache = {"time_s": times, "radius_m": radii,
             "cell_centres_m": model.mesh.centres, "faces_m": model.mesh.faces,
             "weights": model.mesh.weights, "temperature": np.empty((len(times), 21)),
             "moisture": np.empty((len(times), 21))}
    knots = model.environment.times
    cuts = np.unique(np.r_[0., knots[(knots > 0) & (knots < end)], end])
    profile_times, profiles = [0.], [initial.copy()]
    gx, gw = leggauss(8)
    water_out = heat_out = 0.
    water_error = heat_error = semi_heat_error = 0.
    accepted, nfev, njev, nlu = [], 0, 0, 0
    minimum_C, minimum_T = model.C0, model.T0
    b0 = float(model.thermal(np.asarray(model.C0))[0])
    atol = np.tile([atol_temperature, atol_moisture], n)
    for a, b in zip(cuts[:-1], cuts[1:]):
        integrator = InitializedBDF if method == "BDF" else method
        result = solve_ivp(model.rhs, (a, b), state, method=integrator, rtol=rtol,
                           atol=atol, jac=model.jacobian, dense_output=True, max_step=max_step)
        if not result.success or not np.isfinite(result.y).all():
            raise RuntimeError(result.message)
        if np.min(result.y[1::2]) <= 0:
            raise RuntimeError("Nonpositive accepted moisture")
        state = result.y[:, -1].copy()
        minimum_C = min(minimum_C, float(np.min(result.y[1::2])))
        minimum_T = min(minimum_T, float(np.min(result.y[0::2])))
        indices = np.flatnonzero((times >= a) & (times <= b))
        if len(indices):
            fields, extra = reconstruct(model, times[indices], result.sol(times[indices]), radii)
            for name, values in zip(("temperature", "moisture"), fields):
                cache[name][indices] = values
            for key, values in extra.items():
                if key not in cache:
                    cache[key] = np.empty(len(times))
                cache[key][indices] = values
        steps = np.diff(result.t)
        accepted.extend(steps.tolist())
        nfev += result.nfev
        njev += result.njev
        nlu += result.nlu
        # Independent quadrature of the reconstructed boundary-flow history.
        # Variable-b heat balance uses H=mean[b(C)*(T-T0)], minus its b'(C)
        # correction, not an erroneous uncorrected difference of rho*cp*T.
        for first in range(0, len(steps), 64):
            last = min(first+64, len(steps))
            dt = steps[first:last]
            mid = .5*(result.t[first+1:last+1]+result.t[first:last])
            tq = (mid[:, None]+.5*dt[:, None]*gx).ravel()
            yq = result.sol(tq)
            dyq = model.rhs(tq, yq)
            capacity, _, bp, _ = model.thermal(yq[1::2])
            ts, cs = model.surface(tq, yq[-2], yq[-1])
            ta, ce = model.ambient(tq)
            qt, qc = model.h*(ts-ta), model.hm*(cs-ce)
            correction = np.sum(model.mesh.weights[:, None]*bp*(yq[0::2]-model.T0)*dyq[1::2], axis=0)
            heat_rate = 2*qt/model.mesh.radius-correction
            water_rate = 2*qc/model.mesh.radius
            heat_integral = heat_out+np.cumsum(.5*dt*np.sum(heat_rate.reshape(-1, 8)*gw, axis=1))
            water_integral = water_out+np.cumsum(.5*dt*np.sum(water_rate.reshape(-1, 8)*gw, axis=1))
            yend = result.y[:, first+1:last+1]
            bend = model.thermal(yend[1::2])[0]
            heat_content = np.sum(model.mesh.weights[:, None]*bend*(yend[0::2]-model.T0), axis=0)
            water_mean = np.sum(model.mesh.weights[:, None]*yend[1::2], axis=0)
            heat_error = max(heat_error, float(np.max(np.abs(heat_content+heat_integral))))
            water_error = max(water_error, float(np.max(np.abs(water_mean-model.C0+water_integral))))
            heat_out, water_out = float(heat_integral[-1]), float(water_integral[-1])
            semi = np.sum(model.mesh.weights[:, None]*capacity*dyq[0::2], axis=0)+2*qt/model.mesh.radius
            semi_heat_error = max(semi_heat_error, float(np.max(np.abs(semi))))
        profile_times.append(float(b))
        profiles.append(state.copy())
    cache["profile_time_s"] = np.asarray(profile_times)
    cache["profile_state"] = np.asarray(profiles)
    cache["final_state"] = state
    stats = {
        "accepted_steps": len(accepted), "accepted_max_step_s": max(accepted),
        "accepted_min_step_s": min(accepted), "nfev": nfev, "njev": njev, "nlu": nlu,
        "minimum_accepted_temperature_C": minimum_T, "minimum_accepted_moisture": minimum_C,
        "balance_quadrature_order": 8,
        "water_balance_max_abs_dry_basis": water_error,
        "water_balance_normalizer": model.C0, "water_balance_relative": water_error/model.C0,
        "heat_corrected_balance_max_abs_J_m3": heat_error,
        "heat_balance_normalizer_J_m3": b0*22., "heat_corrected_balance_relative": heat_error/(b0*22.),
        "heat_semidiscrete_balance_max_abs_W_m3": semi_heat_error,
        "heat_balance_definition": "mean[b(C)*(T-T0)] + integral(2*qT/R - mean[bprime*(T-T0)*Cdot])",
        "elapsed_s": time.perf_counter()-started,
    }
    return cache, stats
