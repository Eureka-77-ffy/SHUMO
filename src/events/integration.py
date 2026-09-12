"""Continue an unrounded Q23 state and locate its reconstructed global event.

No quasi-isothermal replacement is made. Output is reconstructed in chunks;
only occasional internal profiles are retained. The 4 h discontinuity is an
explicit left/right split, not a rapidly interpolated physical transient.
"""
import time

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.integrate import solve_ivp

from src.boundary import Environment
from src.coupled.integrate import reconstruct
from src.solver import InitializedBDF


class ScenarioEnvironment(Environment):
    def __init__(self, config, scenario="main"):
        super().__init__(config)
        if scenario == "last_observation":
            self.platform = self.data[-1, 1:].copy()
        elif scenario == "mean_9000_14400":
            self.platform = self.data[self.times >= 9000, 1:].mean(axis=0)
        elif scenario != "main":
            raise ValueError("Unknown environment extension")
        self.scenario = scenario
        self.right_at_join = False

    def value(self, t, field):
        t = np.asarray(t)
        value = super().value(t, field)
        if self.right_at_join:
            col = 0 if field == "temperature" else 1
            platform = self.platform[col] * (1 if col == 0 else self.beta)
            value = np.where(t >= self.times[-1], platform, value)
        return value


def global_maximum(model, t, state):
    """Includes even-axis value, every cell and the actual surface.

PCHIP cannot overshoot interval endpoints. The even quadratic on the first
half-cell is monotone in r^2, so its maximum is also at an endpoint.
"""
    C = state[1::2]
    x0, x1 = model.mesh.centres[:2]
    axis = (x1*x1*C[0]-x0*x0*C[1])/(x1*x1-x0*x0)
    surface = model.surface(t, state[-2], state[-1])[1]
    values = np.concatenate((np.asarray(axis).reshape((1,)+C.shape[1:]), C,
                             np.asarray(surface).reshape((1,)+C.shape[1:])), axis=0)
    return np.max(values, axis=0)


def continue_to_event(model, initial_state, start=10800., limit=604800.,
                      rtol=1e-10, atol_temperature=1e-11, atol_moisture=1e-13,
                      max_step=60., method="BDF", output_step=60., chunk_s=1800.):
    if start != 10800. or limit <= start:
        raise ValueError("Continuation must use the 3 h Q2 checkpoint and original clock")
    begun = time.perf_counter()
    n, R = model.mesh.n, model.mesh.radius
    initial_state = np.asarray(initial_state).copy()
    if initial_state.shape != (2*n,) or not np.isfinite(initial_state).all():
        raise ValueError("Checkpoint dimension/state mismatch")
    threshold = model.config["events"]["moisture_threshold"]
    if global_maximum(model, start, initial_state) <= threshold:
        raise ValueError("Checkpoint already meets threshold; earlier event must be searched")
    state = initial_state.copy()
    atol = np.tile([atol_temperature, atol_moisture], n)
    integrator = InitializedBDF if method == "BDF" else method
    radii = np.arange(21)*.001
    times_out, Ts, Cs, extras = [], [], [], {}
    profile_t, profile_y = [start], [state.copy()]
    cuts = np.unique(np.r_[start, model.environment.times[(model.environment.times > start)],
                           np.arange(14400., limit, chunk_s), limit])
    cuts = cuts[(cuts >= start) & (cuts <= limit)]
    gx, gw = leggauss(8)
    heat_out = water_out = 0.
    water_error = heat_error = 0.
    binit = model.thermal(state[1::2])[0]
    Hinit = float(np.sum(model.mesh.weights*binit*(state[::2]-model.T0)))
    Cinit = float(np.sum(model.mesh.weights*state[1::2]))
    initial_capacity = float(model.thermal(np.asarray(model.C0))[0])
    accepted, accepted_maximum, accepted_times = [], [], []
    nfev = nlu = njev = 0
    max_temporal_increase = 0.
    event_time = event_state = bracket = None
    minimum_C = float(np.min(state[1::2]))

    def event(t, y):
        return float(global_maximum(model, t, y)-threshold)
    event.direction, event.terminal = -1, True

    for a, b in zip(cuts[:-1], cuts[1:]):
        model.environment.right_at_join = a >= 14400.
        result = solve_ivp(model.rhs, (a, b), state, method=integrator, rtol=rtol,
                           atol=atol, jac=model.jacobian, dense_output=True,
                           max_step=max_step, events=event)
        if not result.success or not np.isfinite(result.y).all():
            raise RuntimeError(result.message)
        state = result.y[:, -1].copy()
        stop = float(result.t[-1])
        minimum_C = min(minimum_C, float(np.min(result.y[1::2])))
        if minimum_C <= 0:
            raise RuntimeError("Nonpositive accepted moisture")
        steps = np.diff(result.t)
        accepted.extend(steps.tolist())
        nfev += result.nfev; nlu += result.nlu; njev += result.njev
        # Accepted nodes plus midpoints: check spatial maxima, no centre-only surrogate.
        check_times = np.sort(np.r_[result.t, .5*(result.t[:-1]+result.t[1:])])
        check_max = global_maximum(model, check_times, result.sol(check_times))
        max_temporal_increase = max(max_temporal_increase, float(np.max(np.diff(check_max))))
        accepted_times.extend(check_times.tolist())
        accepted_maximum.extend(check_max.tolist())
        tsample = np.arange(np.floor(a/output_step)*output_step+output_step, stop+1e-9, output_step)
        tsample = tsample[(tsample > a) & (tsample <= stop)]
        if a == start:
            tsample = np.r_[start, tsample]
        if result.t_events[0].size:
            event_time = float(result.t_events[0][0])
            event_state = result.y_events[0][0].copy()
            if not len(tsample) or abs(tsample[-1]-event_time) > 1e-9:
                tsample = np.r_[tsample, event_time]
            else:
                tsample[-1] = event_time
            bracket = [float(result.t[-2]), event_time]
        # Bounded-memory spatial reconstruction of each output chunk.
        if len(tsample):
            fields, diag = reconstruct(model, tsample, result.sol(tsample), radii)
            times_out.append(tsample); Ts.append(fields[0]); Cs.append(fields[1])
            for key, value in diag.items():
                extras.setdefault(key, []).append(value)
        for first in range(0, len(steps), 32):
            last = min(first+32, len(steps))
            dt = steps[first:last]
            mid = .5*(result.t[first+1:last+1]+result.t[first:last])
            tq = (mid[:, None]+.5*dt[:, None]*gx).ravel()
            yq = result.sol(tq)
            dyq = model.rhs(tq, yq)
            bp = model.thermal(yq[1::2])[2]
            st, sc = model.surface(tq, yq[-2], yq[-1])
            ta, ce = model.ambient(tq)
            correction = np.sum(model.mesh.weights[:, None]*bp*(yq[::2]-model.T0)*dyq[1::2], axis=0)
            hrate = 2*model.h/R*(st-ta)-correction
            crate = 2*model.hm/R*(sc-ce)
            hi = heat_out+np.cumsum(.5*dt*np.sum(hrate.reshape(-1, 8)*gw, axis=1))
            ci = water_out+np.cumsum(.5*dt*np.sum(crate.reshape(-1, 8)*gw, axis=1))
            ye = result.y[:, first+1:last+1]
            H = np.sum(model.mesh.weights[:, None]*model.thermal(ye[1::2])[0]*(ye[::2]-model.T0), axis=0)
            mean = np.sum(model.mesh.weights[:, None]*ye[1::2], axis=0)
            heat_error = max(heat_error, float(np.max(np.abs(H-Hinit+hi))))
            water_error = max(water_error, float(np.max(np.abs(mean-Cinit+ci))))
            heat_out, water_out = float(hi[-1]), float(ci[-1])
        profile_t.append(stop); profile_y.append(state.copy())
        if result.t_events[0].size:
            break
        if b % 21600 == 0:
            print(f"t={b/3600:.1f} h, full-domain Cmax={check_max[-1]:.7f}", flush=True)

    cache = {"time_s": np.concatenate(times_out), "radius_m": radii,
             "temperature": np.concatenate(Ts), "moisture": np.concatenate(Cs),
             "profile_time_s": np.asarray(profile_t), "profile_state": np.asarray(profile_y),
             "final_state": state, "cell_centres_m": model.mesh.centres,
             "faces_m": model.mesh.faces, "weights": model.mesh.weights,
             "accepted_check_time_s": np.asarray(accepted_times),
             "accepted_check_global_max": np.asarray(accepted_maximum)}
    cache.update({key: np.concatenate(value) for key, value in extras.items()})
    stats = {"accepted_steps": len(accepted), "accepted_max_step_s": max(accepted),
             "accepted_min_step_s": min(accepted), "nfev": nfev, "nlu": nlu, "njev": njev,
             "minimum_accepted_moisture": minimum_C,
             "global_max_temporal_increase_at_checked_points": max_temporal_increase,
             "water_balance_tail_max_abs": water_error, "water_balance_relative": water_error/model.C0,
             "heat_corrected_balance_tail_max_abs_J_m3": heat_error,
             "heat_corrected_balance_relative": heat_error/(initial_capacity*22),
             "event_found": event_time is not None, "event_time_s": event_time,
             "boundary_extension": model.environment.scenario,
             "boundary_platform_T_Ceq": [float(model.environment.platform[0]), float(model.environment.platform[1]*model.environment.beta*model.beta)],
             "event_search": "max(all_cells,even_axis,true_surface), PCHIP bounded by nodes",
             "event_bracket_s": bracket, "balance_quadrature_order": 8}
    if event_time is not None:
        # Independently continue beyond the zero instead of calling it strictly dry.
        later = float(np.floor(event_time)+1)
        after = max(event_time+1., later)
        post = solve_ivp(model.rhs, (event_time, after), event_state, method=integrator,
                         rtol=rtol, atol=atol, jac=model.jacobian, dense_output=True,
                         max_step=min(1., max_step))
        if not post.success:
            raise RuntimeError("Post-event verification failed")
        pre_time = max(result.t[0], event_time-1.)
        check_t = np.asarray([pre_time, event_time, later, after])
        values = [event(pre_time, result.sol(pre_time)), event(event_time, event_state),
                  event(later, post.sol(later)), event(after, post.y[:, -1])]
        assert values[0] > 0 and values[2] < 0 and values[3] < 0
        slope = (values[-1]-values[0])/(after-pre_time)
        stats.update({"post_event_times_s": check_t.tolist(), "post_event_g": values,
                      "event_slope_per_s": float(slope), "first_checked_integer_second": later,
                      "strict_inequality_holds_at_critical_time": False})
        cache["post_event_state"] = post.y[:, -1]
        cache["event_state"] = event_state
    stats["elapsed_s"] = time.perf_counter()-begun
    return cache, stats
