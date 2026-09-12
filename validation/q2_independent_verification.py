#!/usr/bin/env python3
"""Refine the independent FE implementation and compare unrounded fields."""
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np

from scripts.question2 import read_case, fingerprint
from src.boundary import Environment
from src.coupled.model import CoupledFVM
from src.coupled.integrate import reconstruct
from src.data_io import ROOT, digest, load_config, save_json
from src.fvm import RadialMesh
from validation.q2_independent_fem import run, IndependentFEM


def read_or_run(n, dt, method):
    label = f"fem_{method}_n{n}_dt{dt:g}"
    folder = ROOT/"results/cache/q2_independent"
    path, metapath = folder/(label+".npz"), folder/(label+".json")
    if path.exists() and metapath.exists():
        meta = json.loads(metapath.read_text())
        if (meta["script_sha256"] == digest(ROOT/"validation/q2_independent_fem.py")
                and meta["config_sha256"] == digest(ROOT/"config/model_config.json")
                and meta["cache_sha256"] == digest(path)):
            with np.load(path) as raw:
                return {k: raw[k] for k in raw.files}, meta
    return run(n, dt, method=method)


def main():
    config = load_config()
    test = IndependentFEM(11, config)
    assert abs(np.sum(test.multiply(test.mass, np.ones(12)))-.02**2/2) < 1e-17
    # The variable-step BDF2 expression differentiates quadratic polynomials.
    for ratio in (.1, .5, 1., 1.8):
        h = .2*ratio
        now = 2.
        older = now-h-.2
        value = ((1+2*ratio)/(1+ratio)*now**2-(1+ratio)*(now-h)**2
                 + ratio**2/(1+ratio)*older**2)/h
        assert abs(value-2*now) < 1e-12
    primary, meta = read_case("q2_final")
    model = CoupledFVM(RadialMesh(meta["n"], .02, meta["grading"]), config, Environment(config))
    report = {"source_fingerprint": fingerprint(), "script_sha256": digest(__file__),
              "independent_script_sha256": digest(ROOT/"validation/q2_independent_fem.py"),
              "independence": "own_nodal_mesh_consistent_mass_Gauss_assembly_direct_Robin_and_time_discretization",
              "time_s": [], "backward_euler": [], "bdf2": [], "complete": False}
    cases = [("backward_euler", 160, 30.), ("backward_euler", 320, 15.), ("backward_euler", 640, 7.5),
             ("bdf2", 320, 15.), ("bdf2", 640, 7.5), ("bdf2", 1280, 3.75),
             ("bdf2", 1280, 1.875), ("bdf2", 1280, .9375), ("bdf2", 640, .9375)]
    outputs = {}
    for method, n, dt in cases:
        candidate, cm = read_or_run(n, dt, method)
        outputs[(method, n, dt)] = candidate
        indices = candidate["time_s"].astype(int)
        report["time_s"] = candidate["time_s"].tolist()
        row = {"label": cm["label"], "elements": n, "max_step_s": dt, "comparisons_to_primary": {}}
        for field in ("temperature", "moisture"):
            diff = np.abs(candidate[field][1:]-primary[field][indices[1:]])
            where = np.unravel_index(np.argmax(diff), diff.shape)
            row["comparisons_to_primary"][field] = {"max_abs": float(np.max(diff)),
                "time_s": float(candidate["time_s"][where[0]+1]), "radius_cm": float(candidate["radius_m"][where[1]]*100)}
        report[method].append(row)
        save_json(ROOT/"reports/q2_independent_verification.json", report)
    final = outputs[("bdf2", 1280, .9375)]
    report["separate_refinements"] = {}
    for key, comparator in (("time_3p75_to_1p875", ("bdf2", 1280, 3.75)),
                             ("time_1p875_to_0p9375", ("bdf2", 1280, 1.875)),
                             ("space_640_to_1280_at_0p9375", ("bdf2", 640, .9375))):
        fine = outputs[("bdf2", 1280, 1.875)] if key == "time_3p75_to_1p875" else final
        coarse = outputs[comparator]
        report["separate_refinements"][key] = {f: float(np.max(np.abs(fine[f]-coarse[f]))) for f in ("temperature", "moisture")}
    paper = np.arange(1800., 10801., 1800.)
    pidx = np.searchsorted(primary["profile_time_s"], paper)
    fidx = np.searchsorted(final["time_s"], paper)
    assert np.array_equal(primary["profile_time_s"][pidx], paper)
    sampled, _ = reconstruct(model, paper, primary["profile_state"][pidx].T, final["nodes_m"])
    dense = {f: float(np.max(np.abs(sampled[j]-final[f+"_nodes"][fidx])))
             for j, f in enumerate(("temperature", "moisture"))}
    report["dense_paper_profile_comparison"] = {"times_s": paper.tolist(), "radial_points": len(final["nodes_m"]), "max_abs": dense}
    target = config["diagnostics"]["numerical_absolute_target_for_four_decimal_claim"]
    finest_row = next(r for r in report["bdf2"] if r["elements"] == 1280 and r["max_step_s"] == .9375)
    assert all(v["max_abs"] < target for v in finest_row["comparisons_to_primary"].values())
    assert all(v < target for v in dense.values())
    report["status"] = "passed_Q2_independent_spatial_and_temporal_sample_checks_not_continuous_error_bound"
    report["complete"] = True
    save_json(ROOT/"reports/q2_independent_verification.json", report)
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
