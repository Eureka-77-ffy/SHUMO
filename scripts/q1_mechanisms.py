#!/usr/bin/env python3
"""Derive Q1 mechanism diagnostics from the audited, unrounded trajectory.

No new PDE solve, fitted parameter, or statistical uncertainty is introduced.
"""
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np

from scripts.export_question1 import read_case
from src.boundary import Environment
from src.data_io import ROOT, digest, load_config, save_json, verify_sources


def main():
    verify_sources()
    audit = json.loads((ROOT / "reports/q1_export_audit.json").read_text())
    if audit["status"] != "passed_Q1_local_acceptance_and_output_audit":
        raise ValueError("Run the Q1 export acceptance first")
    for relative, expected in audit["artifacts_sha256"].items():
        if digest(ROOT / relative) != expected:
            raise ValueError("Changed audited artifact: " + relative)
    cache, metadata = read_case("q1_final")
    config = load_config()
    env = Environment(config)
    radius = config["geometry"]["reference_radius_m"]
    initial = config["initial"]["moisture_dry_basis"]
    properties = config["property_sets"]["q1"]
    alpha = properties["k_W_m_K"]["offset"] / (
        properties["rho_kg_m3"]["offset"] * properties["cp_J_kg_K"]["offset"]
    )
    prefactor = properties["D_m2_s"]["prefactor"]
    exponent = properties["D_m2_s"]["moisture_exponent"]
    diffusion = prefactor * np.exp(-exponent / cache["moisture"])
    initial_diffusion = prefactor * np.exp(-exponent / initial)
    times = cache["time_s"]
    temp, moisture = cache["temperature"], cache["moisture"]
    air_t, air_c = env.value(times, "temperature"), env.value(times, "moisture")
    h = config["boundary"]["heat_transfer_W_m2_K"]
    hm = config["boundary"]["normalized_mass_transfer_m_s"]
    loss = initial - cache["moisture_mean"][-1]
    shells = []
    for inner in (0.010, 0.015):
        # Exact annular overlap weights, using the finite-volume cell averages.
        left = np.maximum(cache["faces_m"][:-1], inner)
        right = cache["faces_m"][1:]
        weights = np.maximum(right**2 - left**2, 0.0) / radius**2
        contribution = float(np.dot(weights, initial - cache["moisture_cells"][-1]))
        shells.append({
            "inner_radius_cm": inner * 100,
            "outer_radius_cm": radius * 100,
            "initial_dry_mass_fraction": float(weights.sum()),
            "loss_per_total_dry_mass_kg_kg": contribution,
            "fraction_of_total_water_loss": contribution / loss,
            "quadrature": "piecewise_constant_cell_average_with_exact_annular_overlap",
        })
    accumulated = np.trapezoid(diffusion, times, axis=0) / radius**2
    coarse = np.trapezoid(diffusion[::2], times[::2], axis=0) / radius**2
    report = {
        "scope": "Q1 model-conditional mechanism diagnostics; not experimental causality",
        "time_s": 1800,
        "source_fingerprint": metadata["source_fingerprint"],
        "cache_sha256": metadata["cache_sha256"],
        "script_sha256": digest(__file__),
        "effective_air_to_material_boundary_mapping": "C_eq = beta * Y; beta = 1",
        "scales": {
            "alpha_m2_s": alpha,
            "D_initial_m2_s": initial_diffusion,
            "alpha_over_D_initial": alpha / initial_diffusion,
            "radial_heat_scale_s": radius**2 / alpha,
            "radial_water_scale_s": radius**2 / initial_diffusion,
            "Fo_heat_1800": alpha * 1800 / radius**2,
            "Fo_water_initial_D_1800": initial_diffusion * 1800 / radius**2,
            "accumulated_local_Fo_water_centre": float(accumulated[0]),
            "accumulated_local_Fo_water_surface": float(accumulated[-1]),
            "accumulated_Fo_1s_vs_2s_quadrature_max_difference": float(np.max(np.abs(accumulated - coarse))),
            "accumulated_Fo_is_diagnostic_not_constant_coefficient_solution": True,
        },
        "temperature": {
            "centre_C": float(temp[-1, 0]),
            "surface_C": float(temp[-1, -1]),
            "mean_C": float(cache["temperature_mean"][-1]),
            "air_C": float(air_t[-1]),
            "air_minus_centre_C": float(air_t[-1] - temp[-1, 0]),
            "surface_minus_centre_C": float(temp[-1, -1] - temp[-1, 0]),
            "surface_outward_heat_flux_W_m2": float(h * (temp[-1, -1] - air_t[-1])),
        },
        "moisture": {
            "centre_dry_basis": float(moisture[-1, 0]),
            "surface_dry_basis": float(moisture[-1, -1]),
            "mean_dry_basis": float(cache["moisture_mean"][-1]),
            "centre_drop_from_initial": float(initial - moisture[-1, 0]),
            "loss_per_dry_mass_kg_kg": float(loss),
            "fraction_of_initial_water_removed": float(loss / initial),
            "surface_normalized_outward_flux_m_s": float(hm * (moisture[-1, -1] - air_c[-1])),
            "instantaneous_mean_rate_dry_basis_per_s": float(-2 * hm / radius * (moisture[-1, -1] - air_c[-1])),
            "surface_D_over_initial_D": float(diffusion[-1, -1] / initial_diffusion),
            "shell_loss_contributions": shells,
        },
    }
    assert 0 < loss < initial
    assert all(0 <= shell["fraction_of_total_water_loss"] <= 1 + 1e-10 for shell in shells)
    assert report["temperature"]["surface_outward_heat_flux_W_m2"] < 0
    assert report["moisture"]["instantaneous_mean_rate_dry_basis_per_s"] < 0
    save_json(ROOT / "reports/q1_mechanisms.json", report)
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
