#!/usr/bin/env python3
"""Generate the input-holdout and end-face diagnostic figures and register the full bundle."""

from __future__ import annotations

import csv
import hashlib
import json
import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image, ImageDraw
from scipy.interpolate import PchipInterpolator

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from make_all_figures import (  # noqa: E402
    COLORS,
    FIG_DIR,
    MOIST_CMAP,
    ascii_log_axis,
    grid,
    load_measurements,
    panel,
    save_figure,
    style,
)
from src.data_io import digest, load_config, verify_sources  # noqa: E402


CONTRACTS = {
    "fig00a_research_roadmap": {
        "caption": "整体技术路线",
        "section": "问题分析",
        "core_conclusion": "四问按物性、耦合、终点与移动域逐层扩展，并由同一套验证链约束。",
        "archetype": "schematic-led composite",
        "panel_map": {"main": "输入口径、四问递进、统一验证与结果交付"},
        "reviewer_risk": "箭头表示模型复杂度递进；Q2和Q4均从原始初态独立重算。",
        "source_data": ["reports/model_contract.md", "reports/final_result_validation.json"],
    },
    "fig00b_model_geometry": {
        "caption": "几何边界与收缩坐标",
        "section": "模型假设与符号说明",
        "core_conclusion": "一维主模型对应中截面径向有效情景，Q4通过材料坐标处理径向收缩，端面暴露仅作条件对照。",
        "archetype": "schematic-led composite",
        "panel_map": {"a": "圆柱几何、侧面Robin边界与端面口径", "b": "固定域和移动域坐标映射"},
        "reviewer_risk": "端面状态并非题面已知事实；二维暴露情景不替换正式一维结果。",
        "source_data": ["config/model_config.json", "reports/model_contract.md", "reports/endface_assessment.md"],
    },
    "figS01_input_holdout": {
        "caption": "输入插值留出检验",
        "section": "附录：输入数据检验",
        "core_conclusion": "线性插值在三组输入的隔点留出中均优于持续值基线，PCHIP没有稳定优势。",
        "archetype": "quantitative grid",
        "panel_map": {
            "a-c": "三组观测的隔点留出预测",
            "d-f": "线性、PCHIP与持续值基线的留出残差",
        },
        "reviewer_risk": "只验证观测结点间的输入插值，不验证药材内部温湿场，也不验证4小时后的环境延拓。",
        "source_data": ["original_appendix_1", "original_appendix_2", "reports/input_holdout_validation.json"],
    },
    "figS02_endface_assessment": {
        "caption": "端面假设条件对照",
        "section": "附录：端面假设审查",
        "core_conclusion": "暴露端面的影响集中在端部；中截面场差较小，但整根平均含水率与Q3终点仍可出现更明显差异。",
        "archetype": "asymmetric mixed-modality figure",
        "panel_map": {
            "a-b": "3小时端面暴露相对封闭情景的二维含水率差",
            "c": "整根平均含水率的情景差随时间变化",
            "d": "中截面、整根平均量和终点相对显示半单位的影响量级",
        },
        "reviewer_risk": "这是端面与侧面采用同一系数的条件情景，不是未知真实端面条件的误差上界。",
        "source_data": ["reports/endface_audit.json", "reports/endface_diagnostics.json", "results/cache/endface"],
    },
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def holdout(times: np.ndarray, values: np.ndarray):
    held = np.arange(5, len(times) - 1, 5)
    keep = np.ones(len(times), dtype=bool)
    keep[held] = False
    observed = values[held]
    linear = np.interp(times[held], times[keep], values[keep])
    pchip = PchipInterpolator(times[keep], values[keep], extrapolate=False)(times[held])
    persistence = values[held - 1]
    return held, observed, linear, pchip, persistence


def rmse(observed, predicted) -> float:
    return float(np.sqrt(np.mean((np.asarray(predicted) - np.asarray(observed)) ** 2)))


def fig_input_holdout(cfg):
    env, radius = load_measurements(cfg)
    series = [
        ("环境温度", "℃", env[:, 0], env[:, 1], COLORS["temperature"]),
        ("环境含湿量", "kg/kg", env[:, 0], env[:, 2], COLORS["moisture"]),
        ("观测半径", "cm", radius[:, 0], radius[:, 1], COLORS["radius"]),
    ]
    fig, axes = plt.subplots(2, 3, figsize=(7.2, 5.25))
    fig.subplots_adjust(left=0.085, right=0.985, top=0.94, bottom=0.105, wspace=0.38, hspace=0.48)
    source_rows = []
    for j, (name, unit, times, values, color) in enumerate(series):
        held, observed, linear, pchip, persistence = holdout(times, values)
        time_h = times / 3600
        held_h = times[held] / 3600
        top = axes[0, j]
        top.plot(time_h, values, color="#B9B9B9", lw=0.85, label="全部观测")
        top.scatter(held_h, observed, s=18, color=color, edgecolor="white", lw=0.35, zorder=3, label="留出观测")
        top.scatter(held_h, linear, s=18, marker="x", color="#E58C46", lw=0.9, zorder=4, label="线性预测")
        top.set(xlabel="时间 / h", ylabel=f"{name} / {unit}")
        grid(top)
        panel(top, "abc"[j], f"{name}的隔点留出")

        bottom = axes[1, j]
        residual_sets = [
            ("线性", linear - observed, COLORS["moisture"]),
            ("PCHIP", pchip - observed, COLORS["accent"]),
            ("持续值", persistence - observed, COLORS["neutral"]),
        ]
        for label, residual, c in residual_sets:
            bottom.plot(held_h, residual, lw=0.8, marker="o", ms=2.3, color=c, alpha=0.85, label=label)
        bottom.axhline(0, color="#2C2C2C", lw=0.75, ls="--")
        bottom.set(xlabel="留出时刻 / h", ylabel=f"预测残差 / {unit}")
        grid(bottom)
        panel(bottom, "def"[j], "三种插值/基线的残差")
        r_lin, r_pch, r_per = rmse(observed, linear), rmse(observed, pchip), rmse(observed, persistence)
        rmse_note = bottom.text(
            0.02,
            0.97,
            f"RMSE\n线性 {r_lin:.4g}\nPCHIP {r_pch:.4g}\n持续值 {r_per:.4g}",
            transform=bottom.transAxes,
            va="top",
            fontsize=5.9,
            color="#333333",
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.82, "pad": 1.5},
        )
        rmse_note.set_in_layout(False)
        if j == 0:
            top.legend(loc="lower right", fontsize=6.0)
            bottom.legend(loc="lower right", ncol=3, fontsize=5.9)
        for i, idx in enumerate(held):
            source_rows.append(
                {
                    "series": name,
                    "unit": unit,
                    "held_index": int(idx),
                    "time_s": float(times[idx]),
                    "observed": float(observed[i]),
                    "linear": float(linear[i]),
                    "pchip": float(pchip[i]),
                    "persistence": float(persistence[i]),
                    "linear_residual": float(linear[i] - observed[i]),
                    "pchip_residual": float(pchip[i] - observed[i]),
                    "persistence_residual": float(persistence[i] - observed[i]),
                }
            )
    fig.text(0.5, 0.018, "注：留出对象仅为环境与半径输入结点，不是药材内部温湿响应。", ha="center", fontsize=6.6, color="#555555")
    source_path = FIG_DIR / "source_data" / "figS01_input_holdout.csv"
    source_path.parent.mkdir(exist_ok=True)
    with source_path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(source_rows[0]))
        writer.writeheader()
        writer.writerows(source_rows)
    return save_figure(fig, "figS01_input_holdout"), source_path


def load_pair(group: str, exposed_label: str, insulated_label: str, nz: int, nr: int):
    exposed = dict(np.load(ROOT / "results/cache/endface" / f"{exposed_label}.npz"))
    insulated = dict(np.load(ROOT / "results/cache/endface" / f"{insulated_label}.npz"))
    return {
        "group": group,
        "e": exposed,
        "i": insulated,
        "nz": nz,
        "nr": nr,
    }


def moisture_cells(cache, nz, nr):
    return cache["state"][:, 1::2].reshape(len(cache["time_s"]), nz, nr)


def common_series(pair):
    e, i = pair["e"], pair["i"]
    times, ie, ii = np.intersect1d(e["time_s"], i["time_s"], return_indices=True)
    return times, ie, ii


def field_difference_at(pair, target_s: float):
    e, i = pair["e"], pair["i"]
    ie = int(np.argmin(np.abs(e["time_s"] - target_s)))
    ii = int(np.argmin(np.abs(i["time_s"] - target_s)))
    if abs(e["time_s"][ie] - target_s) > 1e-6 or abs(i["time_s"][ii] - target_s) > 1e-6:
        raise RuntimeError(f"Missing target time {target_s:g} s for {pair['group']}")
    ce = moisture_cells(e, pair["nz"], pair["nr"])[ie]
    ci = moisture_cells(i, 3, pair["nr"])[ii, 0]
    radius_cm = e["xi_centres"] * e["radius_m"][ie] * 100
    z_cm = e["z_centres_m"] * 100
    return radius_cm, z_cm, ce - ci[None, :]


def fig_endface_assessment():
    report = json.loads((ROOT / "reports/endface_audit.json").read_text(encoding="utf-8"))
    q23_case = next(
        x for x in report["cases"]
        if x["group"] == "q23" and x["nr"] == 64 and x["nz"] == 64 and x["rtol"] == 1e-9 and x["max_step_s"] == 60.0
    )
    q4_case = next(
        x for x in report["cases"]
        if x["group"] == "q4" and x["nr"] == 128 and x["nz"] == 64 and x["rtol"] == 1e-10 and x["max_step_s"] == 30.0
    )
    pairs = [
        load_pair("Q2/Q3", q23_case["exposed_label"], q23_case["insulated_label"], 64, 64),
        load_pair("Q4", q4_case["exposed_label"], q4_case["insulated_label"], 64, 128),
    ]
    fig = plt.figure(figsize=(7.2, 5.55), layout="constrained")
    gs = fig.add_gridspec(2, 2, height_ratios=[1.05, 1])
    axes = [fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1]), fig.add_subplot(gs[1, 0]), fig.add_subplot(gs[1, 1])]
    source_rows = []
    cmap = mpl.colors.LinearSegmentedColormap.from_list("endface_delta", ["#08306B", "#6BAED6", "#F7FBFF"])
    for ax, pair, letter in zip(axes[:2], pairs, "ab"):
        r_cm, z_cm, diff = field_difference_at(pair, 10800.0)
        vmax = 0.0
        vmin = min(-0.05, float(np.min(diff)))
        mesh = ax.pcolormesh(r_cm, z_cm, diff, shading="nearest", cmap=cmap, vmin=vmin, vmax=vmax, rasterized=True)
        cbar = fig.colorbar(mesh, ax=ax, pad=0.018, fraction=0.047)
        cbar.set_label(r"$\Delta C=C_{\rm exposed}-C_{\rm insulated}$")
        cbar.outline.set_linewidth(0.6)
        ax.set(xlabel="距轴线 / cm", ylabel="距中截面 / cm")
        ax.axhline(0, color="#222222", lw=0.7)
        panel(ax, letter, f"{pair['group']}在3 h的端面影响")
        for iz, z in enumerate(z_cm):
            for ir, r in enumerate(r_cm):
                source_rows.append({"record": "field", "group": pair["group"], "time_s": 10800.0, "r_cm": float(r), "z_cm": float(z), "value": float(diff[iz, ir]), "metric": "delta_C"})

    ax = axes[2]
    for pair, color in zip(pairs, [COLORS["moisture"], COLORS["temperature"]]):
        times, ie, ii = common_series(pair)
        delta = pair["e"]["mean_C"][ie] - pair["i"]["mean_C"][ii]
        ax.plot(times / 3600, delta, color=color, lw=1.4, marker="o", ms=2.3, label=pair["group"])
        for t, value in zip(times, delta):
            source_rows.append({"record": "mean_series", "group": pair["group"], "time_s": float(t), "r_cm": "", "z_cm": "", "value": float(value), "metric": "delta_mean_C"})
    ax.axhline(0, color="#222222", lw=0.75, ls="--")
    ax.set(xlabel="时间 / h", ylabel=r"整根平均含水率差 $\Delta\bar C$ / (kg/kg)")
    ax.legend()
    grid(ax)
    panel(ax, "c", "整根平均失水对端面更敏感")

    ax = axes[3]
    categories = ["中截面温度", "中截面含水率", "3 h整根平均", "临界时间"]
    scales = np.array([5e-5, 5e-5, 5e-5, 0.18])
    q23_values = np.array([
        q23_case["fields"]["temperature"]["midplane_max_abs"],
        q23_case["fields"]["moisture"]["midplane_max_abs"],
        abs(next(x["mean_C_difference"] for x in q23_case["snapshots"] if x["time_s"] == 10800.0)),
        abs(q23_case["paired_event_difference_s"]),
    ]) / scales
    q4_values = np.array([
        q4_case["fields"]["temperature"]["midplane_max_abs"],
        q4_case["fields"]["moisture"]["midplane_max_abs"],
        abs(next(x["mean_C_difference"] for x in q4_case["snapshots"] if x["time_s"] == 10800.0)),
        abs(q4_case["paired_event_difference_s"]),
    ]) / scales
    x = np.arange(len(categories))
    ax.scatter(x - 0.08, q23_values, color=COLORS["moisture"], s=28, label="Q2/Q3", zorder=3)
    ax.scatter(x + 0.08, q4_values, color=COLORS["temperature"], marker="s", s=26, label="Q4", zorder=3)
    for xi, a, b in zip(x, q23_values, q4_values):
        ax.plot([xi - 0.08, xi + 0.08], [a, b], color="#D1D1D1", lw=0.8, zorder=1)
    ax.axhline(1, color=COLORS["threshold"], lw=0.9, ls="--", label="显示半单位")
    ax.set_yscale("log")
    ascii_log_axis(ax)
    ax.set_xticks(x, categories, rotation=15, ha="right")
    ax.set_ylabel("情景差绝对值 / 显示半单位")
    ax.legend(ncol=2, fontsize=6.0)
    grid(ax)
    panel(ax, "d", "影响随输出指标而异")
    ax.text(0.98, 0.04, "临界时间差均为提前", transform=ax.transAxes, ha="right", va="bottom", fontsize=6.3, color="#555555")
    for group, values in (("Q2/Q3", q23_values), ("Q4", q4_values)):
        for category, value in zip(categories, values):
            source_rows.append({"record": "normalized_magnitude", "group": group, "time_s": "", "r_cm": "", "z_cm": "", "value": float(value), "metric": category})

    source_path = FIG_DIR / "source_data" / "figS02_endface_assessment.csv"
    source_path.parent.mkdir(exist_ok=True)
    with source_path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["record", "group", "time_s", "r_cm", "z_cm", "value", "metric"])
        writer.writeheader()
        writer.writerows(source_rows)
    return save_figure(fig, "figS02_endface_assessment"), source_path


def build_contact_sheet():
    ordered = [
        ("fig00a_research_roadmap", "Fig. 1"),
        ("fig00b_model_geometry", "Fig. 2"),
        ("fig01_inputs_scales", "Fig. 3"),
        ("fig02_q1_fields", "Fig. 4"),
        ("fig03_q2_coupling", "Fig. 5"),
        ("fig04_q3_endpoint", "Fig. 6"),
        ("fig05_q4_shrinkage", "Fig. 7"),
        ("fig06_numerical_validation", "Fig. 8"),
        ("fig07_sensitivity_risk", "Fig. 9"),
        ("figS01_input_holdout", "Fig. S1"),
        ("figS02_endface_assessment", "Fig. S2"),
    ]
    width, thumb_w, thumb_h = 1500, 720, 520
    rows = (len(ordered) + 1) // 2
    canvas = Image.new("RGB", (width, rows * (thumb_h + 55) + 20), "white")
    draw = ImageDraw.Draw(canvas)
    for index, (stem, label) in enumerate(ordered):
        image = Image.open(FIG_DIR / f"{stem}.png").convert("RGB")
        image.thumbnail((thumb_w, thumb_h))
        x = 20 + (index % 2) * 740
        y = 35 + (index // 2) * (thumb_h + 55)
        canvas.paste(image, (x + (thumb_w - image.width) // 2, y))
        draw.text((x, 10 + (index // 2) * (thumb_h + 55)), label, fill="black")
    path = FIG_DIR / "contact_sheet_complete.png"
    canvas.save(path, dpi=(180, 180))
    return path


def register_bundle(data_outputs, source_paths, input_count):
    manifest_path = FIG_DIR / "figure_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    new_stems = ["fig00a_research_roadmap", "fig00b_model_geometry", "figS01_input_holdout", "figS02_endface_assessment"]
    source_map = {
        "fig00a_research_roadmap": [FIG_DIR / "fig00a_research_roadmap.drawio", FIG_DIR / "fig00a_research_roadmap.json"],
        "fig00b_model_geometry": [FIG_DIR / "fig00b_model_geometry.tex"],
        "figS01_input_holdout": [source_paths[0], Path(__file__)],
        "figS02_endface_assessment": [source_paths[1], Path(__file__)],
    }
    for stem in new_stems:
        candidates = [FIG_DIR / f"{stem}.{ext}" for ext in ("svg", "pdf", "png", "tiff")]
        files = {p.suffix.lstrip("."): {"bytes": p.stat().st_size, "sha256": sha256(p)} for p in candidates if p.exists()}
        if "pdf" not in files or "png" not in files:
            raise RuntimeError(f"Incomplete rendered figure: {stem}")
        manifest["figures"][stem] = {
            "files": files,
            "source_files": {str(p.relative_to(ROOT)): {"bytes": p.stat().st_size, "sha256": sha256(p)} for p in source_map[stem]},
            "contract": CONTRACTS[stem],
        }
    manifest["status"] = "generated_pending_complete_visual_qa"
    manifest["backend"] = "Python/matplotlib for quantitative figures; draw.io for roadmap; TikZ/XeLaTeX for geometry"
    manifest["paper_order"] = [
        "fig00a_research_roadmap", "fig00b_model_geometry", "fig01_inputs_scales", "fig02_q1_fields",
        "fig03_q2_coupling", "fig04_q3_endpoint", "fig05_q4_shrinkage", "fig06_numerical_validation",
        "fig07_sensitivity_risk", "figS01_input_holdout", "figS02_endface_assessment",
    ]
    manifest["input_manifest_items_verified"] = input_count
    contact = build_contact_sheet()
    manifest["contact_sheet_complete"] = {"path": str(contact.relative_to(ROOT)), "sha256": sha256(contact)}
    manifest["supplement_script_sha256"] = digest(__file__)
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return contact


def main():
    style()
    inputs = verify_sources()
    cfg = load_config()
    holdout_outputs, holdout_source = fig_input_holdout(cfg)
    endface_outputs, endface_source = fig_endface_assessment()
    contact = register_bundle([holdout_outputs, endface_outputs], [holdout_source, endface_source], len(inputs["files"]))
    print(json.dumps({"new_figures": 4, "new_quantitative_files": len(holdout_outputs) + len(endface_outputs), "contact_sheet": str(contact)}, ensure_ascii=False))


if __name__ == "__main__":
    main()
