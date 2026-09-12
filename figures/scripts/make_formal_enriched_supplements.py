#!/usr/bin/env python3
"""Generate the information-enriched diagnostic figures and register all figures."""

from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from figures.scripts import make_formal_enriched_figures as formal  # noqa: E402
from figures.scripts import make_supplement_figures as legacy  # noqa: E402
from src.data_io import digest, load_config, verify_sources  # noqa: E402


FIG_DIR = ROOT / "figures"
COLORS = formal.COLORS


CONTRACTS = {
    "fig00a_research_roadmap": {
            "caption": "四问统一建模与验证框架",
        "section": "问题重述",
        "core_conclusion": "四类输入汇入共享守恒内核，四问按真实初态关系分流，并由三类证据约束。",
        "archetype": "schematic-led composite",
        "panel_map": {"top": "输入", "middle": "共享内核与四问关系", "bottom": "验证证据"},
        "reviewer_risk": "只有Q2到Q3的金色箭头表示状态继续；Q1、Q2和Q4均由题给初态独立计算。",
        "source_data": ["reports/framework_reference_audit.md", "reports/model_contract.md", "reports/final_result_validation.json"],
    },
    "fig00b_model_geometry": {
        "caption": "几何边界与收缩坐标",
        "section": "模型假设与符号说明",
        "core_conclusion": "薄环控制体、侧面Robin边界和材料坐标共同闭合固定域与收缩域离散。",
        "archetype": "schematic-led composite",
        "panel_map": {"a": "圆柱几何与端面条件", "b": "环形控制体、通量和材料坐标映射"},
        "reviewer_risk": "端面封闭是主情景假设；暴露端面只作条件对照。",
        "source_data": ["config/model_config.json", "reports/model_contract.md", "reports/endface_assessment.md"],
    },
    "figS01_input_holdout": {
        "caption": "输入插值留出检验",
        "section": "附录：输入数据检验",
        "core_conclusion": "线性插值在三类输入的隔点留出中均优于持续值基线，PCHIP没有稳定优势。",
        "archetype": "asymmetric quantitative composite",
        "panel_map": {"a-c": "留出观测与插值重建", "d-f": "逐时刻残差、RMSE/MAE和最大残差", "g": "跨量纲归一化RMSE与相对基线改善"},
        "evidence_hierarchy": {"hero": "a-c", "diagnostic": "d-f", "decision": "g"},
        "reviewer_risk": "检验只支持观测结点之间的插值方案，不验证内部场或观测区间之外的延拓。",
        "source_data": ["original_appendix_1", "original_appendix_2", "reports/input_holdout_validation.json"],
    },
    "figS02_endface_assessment": {
        "caption": "端面假设条件对照",
        "section": "附录：端面假设审查",
        "core_conclusion": "端面暴露的局部影响随距端面距离快速衰减；中截面较稳健，但整根平均量和终点仍需条件化表述。",
        "archetype": "asymmetric mixed-modality figure",
        "panel_map": {"a-b": "3小时二维场差", "c": "端面影响的轴向衰减", "d": "整根平均含水率差", "e": "各输出相对显示单位的影响", "f": "配对端面效应在细化情景下的稳定性"},
        "evidence_hierarchy": {"hero": "a-c", "global_effect": "d-e", "robustness": "f"},
        "reviewer_risk": "暴露端面采用与侧面相同系数，是条件情景而非未知真实端面条件的误差上界。",
        "source_data": ["reports/endface_audit.json", "reports/endface_diagnostics.json", "results/cache/endface"],
    },
}


def _mae(observed, predicted) -> float:
    return float(np.mean(np.abs(np.asarray(predicted) - np.asarray(observed))))


def fig_input_holdout(cfg):
    env, radius = formal.load_measurements(cfg)
    series = [
        ("环境温度", "℃", env[:, 0], env[:, 1], COLORS["temperature"]),
        ("环境含湿量", "kg/kg", env[:, 0], env[:, 2], COLORS["moisture"]),
        ("观测半径", "cm", radius[:, 0], radius[:, 1], COLORS["radius"]),
    ]
    fig = plt.figure(figsize=(7.2, 7.2))
    gs = fig.add_gridspec(3, 6, height_ratios=[0.92, 0.92, 0.78], left=0.08, right=0.985, bottom=0.075, top=0.965, hspace=0.64, wspace=0.9)
    top = [fig.add_subplot(gs[0, 0:2]), fig.add_subplot(gs[0, 2:4]), fig.add_subplot(gs[0, 4:6])]
    bottom = [fig.add_subplot(gs[1, 0:2]), fig.add_subplot(gs[1, 2:4]), fig.add_subplot(gs[1, 4:6])]
    summary = fig.add_subplot(gs[2, :])
    source_rows = []
    rmse_matrix = []
    for j, (name, unit, times, values, color) in enumerate(series):
        held, observed, linear, pchip, persistence = legacy.holdout(times, values)
        time_h = times / 3600.0
        held_h = times[held] / 3600.0
        ax = top[j]
        ax.plot(time_h, values, color=COLORS["neutral_light"], lw=1.0, label="全部观测")
        ax.scatter(held_h, observed, s=18, color=color, edgecolor="white", lw=0.35, zorder=4, label="留出观测")
        ax.plot(held_h, linear, color=COLORS["gold"], lw=0.8, marker="x", ms=3.0, label="线性预测")
        ax.plot(held_h, pchip, color=COLORS["accent"], lw=0.7, ls="--", marker="+", ms=2.8, label="PCHIP预测")
        ax.set(xlabel="时间 / h", ylabel=f"{name} / {unit}")
        if j == 0:
            ax.legend(ncol=2, loc="lower right", fontsize=5.7)
        formal.grid(ax); formal.panel(ax, "abc"[j], f"{name}｜隔点留出重建")

        residual_sets = [("线性", linear - observed, COLORS["moisture"]), ("PCHIP", pchip - observed, COLORS["accent"]), ("持续值", persistence - observed, COLORS["neutral"])]
        ax = bottom[j]
        for label, residual, c in residual_sets:
            ax.plot(held_h, residual, lw=0.8, marker="o", ms=2.1, color=c, alpha=0.88, label=label)
        ax.axhline(0, color=COLORS["threshold"], lw=0.75, ls="--")
        lin_res = linear - observed
        worst = int(np.argmax(np.abs(lin_res)))
        ax.scatter(held_h[worst], lin_res[worst], s=34, marker="*", color=COLORS["gold"], edgecolor="white", lw=0.4, zorder=5)
        r_lin, r_pch, r_per = legacy.rmse(observed, linear), legacy.rmse(observed, pchip), legacy.rmse(observed, persistence)
        rmse_matrix.append((r_lin, r_pch, r_per))
        ax.text(0.02, 0.97, f"线性 RMSE={r_lin:.3g}\n线性 MAE={_mae(observed, linear):.3g}\n最大残差={abs(lin_res[worst]):.3g}", transform=ax.transAxes, va="top", fontsize=5.7, color=COLORS["neutral"], bbox={"fc": "white", "ec": "none", "alpha": 0.82, "pad": 1.2})
        ax.set(xlabel="留出时刻 / h", ylabel=f"预测残差 / {unit}")
        if j == 0:
            ax.legend(ncol=3, loc="lower right", fontsize=5.6)
        formal.grid(ax); formal.panel(ax, "def"[j], "残差轨迹与最不利留出点")
        for i, idx in enumerate(held):
            source_rows.append({"series": name, "unit": unit, "held_index": int(idx), "time_s": float(times[idx]), "observed": float(observed[i]), "linear": float(linear[i]), "pchip": float(pchip[i]), "persistence": float(persistence[i]), "linear_residual": float(linear[i] - observed[i]), "pchip_residual": float(pchip[i] - observed[i]), "persistence_residual": float(persistence[i] - observed[i])})

    ratios = np.asarray([[row[0] / row[2], row[1] / row[2], 1.0] for row in rmse_matrix])
    x = np.arange(3); w = 0.23
    method_colors = [COLORS["moisture"], COLORS["accent"], COLORS["neutral_light"]]
    method_edges = [COLORS["moisture"], COLORS["accent"], COLORS["neutral"]]
    for k, label in enumerate(["线性", "PCHIP", "持续值"]):
        bars = summary.bar(x + (k - 1) * w, ratios[:, k], w, color=method_colors[k], edgecolor=method_edges[k], lw=0.65, hatch="///" if k == 2 else "", label=label)
        for bar, value in zip(bars, ratios[:, k]):
            summary.text(bar.get_x() + bar.get_width() / 2, value + 0.025, f"{value:.2f}", ha="center", va="bottom", fontsize=5.7)
    summary.axhline(1, color=COLORS["threshold"], lw=0.85, ls="--", label="持续值基线")
    summary.set_xticks(x, ["环境温度", "环境含湿量", "观测半径"])
    summary.set(ylabel="RMSE / 持续值RMSE", ylim=(0, 1.16))
    summary.legend(ncol=4, loc="upper center")
    gains = 100 * (1 - ratios[:, 0])
    for xi, gain, pair in zip(x, gains, ratios[:, :2]):
        summary.text(xi, min(0.93, float(np.max(pair)) + 0.08), f"线性改善 {gain:.1f}%", ha="center", va="bottom", fontsize=5.9, color=COLORS["gold"])
    formal.grid(summary); formal.panel(summary, "g", "跨量纲归一化误差支持线性插值")

    source_path = FIG_DIR / "source_data" / "figS01_input_holdout.csv"
    source_path.parent.mkdir(exist_ok=True)
    with source_path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(source_rows[0]))
        writer.writeheader(); writer.writerows(source_rows)
    return formal.save_figure(fig, "figS01_input_holdout"), source_path


def _paired_effects(report: dict, group: str):
    order = ["both_space_axes", "radial", "axial", "time_and_tolerance"]
    labels = ["两轴细化", "径向细化", "轴向细化", "时步/容差"]
    values = []
    for kind in order:
        rows = [x for x in report["paired_refinements"] if x["group"] == group and x["kind"] == kind]
        if group == "q4" and kind == "radial":
            tight = [x for x in report["paired_refinements"] if x["group"] == group and x["kind"] == "radial_at_tight_time"]
            if tight:
                rows = tight
        values.append(abs(float(rows[-1]["fine_paired_effect_s"])) if rows else np.nan)
    return labels, np.asarray(values)


def fig_endface_assessment():
    report = json.loads((ROOT / "reports/endface_audit.json").read_text(encoding="utf-8"))
    diagnostics = json.loads((ROOT / "reports/endface_diagnostics.json").read_text(encoding="utf-8"))
    q23_case = next(x for x in report["cases"] if x["group"] == "q23" and x["nr"] == 64 and x["nz"] == 64 and x["rtol"] == 1e-9 and x["max_step_s"] == 60.0)
    q4_case = next(x for x in report["cases"] if x["group"] == "q4" and x["nr"] == 128 and x["nz"] == 64 and x["rtol"] == 1e-10 and x["max_step_s"] == 30.0)
    pairs = [legacy.load_pair("Q2/Q3", q23_case["exposed_label"], q23_case["insulated_label"], 64, 64), legacy.load_pair("Q4", q4_case["exposed_label"], q4_case["insulated_label"], 64, 128)]
    fig = plt.figure(figsize=(7.2, 7.0))
    gs = fig.add_gridspec(3, 4, height_ratios=[1.08, 0.9, 0.95], left=0.09, right=0.985, bottom=0.075, top=0.965, hspace=0.58, wspace=0.9)
    axes = [fig.add_subplot(gs[0, :2]), fig.add_subplot(gs[0, 2:]), fig.add_subplot(gs[1, :2]), fig.add_subplot(gs[1, 2:]), fig.add_subplot(gs[2, :2]), fig.add_subplot(gs[2, 2:])]
    source_rows = []
    cmap = mpl.colors.LinearSegmentedColormap.from_list("endface_delta_b", ["#30356F", "#7A82C2", "#D2D1EE", "#F7F6FC"])
    spatial = []
    for ax, pair, letter in zip(axes[:2], pairs, "ab"):
        r_cm, z_cm, diff = legacy.field_difference_at(pair, 10800.0)
        vmin = min(-0.05, float(np.min(diff)))
        mesh = ax.pcolormesh(r_cm, z_cm, diff, shading="nearest", cmap=cmap, vmin=vmin, vmax=0.0, rasterized=True)
        cbar = fig.colorbar(mesh, ax=ax, pad=0.014, fraction=0.043)
        cbar.set_label("ΔC = 暴露端面 - 封闭端面")
        cbar.outline.set_linewidth(0.55)
        ax.set(xlabel="距轴线 / cm", ylabel="距中截面 / cm")
        ax.axhline(0, color=COLORS["threshold"], lw=0.7)
        formal.panel(ax, letter, f"{pair['group']}在3 h的端面局部影响")
        spatial.append((pair["group"], z_cm, np.max(np.abs(diff), axis=1)))
        for iz, z in enumerate(z_cm):
            for ir, r in enumerate(r_cm):
                source_rows.append({"record": "field", "group": pair["group"], "time_s": 10800.0, "r_cm": float(r), "z_cm": float(z), "value": float(diff[iz, ir]), "metric": "delta_C"})

    ax = axes[2]
    half_length_cm = 12.5
    for (name, z_cm, envelope), color, marker in zip(spatial, [COLORS["moisture"], COLORS["temperature"]], ["o", "s"]):
        distance_from_end = half_length_cm - z_cm
        ax.semilogy(distance_from_end, envelope, color=color, lw=1.4, marker=marker, markevery=max(1, len(z_cm)//8), ms=2.8, label=name)
    ax.set(xlabel="距暴露端面距离 / cm", ylabel="径向最大 |ΔC|", xlim=(0, 12.5))
    formal.ascii_log_axis(ax)
    ax.yaxis.set_major_locator(mpl.ticker.LogLocator(base=10, subs=(1.0,), numticks=10))
    ax.legend(); formal.grid(ax, "both"); formal.panel(ax, "c", "局部影响由端面向中截面衰减")

    ax = axes[3]
    for pair, color, marker in zip(pairs, [COLORS["moisture"], COLORS["temperature"]], ["o", "s"]):
        times, ie, ii = legacy.common_series(pair)
        delta = pair["e"]["mean_C"][ie] - pair["i"]["mean_C"][ii]
        ax.plot(times / 3600.0, delta, color=color, lw=1.4, marker=marker, ms=2.2, label=pair["group"])
        for t, value in zip(times, delta):
            source_rows.append({"record": "mean_series", "group": pair["group"], "time_s": float(t), "r_cm": "", "z_cm": "", "value": float(value), "metric": "delta_mean_C"})
    ax.axhline(0, color=COLORS["threshold"], lw=0.75, ls="--")
    ax.set(xlabel="时间 / h", ylabel="整根平均含水率差 ΔCbar")
    ax.legend(); formal.grid(ax); formal.panel(ax, "d", "整根平均量对端面更敏感")

    categories = ["中截面温度", "中截面含水率", "3 h整根平均", "临界时间"]
    scales = np.asarray([5e-5, 5e-5, 5e-5, 0.18])
    q23_values = np.asarray([q23_case["fields"]["temperature"]["midplane_max_abs"], q23_case["fields"]["moisture"]["midplane_max_abs"], abs(next(x["mean_C_difference"] for x in q23_case["snapshots"] if x["time_s"] == 10800.0)), abs(q23_case["paired_event_difference_s"])]) / scales
    q4_values = np.asarray([q4_case["fields"]["temperature"]["midplane_max_abs"], q4_case["fields"]["moisture"]["midplane_max_abs"], abs(next(x["mean_C_difference"] for x in q4_case["snapshots"] if x["time_s"] == 10800.0)), abs(q4_case["paired_event_difference_s"])]) / scales
    ax = axes[4]; x = np.arange(len(categories))
    ax.scatter(x - 0.08, q23_values, color=COLORS["moisture"], s=30, label="Q2/Q3", edgecolor="white", lw=0.4, zorder=3)
    ax.scatter(x + 0.08, q4_values, color=COLORS["temperature"], marker="s", s=29, label="Q4", edgecolor="white", lw=0.4, zorder=3)
    for xi, a, b in zip(x, q23_values, q4_values):
        ax.plot([xi - 0.08, xi + 0.08], [a, b], color=COLORS["neutral_light"], lw=1.0, zorder=1)
    ax.axhline(1, color=COLORS["gold"], lw=0.9, ls="--", label="显示半单位")
    ax.set_yscale("log"); formal.ascii_log_axis(ax)
    ax.yaxis.set_major_locator(mpl.ticker.LogLocator(base=10, subs=(1.0,), numticks=8))
    ax.set_xticks(x, categories, rotation=18, ha="right")
    ax.set_ylabel("情景差绝对值 / 显示半单位")
    ax.legend(ncol=2, fontsize=5.8); formal.grid(ax); formal.panel(ax, "e", "端面影响取决于所考察输出")

    labels, v23 = _paired_effects(diagnostics, "q23")
    _, v4 = _paired_effects(diagnostics, "q4")
    ax = axes[5]; x = np.arange(len(labels))
    ax.plot(x, v23, color=COLORS["moisture"], marker="o", lw=1.4, ms=4, label="Q2/Q3")
    ax.plot(x, v4, color=COLORS["temperature"], marker="s", lw=1.4, ms=4, label="Q4")
    ax.set_yscale("log"); formal.ascii_log_axis(ax)
    ax.yaxis.set_major_locator(mpl.ticker.LogLocator(base=10, subs=(1.0,), numticks=8))
    ax.set_xticks(x, labels, rotation=18, ha="right")
    ax.set(ylabel="配对终点效应绝对值 / s", ylim=(8e-3, 15))
    ax.legend(ncol=1, loc="center left", bbox_to_anchor=(0.03, 0.45))
    ax.text(0.98, 0.60, "细化后仍稳定在各自数量级", transform=ax.transAxes, ha="right", fontsize=6.0, color=COLORS["gold"])
    formal.grid(ax); formal.panel(ax, "f", "配对端面效应的离散稳健性")

    for group, values in (("Q2/Q3", q23_values), ("Q4", q4_values)):
        for category, value in zip(categories, values):
            source_rows.append({"record": "normalized_magnitude", "group": group, "time_s": "", "r_cm": "", "z_cm": "", "value": float(value), "metric": category})
    source_path = FIG_DIR / "source_data" / "figS02_endface_assessment.csv"
    source_path.parent.mkdir(exist_ok=True)
    with source_path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["record", "group", "time_s", "r_cm", "z_cm", "value", "metric"])
        writer.writeheader(); writer.writerows(source_rows)
    return formal.save_figure(fig, "figS02_endface_assessment"), source_path


def register_bundle(source_paths, input_count):
    manifest_path = FIG_DIR / "figure_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    sources = {
        "fig00a_research_roadmap": [FIG_DIR / "fig00a_research_roadmap.drawio", FIG_DIR / "fig00a_research_roadmap.json"],
        "fig00b_model_geometry": [FIG_DIR / "fig00b_model_geometry.tex"],
        "figS01_input_holdout": [source_paths[0], Path(__file__)],
        "figS02_endface_assessment": [source_paths[1], Path(__file__)],
    }
    for stem in CONTRACTS:
        candidates = [FIG_DIR / f"{stem}.{ext}" for ext in ("svg", "pdf", "png", "tiff")]
        files = {p.suffix.lstrip("."): {"bytes": p.stat().st_size, "sha256": formal.sha256(p)} for p in candidates if p.exists()}
        if "pdf" not in files or "png" not in files:
            raise RuntimeError(f"Incomplete rendered figure: {stem}")
        manifest["figures"][stem] = {"files": files, "source_files": {str(p.relative_to(ROOT)): {"bytes": p.stat().st_size, "sha256": formal.sha256(p)} for p in sources[stem]}, "contract": CONTRACTS[stem]}
    manifest["paper_order"] = ["fig00a_research_roadmap", "fig00b_model_geometry", "fig01_inputs_scales", "fig02_q1_fields", "fig03_q2_coupling", "fig04_q3_endpoint", "fig05_q4_shrinkage", "fig06_numerical_validation", "fig07_sensitivity_risk", "figS01_input_holdout", "figS02_endface_assessment"]
    manifest["input_manifest_items_verified"] = input_count
    contact = legacy.build_contact_sheet()
    manifest["contact_sheet_complete"] = {"path": str(contact.relative_to(ROOT)), "sha256": formal.sha256(contact)}
    manifest["supplement_script_sha256"] = digest(__file__)
    manifest["status"] = "generated_information_enriched_pending_complete_visual_qa"
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return contact


def main():
    formal.style()
    inputs = verify_sources()
    cfg = load_config()
    holdout_outputs, holdout_source = fig_input_holdout(cfg)
    endface_outputs, endface_source = fig_endface_assessment()
    contact = register_bundle([holdout_source, endface_source], len(inputs["files"]))
    print(json.dumps({"diagnostic_figures": 2, "files": len(holdout_outputs) + len(endface_outputs), "contact_sheet": str(contact)}, ensure_ascii=False))


if __name__ == "__main__":
    main()
