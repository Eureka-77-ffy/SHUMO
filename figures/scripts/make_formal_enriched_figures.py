#!/usr/bin/env python3
"""Generate the information-enriched B-palette figures used by the paper.

Figure contract
---------------
Core conclusion:
    The four tasks form one traceable thermal-moisture argument: the boundary
    input stabilises early, internal moisture transport controls the long tail,
    the stopping event must be defined over the full domain, and shrinkage
    changes both transport length and surface exchange.
Archetype:
    Asymmetric quantitative composites with one or two hero panels.
Target/output:
    CUMCM A4 manuscript, 182.9 mm wide; editable SVG/PDF plus 600 dpi PNG/TIFF.
Backend:
    Python/matplotlib only for every quantitative render and preview.
Statistics:
    Deterministic PDE scenarios; no replicates, confidence intervals or
    probability interpretation are claimed.
Image integrity:
    Heatmaps use the accepted full-precision caches without smoothing. Display
    sampling is temporal only. The Q4 physical-domain exterior remains white.
Reviewer risk:
    Scale, sensitivity, ablation and factorial panels are model diagnostics,
    not experimental causal identification or physical error bounds.
"""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import LogLocator

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from figures.scripts import make_all_figures as legacy  # noqa: E402
from figures.scripts import make_enriched_layout_previews as accepted  # noqa: E402
from src.data_io import digest, load_config, verify_sources  # noqa: E402


FIG_DIR = ROOT / "figures"
FORMATS = ("svg", "pdf", "png", "tiff")
COLORS = accepted.COLORS
TEMP_CMAP = accepted.TEMP_CMAP
MOIST_CMAP = accepted.MOIST_CMAP


CONTRACTS = {
    "fig01_inputs_scales": {
        "caption": "环境驱动与特征尺度",
        "section": "总体分析",
        "core_conclusion": "环境输入在4小时内趋稳，但水分扩散、收缩和终点搜索持续至50小时以上，分布参数模型与长期边界闭合均不可省略。",
        "archetype": "asymmetric mixed-modality figure",
        "panel_map": {
            "a": "环境温度和含湿量观测及平台区间",
            "b": "半径观测、Q4临界半径和观测支持范围",
            "c": "四问计算时域与输入观测覆盖关系",
            "d": "三套物性的初始热、质扩散尺度",
            "e": "初始热、质Biot数与集总近似判据",
        },
        "evidence_hierarchy": {"hero": "a-c", "model_choice": "d-e"},
        "reviewer_risk": "特征时间与Biot数仅用于量级判断和模型选型，不替代数值终点。",
        "source_data": ["original_appendix_1", "original_appendix_2", "reports/scales_snapshot.csv"],
    },
    "fig02_q1_fields": {
        "caption": "问题一温湿场演化",
        "section": "问题一",
        "core_conclusion": "前30分钟热响应已深入内部，而水分损失仍受表层梯度控制，热、质过程具有清晰的时空分离。",
        "archetype": "asymmetric quantitative composite",
        "panel_map": {
            "a": "温度时空场（主证据）",
            "b": "含水率时空场（主证据）",
            "c": "典型温度径向剖面",
            "d": "典型含水率径向剖面",
            "e": "体积平均温度与含水率响应",
            "f": "中心—表面热湿梯度的同步演化",
        },
        "evidence_hierarchy": {"hero": "a-b", "spatial": "c-d", "summary": "e-f"},
        "reviewer_risk": "热图、剖面和统计量来自同一未舍入解；平均量不替代局部梯度。",
        "source_data": ["results/final/q1_full_precision.npz", "original_appendix_1"],
    },
    "fig03_q2_coupling": {
        **accepted.CONTRACTS["fig03_q2_enriched"],
        "caption": "问题二耦合响应",
        "section": "问题二",
    },
    "fig04_q3_endpoint": {
        **accepted.CONTRACTS["fig04_q3_enriched"],
        "caption": "问题三全域干燥终点",
        "section": "问题三",
    },
    "fig05_q4_shrinkage": {
        **accepted.CONTRACTS["fig05_q4_enriched"],
        "caption": "问题四收缩域效应",
        "section": "问题四",
    },
    "fig06_numerical_validation": {
        "caption": "数值收敛与独立验证",
        "section": "模型检验",
        "core_conclusion": "空间加密、时间设置、独立有限元和累计守恒从四条独立证据链支持题定数值，而这些差异远小于结构情景差异。",
        "archetype": "quantitative grid",
        "panel_map": {
            "a-b": "Q2温度和水分的空间收敛",
            "c": "Q3/Q4临界时间空间收敛",
            "d": "时间步、容差和积分器对终点的影响",
            "e": "四问与独立有限元的场值差异",
            "f": "四问累计守恒残差",
        },
        "evidence_hierarchy": {"hero": "a-d", "independent_check": "e", "identity_check": "f"},
        "reviewer_risk": "图示差异是观测到的数值差，不是严格误差上界或实验精度。",
        "source_data": ["reports/q2_convergence.json", "reports/q3_convergence.json", "reports/q4_convergence.json", "reports/final_result_validation.json"],
    },
    "fig07_sensitivity_risk": {
        "caption": "终点灵敏度与风险层级",
        "section": "灵敏度分析",
        "core_conclusion": "两问终点都由扩散率尺度主导，传质次之、换热几乎不敏感；边界与结构情景远大于纯数值误差。",
        "archetype": "asymmetric quantitative composite",
        "panel_map": {
            "a-b": "Q3/Q4正负10%单因素情景及非对称响应",
            "c": "有符号无量纲灵敏度及统一排序",
            "d": "数值、输入、参数和结构风险的跨数量级比较",
        },
        "evidence_hierarchy": {"hero": "a-c", "risk_context": "d"},
        "reviewer_risk": "局部单因素扰动和情景包络均未赋予概率分布，不解释为全局敏感性或置信区间。",
        "source_data": ["reports/q3_sensitivity.json", "reports/q4_sensitivity.json", "reports/endface_diagnostics.json", "reports/final_result_validation.json"],
    },
}


def style() -> None:
    accepted.configure_style()


def panel(ax, label: str, title: str) -> None:
    accepted.panel(ax, label, title)


def grid(ax, axis: str = "y") -> None:
    accepted.subtle_grid(ax, axis)


def ascii_log_axis(ax, axis: str = "y") -> None:
    legacy.ascii_log_axis(ax, axis)


def heatmap(*args, **kwargs):
    return accepted.heatmap(*args, **kwargs)


def save_figure(fig, stem: str):
    accepted.OUT = FIG_DIR
    return accepted.save_figure(fig, stem)


def load_measurements(cfg):
    return legacy.load_measurements(cfg)


def json_load(relative: str):
    return json.loads((ROOT / relative).read_text(encoding="utf-8"))


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def fig01(cfg, env: np.ndarray, radius: np.ndarray, final_validation: dict):
    fig = plt.figure(figsize=(7.2, 5.7))
    gs = fig.add_gridspec(2, 6, height_ratios=[1.05, 1.0], left=0.078, right=0.985, bottom=0.09, top=0.965, hspace=0.55, wspace=1.05)
    ax_a = fig.add_subplot(gs[0, :3])
    ax_b = fig.add_subplot(gs[0, 3:])
    ax_c = fig.add_subplot(gs[1, :2])
    ax_d = fig.add_subplot(gs[1, 2:4])
    ax_e = fig.add_subplot(gs[1, 4:])

    t_h = env[:, 0] / 3600.0
    ax_a.plot(t_h, env[:, 1], color=COLORS["temperature"], lw=1.55, label="环境温度")
    ax_a.scatter(t_h[::12], env[::12, 1], s=9, color=COLORS["temperature"], edgecolor="white", lw=0.35, zorder=3)
    ax_a.set(xlabel="时间 / h", ylabel="温度 / ℃", xlim=(0, 4), ylim=(27, 52))
    ax_a2 = ax_a.twinx()
    ax_a2.plot(t_h, env[:, 2], color=COLORS["moisture"], lw=1.4, label="环境含湿量")
    ax_a2.scatter(t_h[::12], env[::12, 2], s=8, color=COLORS["moisture"], edgecolor="white", lw=0.3, zorder=3)
    ax_a2.set_ylabel("环境含湿量 / (kg/kg)", color=COLORS["moisture"])
    ax_a2.tick_params(axis="y", colors=COLORS["moisture"])
    accepted.phase_bands(ax_a, [(0, 0.5, "快速升温"), (0.5, 1.5, "渐近追赶"), (1.5, 4, "平台波动")])
    handles = [mpl.lines.Line2D([], [], color=COLORS["temperature"], lw=1.5), mpl.lines.Line2D([], [], color=COLORS["moisture"], lw=1.4)]
    ax_a.legend(handles, ["环境温度", "环境含湿量"], ncol=2, loc="lower right")
    ax_a.text(0.98, 0.20, "3—4 h均值用于长期平台闭合", transform=ax_a.transAxes, ha="right", va="bottom", fontsize=6.2, color=COLORS["neutral"])
    grid(ax_a)
    panel(ax_a, "a", "环境边界由快速变化转入平台区")

    rt_h, r_cm = radius[:, 0] / 3600.0, radius[:, 1]
    ax_b.plot(rt_h, r_cm, color=COLORS["radius"], lw=1.5)
    ax_b.scatter(rt_h[::4], r_cm[::4], s=9, color=COLORS["radius"], edgecolor="white", lw=0.35, zorder=3, label="半径观测")
    q4_h = float(final_validation["critical_results"]["Q4"]["critical_time_h"])
    q4_r = float(np.interp(q4_h * 3600.0, radius[:, 0], radius[:, 1]))
    ax_b.scatter(q4_h, q4_r, s=48, marker="*", color=COLORS["gold"], edgecolor="white", lw=0.45, zorder=4)
    ax_b.axvline(q4_h, color=COLORS["gold"], lw=0.85, ls=":")
    ax_b.annotate(f"Q4终点\n{q4_h:.2f} h，{q4_r:.3f} cm", (q4_h, q4_r), xytext=(10, 18), textcoords="offset points", color=COLORS["gold"], arrowprops={"arrowstyle": "-", "lw": 0.65, "color": COLORS["gold"]})
    ax_b.set(xlabel="时间 / h", ylabel="", xlim=(0, 72), ylim=(1.14, 2.05))
    ax_b.text(0.015, 0.97, "$R$ / cm", transform=ax_b.transAxes, ha="left", va="top", fontsize=6.4, color=COLORS["radius"])
    ax_b.text(0.98, 0.94, "半径观测覆盖0—72 h，无需外推", transform=ax_b.transAxes, ha="right", va="top", fontsize=6.2, color=COLORS["neutral"])
    grid(ax_b)
    panel(ax_b, "b", "观测收缩覆盖问题四计算时域")

    q3_h = float(final_validation["critical_results"]["Q3"]["critical_time_h"])
    labels = ["环境观测", "半径观测", "Q1", "Q2", "Q4", "Q3"]
    ends = [4.0, 72.0, 0.5, 3.0, q4_h, q3_h]
    colors = [COLORS["temperature"], COLORS["radius"], COLORS["temperature_light"], COLORS["moisture_light"], COLORS["radius"], COLORS["moisture"]]
    y = np.arange(len(labels))[::-1]
    for yi, end, color in zip(y, ends, colors):
        ax_c.hlines(yi, 0, end, color=color, lw=5.0, alpha=0.9)
        ax_c.scatter(end, yi, s=18, color=color, edgecolor="white", lw=0.4, zorder=3)
        ax_c.text(end + 1.1, yi, f"{end:.2f} h" if end >= 10 else f"{end:g} h", va="center", fontsize=5.8, color=color)
    ax_c.set_yticks(y, labels)
    ax_c.set(xlabel="观测或计算时域 / h", xlim=(0, 81), ylim=(-0.6, 5.6))
    grid(ax_c, "x")
    panel(ax_c, "c", "输入支撑与四问时域")

    scale = np.genfromtxt(ROOT / "reports/scales_snapshot.csv", delimiter=",", names=True, dtype=None, encoding="utf-8")
    rows = [row for row in scale if float(row["C"]) == 2.55 and float(row["T_C"]) == 28.0]
    groups = ["Q1", "Q2/Q3", "Q4"]
    heat_h = np.asarray([float(row["tau_T_min"]) / 60.0 for row in rows])
    water_h = np.asarray([float(row["tau_C_h"]) for row in rows])
    x = np.arange(3)
    w = 0.33
    b1 = ax_d.bar(x - w / 2, heat_h, w, color=COLORS["temperature_light"], edgecolor=COLORS["temperature"], lw=0.65, hatch="///", label="热扩散")
    b2 = ax_d.bar(x + w / 2, water_h, w, color=COLORS["moisture_light"], edgecolor=COLORS["moisture"], lw=0.65, label="水分扩散")
    for bars, vals in ((b1, heat_h), (b2, water_h)):
        for bar, value in zip(bars, vals):
            ax_d.text(bar.get_x() + bar.get_width() / 2, value * 1.12, f"{value:.1f}", ha="center", va="bottom", fontsize=5.7)
    ax_d.set_yscale("log")
    ascii_log_axis(ax_d)
    ax_d.set_xticks(x, groups)
    ax_d.set(xlabel="物性体系", ylabel="初态径向特征时间 / h", ylim=(0.35, 190))
    ax_d.legend(ncol=2, loc="upper left")
    grid(ax_d)
    panel(ax_d, "d", "热—质特征时间分离")

    bi_t = np.asarray([float(row["Bi_T_R"]) for row in rows])
    bi_c = np.asarray([float(row["Bi_C_R"]) for row in rows])
    b1 = ax_e.bar(x - w / 2, bi_t, w, color=COLORS["temperature_light"], edgecolor=COLORS["temperature"], lw=0.65, hatch="///", label="BiT")
    b2 = ax_e.bar(x + w / 2, bi_c, w, color=COLORS["moisture_light"], edgecolor=COLORS["moisture"], lw=0.65, label="BiC")
    ax_e.axhspan(0.06, 0.1, color="#EEECE7", alpha=0.8, lw=0)
    ax_e.axhline(0.1, color=COLORS["threshold"], lw=0.8, ls="--")
    ax_e.text(0.03, 0.06, "集总近似常用区", transform=ax_e.transAxes, fontsize=5.7, color=COLORS["neutral"])
    ax_e.set_yscale("log")
    ascii_log_axis(ax_e)
    ax_e.set_xticks(x, groups)
    ax_e.set(xlabel="物性体系", ylabel="初始Biot数", ylim=(0.06, 28))
    ax_e.legend(ncol=2, loc="upper left")
    grid(ax_e)
    panel(ax_e, "e", "径向分布不可忽略")
    return save_figure(fig, "fig01_inputs_scales")


def fig02(q1, env: np.ndarray):
    fig = plt.figure(figsize=(7.2, 6.85))
    gs = fig.add_gridspec(3, 4, height_ratios=[1.15, 0.95, 0.9], left=0.08, right=0.985, bottom=0.075, top=0.965, hspace=0.56, wspace=0.86)
    axes = [fig.add_subplot(gs[0, :2]), fig.add_subplot(gs[0, 2:]), fig.add_subplot(gs[1, :2]), fig.add_subplot(gs[1, 2:]), fig.add_subplot(gs[2, :2]), fig.add_subplot(gs[2, 2:])]
    time_min = q1["time_s"] / 60.0
    radius_cm = q1["radius_m"] * 100.0
    heatmap(axes[0], time_min, radius_cm, q1["temperature"], TEMP_CMAP, "温度 / ℃")
    axes[0].set(xlabel="时间 / min", ylabel="距中心距离 / cm", xlim=(0, 30))
    axes[0].axvline(15, color="white", lw=0.75, ls="--")
    panel(axes[0], "a", "温度时空场｜热影响持续向内部推进")
    heatmap(axes[1], time_min, radius_cm, q1["moisture"], MOIST_CMAP, "干基含水率 / (kg/kg)")
    axes[1].set(xlabel="时间 / min", ylabel="距中心距离 / cm", xlim=(0, 30))
    axes[1].axvline(15, color="white", lw=0.75, ls="--")
    panel(axes[1], "b", "含水率时空场｜失水集中于表层")

    selected = [0, 5, 15, 30]
    cmap_t = mpl.colors.LinearSegmentedColormap.from_list("q1_t", ["#F3D7B5", COLORS["temperature"]])
    cmap_c = mpl.colors.LinearSegmentedColormap.from_list("q1_c", ["#DAD9EE", COLORS["moisture"]])
    for minute, ct, cc in zip(selected, cmap_t(np.linspace(0.12, 0.95, len(selected))), cmap_c(np.linspace(0.12, 0.95, len(selected)))):
        i = int(np.argmin(np.abs(time_min - minute)))
        axes[2].plot(radius_cm, q1["temperature"][i], color=ct, lw=1.55 if minute == 30 else 1.15, label=f"{minute} min")
        axes[3].plot(radius_cm, q1["moisture"][i], color=cc, lw=1.55 if minute == 30 else 1.15, label=f"{minute} min")
    axes[2].set(xlabel="距中心距离 / cm", ylabel="温度 / ℃", xlim=(0, 2))
    axes[3].set(xlabel="距中心距离 / cm", ylabel="干基含水率 / (kg/kg)", xlim=(0, 2))
    axes[2].legend(ncol=2, loc="upper left")
    axes[3].legend(ncol=2, loc="lower left")
    grid(axes[2]); grid(axes[3])
    panel(axes[2], "c", "温度剖面｜中心响应逐步建立")
    panel(axes[3], "d", "含水率剖面｜表面梯度保持陡峭")

    t_plot = time_min[::10]
    ax_e = axes[4]
    ax_e.plot(t_plot, q1["temperature_mean"][::10], color=COLORS["temperature"], lw=1.5, label="平均温度")
    ax_e.set(xlabel="时间 / min", ylabel="体积平均温度 / ℃", xlim=(0, 30))
    ax_e2 = ax_e.twinx()
    ax_e2.plot(t_plot, q1["moisture_mean"][::10], color=COLORS["moisture"], lw=1.45, label="平均含水率")
    ax_e2.set_ylabel("体积平均含水率 / (kg/kg)", color=COLORS["moisture"])
    ax_e2.tick_params(axis="y", colors=COLORS["moisture"])
    handles = [mpl.lines.Line2D([], [], color=COLORS["temperature"], lw=1.5), mpl.lines.Line2D([], [], color=COLORS["moisture"], lw=1.45)]
    ax_e.legend(handles, ["平均温度", "平均含水率"], loc="center right")
    ax_e.text(0.03, 0.94, f"30 min：{q1['temperature_mean'][-1]:.2f} ℃，{q1['moisture_mean'][-1]:.3f} kg/kg", transform=ax_e.transAxes, va="top", fontsize=6.2, color=COLORS["neutral"])
    grid(ax_e)
    panel(ax_e, "e", "全体平均响应仍未达到稳态")

    ax_f = axes[5]
    delta_t = q1["temperature"][::10, -1] - q1["temperature"][::10, 0]
    delta_c = q1["moisture"][::10, 0] - q1["moisture"][::10, -1]
    ax_f.plot(t_plot, delta_t, color=COLORS["temperature"], lw=1.5, label="表面-中心温差")
    ax_f.fill_between(t_plot, 0, delta_t, color=COLORS["temperature_light"], alpha=0.24, lw=0)
    ax_f.set(xlabel="时间 / min", ylabel="表面-中心温差 / ℃", xlim=(0, 30))
    ax_f2 = ax_f.twinx()
    ax_f2.plot(t_plot, delta_c, color=COLORS["moisture"], lw=1.45, label="中心-表面含水率差")
    ax_f2.set_ylabel("中心-表面含水率差 / (kg/kg)", color=COLORS["moisture"])
    ax_f2.tick_params(axis="y", colors=COLORS["moisture"])
    handles = [mpl.lines.Line2D([], [], color=COLORS["temperature"], lw=1.5), mpl.lines.Line2D([], [], color=COLORS["moisture"], lw=1.45)]
    ax_f.legend(handles, ["温差", "含水率差"], loc="upper left")
    ax_f.text(0.98, 0.08, f"30 min：ΔT={delta_t[-1]:.2f} ℃，ΔC={delta_c[-1]:.3f}", transform=ax_f.transAxes, ha="right", fontsize=6.2, color=COLORS["neutral"])
    grid(ax_f)
    panel(ax_f, "f", "中心—表面梯度揭示热质分离")
    return save_figure(fig, "fig02_q1_fields")


def fig03(cfg, env, q23):
    accepted.OUT = FIG_DIR
    return accepted.render_q2(cfg, env, q23, stem="fig03_q2_coupling")


def fig04(q23, final_validation):
    accepted.OUT = FIG_DIR
    return accepted.render_q3(q23, final_validation, stem="fig04_q3_endpoint")


def fig05(q4, radius, final_validation):
    accepted.OUT = FIG_DIR
    return accepted.render_q4(q4, radius, final_validation, stem="fig05_q4_shrinkage")


def _event_temporal(report: dict) -> tuple[list[str], np.ndarray]:
    labels, values = [], []
    rename = {"cap60": "步长60 s", "cap15": "步长15 s", "radau": "Radau", "tol6": "rtol 1e-6", "tol7": "rtol 1e-7", "tol8": "rtol 1e-8"}
    for row in report["temporal"]:
        raw = row["label"].split("_", 1)[1]
        labels.append(rename.get(raw, raw))
        delta = row.get("event_delta_s")
        if delta is None:
            delta = row.get("difference", {}).get("event_delta_s")
        values.append(max(abs(float(delta)), 1e-6))
    return labels, np.asarray(values)


def fig06(final_validation):
    q2c = json_load("reports/q2_convergence.json")
    q3c = json_load("reports/q3_convergence.json")
    q4c = json_load("reports/q4_convergence.json")
    fig = plt.figure(figsize=(7.2, 6.45))
    gs = fig.add_gridspec(2, 6, height_ratios=[1.0, 1.08], left=0.083, right=0.985, bottom=0.09, top=0.965, hspace=0.58, wspace=0.95)
    axes = [fig.add_subplot(gs[0, 0:2]), fig.add_subplot(gs[0, 2:4]), fig.add_subplot(gs[0, 4:6]), fig.add_subplot(gs[1, 0:2]), fig.add_subplot(gs[1, 2:4]), fig.add_subplot(gs[1, 4:6])]

    for ax, quantity, title, ylabel, color, letter in [
        (axes[0], "temperature", "Q2温度空间收敛", "相邻网格最大差 / ℃", COLORS["temperature"], "a"),
        (axes[1], "moisture", "Q2水分空间收敛", "相邻网格最大差 / (kg/kg)", COLORS["moisture"], "b"),
    ]:
        x, y = legacy.mesh_differences(q2c, quantity)
        ax.loglog(x, y, "o-", color=color, lw=1.45, ms=4.0, label="观测差")
        ref = y[-1] * (x[-1] / x) ** 2
        ax.loglog(x, ref, ls="--", color=COLORS["neutral"], lw=0.9, label="二阶参考")
        ax.axhline(5e-5, color=COLORS["gold"], ls=":", lw=0.9, label="5e-5检查尺度")
        order = np.log(y[-2] / y[-1]) / np.log(x[-1] / x[-2])
        ax.text(0.97, 0.08, f"末两级观测阶 p={order:.2f}", transform=ax.transAxes, ha="right", fontsize=6.0, color=color)
        ax.set(xlabel="径向控制体数 N", ylabel=ylabel)
        ascii_log_axis(ax, "x"); ascii_log_axis(ax, "y")
        ax.legend(loc="upper right")
        grid(ax, "both"); panel(ax, letter, title)

    ax = axes[2]
    for report, label, color, marker in [(q3c, "Q3", COLORS["moisture"], "o"), (q4c, "Q4", COLORS["temperature"], "s")]:
        xs, ys = [], []
        rows = report["mesh_series"]
        for i in range(1, len(rows)):
            delta = rows[i].get("event_delta_s")
            if delta is None:
                delta = rows[i]["event_time_s"] - rows[i - 1]["event_time_s"]
            xs.append(rows[i]["n"]); ys.append(abs(float(delta)))
        ax.loglog(xs, ys, marker=marker, color=color, lw=1.4, ms=4.0, label=label)
    ax.axhline(0.18, color=COLORS["gold"], ls=":", lw=0.9, label="四位小时检查尺度")
    ax.fill_between([100, 3000], 1e-4, 0.18, color="#F4EBD9", alpha=0.55, lw=0)
    ax.set(xlabel="径向控制体数 N", ylabel="相邻网格临界时间差 / s", ylim=(3e-2, 2e3))
    ascii_log_axis(ax, "x"); ascii_log_axis(ax, "y")
    ax.legend(ncol=2, loc="upper right")
    grid(ax, "both"); panel(ax, "c", "终点时间空间收敛")

    ax = axes[3]
    labels, v3 = _event_temporal(q3c)
    _, v4 = _event_temporal(q4c)
    x = np.arange(len(labels))
    ax.plot(x, v3, color=COLORS["moisture"], marker="o", lw=1.35, ms=3.8, label="Q3")
    ax.plot(x, v4, color=COLORS["temperature"], marker="s", lw=1.35, ms=3.8, label="Q4")
    ax.axhline(0.18, color=COLORS["gold"], lw=0.85, ls=":")
    ax.set_yscale("log"); ascii_log_axis(ax)
    ax.set_xticks(x, labels, rotation=27, ha="right")
    ax.set(ylabel="临界时间绝对变化 / s", ylim=(3e-5, 10))
    ax.legend(ncol=2, loc="upper right")
    grid(ax); panel(ax, "d", "时间步、容差与积分器复核")

    refs = final_validation["independent_references"]
    qlabels = ["Q1", "Q2", "Q3", "Q4"]
    tvals = [x["comparisons"]["temperature"]["primary_unrounded_vs_reference_max_abs"] for x in refs]
    cvals = [x["comparisons"]["moisture"]["primary_unrounded_vs_reference_max_abs"] for x in refs]
    x = np.arange(4); w = 0.34
    ax = axes[4]
    ax.bar(x - w / 2, tvals, w, color=COLORS["temperature_light"], edgecolor=COLORS["temperature"], lw=0.65, hatch="///", label="温度")
    ax.bar(x + w / 2, cvals, w, color=COLORS["moisture_light"], edgecolor=COLORS["moisture"], lw=0.65, label="含水率")
    ax.axhline(5e-5, color=COLORS["gold"], ls=":", lw=0.9, label="检查尺度")
    ax.set_yscale("log"); ascii_log_axis(ax)
    ax.set_xticks(x, qlabels)
    ax.set(ylabel="与独立有限元最大绝对差", ylim=(1e-7, 1e-4))
    ax.legend(ncol=2, loc="upper left")
    grid(ax); panel(ax, "e", "独立实现交叉核验")

    rows = final_validation["conservation"]["rows"]
    water = [x["water_relative"] for x in rows]
    heat = [x["heat_relative"] for x in rows]
    ax = axes[5]
    ax.bar(x - w / 2, water, w, color=COLORS["moisture_light"], edgecolor=COLORS["moisture"], lw=0.65, label="水分")
    ax.bar(x + w / 2, heat, w, color=COLORS["temperature_light"], edgecolor=COLORS["temperature"], lw=0.65, hatch="///", label="有效显热")
    ax.axhline(1e-6, color=COLORS["gold"], ls=":", lw=0.9, label="验收阈值")
    ax.set_yscale("log"); ascii_log_axis(ax)
    ax.set_xticks(x, qlabels)
    ax.set(ylabel="归一化累计守恒残差", ylim=(1e-12, 3e-6))
    ax.legend(ncol=2, loc="upper left")
    grid(ax); panel(ax, "f", "离散守恒恒等式检查")
    return save_figure(fig, "fig06_numerical_validation")


def _tornado(ax, report, prefix, baseline_s, letter, title):
    values = legacy.sensitivity_map(report, prefix)
    params = ["D", "hm", "beta", "h"]
    labels = ["扩散率 $D$", "传质系数 $h_m$", "边界映射 $\\beta$", "换热系数 $h$"]
    y = np.arange(len(params))[::-1]
    for yi, p in zip(y, params):
        lo = 100.0 * (values[p][0.9] / baseline_s - 1.0)
        hi = 100.0 * (values[p][1.1] / baseline_s - 1.0)
        ax.plot([lo, hi], [yi, yi], color=COLORS["neutral_light"], lw=5.5, solid_capstyle="round", zorder=1)
        ax.scatter(lo, yi, color=COLORS["gold"], marker="v", s=27, edgecolor="white", lw=0.4, zorder=3, label="参数×0.9" if yi == y[0] else None)
        ax.scatter(hi, yi, color=COLORS["moisture"], marker="^", s=27, edgecolor="white", lw=0.4, zorder=3, label="参数×1.1" if yi == y[0] else None)
        ax.text(lo, yi + 0.18, f"{lo:+.2f}%", ha="center", va="bottom", fontsize=5.6, color=COLORS["gold"])
        ax.text(hi, yi - 0.18, f"{hi:+.2f}%", ha="center", va="top", fontsize=5.6, color=COLORS["moisture"])
    ax.axvline(0, color=COLORS["threshold"], lw=0.8)
    ax.set_yticks(y, labels)
    ax.set_xlabel("临界时间相对变化 / %")
    ax.legend(ncol=2, loc="lower right")
    grid(ax, "x"); panel(ax, letter, title)


def fig07(final_validation):
    q3s = json_load("reports/q3_sensitivity.json")
    q4s = json_load("reports/q4_sensitivity.json")
    endface = json_load("reports/endface_diagnostics.json")
    base3 = float(final_validation["critical_results"]["Q3"]["critical_time_s"])
    base4 = float(final_validation["critical_results"]["Q4"]["critical_time_s"])
    fig = plt.figure(figsize=(7.2, 6.45))
    gs = fig.add_gridspec(2, 4, height_ratios=[1.0, 1.08], left=0.105, right=0.985, bottom=0.09, top=0.965, hspace=0.5, wspace=0.9)
    ax_a = fig.add_subplot(gs[0, :2]); ax_b = fig.add_subplot(gs[0, 2:]); ax_c = fig.add_subplot(gs[1, :2]); ax_d = fig.add_subplot(gs[1, 2:])
    _tornado(ax_a, q3s, "q3", base3, "a", "Q3单因素响应｜扩散率主导")
    _tornado(ax_b, q4s, "q4", base4, "b", "Q4单因素响应｜排序保持稳定")

    q3_values = legacy.sensitivity_map(q3s, "q3")
    q4_values = legacy.sensitivity_map(q4s, "q4")
    params = ["D", "hm", "beta", "h"]
    labels_s = ["扩散率 $D$", "传质系数 $h_m$", "边界映射 $\\beta$", "换热系数 $h$"]
    sens3 = np.asarray([(q3_values[p][1.1] - q3_values[p][0.9]) / (0.2 * base3) for p in params])
    sens4 = np.asarray([(q4_values[p][1.1] - q4_values[p][0.9]) / (0.2 * base4) for p in params])
    y = np.arange(len(params))[::-1]
    for yi, a, b in zip(y, sens3, sens4):
        ax_c.plot([a, b], [yi, yi], color=COLORS["neutral_light"], lw=4.5, zorder=1)
        ax_c.scatter(a, yi + 0.08, s=30, color=COLORS["moisture"], marker="o", edgecolor="white", lw=0.4, label="Q3" if yi == y[0] else None, zorder=3)
        ax_c.scatter(b, yi - 0.08, s=28, color=COLORS["temperature"], marker="s", edgecolor="white", lw=0.4, label="Q4" if yi == y[0] else None, zorder=3)
        if min(a, b) < -0.75:
            ax_c.text(max(a, b) + 0.025, yi, f"{a:.3f}/{b:.3f}", ha="left", va="center", fontsize=5.7, color=COLORS["neutral"])
        else:
            ax_c.text(min(a, b) - 0.025, yi, f"{a:.3f}/{b:.3f}", ha="right", va="center", fontsize=5.7, color=COLORS["neutral"])
    ax_c.axvline(0, color=COLORS["threshold"], lw=0.8)
    ax_c.axvspan(-1.02, -0.5, color="#EEEAF8", alpha=0.65, lw=0)
    ax_c.set_yticks(y, labels_s)
    ax_c.set(xlabel="有符号无量纲局部灵敏度 $S_p$", xlim=(-1.04, 0.12))
    ax_c.legend(ncol=2, loc="lower right")
    ax_c.text(0.03, 0.06, r"统一排序：$|S_D|>|S_{h_m}|>|S_{\beta}|>|S_h|$", transform=ax_c.transAxes, fontsize=6.0, color=COLORS["gold"])
    grid(ax_c, "x"); panel(ax_c, "c", "方向、强度与稳健排序")

    selected = final_validation["numerical_error"]["selected"]
    numerical3 = max(x["value"] for x in selected["Q3"] if x["quantity"] == "critical_time")
    numerical4 = max(x["value"] for x in selected["Q4"] if x["quantity"] == "critical_time")
    uncertainty = final_validation["scenario_uncertainty"]
    boundary3 = max(abs(uncertainty["Q3_boundary_extension_plus_baseline"]["minimum_event_time_s"] - base3), abs(uncertainty["Q3_boundary_extension_plus_baseline"]["maximum_event_time_s"] - base3))
    boundary4 = max(abs(uncertainty["Q4_boundary_extension_plus_baseline"]["minimum_event_time_s"] - base4), abs(uncertainty["Q4_boundary_extension_plus_baseline"]["maximum_event_time_s"] - base4))
    param3 = max(abs(x["event_time_s"] - base3) for x in q3s["cases"] if any(f"_{p}_multiplier_" in x["label"] for p in ("h", "hm", "D", "beta")))
    param4 = max(abs(x["event_time_s"] - base4) for x in q4s["cases"] if any(f"_{p}_multiplier_" in x["label"] for p in ("h", "hm", "D", "beta")))
    structural3 = max(abs(x["event_time_s"] - base3) for x in uncertainty["Q3_mechanism_ablations"])
    structural4 = max(abs(x["event_time_s"] - base4) for x in uncertainty["Q4_factorial_counterfactuals"])
    radius4 = abs(uncertainty["Q4_radius_interpolation"][0]["delta_from_linear_s"])
    ef3 = abs(next(x["fine_paired_effect_s"] for x in endface["paired_refinements"] if x["group"] == "q23" and x["kind"] == "axial"))
    ef4 = abs(next(x["fine_paired_effect_s"] for x in endface["paired_refinements"] if x["group"] == "q4" and x["kind"] == "radial_at_tight_time"))
    labels = ["数值离散", "半径插值", "端面情景", "边界延拓", "±10%参数", "机制/因素对照"]
    q3 = np.asarray([numerical3, np.nan, ef3, boundary3, param3, structural3])
    q4 = np.asarray([numerical4, radius4, ef4, boundary4, param4, structural4])
    y = np.arange(len(labels))[::-1]
    ax_d.axvspan(1e-2, 1, color="#E9F1EF", alpha=0.7, lw=0)
    ax_d.axvspan(1, 1e3, color="#F4EBD9", alpha=0.55, lw=0)
    ax_d.axvspan(1e3, 1e6, color="#EEEAF8", alpha=0.55, lw=0)
    ax_d.scatter(q3, y + 0.13, color=COLORS["moisture"], s=28, marker="o", label="Q3", edgecolor="white", lw=0.4, zorder=3)
    ax_d.scatter(q4, y - 0.13, color=COLORS["temperature"], s=27, marker="s", label="Q4", edgecolor="white", lw=0.4, zorder=3)
    for yi, a, b in zip(y, q3, q4):
        finite = [v for v in (a, b) if np.isfinite(v)]
        ax_d.plot([min(finite), max(finite)], [yi, yi], color=COLORS["neutral_light"], lw=1.0, zorder=1)
    ax_d.set_xscale("log"); ascii_log_axis(ax_d, "x")
    ax_d.xaxis.set_major_locator(LogLocator(base=10, subs=(1.0,), numticks=9))
    ax_d.set_yticks(y, labels)
    ax_d.set(xlabel="临界时间最大绝对情景差 / s", xlim=(1e-2, 1e6))
    ax_d.legend(ncol=2, loc="lower right")
    ax_d.text(0.02, 0.97, "数值", transform=ax_d.transAxes, va="top", fontsize=5.8, color=COLORS["radius"])
    ax_d.text(0.30, 0.97, "输入/边界", transform=ax_d.transAxes, va="top", fontsize=5.8, color=COLORS["gold"])
    ax_d.text(0.70, 0.97, "参数/结构", transform=ax_d.transAxes, va="top", fontsize=5.8, color=COLORS["accent"])
    grid(ax_d, "x"); panel(ax_d, "d", "风险量级分层｜非置信区间")
    return save_figure(fig, "fig07_sensitivity_risk")


def build_contact_sheet(paths):
    return legacy.build_contact_sheet(paths)


def main() -> None:
    style()
    inputs = verify_sources()
    cfg = load_config()
    env, radius = load_measurements(cfg)
    final_validation = json_load("reports/final_result_validation.json")
    with np.load(ROOT / "results/final/q1_full_precision.npz") as q1, np.load(ROOT / "results/final/q23_full_precision.npz") as q23, np.load(ROOT / "results/final/q4_full_precision.npz") as q4:
        generated = [
            fig01(cfg, env, radius, final_validation),
            fig02(q1, env),
            fig03(cfg, env, q23),
            fig04(q23, final_validation),
            fig05(q4, radius, final_validation),
            fig06(final_validation),
            fig07(final_validation),
        ]
    contact = build_contact_sheet(generated)
    output_records = {}
    for stem, group in zip(CONTRACTS, generated):
        output_records[stem] = {
            "files": {p.suffix.lstrip("."): {"bytes": p.stat().st_size, "sha256": sha256(p)} for p in group},
            "source_files": {str(Path(__file__).relative_to(ROOT)): {"bytes": Path(__file__).stat().st_size, "sha256": sha256(Path(__file__))}},
            "contract": CONTRACTS[stem],
        }
    manifest = {
        "status": "generated_information_enriched_pending_complete_visual_qa",
        "backend": f"Python {sys.version_info.major}.{sys.version_info.minor} / matplotlib {mpl.__version__}",
        "palette": "B 靛蓝鎏金",
        "figure_width_mm": 182.88,
        "raster_dpi": 600,
        "svg_text_editable": True,
        "font_policy": "Songti SC (SimSun-compatible) for Chinese; editable TrueType text in SVG/PDF",
        "heatmap_integrity": "No smoothing or invented values; display-only temporal sampling is declared; Q4 physical-domain exterior is masked white.",
        "input_manifest_items_verified": len(inputs["files"]),
        "source_artifacts": {
            "results/final/q1_full_precision.npz": digest(ROOT / "results/final/q1_full_precision.npz"),
            "results/final/q23_full_precision.npz": digest(ROOT / "results/final/q23_full_precision.npz"),
            "results/final/q4_full_precision.npz": digest(ROOT / "results/final/q4_full_precision.npz"),
            "reports/final_result_validation.json": digest(ROOT / "reports/final_result_validation.json"),
        },
        "figures": output_records,
        "contact_sheet": {"path": str(contact.relative_to(ROOT)), "sha256": sha256(contact)},
        "script_sha256": digest(__file__),
    }
    (FIG_DIR / "figure_manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"figures": len(generated), "files": sum(len(x) for x in generated), "contact_sheet": str(contact)}, ensure_ascii=False))


if __name__ == "__main__":
    main()
