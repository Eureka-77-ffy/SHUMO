#!/usr/bin/env python3
"""Render enriched B-palette layout prototypes for Q2, Q3 and Q4.

Figure contracts
----------------
Q2 core conclusion:
    The thermal field nearly equilibrates within 3 h, while moisture retains a
    radial gradient; warming raises diffusivity, but dehydration offsets part
    of that gain and separates the centre and surface state paths.
Q3 core conclusion:
    The full-domain stopping event occurs at 57.4740 h and is controlled by the
    last wet point; surface and volume-average criteria would stop much earlier.
Q4 core conclusion:
    Shrinkage shortens the diffusion path and raises surface-to-volume exchange;
    the geometry and property effects interact rather than add independently.

Archetype: asymmetric quantitative composites with a dominant spatial panel.
Target/output: CUMCM double-column figures, 182.9 mm wide; editable SVG/PDF and
600 dpi PNG/TIFF previews.
Backend: Python/matplotlib only.
Statistics: deterministic PDE scenarios; no replicates, confidence intervals,
or probability claims.
Image integrity: heatmaps use accepted full-precision caches without smoothing;
Q4 physical-domain exterior remains masked white.
Reviewer risk: diffusivity decomposition and factorial contrasts are model
diagnostics, not experimental causal identification.

These files are written to a preview-only directory and never replace accepted
manuscript figures.
"""

from __future__ import annotations

import csv
import hashlib
import json
import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from PIL import Image, ImageDraw, ImageFont
from pypdf import PdfReader

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from figures.scripts import make_all_figures as manuscript
from src.data_io import load_config, verify_sources
from src.properties import PropertyLaw


OUT = ROOT / "figures" / "layout_previews" / "b_indigo_gold"
FORMATS = ("svg", "pdf", "png", "tiff")

COLORS = {
    "temperature": "#B54B36",
    "temperature_light": "#E7A36A",
    "moisture": "#3F4A9A",
    "moisture_light": "#AEB8E8",
    "radius": "#26877C",
    "mean": "#7376A7",
    "surface": "#C49A38",
    "threshold": "#25242B",
    "neutral": "#70717A",
    "neutral_light": "#E3E0D8",
    "accent": "#8C5A86",
    "blue_accent": "#356FA1",
    "gold": "#D98B28",
    "paper": "#FBFAF7",
}

TEMP_CMAP = mpl.colors.LinearSegmentedColormap.from_list(
    "indigo_gold_temperature", ["#FFF9E8", "#F1D18A", "#D77A32", "#8F2D2D"]
)
MOIST_CMAP = mpl.colors.LinearSegmentedColormap.from_list(
    "indigo_gold_moisture", ["#F7F6FC", "#D2D1EE", "#7A82C2", "#30356F"]
)
MOIST_CMAP.set_bad("white")

CONTRACTS = {
    "fig03_q2_enriched": {
        "core_conclusion": "3小时内热场趋稳而水分梯度仍显著；温升提高扩散率，但失水抵消部分增益并造成中心—表面状态路径分离。",
        "archetype": "asymmetric mixed-modality figure",
        "panel_map": {
            "a": "温度时空场（主证据）",
            "b": "含水率时空场（主证据）",
            "c": "中心、表面与环境温度的分阶段响应",
            "d": "中心与表面的温度—含水率状态路径",
            "e": "中心与表面局部扩散率轨迹",
            "f": "3小时扩散率对数变化的温度项、含水项及净效应分解",
        },
        "evidence_hierarchy": {
            "hero": "a-b",
            "mechanism": "d-f",
            "support": "c-e",
        },
        "reviewer_risk": "扩散率与分解项来自题定经验物性式，只作模型机制诊断。",
        "source_data": ["results/final/q23_full_precision.npz", "config/model_config.json", "original_appendix_1"],
    },
    "fig04_q3_enriched": {
        "core_conclusion": "全域终点由最后湿点控制，57.4740小时才达到阈值；表面值或体积平均值会分别提前44.92和21.45小时停机。",
        "archetype": "asymmetric mixed-modality figure",
        "panel_map": {
            "a": "全时域含水率场（主证据）",
            "b": "三种停机判据的阈值到达时刻",
            "c": "全域最大值、平均值与表面值轨迹",
            "d": "截面达标比例随时间增长",
            "e": "典型时刻径向剖面",
            "f": "临界点邻域与首个严格达标整数秒",
        },
        "evidence_hierarchy": {
            "hero": "a",
            "event_definition": "b-c-f",
            "spatial_mechanism": "d-e",
        },
        "reviewer_risk": "临界时刻对应等号，严格小于阈值从其后时刻开始；截面比例是模型场的几何诊断量。",
        "source_data": ["results/final/q23_full_precision.npz", "reports/final_result_validation.json"],
    },
    "fig05_q4_enriched": {
        "core_conclusion": "收缩将扩散路径倍率提高到2.78、交换尺度提高到1.67，并与物性形成显著非加性交互。",
        "archetype": "asymmetric mixed-modality figure",
        "panel_map": {
            "a": "材料坐标含水率场（主证据）",
            "b": "物理坐标含水率场及移动边界（主证据）",
            "c": "半径、路径与表面积体积比的同步几何变化",
            "d": "收缩域内全域判据与统计量",
            "e": "物性×几何的2×2交互图",
            "f": "不同时间物理半径上的剖面与收缩端点",
        },
        "evidence_hierarchy": {
            "hero": "a-b",
            "mechanism": "c-e",
            "support": "d-f",
        },
        "reviewer_risk": "四组因子对照是受控模型情景，不是实验因果识别；物理域外必须保持空白。",
        "source_data": ["results/final/q4_full_precision.npz", "tables/q4_factorial.csv", "original_appendix_2"],
    },
}


def configure_style() -> None:
    manuscript.style()
    mpl.rcParams.update(
        {
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "axes.titleweight": "semibold",
            "axes.edgecolor": "#3B3A40",
            "axes.labelcolor": "#34333A",
            "xtick.color": "#4F4E56",
            "ytick.color": "#4F4E56",
            "legend.handlelength": 1.7,
            "legend.columnspacing": 0.9,
            "legend.handletextpad": 0.45,
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
        }
    )


def panel(ax, label: str, title: str) -> None:
    ax.text(
        -0.115,
        1.035,
        label,
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=9.2,
        fontweight="bold",
        color=COLORS["threshold"],
    )
    ax.set_title(title, loc="left", pad=4.5)


def subtle_grid(ax, axis: str = "y") -> None:
    ax.grid(axis=axis, color="#E9E6E0", lw=0.48, zorder=0)
    ax.set_axisbelow(True)


def phase_bands(ax, spans: list[tuple[float, float, str]]) -> None:
    fills = ["#F6E8DA", "#EEEAF8", "#E5F0EE"]
    y0, y1 = ax.get_ylim()
    for (lo, hi, label), fill in zip(spans, fills):
        ax.axvspan(lo, hi, color=fill, alpha=0.48, lw=0, zorder=0)
        ax.text((lo + hi) / 2, y1 - 0.035 * (y1 - y0), label, ha="center", va="top", fontsize=6.0, color=COLORS["neutral"])


def heatmap(
    ax,
    time,
    space,
    values,
    cmap,
    label: str | None,
    *,
    vmin=None,
    vmax=None,
    add_colorbar: bool = True,
):
    x = manuscript.centers_to_edges(np.asarray(time), lower=float(np.min(time)), upper=float(np.max(time)))
    y = manuscript.centers_to_edges(np.asarray(space), lower=float(np.min(space)), upper=float(np.max(space)))
    mesh = ax.pcolormesh(
        x,
        y,
        np.asarray(values).T,
        cmap=cmap,
        shading="auto",
        vmin=vmin,
        vmax=vmax,
        rasterized=True,
    )
    if add_colorbar:
        cbar = ax.figure.colorbar(mesh, ax=ax, pad=0.014, fraction=0.043)
        cbar.set_label(label)
        cbar.outline.set_linewidth(0.55)
        cbar.ax.tick_params(length=2.2, width=0.55)
    return mesh


def save_figure(fig, stem: str) -> list[Path]:
    outputs = []
    for fmt in FORMATS:
        path = OUT / f"{stem}.{fmt}"
        kwargs = {"bbox_inches": "tight", "facecolor": "white"}
        if fmt in {"png", "tiff"}:
            kwargs["dpi"] = 600
        if fmt == "tiff":
            kwargs["pil_kwargs"] = {"compression": "tiff_lzw"}
        fig.savefig(path, **kwargs)
        outputs.append(path)
    plt.close(fig)
    return outputs


def crossing_time_hours(time_h: np.ndarray, series: np.ndarray, threshold: float = 0.15) -> float:
    indices = np.flatnonzero(series <= threshold)
    if not len(indices):
        return float("nan")
    i = int(indices[0])
    if i == 0:
        return float(time_h[0])
    return float(np.interp(threshold, [series[i], series[i - 1]], [time_h[i], time_h[i - 1]]))


def dry_cross_section_fraction(radius_m: np.ndarray, moisture: np.ndarray, threshold: float = 0.15) -> np.ndarray:
    out = np.zeros(len(moisture), dtype=float)
    R = float(radius_m[-1])
    for i, row in enumerate(moisture):
        if row[0] <= threshold:
            out[i] = 1.0
            continue
        if row[-1] > threshold:
            continue
        j = int(np.flatnonzero(row <= threshold)[0])
        r_cross = float(np.interp(threshold, [row[j], row[j - 1]], [radius_m[j], radius_m[j - 1]]))
        out[i] = 1.0 - (r_cross / R) ** 2
    return out


def render_q2(cfg, env: np.ndarray, q23, stem: str = "fig03_q2_enriched") -> list[Path]:
    end = int(np.searchsorted(q23["time_s"], 10800.0, side="right"))
    sample = np.unique(np.r_[np.arange(0, end, 12), end - 1])
    stride = 24
    time_h = q23["time_s"][:end] / 3600.0
    t = time_h[::stride]
    radius_cm = q23["radius_m"] * 100.0

    fig = plt.figure(figsize=(7.2, 7.05), layout="constrained")
    gs = fig.add_gridspec(3, 4, height_ratios=[1.12, 0.92, 0.92], hspace=0.10, wspace=0.18)
    ax_a = fig.add_subplot(gs[0, :2])
    ax_b = fig.add_subplot(gs[0, 2:])
    ax_c = fig.add_subplot(gs[1, :2])
    ax_d = fig.add_subplot(gs[1, 2:])
    ax_e = fig.add_subplot(gs[2, :2])
    ax_f = fig.add_subplot(gs[2, 2:])

    heatmap(ax_a, time_h[sample], radius_cm, q23["temperature"][:end][sample], TEMP_CMAP, "温度 / ℃")
    ax_a.set(xlabel="时间 / h", ylabel="距中心距离 / cm", xlim=(0, 3))
    ax_a.axvline(0.5, color="white", lw=0.8, ls="--", alpha=0.9)
    panel(ax_a, "a", "温度时空场｜热响应快速贯穿")

    heatmap(ax_b, time_h[sample], radius_cm, q23["moisture"][:end][sample], MOIST_CMAP, "干基含水率 / (kg/kg)")
    ax_b.set(xlabel="时间 / h", ylabel="距中心距离 / cm", xlim=(0, 3))
    ax_b.axvline(0.5, color="white", lw=0.8, ls="--", alpha=0.9)
    panel(ax_b, "b", "含水率时空场｜径向梯度持续存在")

    Tc = q23["temperature"][:end:stride, 0]
    Ts = q23["temperature"][:end:stride, -1]
    Cc = q23["moisture"][:end:stride, 0]
    Cs = q23["moisture"][:end:stride, -1]
    env_t = np.interp(q23["time_s"][:end:stride], env[:, 0], env[:, 1])
    ax_c.plot(t, Tc, color=COLORS["temperature"], lw=1.55, label="中心")
    ax_c.plot(t, Ts, color=COLORS["temperature_light"], lw=1.35, label="表面")
    ax_c.plot(t, env_t, color=COLORS["neutral"], lw=1.15, ls="--", label="环境")
    ax_c.set(xlabel="时间 / h", ylabel="温度 / ℃", xlim=(0, 3), ylim=(27, 52.5))
    phase_bands(ax_c, [(0, 0.5, "快速预热"), (0.5, 1.5, "内外追赶"), (1.5, 3, "平台趋稳")])
    ax_c.legend(ncol=3, loc="lower right")
    subtle_grid(ax_c)
    ax_c.annotate(
        f"3 h中心 {Tc[-1]:.2f} ℃",
        (t[-1], Tc[-1]),
        xytext=(-8, -18),
        textcoords="offset points",
        ha="right",
        color=COLORS["temperature"],
        arrowprops={"arrowstyle": "-", "lw": 0.6, "color": COLORS["temperature"]},
    )
    panel(ax_c, "c", "温度响应的阶段划分")

    ax_d.plot(Tc, Cc, color=COLORS["moisture"], lw=1.55, label="中心路径")
    ax_d.plot(Ts, Cs, color=COLORS["surface"], lw=1.45, label="表面路径")
    for hour in (0.5, 1.0, 2.0, 3.0):
        i = int(np.argmin(np.abs(t - hour)))
        ax_d.scatter(Tc[i], Cc[i], s=18, color=COLORS["moisture"], edgecolor="white", lw=0.45, zorder=4)
        ax_d.scatter(Ts[i], Cs[i], s=18, color=COLORS["surface"], edgecolor="white", lw=0.45, zorder=4)
        if hour in (0.5, 3.0):
            ax_d.text(Ts[i] + 0.20, Cs[i] - 0.055, f"{hour:g} h", fontsize=6.1, color=COLORS["neutral"])
    ax_d.scatter([Tc[0]], [Cc[0]], marker="D", s=24, facecolor="white", edgecolor=COLORS["threshold"], lw=0.8, zorder=5)
    ax_d.set(xlabel="局部温度 / ℃", ylabel="局部干基含水率 / (kg/kg)")
    ax_d.legend(loc="lower left")
    subtle_grid(ax_d)
    panel(ax_d, "d", "热湿状态路径｜同温度、不同含水率")

    prop = PropertyLaw(cfg, "q23")
    _, _, Dc = prop.evaluate(Cc, Tc)
    _, _, Ds = prop.evaluate(Cs, Ts)
    ax_e.plot(t, Dc, color=COLORS["moisture"], lw=1.55, label="中心")
    ax_e.plot(t, Ds, color=COLORS["surface"], lw=1.4, label="表面")
    ax_e.fill_between(t, Ds, Dc, where=Dc >= Ds, color=COLORS["moisture_light"], alpha=0.28, lw=0)
    ax_e.set_yscale("log")
    manuscript.ascii_log_axis(ax_e)
    ax_e.set(xlabel="时间 / h", ylabel="局部扩散率 D / (m²/s)", xlim=(0, 3))
    ax_e.legend(loc="lower right")
    subtle_grid(ax_e)
    ax_e.text(0.98, 0.94, f"3 h比值  D中心/D表面 = {Dc[-1] / Ds[-1]:.2f}", transform=ax_e.transAxes, ha="right", va="top", fontsize=6.35, color=COLORS["neutral"])
    panel(ax_e, "e", "局部扩散率的中心—表面分化")

    dcfg = cfg["property_sets"]["q23"]["D_m2_s"]
    a = float(dcfg["moisture_exponent"])
    E = float(dcfg["thermal_exponent_K"])
    terms = []
    for j in (0, -1):
        T0 = float(q23["temperature"][0, j])
        C0 = float(q23["moisture"][0, j])
        T3 = float(q23["temperature"][end - 1, j])
        C3 = float(q23["moisture"][end - 1, j])
        thermal = E * (1.0 / (T0 + 273.15) - 1.0 / (T3 + 273.15))
        moisture = a * (1.0 / C0 - 1.0 / C3)
        terms.append((thermal, moisture, thermal + moisture))
    x = np.arange(2)
    width = 0.23
    labels = ["温度贡献", "含水贡献", "净变化"]
    bar_colors = [COLORS["temperature_light"], COLORS["moisture_light"], COLORS["gold"]]
    edges = [COLORS["temperature"], COLORS["moisture"], "#8D641A"]
    hatches = ["///", "\\\\", ""]
    for k in range(3):
        vals = [terms[0][k], terms[1][k]]
        bars = ax_f.bar(x + (k - 1) * width, vals, width, color=bar_colors[k], edgecolor=edges[k], lw=0.65, hatch=hatches[k], label=labels[k], zorder=2)
        for bar, val in zip(bars, vals):
            ax_f.text(bar.get_x() + bar.get_width() / 2, val + (0.025 if val >= 0 else -0.035), f"{val:+.2f}", ha="center", va="bottom" if val >= 0 else "top", fontsize=5.8)
    ax_f.axhline(0, color=COLORS["threshold"], lw=0.75)
    ax_f.set_xticks(x, ["中心", "表面"])
    ax_f.set(ylabel="相对初态的 Δln D", ylim=(-0.38, 1.05))
    ax_f.legend(ncol=3, loc="upper right")
    subtle_grid(ax_f)
    panel(ax_f, "f", "3 h扩散率变化的机制分解")

    return save_figure(fig, stem)


def render_q3(q23, final_validation: dict, stem: str = "fig04_q3_enriched") -> list[Path]:
    event_h = float(final_validation["critical_results"]["Q3"]["critical_time_h"])
    event_s = float(final_validation["critical_results"]["Q3"]["critical_time_s"])
    strict_s = float(final_validation["critical_results"]["Q3"]["first_strict_integer_second"])
    time_h = q23["time_s"] / 3600.0
    radius_cm = q23["radius_m"] * 100.0
    sample = np.unique(np.r_[np.arange(0, len(time_h), 72), len(time_h) - 1])
    stride = 60

    fig = plt.figure(figsize=(7.2, 7.35))
    gs = fig.add_gridspec(
        3,
        6,
        height_ratios=[1.12, 0.92, 0.92],
        left=0.085,
        right=0.985,
        bottom=0.075,
        top=0.965,
        hspace=0.55,
        wspace=0.85,
    )
    ax_a = fig.add_subplot(gs[0, :4])
    ax_b = fig.add_subplot(gs[0, 4:])
    ax_c = fig.add_subplot(gs[1, :3])
    ax_d = fig.add_subplot(gs[1, 3:])
    ax_e = fig.add_subplot(gs[2, :3])
    ax_f = fig.add_subplot(gs[2, 3:])

    mesh = heatmap(ax_a, time_h[sample], radius_cm, q23["moisture"][sample], MOIST_CMAP, None, vmin=0.05, vmax=2.55, add_colorbar=False)
    ax_a.axvline(event_h, color="white", lw=1.15, ls="--")
    ax_a.set(xlabel="时间 / h", ylabel="距中心距离 / cm", xlim=(0, event_h))
    cax = ax_a.inset_axes([0.61, 0.84, 0.34, 0.037])
    cax.set_in_layout(False)
    cbar = fig.colorbar(mesh, cax=cax, orientation="horizontal")
    cbar.set_label("干基含水率 C / (kg/kg)", fontsize=5.8, color=COLORS["neutral"])
    cbar.ax.xaxis.set_label_position("top")
    cbar.ax.xaxis.set_ticks_position("bottom")
    cbar.ax.tick_params(labelsize=5.4, length=1.8, width=0.45, colors=COLORS["neutral"])
    cbar.outline.set_linewidth(0.45)
    panel(ax_a, "a", "全时域水分场｜干燥前沿由表面向中心推进")

    crossing = [
        crossing_time_hours(time_h, q23["moisture"][:, -1]),
        crossing_time_hours(time_h, q23["moisture_mean"]),
        crossing_time_hours(time_h, q23["moisture_global_max"]),
    ]
    names = ["表面", "平均", "全域最大"]
    y = np.arange(3)[::-1]
    ax_b.hlines(y, 0, crossing, color=COLORS["neutral_light"], lw=4.5, zorder=1)
    for yi, value, color, marker in zip(y, crossing, [COLORS["surface"], COLORS["mean"], COLORS["moisture"]], ["s", "D", "o"]):
        ax_b.scatter(value, yi, s=38, color=color, marker=marker, edgecolor="white", lw=0.6, zorder=3)
        ax_b.text(value - 1.0, yi + 0.16, f"{value:.2f} h", ha="right", va="bottom", fontsize=6.3, color=color, fontweight="semibold")
    ax_b.axvline(event_h, color=COLORS["gold"], lw=0.9, ls=":")
    ax_b.set_yticks(y, names)
    ax_b.set(xlabel="首次达到0.15的时刻 / h", xlim=(0, 61), ylim=(-0.65, 2.6))
    ax_b.text(0.98, 0.05, f"平均判据将提前\n{event_h - crossing[1]:.2f} h", transform=ax_b.transAxes, ha="right", va="bottom", fontsize=6.15, color=COLORS["gold"])
    subtle_grid(ax_b, "x")
    panel(ax_b, "b", "三种停机判据的时间偏差")

    ax_c.plot(time_h[::stride], q23["moisture_global_max"][::stride], color=COLORS["moisture"], lw=1.55, label="全域最大值")
    ax_c.plot(time_h[::stride], q23["moisture_mean"][::stride], color=COLORS["mean"], lw=1.3, label="体积平均值")
    ax_c.plot(time_h[::stride], q23["moisture"][::stride, -1], color=COLORS["surface"], lw=1.25, label="表面值")
    ax_c.axhline(0.15, color=COLORS["threshold"], lw=0.9, ls="--", label="阈值")
    ax_c.axvline(event_h, color=COLORS["gold"], lw=0.85, ls=":")
    ax_c.set(xlabel="时间 / h", ylabel="干基含水率 / (kg/kg)", xlim=(0, event_h), ylim=(0.03, 2.65))
    ax_c.legend(ncol=2, loc="upper right")
    subtle_grid(ax_c)
    panel(ax_c, "c", "判据量的长期分离")

    sampled_moisture = q23["moisture"][::stride]
    fractions = dry_cross_section_fraction(q23["radius_m"], sampled_moisture)
    sampled_t = time_h[::stride]
    ax_d.fill_between(sampled_t, 0, 100 * fractions, color=COLORS["moisture_light"], alpha=0.62, lw=0)
    ax_d.plot(sampled_t, 100 * fractions, color=COLORS["moisture"], lw=1.5)
    for hour in (24, 36, 48, event_h):
        value = 100 * float(np.interp(hour, sampled_t, fractions))
        ax_d.scatter(hour, value, s=19, color=COLORS["gold"] if hour == event_h else COLORS["moisture"], edgecolor="white", lw=0.45, zorder=3)
        if hour in (36, event_h):
            ax_d.text(hour - 0.8, value - (7 if hour == event_h else -4), f"{value:.0f}%", ha="right", va="top" if hour == event_h else "bottom", fontsize=6.2)
    ax_d.axvline(event_h, color=COLORS["gold"], lw=0.85, ls=":")
    ax_d.set(xlabel="时间 / h", ylabel="截面达标率 / %", xlim=(0, event_h), ylim=(0, 104))
    subtle_grid(ax_d)
    panel(ax_d, "d", "达标区域的截面占比")

    selected_h = [3, 12, 24, 36, 48, event_h]
    line_cmap = mpl.colors.LinearSegmentedColormap.from_list("q3_profiles", ["#D8D6ED", COLORS["moisture"], "#252452"])
    line_colors = line_cmap(np.linspace(0.08, 0.95, len(selected_h)))
    for hour, color in zip(selected_h, line_colors):
        i = int(np.argmin(np.abs(time_h - hour)))
        label = f"终点 {event_h:.2f} h" if hour == event_h else f"{hour:g} h"
        lw = 1.65 if hour == event_h else 1.18
        ax_e.plot(radius_cm, q23["moisture"][i], color=color, lw=lw, label=label)
    ax_e.axhline(0.15, color=COLORS["threshold"], lw=0.85, ls="--")
    ax_e.set(xlabel="距中心距离 / cm", ylabel="干基含水率 / (kg/kg)", xlim=(0, 2.0), ylim=(0.03, 1.9))
    ax_e.legend(ncol=2, loc="upper right")
    subtle_grid(ax_e)
    panel(ax_e, "e", "最后湿点锁定在中心区域")

    mask = q23["time_s"] >= event_s - 3600.0
    rel_min = (q23["time_s"][mask] - event_s) / 60.0
    residual = 1e3 * (q23["moisture_global_max"][mask] - 0.15)
    ax_f.plot(rel_min, residual, color=COLORS["moisture"], lw=1.55)
    ax_f.fill_between(rel_min, 0, residual, where=residual >= 0, color=COLORS["moisture_light"], alpha=0.32, lw=0)
    ax_f.axhline(0, color=COLORS["threshold"], lw=0.85, ls="--")
    ax_f.axvline(0, color=COLORS["gold"], lw=0.95, ls=":")
    strict_y = float(np.interp(strict_s, q23["time_s"], q23["moisture_global_max"]) * 1e3 - 150.0)
    ax_f.scatter((strict_s - event_s) / 60.0, strict_y, s=28, color=COLORS["blue_accent"], marker="D", edgecolor="white", lw=0.55, zorder=4, label="首个严格达标整数秒")
    ax_f.set(xlabel="相对临界时刻 / min", ylabel="1000(Cmax-0.15)", xlim=(-60, 0.55))
    ax_f.legend(loc="upper right")
    subtle_grid(ax_f)
    ax_f.text(0.03, 0.10, f"等号时刻  {event_h:.4f} h\n严格小于  {strict_s:.0f} s", transform=ax_f.transAxes, ha="left", va="bottom", fontsize=6.25, color=COLORS["neutral"])
    panel(ax_f, "f", "终点根的局部解析")

    return save_figure(fig, stem)


def render_q4(q4, radius_data: np.ndarray, final_validation: dict, stem: str = "fig05_q4_enriched") -> list[Path]:
    event_h = float(final_validation["critical_results"]["Q4"]["critical_time_h"])
    pt, xi, material, physical_r_cm, physical, radii_m = manuscript.interpolate_material_profiles(q4, radius_data)
    ph = pt / 3600.0
    t_h = q4["time_s"] / 3600.0
    R_cm = q4["radius_current_m"] * 100.0
    s = R_cm / 2.0

    fig = plt.figure(figsize=(7.2, 7.15), layout="constrained")
    gs = fig.add_gridspec(3, 4, height_ratios=[1.12, 0.92, 0.95], hspace=0.11, wspace=0.18)
    ax_a = fig.add_subplot(gs[0, :2])
    ax_b = fig.add_subplot(gs[0, 2:])
    ax_c = fig.add_subplot(gs[1, :2])
    ax_d = fig.add_subplot(gs[1, 2:])
    ax_e = fig.add_subplot(gs[2, :2])
    ax_f = fig.add_subplot(gs[2, 2:])

    mesh_a = heatmap(ax_a, ph, xi, material, MOIST_CMAP, None, vmin=0.05, vmax=2.55, add_colorbar=False)
    ax_a.axvline(event_h, color="white", lw=1.05, ls="--")
    ax_a.set(xlabel="时间 / h", ylabel="材料坐标 ξ", xlim=(0, event_h), ylim=(0, 1))
    panel(ax_a, "a", "材料坐标场｜固定计算域")

    heatmap(ax_b, ph, physical_r_cm, np.ma.masked_invalid(physical), MOIST_CMAP, None, vmin=0.05, vmax=2.55, add_colorbar=False)
    ax_b.plot(ph, radii_m * 100.0, color=COLORS["threshold"], lw=1.1, label="移动表面 R(t)")
    ax_b.fill_between(ph, radii_m * 100.0, 2.0, color="white", alpha=0.98, lw=0, zorder=2)
    ax_b.plot(ph, radii_m * 100.0, color=COLORS["threshold"], lw=1.1, zorder=3)
    ax_b.set(xlabel="时间 / h", ylabel="物理半径 r / cm", xlim=(0, event_h), ylim=(0, 2.0))
    ax_b.legend(loc="upper right")
    panel(ax_b, "b", "物理坐标场｜域外留白并显式跟踪边界")
    cbar = fig.colorbar(mesh_a, ax=[ax_a, ax_b], pad=0.012, fraction=0.025)
    cbar.set_label("干基含水率 / (kg/kg)")
    cbar.outline.set_linewidth(0.55)

    ax_c.plot(t_h, s, color=COLORS["radius"], lw=1.5, label="R/R0")
    ax_c.plot(t_h, 1.0 / s, color=COLORS["surface"], lw=1.4, label="R0/R（交换尺度）")
    ax_c.plot(t_h, 1.0 / s**2, color=COLORS["accent"], lw=1.5, label="R0²/R²（扩散路径）")
    ax_c.axvline(event_h, color=COLORS["gold"], lw=0.85, ls=":")
    ax_c.set(xlabel="时间 / h", ylabel="相对初始值 / 倍", xlim=(0, event_h), ylim=(0.52, 2.95))
    ax_c.legend(loc="upper left")
    subtle_grid(ax_c)
    ax_c.text(event_h - 0.7, 1.0 / s[-1] ** 2 + 0.08, f"{1.0 / s[-1] ** 2:.2f}×", ha="right", color=COLORS["accent"], fontsize=6.4, fontweight="semibold")
    ax_c.text(event_h - 0.7, 1.0 / s[-1] - 0.14, f"{1.0 / s[-1]:.2f}×", ha="right", color=COLORS["surface"], fontsize=6.4, fontweight="semibold")
    panel(ax_c, "c", "收缩同时改变两类几何尺度")

    stride = 5
    ax_d.plot(t_h[::stride], q4["moisture_global_max"][::stride], color=COLORS["moisture"], lw=1.55, label="全域最大值")
    ax_d.plot(t_h[::stride], q4["moisture_mean"][::stride], color=COLORS["mean"], lw=1.3, label="体积平均值")
    ax_d.plot(t_h[::stride], q4["moisture"][::stride, -1], color=COLORS["surface"], lw=1.25, label="真实表面值")
    ax_d.axhline(0.15, color=COLORS["threshold"], lw=0.9, ls="--", label="阈值")
    ax_d.axvline(event_h, color=COLORS["gold"], lw=0.85, ls=":")
    ax_d.set(xlabel="时间 / h", ylabel="干基含水率 / (kg/kg)", xlim=(0, event_h), ylim=(0.03, 2.65))
    ax_d.legend(ncol=2, loc="upper right")
    subtle_grid(ax_d)
    panel(ax_d, "d", "收缩域内的全域判据")

    factor = final_validation["scenario_uncertainty"]["Q4_factorial_counterfactuals"]
    lookup = {(row["group"], row["geometry"]): float(row["event_h"]) for row in factor}
    x = np.array([0.0, 1.0])
    q23_vals = [lookup[("q23", "fixed")], lookup[("q23", "linear")]]
    q4_vals = [lookup[("q4", "fixed")], lookup[("q4", "linear")]]
    ax_e.plot(x, q23_vals, color=COLORS["moisture"], marker="o", ms=5.2, lw=1.55, label="附件3物性")
    ax_e.plot(x, q4_vals, color=COLORS["temperature"], marker="s", ms=5.0, lw=1.55, label="附件4物性")
    for xpos, value, label, color in [
        (0, q23_vals[0], "A", COLORS["moisture"]),
        (1, q23_vals[1], "B", COLORS["moisture"]),
        (0, q4_vals[0], "C", COLORS["temperature"]),
        (1, q4_vals[1], "D", COLORS["temperature"]),
    ]:
        dx = -0.04 if xpos == 0 else 0.04
        ax_e.text(xpos + dx, value + 4.2, f"{label}  {value:.2f} h", ha="left" if xpos == 1 else "right", va="bottom", fontsize=6.15, color=color)
    interaction_h = float(json.loads((ROOT / "reports" / "q4_factorial.json").read_text(encoding="utf-8"))["decomposition"]["interaction_s"]) / 3600.0
    ax_e.text(0.98, 0.96, f"非加性交互项  {interaction_h:+.2f} h", transform=ax_e.transAxes, ha="right", va="top", fontsize=6.25, color=COLORS["gold"], bbox={"boxstyle": "round,pad=0.25", "fc": "#FFF7E5", "ec": "#D7B975", "lw": 0.55})
    ax_e.set_xticks(x, ["固定半径", "观测收缩"])
    ax_e.set(ylabel="临界时间 / h", xlim=(-0.18, 1.18), ylim=(0, 145))
    ax_e.legend(loc="lower left")
    subtle_grid(ax_e)
    panel(ax_e, "e", "物性×几何的2×2交互图")

    selected_h = [12, 24, 36, 48, event_h]
    line_cmap = mpl.colors.LinearSegmentedColormap.from_list("q4_profiles", ["#D8D6ED", COLORS["moisture"], "#252452"])
    line_colors = line_cmap(np.linspace(0.12, 0.95, len(selected_h)))
    for hour, color in zip(selected_h, line_colors):
        i = int(np.argmin(np.abs(ph - hour)))
        r_plot = xi * radii_m[i] * 100.0
        label = f"终点 {event_h:.2f} h" if hour == event_h else f"{hour:g} h"
        lw = 1.65 if hour == event_h else 1.18
        ax_f.plot(r_plot, material[i], color=color, lw=lw, label=label)
        ax_f.scatter(r_plot[-1], material[i, -1], s=13, color=color, edgecolor="white", lw=0.35, zorder=3)
    ax_f.axhline(0.15, color=COLORS["threshold"], lw=0.85, ls="--")
    ax_f.set(xlabel="当前物理半径 r / cm", ylabel="干基含水率 / (kg/kg)", xlim=(0, 2.0), ylim=(0.03, 0.82))
    ax_f.legend(ncol=2, loc="upper right")
    subtle_grid(ax_f)
    panel(ax_f, "f", "物理剖面及其收缩端点")

    return save_figure(fig, stem)


def font(size: int):
    candidates = [
        "/System/Library/Fonts/Supplemental/Songti.ttc",
        "/System/Library/Fonts/Supplemental/Arial Unicode.ttf",
    ]
    for candidate in candidates:
        if Path(candidate).is_file():
            return ImageFont.truetype(candidate, size=size)
    return ImageFont.load_default()


def build_board(stems: list[str]) -> Path:
    width, header_h, row_h = 2600, 180, 2100
    canvas = Image.new("RGB", (width, header_h + row_h * len(stems)), "#F3F1EC")
    draw = ImageDraw.Draw(canvas)
    draw.text((72, 34), "B配色｜Q2—Q4信息结构升级样稿", fill="#24232A", font=font(47))
    draw.text((72, 100), "主证据 + 机制分解 + 事件诊断；全部来自已验收未舍入结果", fill="#65636C", font=font(26))
    labels = ["Q2  热湿耦合机制", "Q3  全域事件定义", "Q4  移动域与因子交互"]
    for row, (stem, label) in enumerate(zip(stems, labels)):
        y = header_h + row * row_h
        draw.rounded_rectangle((45, y + 20, width - 45, y + row_h - 25), radius=28, fill="white", outline="#D9D4CA", width=2)
        draw.text((82, y + 48), label, fill="#34323A", font=font(35))
        image = Image.open(OUT / f"{stem}.png").convert("RGB")
        image.thumbnail((2420, 1900), Image.Resampling.LANCZOS)
        x = (width - image.width) // 2
        iy = y + 125 + (1900 - image.height) // 2
        canvas.paste(image, (x, iy))
    path = OUT / "enriched_layout_board.png"
    canvas.save(path, dpi=(220, 220))
    return path


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def verify_formal_figures_unchanged() -> dict[str, bool]:
    manifest = json.loads((ROOT / "figures" / "figure_manifest.json").read_text(encoding="utf-8"))
    checks = {}
    for stem in ("fig03_q2_coupling", "fig04_q3_endpoint", "fig05_q4_shrinkage"):
        for ext, meta in manifest["figures"][stem]["files"].items():
            path = ROOT / "figures" / f"{stem}.{ext}"
            checks[f"{stem}.{ext}"] = path.is_file() and digest(path) == meta["sha256"]
    return checks


def qa(outputs: dict[str, list[Path]]) -> dict:
    checks = {}
    for stem, paths in outputs.items():
        row = {}
        by_ext = {p.suffix.lstrip("."): p for p in paths}
        svg_text = by_ext["svg"].read_text(encoding="utf-8")
        row["svg_editable_text"] = "<text" in svg_text and "Songti SC" in svg_text and "�" not in svg_text
        row["pdf_single_page"] = len(PdfReader(str(by_ext["pdf"])).pages) == 1
        for ext in ("png", "tiff"):
            with Image.open(by_ext[ext]) as image:
                dpi = image.info.get("dpi", (0, 0))
                row[f"{ext}_pixels"] = list(image.size)
                row[f"{ext}_dpi"] = [float(dpi[0]), float(dpi[1])]
                row[f"{ext}_600dpi"] = min(float(dpi[0]), float(dpi[1])) >= 590.0
        row["files"] = {p.suffix.lstrip("."): {"bytes": p.stat().st_size, "sha256": digest(p)} for p in paths}
        checks[stem] = row
    return checks


def export_derived_source_data(q23, final_validation: dict) -> list[str]:
    source_dir = OUT / "source_data"
    source_dir.mkdir(parents=True, exist_ok=True)
    time_h = q23["time_s"] / 3600.0
    rows = [
        ("surface", crossing_time_hours(time_h, q23["moisture"][:, -1])),
        ("volume_mean", crossing_time_hours(time_h, q23["moisture_mean"])),
        ("global_max", crossing_time_hours(time_h, q23["moisture_global_max"])),
    ]
    crossing_path = source_dir / "q3_threshold_crossings.csv"
    with crossing_path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["criterion", "threshold", "first_crossing_h"])
        for name, value in rows:
            writer.writerow([name, 0.15, f"{value:.12f}"])

    factor_path = source_dir / "q4_factorial_used.csv"
    with factor_path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["scenario", "properties", "geometry", "critical_time_h"])
        for row in final_validation["scenario_uncertainty"]["Q4_factorial_counterfactuals"]:
            writer.writerow([row["letter"], row["group"], row["geometry"], f"{float(row['event_h']):.12f}"])
    return [str(crossing_path.relative_to(ROOT)), str(factor_path.relative_to(ROOT))]


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    verify_sources()
    configure_style()
    cfg = load_config()
    env, radius_data = manuscript.load_measurements(cfg)
    final_validation = manuscript.json_load("reports/final_result_validation.json")
    with np.load(ROOT / "results" / "final" / "q23_full_precision.npz") as q23, np.load(ROOT / "results" / "final" / "q4_full_precision.npz") as q4:
        outputs = {
            "fig03_q2_enriched": render_q2(cfg, env, q23),
            "fig04_q3_enriched": render_q3(q23, final_validation),
            "fig05_q4_enriched": render_q4(q4, radius_data, final_validation),
        }
        derived = export_derived_source_data(q23, final_validation)
    board = build_board(list(outputs))
    formal_checks = verify_formal_figures_unchanged()
    report = {
        "status": "generated_for_user_layout_selection",
        "backend": "Python/matplotlib only",
        "palette": "B 靛蓝鎏金",
        "official_figures_overwritten": False,
        "formal_figure_hashes_unchanged": all(formal_checks.values()),
        "formal_figure_checks": formal_checks,
        "comparison_board": str(board.relative_to(ROOT)),
        "contracts": CONTRACTS,
        "derived_source_data": derived,
        "qa": qa(outputs),
    }
    (OUT / "preview_manifest.json").write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": report["status"], "figures": len(outputs), "board": str(board), "formal_unchanged": all(formal_checks.values())}, ensure_ascii=False))


if __name__ == "__main__":
    main()
