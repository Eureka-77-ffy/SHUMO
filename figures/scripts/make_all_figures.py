#!/usr/bin/env python3
"""Generate the seven argument-led manuscript figures from accepted full-precision data."""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.ticker import FuncFormatter, LogLocator, NullFormatter
from openpyxl import load_workbook

# Publication contract: editable SVG text and an explicit sans-serif stack.
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial", "DejaVu Sans", "Liberation Sans"]
plt.rcParams["svg.fonttype"] = "none"

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from src.data_io import digest, load_config, verify_sources  # noqa: E402
from src.properties import PropertyLaw  # noqa: E402


FIG_DIR = ROOT / "figures"
FORMATS = ("svg", "pdf", "png", "tiff")

COLORS = {
    "temperature": "#B64342",
    "temperature_light": "#E9A6A1",
    "moisture": "#0F4D92",
    "moisture_light": "#B4C0E4",
    "radius": "#42949E",
    "mean": "#7884B4",
    "surface": "#8BCF8B",
    "threshold": "#272727",
    "neutral": "#767676",
    "neutral_light": "#D8D8D8",
    "accent": "#9A4D8E",
    "positive": "#2E9E44",
    "negative": "#E53935",
}

TEMP_CMAP = mpl.colors.LinearSegmentedColormap.from_list(
    "paper_temperature", ["#fff7ec", "#fdd49e", "#ef6548", "#990000"]
)
MOIST_CMAP = mpl.colormaps["Blues"].copy()
MOIST_CMAP.set_bad("white")


CONTRACTS = {
    "fig01_inputs_scales": {
        "caption": "环境驱动与特征尺度",
        "section": "数据分析与尺度分析",
        "core_conclusion": "边界输入先趋于平台，而内部水分扩散与收缩持续更久，因此必须保留时变边界并区分热、质与几何尺度。",
        "archetype": "quantitative grid",
        "panel_map": {
            "a": "实测环境温度与环境含湿量",
            "b": "实测半径轨迹及Q4临界时刻",
            "c": "三套物性的初态径向热、质扩散尺度",
        },
        "reviewer_risk": "特征尺度仅用于量级判断，不解释为预测干燥时间。",
        "source_data": ["original_appendix_1", "original_appendix_2", "reports/scales_snapshot.csv"],
    },
    "fig02_q1_fields": {
        "caption": "问题一温湿场演化",
        "section": "问题一",
        "core_conclusion": "前30分钟热影响已进入内部，而显著失水仍集中在表层。",
        "archetype": "quantitative grid",
        "panel_map": {
            "a": "温度时空场",
            "b": "干基含水率时空场",
            "c": "典型温度径向剖面",
            "d": "典型含水率径向剖面",
        },
        "reviewer_risk": "热图与剖面来自同一未舍入解，但分别回答全局传播和局部梯度问题。",
        "source_data": ["results/final/q1_full_precision.npz"],
    },
    "fig03_q2_coupling": {
        "caption": "问题二耦合响应",
        "section": "问题二",
        "core_conclusion": "3小时内温度趋近环境平台，但水分梯度仍显著，温升促进与表层低含水抑制共同控制局部扩散率。",
        "archetype": "quantitative grid",
        "panel_map": {
            "a": "温度时空场",
            "b": "干基含水率时空场",
            "c": "中心、表面与环境温度响应",
            "d": "中心与表面局部扩散率",
        },
        "reviewer_risk": "扩散率是经验模型的局部诊断量，不作实验因果解释。",
        "source_data": ["results/final/q23_full_precision.npz", "config/model_config.json"],
    },
    "fig04_q3_endpoint": {
        "caption": "问题三全域干燥终点",
        "section": "问题三",
        "core_conclusion": "后期干燥由全域最大含水率的缓慢下降控制，临界阈值在57.4740小时首次到达。",
        "archetype": "asymmetric mixed-modality figure",
        "panel_map": {
            "a": "全时域水分场",
            "b": "全域最大值、体积平均值与表面值",
            "c": "典型径向剖面",
            "d": "临界点邻域及严格达标整数秒",
        },
        "reviewer_risk": "临界时刻对应等号，严格小于阈值从其后的时刻开始。",
        "source_data": ["results/final/q23_full_precision.npz", "reports/final_result_validation.json"],
    },
    "fig05_q4_shrinkage": {
        "caption": "问题四收缩域效应",
        "section": "问题四",
        "core_conclusion": "收缩同时缩短扩散路径并提高单位体积交换尺度；材料坐标与物理坐标给出一致场演化，且物性与几何作用不能混为一谈。",
        "archetype": "asymmetric mixed-modality figure",
        "panel_map": {
            "a": "观测半径与几何倍率",
            "b": "材料坐标中的水分场",
            "c": "物理坐标中的水分场与域外遮罩",
            "d": "物性和几何的四组受控对照",
        },
        "reviewer_risk": "物理域外保持空白；四组对照是模型情景，不是实验因果识别。",
        "source_data": ["results/final/q4_full_precision.npz", "tables/q4_factorial.csv", "original_appendix_2"],
    },
    "fig06_numerical_validation": {
        "caption": "数值收敛与独立验证",
        "section": "模型检验",
        "core_conclusion": "空间、时间和独立实现差异均收敛到远小于主要情景差异的水平，离散守恒残差满足预设阈值。",
        "archetype": "quantitative grid",
        "panel_map": {
            "a": "Q2温度空间收敛",
            "b": "Q2水分空间收敛",
            "c": "Q3/Q4临界时间空间收敛",
            "d": "四问独立实现的场值差异",
            "e": "四问归一化守恒残差",
        },
        "reviewer_risk": "所有差异均为观测到的数值差，而非严格误差上界或物理精度证明。",
        "source_data": ["reports/q2_convergence.json", "reports/q3_convergence.json", "reports/q4_convergence.json", "reports/final_result_validation.json"],
    },
    "fig07_sensitivity_risk": {
        "caption": "终点灵敏度与风险层级",
        "section": "灵敏度分析",
        "core_conclusion": "Q3和Q4终点均对扩散率尺度最敏感、对换热系数近乎不敏感；边界和结构性情景差异显著大于纯数值误差。",
        "archetype": "quantitative grid",
        "panel_map": {
            "a": "Q3正负10%单因素情景",
            "b": "Q4正负10%单因素情景",
            "c": "Q3和Q4的有符号无量纲局部灵敏度系数",
            "d": "数值、边界、端面、参数及结构情景的影响量级",
        },
        "reviewer_risk": "单因素局部扰动与情景包络均未赋予概率分布，不应解释为全局敏感性或置信区间。",
        "source_data": ["reports/q3_sensitivity.json", "reports/q4_sensitivity.json", "reports/endface_diagnostics.json", "reports/final_result_validation.json"],
    },
}


def style() -> None:
    mpl.rcParams.update(
        {
            # Songti SC is the installed SimSun-compatible Chinese typeface.
            # It is selected globally so Matplotlib does not silently keep Arial
            # for strings that contain CJK glyphs.
            "font.family": "Songti SC",
            "font.size": 7.3,
            "axes.titlesize": 8.1,
            "axes.labelsize": 7.5,
            "xtick.labelsize": 6.7,
            "ytick.labelsize": 6.7,
            "legend.fontsize": 6.6,
            "axes.spines.right": False,
            "axes.spines.top": False,
            "axes.linewidth": 0.8,
            "xtick.major.width": 0.7,
            "ytick.major.width": 0.7,
            "xtick.major.size": 3.0,
            "ytick.major.size": 3.0,
            "legend.frameon": False,
            "axes.unicode_minus": False,
            "mathtext.fontset": "stix",
            "pdf.fonttype": 42,
            "savefig.facecolor": "white",
        }
    )


def json_load(relative: str):
    return json.loads((ROOT / relative).read_text(encoding="utf-8"))


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def panel(ax, label: str, title: str) -> None:
    ax.text(-0.13, 1.04, label, transform=ax.transAxes, ha="left", va="bottom", fontweight="bold", fontsize=9)
    ax.set_title(title, loc="left", pad=4)


def grid(ax, axis: str = "y") -> None:
    ax.grid(axis=axis, color="#E6E6E6", lw=0.55, zorder=0)
    ax.set_axisbelow(True)


def ascii_log_axis(ax, axis: str = "y") -> None:
    """Use unambiguous ASCII labels; the Chinese font lacks U+2212."""
    def formatter(value, _position):
        if value <= 0 or not np.isfinite(value):
            return ""
        exponent = int(np.floor(np.log10(value) + 1e-12))
        coefficient = value / (10.0 ** exponent)
        if 0.01 <= value < 10000:
            return f"{value:g}"
        if abs(coefficient - 1.0) < 1e-8:
            return f"1e{exponent:d}"
        return f"{coefficient:g}e{exponent:d}"

    target = ax.yaxis if axis == "y" else ax.xaxis
    target.set_major_locator(LogLocator(base=10, subs=(1.0, 2.0, 5.0), numticks=12))
    target.set_major_formatter(FuncFormatter(formatter))
    target.set_minor_formatter(NullFormatter())


def centers_to_edges(values: np.ndarray, lower=None, upper=None) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    if values.size == 1:
        edges = np.array([values[0] - 0.5, values[0] + 0.5])
    else:
        edges = np.r_[values[0] - 0.5 * (values[1] - values[0]), 0.5 * (values[:-1] + values[1:]), values[-1] + 0.5 * (values[-1] - values[-2])]
    if lower is not None:
        edges[0] = lower
    if upper is not None:
        edges[-1] = upper
    return edges


def heatmap(ax, time, space, values, cmap, label, *, vmin=None, vmax=None):
    x = centers_to_edges(np.asarray(time), lower=float(np.min(time)), upper=float(np.max(time)))
    y = centers_to_edges(np.asarray(space), lower=float(np.min(space)), upper=float(np.max(space)))
    mesh = ax.pcolormesh(x, y, np.asarray(values).T, cmap=cmap, shading="auto", vmin=vmin, vmax=vmax, rasterized=True)
    cbar = ax.figure.colorbar(mesh, ax=ax, pad=0.018, fraction=0.047)
    cbar.set_label(label)
    cbar.outline.set_linewidth(0.6)
    return mesh


def save_figure(fig, stem: str):
    outputs = []
    for fmt in FORMATS:
        path = FIG_DIR / f"{stem}.{fmt}"
        kwargs = {"bbox_inches": "tight", "facecolor": "white"}
        if fmt in {"png", "tiff"}:
            kwargs["dpi"] = 600
        if fmt == "tiff":
            kwargs["pil_kwargs"] = {"compression": "tiff_lzw"}
        fig.savefig(path, **kwargs)
        outputs.append(path)
    plt.close(fig)
    return outputs


def load_measurements(cfg):
    source = Path(cfg["source_directory"])
    wb = load_workbook(source / "附件1.xlsx", read_only=True, data_only=True)
    env = np.asarray(list(wb.active.values)[1:], dtype=float)
    wb.close()
    wb = load_workbook(source / "附件2.xlsx", read_only=True, data_only=True)
    radius = np.asarray(list(wb.active.values)[1:], dtype=float)
    wb.close()
    return env, radius


def fig01(cfg, env, radius):
    fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.45), layout="constrained")

    ax = axes[0]
    t_h = env[:, 0] / 3600
    l1 = ax.plot(t_h, env[:, 1], color=COLORS["temperature"], lw=1.5, label="环境温度")[0]
    ax.set(xlabel="时间 / h", ylabel="温度 / ℃", xlim=(0, 4))
    ax2 = ax.twinx()
    l2 = ax2.plot(t_h, env[:, 2], color=COLORS["moisture"], lw=1.35, label="环境含湿量")[0]
    ax2.set_ylabel("环境含湿量 / (kg/kg)", color=COLORS["moisture"])
    ax2.tick_params(axis="y", colors=COLORS["moisture"])
    ax.legend([l1, l2], ["环境温度", "环境含湿量"], loc="lower right")
    panel(ax, "a", "边界输入")

    ax = axes[1]
    rt_h, r_cm = radius[:, 0] / 3600, radius[:, 1]
    ax.plot(rt_h, r_cm, color=COLORS["radius"], lw=1.3)
    ax.scatter(rt_h, r_cm, s=7, color=COLORS["radius"], edgecolor="white", lw=0.25, zorder=3)
    q4_h = 51.09202464437358
    q4_r = float(np.interp(q4_h * 3600, radius[:, 0], radius[:, 1]))
    ax.scatter([q4_h], [q4_r], s=30, marker="*", color=COLORS["negative"], zorder=4)
    ax.annotate(f"Q4临界点\n{q4_h:.4f} h, {q4_r:.3f} cm", (q4_h, q4_r), xytext=(-44, 16), textcoords="offset points", arrowprops={"arrowstyle": "-", "lw": 0.7, "color": COLORS["neutral"]})
    ax.set(xlabel="时间 / h", ylabel="半径 / cm", xlim=(0, 72), ylim=(1.15, 2.05))
    grid(ax)
    panel(ax, "b", "观测收缩")

    scale = np.genfromtxt(ROOT / "reports/scales_snapshot.csv", delimiter=",", names=True, dtype=None, encoding="utf-8")
    rows = [row for row in scale if row["C"] == 2.55 and row["T_C"] == 28]
    groups = ["Q1", "Q2/Q3", "Q4"]
    heat_h = np.array([float(row["tau_T_min"]) / 60 for row in rows])
    water_h = np.array([float(row["tau_C_h"]) for row in rows])
    x = np.arange(3)
    w = 0.34
    ax = axes[2]
    ax.bar(x - w / 2, heat_h, w, color=COLORS["temperature_light"], edgecolor=COLORS["temperature"], lw=0.7, label="热扩散尺度")
    ax.bar(x + w / 2, water_h, w, color=COLORS["moisture_light"], edgecolor=COLORS["moisture"], lw=0.7, label="水分扩散尺度")
    ax.set_yscale("log")
    ascii_log_axis(ax)
    ax.set_xticks(x, groups)
    ax.set(xlabel="物性体系", ylabel="初态径向尺度 / h", ylim=(0.3, 180))
    ax.legend(loc="upper left")
    grid(ax)
    panel(ax, "c", "热—质时间尺度分离")
    return save_figure(fig, "fig01_inputs_scales")


def fig02(q1):
    fig, axes = plt.subplots(2, 2, figsize=(7.2, 5.3), layout="constrained")
    time_min = q1["time_s"] / 60
    radius_cm = q1["radius_m"] * 100
    heatmap(axes[0, 0], time_min, radius_cm, q1["temperature"], TEMP_CMAP, "温度 / ℃")
    axes[0, 0].set(xlabel="时间 / min", ylabel="距中心距离 / cm")
    panel(axes[0, 0], "a", "温度时空传播")
    heatmap(axes[0, 1], time_min, radius_cm, q1["moisture"], MOIST_CMAP, "干基含水率 / (kg/kg)")
    axes[0, 1].set(xlabel="时间 / min", ylabel="距中心距离 / cm")
    panel(axes[0, 1], "b", "水分时空迁移")

    selected = [0, 5, 15, 30]
    line_colors = mpl.colormaps["viridis"](np.linspace(0.18, 0.88, len(selected)))
    for minute, color in zip(selected, line_colors):
        i = int(minute * 60)
        axes[1, 0].plot(radius_cm, q1["temperature"][i], color=color, lw=1.45, label=f"{minute} min")
        axes[1, 1].plot(radius_cm, q1["moisture"][i], color=color, lw=1.45, label=f"{minute} min")
    axes[1, 0].set(xlabel="距中心距离 / cm", ylabel="温度 / ℃", xlim=(0, 2))
    axes[1, 1].set(xlabel="距中心距离 / cm", ylabel="干基含水率 / (kg/kg)", xlim=(0, 2))
    axes[1, 0].legend(ncol=2)
    axes[1, 1].legend(ncol=2)
    grid(axes[1, 0])
    grid(axes[1, 1])
    panel(axes[1, 0], "c", "温度径向剖面")
    panel(axes[1, 1], "d", "含水率径向剖面")
    return save_figure(fig, "fig02_q1_fields")


def fig03(cfg, env, q23):
    end = int(np.searchsorted(q23["time_s"], 10800, side="right"))
    sample = np.arange(0, end, 10)
    if sample[-1] != end - 1:
        sample = np.r_[sample, end - 1]
    time_h = q23["time_s"][:end] / 3600
    radius_cm = q23["radius_m"] * 100
    fig, axes = plt.subplots(2, 2, figsize=(7.2, 5.3), layout="constrained")
    heatmap(axes[0, 0], time_h[sample], radius_cm, q23["temperature"][:end][sample], TEMP_CMAP, "温度 / ℃")
    axes[0, 0].set(xlabel="时间 / h", ylabel="距中心距离 / cm")
    panel(axes[0, 0], "a", "温度场")
    heatmap(axes[0, 1], time_h[sample], radius_cm, q23["moisture"][:end][sample], MOIST_CMAP, "干基含水率 / (kg/kg)")
    axes[0, 1].set(xlabel="时间 / h", ylabel="距中心距离 / cm")
    panel(axes[0, 1], "b", "含水率场")

    ax = axes[1, 0]
    stride = 30
    t = time_h[::stride]
    ax.plot(t, q23["temperature"][:end:stride, 0], color=COLORS["temperature"], lw=1.45, label="中心")
    ax.plot(t, q23["temperature"][:end:stride, -1], color=COLORS["temperature_light"], lw=1.35, label="表面")
    env_t = np.interp(q23["time_s"][:end:stride], env[:, 0], env[:, 1])
    ax.plot(t, env_t, color=COLORS["neutral"], lw=1.15, ls="--", label="环境")
    ax.set(xlabel="时间 / h", ylabel="温度 / ℃", xlim=(0, 3))
    ax.legend(ncol=3)
    grid(ax)
    panel(ax, "c", "内外温度响应")

    prop = PropertyLaw(cfg, "q23")
    _, _, d_center = prop.evaluate(q23["moisture"][:end:stride, 0], q23["temperature"][:end:stride, 0])
    _, _, d_surface = prop.evaluate(q23["moisture"][:end:stride, -1], q23["temperature"][:end:stride, -1])
    ax = axes[1, 1]
    ax.plot(t, d_center, color=COLORS["moisture"], lw=1.5, label="中心")
    ax.plot(t, d_surface, color=COLORS["surface"], lw=1.4, label="表面")
    ax.set_yscale("log")
    ascii_log_axis(ax)
    ax.set(xlabel="时间 / h", ylabel="局部扩散率 D / (m²/s)", xlim=(0, 3))
    ax.legend()
    grid(ax)
    panel(ax, "d", "温湿耦合后的扩散率")
    return save_figure(fig, "fig03_q2_coupling")


def fig04(q23, final_validation):
    fig = plt.figure(figsize=(7.2, 5.6), layout="constrained")
    gs = fig.add_gridspec(2, 4, width_ratios=[1.25, 1.25, 1, 1])
    ax_a = fig.add_subplot(gs[0, :2])
    ax_b = fig.add_subplot(gs[0, 2:])
    ax_c = fig.add_subplot(gs[1, :2])
    ax_d = fig.add_subplot(gs[1, 2:])
    event_h = float(final_validation["critical_results"]["Q3"]["critical_time_h"])
    time_h = q23["time_s"] / 3600
    radius_cm = q23["radius_m"] * 100
    sample = np.unique(np.r_[np.arange(0, len(time_h), 60), len(time_h) - 1])
    heatmap(ax_a, time_h[sample], radius_cm, q23["moisture"][sample], MOIST_CMAP, "干基含水率 / (kg/kg)", vmin=0.05, vmax=2.55)
    ax_a.axvline(event_h, color="white", lw=1.0, ls="--")
    ax_a.set(xlabel="时间 / h", ylabel="距中心距离 / cm")
    panel(ax_a, "a", "全时域水分场")

    stride = 60
    ax_b.plot(time_h[::stride], q23["moisture_global_max"][::stride], color=COLORS["moisture"], lw=1.5, label="全域最大值")
    ax_b.plot(time_h[::stride], q23["moisture_mean"][::stride], color=COLORS["mean"], lw=1.25, label="体积平均值")
    ax_b.plot(time_h[::stride], q23["moisture"][::stride, -1], color=COLORS["surface"], lw=1.2, label="表面值")
    ax_b.axhline(0.15, color=COLORS["threshold"], lw=1.0, ls="--", label="阈值")
    ax_b.axvline(event_h, color=COLORS["negative"], lw=0.9, ls=":")
    ax_b.set(xlabel="时间 / h", ylabel="干基含水率 / (kg/kg)", xlim=(0, event_h))
    ax_b.legend(ncol=2)
    grid(ax_b)
    panel(ax_b, "b", "全域判据与统计量")

    selected_h = [3, 12, 24, 36, 48, event_h]
    line_colors = mpl.colormaps["Blues"](np.linspace(0.28, 0.92, len(selected_h)))
    for hour, color in zip(selected_h, line_colors):
        i = int(np.argmin(np.abs(time_h - hour)))
        label = "临界" if hour == event_h else f"{hour:g} h"
        ax_c.plot(radius_cm, q23["moisture"][i], color=color, lw=1.35, label=label)
    ax_c.axhline(0.15, color=COLORS["threshold"], lw=0.9, ls="--")
    ax_c.set(xlabel="距中心距离 / cm", ylabel="干基含水率 / (kg/kg)", xlim=(0, 2))
    ax_c.legend(ncol=3)
    grid(ax_c)
    panel(ax_c, "c", "径向梯度的长期演化")

    event_s = float(final_validation["critical_results"]["Q3"]["critical_time_s"])
    strict_s = float(final_validation["critical_results"]["Q3"]["first_strict_integer_second"])
    mask = q23["time_s"] >= event_s - 3600
    ax_d.plot((q23["time_s"][mask] - event_s) / 60, 1e3 * (q23["moisture_global_max"][mask] - 0.15), color=COLORS["moisture"], lw=1.5)
    ax_d.axhline(0, color=COLORS["threshold"], lw=0.9, ls="--")
    ax_d.axvline(0, color=COLORS["negative"], lw=0.9, ls=":")
    ax_d.scatter([(strict_s - event_s) / 60], [np.interp(strict_s, q23["time_s"], q23["moisture_global_max"]) * 1e3 - 150], s=20, color=COLORS["positive"], zorder=3, label="首个严格达标整数秒")
    ax_d.set(xlabel="相对临界时刻 / min", ylabel="1000(Cmax-0.15)", xlim=(-60, 0.5))
    ax_d.legend(loc="upper right")
    grid(ax_d)
    panel(ax_d, "d", f"临界邻域（{event_h:.4f} h）")
    return save_figure(fig, "fig04_q3_endpoint")


def interpolate_material_profiles(q4, radius_data):
    times = q4["profile_time_s"]
    states = q4["profile_state"]
    moisture = states[:, 1::2]
    xi_cells = q4["cell_centres_reference_m"] / 0.02
    xi_grid = np.linspace(0, 1, 241)
    material = np.empty((len(times), len(xi_grid)))
    physical_r_cm = np.linspace(0, 2, 241)
    physical = np.full((len(times), len(physical_r_cm)), np.nan)
    radii_m = np.interp(times, radius_data[:, 0], radius_data[:, 1] * 0.01)
    true_surface = np.interp(times, q4["time_s"], q4["moisture"][:, -1])
    for i, cells in enumerate(moisture):
        x0, x1 = xi_cells[:2]
        axis_value = (x1 * x1 * cells[0] - x0 * x0 * cells[1]) / (x1 * x1 - x0 * x0)
        nodes_x = np.r_[0.0, xi_cells, 1.0]
        nodes_c = np.r_[axis_value, cells, true_surface[i]]
        material[i] = np.interp(xi_grid, nodes_x, nodes_c)
        valid = physical_r_cm <= radii_m[i] * 100 + 1e-12
        xi = physical_r_cm[valid] / (radii_m[i] * 100)
        physical[i, valid] = np.interp(xi, nodes_x, nodes_c)
    return times, xi_grid, material, physical_r_cm, physical, radii_m


def fig05(q4, radius_data, final_validation):
    fig = plt.figure(figsize=(7.2, 5.7), layout="constrained")
    gs = fig.add_gridspec(2, 4, width_ratios=[1, 1, 1.12, 1.12])
    ax_a = fig.add_subplot(gs[0, :2])
    ax_b = fig.add_subplot(gs[0, 2:])
    ax_c = fig.add_subplot(gs[1, :2])
    ax_d = fig.add_subplot(gs[1, 2:])

    event_h = float(final_validation["critical_results"]["Q4"]["critical_time_h"])
    t_h = q4["time_s"] / 3600
    radius_cm = q4["radius_current_m"] * 100
    ax_a.plot(t_h, radius_cm, color=COLORS["radius"], lw=1.5, label="半径")
    ax_a.set(xlabel="时间 / h", ylabel="半径 / cm", xlim=(0, event_h), ylim=(1.15, 2.05))
    ax_a2 = ax_a.twinx()
    multiplier = (2.0 / radius_cm) ** 2
    ax_a2.plot(t_h, multiplier, color=COLORS["accent"], lw=1.2, ls="--", label="R0²/R²")
    ax_a2.set_ylabel("扩散几何倍率 R0²/R²", color=COLORS["accent"])
    ax_a2.tick_params(axis="y", colors=COLORS["accent"])
    ax_a.axvline(event_h, color=COLORS["negative"], lw=0.9, ls=":")
    handles = [Line2D([0], [0], color=COLORS["radius"], lw=1.5), Line2D([0], [0], color=COLORS["accent"], lw=1.2, ls="--")]
    ax_a.legend(handles, ["半径", "扩散几何倍率"], loc="center right")
    grid(ax_a)
    panel(ax_a, "a", "收缩驱动与几何放大")

    pt, xi, material, pr, physical, radii_m = interpolate_material_profiles(q4, radius_data)
    ph = pt / 3600
    heatmap(ax_b, ph, xi, material, MOIST_CMAP, "干基含水率 / (kg/kg)", vmin=0.05, vmax=2.55)
    ax_b.set(xlabel="时间 / h", ylabel="材料坐标 ξ", xlim=(0, event_h))
    panel(ax_b, "b", "材料坐标场")

    heatmap(ax_c, ph, pr, np.ma.masked_invalid(physical), MOIST_CMAP, "干基含水率 / (kg/kg)", vmin=0.05, vmax=2.55)
    ax_c.plot(ph, radii_m * 100, color=COLORS["threshold"], lw=1.0, label="当前表面")
    ax_c.set(xlabel="时间 / h", ylabel="物理半径 r / cm", xlim=(0, event_h), ylim=(0, 2))
    ax_c.legend(loc="upper right")
    panel(ax_c, "c", "物理坐标场（域外留白）")

    factor = final_validation["scenario_uncertainty"]["Q4_factorial_counterfactuals"]
    labels = [f"{x['letter']}\n{x['group'].upper()}\n{'收缩' if x['geometry']=='linear' else '定径'}" for x in factor]
    values = [x["event_h"] for x in factor]
    colors = [COLORS["mean"], COLORS["moisture"], COLORS["temperature_light"], COLORS["temperature"]]
    bars = ax_d.bar(np.arange(4), values, color=colors, edgecolor="#4D4D4D", lw=0.6, zorder=2)
    for bar, value in zip(bars, values):
        ax_d.text(bar.get_x() + bar.get_width() / 2, value + 3, f"{value:.2f}", ha="center", va="bottom", fontsize=6.5)
    ax_d.set_xticks(np.arange(4), labels)
    ax_d.set(xlabel="受控模型情景", ylabel="临界时间 / h", ylim=(0, 145))
    grid(ax_d)
    panel(ax_d, "d", "物性×几何四组对照")
    return save_figure(fig, "fig05_q4_shrinkage")


def mesh_differences(report, quantity):
    x, y = [], []
    for row in report["mesh_series"]:
        diff = row.get("difference_from_previous", {}).get(quantity)
        if diff is None:
            diff = row.get("difference_from_previous", {}).get("fields", {}).get(quantity)
        if diff:
            x.append(row["n"])
            y.append(abs(float(diff["max_abs"])))
    return np.asarray(x), np.asarray(y)


def fig06(final_validation):
    q2c = json_load("reports/q2_convergence.json")
    q3c = json_load("reports/q3_convergence.json")
    q4c = json_load("reports/q4_convergence.json")
    fig = plt.figure(figsize=(7.2, 5.5), layout="constrained")
    gs = fig.add_gridspec(2, 6)
    axes = [fig.add_subplot(gs[0, 0:2]), fig.add_subplot(gs[0, 2:4]), fig.add_subplot(gs[0, 4:6]), fig.add_subplot(gs[1, 0:3]), fig.add_subplot(gs[1, 3:6])]

    for ax, quantity, title, ylabel, color in [
        (axes[0], "temperature", "温度空间收敛", "相邻网格最大差 / ℃", COLORS["temperature"]),
        (axes[1], "moisture", "水分空间收敛", "相邻网格最大差 / (kg/kg)", COLORS["moisture"]),
    ]:
        x, y = mesh_differences(q2c, quantity)
        ax.loglog(x, y, "o-", color=color, lw=1.35, ms=3.8, label="Q2")
        if len(x) >= 2:
            ref = y[-1] * (x[-1] / x) ** 2
            ax.loglog(x, ref, ls="--", color=COLORS["neutral"], lw=0.9, label="二阶参考")
        ax.axhline(5e-5, color=COLORS["threshold"], ls=":", lw=0.8, label="5e-5检查尺度")
        ax.set(xlabel="径向控制体数 N", ylabel=ylabel)
        ascii_log_axis(ax, "x")
        ascii_log_axis(ax, "y")
        ax.legend()
        grid(ax, "both")
        panel(ax, "a" if quantity == "temperature" else "b", title)

    ax = axes[2]
    for report, label, color, marker in [(q3c, "Q3", COLORS["moisture"], "o"), (q4c, "Q4", COLORS["temperature"], "s")]:
        xs, ys = [], []
        rows = report["mesh_series"]
        for i in range(1, len(rows)):
            delta = rows[i].get("event_delta_s")
            if delta is None:
                delta = rows[i]["event_time_s"] - rows[i - 1]["event_time_s"]
            xs.append(rows[i]["n"])
            ys.append(abs(float(delta)))
        ax.loglog(xs, ys, marker=marker, color=color, lw=1.3, ms=3.8, label=label)
    ax.axhline(0.18, color=COLORS["threshold"], ls=":", lw=0.8, label="四位小时检查尺度")
    ax.set(xlabel="径向控制体数 N", ylabel="相邻网格临界时间差 / s")
    ascii_log_axis(ax, "x")
    ascii_log_axis(ax, "y")
    ax.legend(ncol=2)
    grid(ax, "both")
    panel(ax, "c", "终点时间空间收敛")

    refs = final_validation["independent_references"]
    qlabels = ["Q1", "Q2", "Q3", "Q4"]
    tvals = [x["comparisons"]["temperature"]["primary_unrounded_vs_reference_max_abs"] for x in refs]
    cvals = [x["comparisons"]["moisture"]["primary_unrounded_vs_reference_max_abs"] for x in refs]
    x = np.arange(4)
    w = 0.34
    ax = axes[3]
    ax.bar(x - w / 2, tvals, w, color=COLORS["temperature_light"], edgecolor=COLORS["temperature"], lw=0.6, label="温度")
    ax.bar(x + w / 2, cvals, w, color=COLORS["moisture_light"], edgecolor=COLORS["moisture"], lw=0.6, label="含水率")
    ax.axhline(5e-5, color=COLORS["threshold"], ls=":", lw=0.8, label="5e-5检查尺度")
    ax.set_yscale("log")
    ascii_log_axis(ax)
    ax.set_xticks(x, qlabels)
    ax.set(xlabel="问题", ylabel="与独立有限元的最大绝对差 / 相应场单位", ylim=(1e-7, 1e-4))
    ax.legend(ncol=3)
    grid(ax)
    panel(ax, "d", "独立实现交叉核验")

    rows = final_validation["conservation"]["rows"]
    qlabels = ["Q1", "Q2", "Q3", "Q4"]
    water = [x["water_relative"] for x in rows]
    heat = [x["heat_relative"] for x in rows]
    ax = axes[4]
    ax.bar(x - w / 2, water, w, color=COLORS["moisture_light"], edgecolor=COLORS["moisture"], lw=0.6, label="水分")
    ax.bar(x + w / 2, heat, w, color=COLORS["temperature_light"], edgecolor=COLORS["temperature"], lw=0.6, label="有效显热")
    ax.axhline(1e-6, color=COLORS["threshold"], ls=":", lw=0.8, label="验收阈值")
    ax.set_yscale("log")
    ascii_log_axis(ax)
    ax.set_xticks(x, qlabels)
    ax.set(xlabel="问题", ylabel="归一化累计守恒残差", ylim=(1e-12, 3e-6))
    ax.legend(ncol=3)
    grid(ax)
    panel(ax, "e", "离散守恒检查")
    return save_figure(fig, "fig06_numerical_validation")


def sensitivity_map(report, prefix):
    out = {}
    for row in report["cases"]:
        label = row["label"]
        if not label.startswith(prefix + "_"):
            continue
        for param in ("h", "hm", "D", "beta"):
            marker = f"_{param}_multiplier_"
            if marker in label:
                mult = float(label.rsplit("_", 1)[-1])
                out.setdefault(param, {})[mult] = float(row["event_time_s"])
    return out


def tornado(ax, report, prefix, baseline_s, letter, title):
    values = sensitivity_map(report, prefix)
    params = ["D", "hm", "beta", "h"]
    labels = ["扩散率 $D$", "传质系数 $h_m$", "边界映射 $\\beta$", "换热系数 $h$"]
    y = np.arange(len(params))[::-1]
    for yi, p in zip(y, params):
        lo = 100 * (values[p][0.9] / baseline_s - 1)
        hi = 100 * (values[p][1.1] / baseline_s - 1)
        ax.plot([lo, hi], [yi, yi], color=COLORS["neutral_light"], lw=5, solid_capstyle="round", zorder=1)
        ax.scatter([lo], [yi], color=COLORS["negative"], marker="v", s=22, zorder=2, label="参数×0.9" if yi == y[0] else None)
        ax.scatter([hi], [yi], color=COLORS["positive"], marker="^", s=22, zorder=2, label="参数×1.1" if yi == y[0] else None)
    ax.axvline(0, color=COLORS["threshold"], lw=0.8)
    ax.set_yticks(y, labels)
    ax.set_xlabel("临界时间相对变化 / %")
    ax.legend(ncol=2, loc="lower right")
    grid(ax, "x")
    panel(ax, letter, title)


def fig07(final_validation):
    q3s = json_load("reports/q3_sensitivity.json")
    q4s = json_load("reports/q4_sensitivity.json")
    endface = json_load("reports/endface_diagnostics.json")
    base3 = float(final_validation["critical_results"]["Q3"]["critical_time_s"])
    base4 = float(final_validation["critical_results"]["Q4"]["critical_time_s"])
    fig = plt.figure(figsize=(7.2, 6.25), layout="constrained")
    gs = fig.add_gridspec(2, 2, height_ratios=[1, 1.12])
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[1, 0])
    ax_d = fig.add_subplot(gs[1, 1])
    tornado(ax_a, q3s, "q3", base3, "a", "Q3单因素敏感性")
    tornado(ax_b, q4s, "q4", base4, "b", "Q4单因素敏感性")

    q3_values = sensitivity_map(q3s, "q3")
    q4_values = sensitivity_map(q4s, "q4")
    params = ["D", "hm", "beta", "h"]
    labels_s = ["扩散率 $D$", "传质系数 $h_m$", "边界映射 $\\beta$", "换热系数 $h$"]
    sens3 = [(q3_values[p][1.1] - q3_values[p][0.9]) / (0.2 * base3) for p in params]
    sens4 = [(q4_values[p][1.1] - q4_values[p][0.9]) / (0.2 * base4) for p in params]
    y_s = np.arange(len(params))[::-1]
    offset = 0.16
    ax_c.barh(y_s + offset, sens3, height=0.28, color=COLORS["moisture"], label="Q3")
    ax_c.barh(y_s - offset, sens4, height=0.28, color=COLORS["temperature"], label="Q4")
    ax_c.axvline(0, color=COLORS["threshold"], lw=0.8)
    ax_c.set_yticks(y_s, labels_s)
    ax_c.set(xlabel="无量纲局部灵敏度系数 $S_p$", xlim=(-1.02, 0.12))
    ax_c.legend(ncol=2, loc="lower right")
    grid(ax_c, "x")
    panel(ax_c, "c", "扰动方向与相对强度")

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
    q3 = [numerical3, np.nan, ef3, boundary3, param3, structural3]
    q4 = [numerical4, radius4, ef4, boundary4, param4, structural4]
    y = np.arange(len(labels))[::-1]
    ax_d.scatter(q3, y + 0.13, color=COLORS["moisture"], s=25, marker="o", label="Q3", zorder=3)
    ax_d.scatter(q4, y - 0.13, color=COLORS["temperature"], s=25, marker="s", label="Q4", zorder=3)
    for yi, v3, v4 in zip(y, q3, q4):
        finite = [v for v in (v3, v4) if np.isfinite(v)]
        ax_d.plot([min(finite), max(finite)], [yi, yi], color=COLORS["neutral_light"], lw=1.0, zorder=1)
    ax_d.set_xscale("log")
    ascii_log_axis(ax_d, "x")
    ax_d.xaxis.set_major_locator(LogLocator(base=10, subs=(1.0,), numticks=9))
    ax_d.set_yticks(y, labels)
    ax_d.set(xlabel="临界时间最大绝对情景差 / s", xlim=(1e-2, 1e6))
    ax_d.legend(ncol=2)
    grid(ax_d, "x")
    panel(ax_d, "d", "影响量级分层（非置信区间）")
    return save_figure(fig, "fig07_sensitivity_risk")


def build_contact_sheet(paths):
    from PIL import Image, ImageDraw

    previews = [Image.open(next(p for p in group if p.suffix == ".png")).convert("RGB") for group in paths]
    width = 1500
    thumb_w = 720
    thumb_h = 560
    rows = (len(previews) + 1) // 2
    canvas = Image.new("RGB", (width, rows * (thumb_h + 55) + 20), "white")
    draw = ImageDraw.Draw(canvas)
    for i, image in enumerate(previews):
        image.thumbnail((thumb_w, thumb_h))
        x = 20 + (i % 2) * 740
        y = 35 + (i // 2) * (thumb_h + 55)
        canvas.paste(image, (x + (thumb_w - image.width) // 2, y))
        draw.text((x, 10 + (i // 2) * (thumb_h + 55)), f"Fig. {i + 1}", fill="black")
    path = FIG_DIR / "contact_sheet.png"
    canvas.save(path, dpi=(180, 180))
    return path


def main():
    style()
    inputs = verify_sources()
    cfg = load_config()
    env, radius = load_measurements(cfg)
    q1 = np.load(ROOT / "results/final/q1_full_precision.npz")
    q23 = np.load(ROOT / "results/final/q23_full_precision.npz")
    q4 = np.load(ROOT / "results/final/q4_full_precision.npz")
    final_validation = json_load("reports/final_result_validation.json")

    generated = [
        fig01(cfg, env, radius),
        fig02(q1),
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
            "contract": CONTRACTS[stem],
        }
    manifest = {
        "status": "generated_pending_visual_qa",
        "backend": "Python 3.9 / matplotlib 3.9.4",
        "figure_width_mm": 182.88,
        "raster_dpi": 600,
        "svg_text_editable": True,
        "font_policy": "Arial for Latin; Songti SC as the installed SimSun-compatible Chinese font",
        "heatmap_integrity": "No smoothing or invented values; display-only temporal sampling is declared in the script; Q4 physical-domain exterior is masked white.",
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
    # Keep this long-standing entry point stable while making the accepted
    # information-enriched B-palette renderer the canonical implementation.
    from figures.scripts.make_formal_enriched_figures import main as enriched_main

    enriched_main()
