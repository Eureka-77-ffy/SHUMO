#!/usr/bin/env python3
"""Render three palette candidates from accepted manuscript data.

Figure contract
---------------
Core conclusion: richer semantic color can strengthen hierarchy without
changing any data, axes, statistics, or manuscript claim.
Archetype: quantitative comparison grid.
Backend: Python/matplotlib only.
Output: separate SVG/PDF/600 dpi PNG/TIFF files plus a selection board.
Panel evidence: Q2 tests sequential heatmaps and trajectories; sensitivity
tests signed markers, bars, neutrals, and log-scale risk levels.
Reviewer risk: palette candidates must remain color-vision and print safe;
these previews never overwrite the accepted manuscript figures.
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image, ImageDraw, ImageFont

from figures.scripts import make_all_figures as manuscript
from src.data_io import load_config, verify_sources


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "figures" / "style_previews"


PALETTES = {
    "a_ocean_coral": {
        "title": "A  深海珊瑚｜现代期刊",
        "description": "冷暖对照清楚，热湿语义最直观；推荐作为全文基准。",
        "colors": {
            "temperature": "#C94C4C",
            "temperature_light": "#F09A8E",
            "moisture": "#2457A6",
            "moisture_light": "#A9C4EB",
            "radius": "#238B7E",
            "mean": "#6D75A8",
            "surface": "#65A98F",
            "threshold": "#27313A",
            "neutral": "#6F7782",
            "neutral_light": "#DDE3E9",
            "accent": "#8A5A9B",
            "positive": "#2A9D8F",
            "negative": "#E76F51",
        },
        "temperature_cmap": ["#FFF8F2", "#F7C6AA", "#E76F51", "#9B2226"],
        "moisture_cmap": ["#F5F9FD", "#C9DCF0", "#5B8FC8", "#173F7A"],
    },
    "b_indigo_gold": {
        "title": "B  靛蓝鎏金｜竞赛醒目",
        "description": "靛蓝稳重、金橙抓重点，打印对比最强，视觉更有冲击力。",
        "colors": {
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
            "positive": "#356FA1",
            "negative": "#D98B28",
        },
        "temperature_cmap": ["#FFF9E8", "#F1D18A", "#D77A32", "#8F2D2D"],
        "moisture_cmap": ["#F7F6FC", "#D2D1EE", "#7A82C2", "#30356F"],
    },
    "c_ink_crimson": {
        "title": "C  青黛绛红｜中式克制",
        "description": "青黛与绛红低饱和统一，论文气质最沉稳，长文阅读最舒适。",
        "colors": {
            "temperature": "#A63D4A",
            "temperature_light": "#D98B84",
            "moisture": "#315D7D",
            "moisture_light": "#A9C1D0",
            "radius": "#397D73",
            "mean": "#66788D",
            "surface": "#9A8B62",
            "threshold": "#292C2F",
            "neutral": "#74797C",
            "neutral_light": "#DCE1E0",
            "accent": "#725A78",
            "positive": "#397D73",
            "negative": "#C46A45",
        },
        "temperature_cmap": ["#FBF7F1", "#E8C0B5", "#C86C6A", "#742F3B"],
        "moisture_cmap": ["#F3F7F6", "#C5D8D5", "#6E9BA4", "#294F68"],
    },
}


def set_palette(spec: dict) -> None:
    manuscript.COLORS.clear()
    manuscript.COLORS.update(spec["colors"])
    manuscript.TEMP_CMAP = mpl.colors.LinearSegmentedColormap.from_list(
        "preview_temperature", spec["temperature_cmap"]
    )
    manuscript.MOIST_CMAP = mpl.colors.LinearSegmentedColormap.from_list(
        "preview_moisture", spec["moisture_cmap"]
    )
    manuscript.MOIST_CMAP.set_bad("white")


def render_candidates() -> list[dict]:
    verify_sources()
    cfg = load_config()
    env, _ = manuscript.load_measurements(cfg)
    q23 = np.load(ROOT / "results" / "final" / "q23_full_precision.npz")
    final_validation = manuscript.json_load("reports/final_result_validation.json")
    records = []
    for key, spec in PALETTES.items():
        folder = OUT / key
        folder.mkdir(parents=True, exist_ok=True)
        manuscript.FIG_DIR = folder
        manuscript.style()
        set_palette(spec)
        q2_files = manuscript.fig03(cfg, env, q23)
        sensitivity_files = manuscript.fig07(final_validation)
        records.append(
            {
                "id": key,
                "title": spec["title"],
                "description": spec["description"],
                "files": [str(path.relative_to(ROOT)) for path in q2_files + sensitivity_files],
            }
        )
    return records


def font(size: int):
    candidates = [
        "/System/Library/Fonts/Supplemental/Songti.ttc",
        "/System/Library/Fonts/Supplemental/Arial Unicode.ttf",
    ]
    for candidate in candidates:
        if Path(candidate).is_file():
            return ImageFont.truetype(candidate, size=size)
    return ImageFont.load_default()


def build_board() -> Path:
    width, row_height = 2600, 1020
    canvas = Image.new("RGB", (width, 170 + len(PALETTES) * row_height), "#F7F8FA")
    draw = ImageDraw.Draw(canvas)
    draw.text((80, 42), "论文配色样稿｜真实数据对照", fill="#20242A", font=font(48))
    draw.text((80, 104), "左：Q2热湿耦合　右：终点灵敏度与风险层级", fill="#606872", font=font(27))
    for row, (key, spec) in enumerate(PALETTES.items()):
        y = 170 + row * row_height
        draw.rounded_rectangle((45, y + 18, width - 45, y + row_height - 22), radius=28, fill="white", outline="#D9DEE5", width=2)
        draw.text((85, y + 48), spec["title"], fill="#20242A", font=font(38))
        draw.text((85, y + 102), spec["description"], fill="#606872", font=font(25))
        for col, stem in enumerate(("fig03_q2_coupling", "fig07_sensitivity_risk")):
            image = Image.open(OUT / key / f"{stem}.png").convert("RGB")
            image.thumbnail((1160, 790), Image.Resampling.LANCZOS)
            x = 85 + col * 1220 + (1160 - image.width) // 2
            image_y = y + 165 + (790 - image.height) // 2
            canvas.paste(image, (x, image_y))
    path = OUT / "palette_comparison.png"
    canvas.save(path, dpi=(220, 220))
    return path


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    records = render_candidates()
    board = build_board()
    manifest = {
        "status": "palette_candidates_for_user_selection",
        "backend": "Python/matplotlib",
        "official_figures_overwritten": False,
        "comparison_board": str(board.relative_to(ROOT)),
        "candidates": records,
    }
    (OUT / "preview_manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps({"candidates": len(records), "board": str(board)}, ensure_ascii=False))


if __name__ == "__main__":
    main()
