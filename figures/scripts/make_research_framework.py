#!/usr/bin/env python3
"""Generate a restrained, editable conference-style research framework."""

from __future__ import annotations

import json
from pathlib import Path
import xml.etree.ElementTree as ET


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "figures" / "fig00a_research_roadmap.drawio"
META = ROOT / "figures" / "fig00a_research_roadmap.json"

FONT = "Microsoft YaHei,PingFang SC,Songti SC,SimSun,Helvetica"
NAVY = "#314E73"
NAVY_DARK = "#243A57"
NAVY_LIGHT = "#EEF3F8"
GOLD = "#B88A3B"
INK = "#252A32"
MID = "#697482"
LINE = "#AEB7C2"
PALE = "#F6F7F9"


def vertex(root, cid, value, x, y, w, h, style):
    cell = ET.SubElement(
        root, "mxCell", id=cid, value=value, style=style, vertex="1", parent="1"
    )
    ET.SubElement(
        cell,
        "mxGeometry",
        x=str(x),
        y=str(y),
        width=str(w),
        height=str(h),
        **{"as": "geometry"},
    )
    return cell


def edge(root, cid, points, color=LINE, width=1.6, dashed=False, arrow=True):
    style = (
        "edgeStyle=orthogonalEdgeStyle;rounded=1;html=1;"
        f"strokeColor={color};strokeWidth={width};"
        f"dashed={1 if dashed else 0};dashPattern=6 5;"
        f"endArrow={'block' if arrow else 'none'};endFill=1;endSize=6;"
    )
    cell = ET.SubElement(
        root, "mxCell", id=cid, value="", style=style, edge="1", parent="1"
    )
    geo = ET.SubElement(cell, "mxGeometry", relative="1", **{"as": "geometry"})
    ET.SubElement(
        geo,
        "mxPoint",
        x=str(points[0][0]),
        y=str(points[0][1]),
        **{"as": "sourcePoint"},
    )
    ET.SubElement(
        geo,
        "mxPoint",
        x=str(points[-1][0]),
        y=str(points[-1][1]),
        **{"as": "targetPoint"},
    )
    if len(points) > 2:
        arr = ET.SubElement(geo, "Array", **{"as": "points"})
        for x, y in points[1:-1]:
            ET.SubElement(arr, "mxPoint", x=str(x), y=str(y))
    return cell


def text_style(size=22, color=INK, bold=False, align="center"):
    return (
        "text;html=1;strokeColor=none;fillColor=none;whiteSpace=wrap;"
        f"fontSize={size};fontStyle={1 if bold else 0};fontColor={color};"
        f"fontFamily={FONT};align={align};verticalAlign=middle;"
    )


def box_style(fill="#FFFFFF", stroke=LINE, size=22, bold=False, radius=10):
    return (
        f"rounded=1;arcSize={radius};whiteSpace=wrap;html=1;fillColor={fill};"
        f"strokeColor={stroke};strokeWidth=1.5;fontSize={size};"
        f"fontStyle={1 if bold else 0};fontColor={INK};fontFamily={FONT};"
        "align=center;verticalAlign=middle;spacingLeft=10;spacingRight=10;"
    )


def section_label(root, cid, value, y):
    vertex(
        root,
        f"{cid}_bar",
        "",
        28,
        y + 5,
        5,
        22,
        "rounded=1;arcSize=50;html=1;fillColor=#314E73;strokeColor=none;",
    )
    vertex(
        root,
        cid,
        value,
        42,
        y,
        100,
        32,
        text_style(20, NAVY_DARK, True, "left"),
    )


def make_graph():
    mxfile = ET.Element(
        "mxfile",
        host="app.diagrams.net",
        agent="mathmodel",
        version="24.7.17",
        pages="1",
    )
    diagram = ET.SubElement(
        mxfile, "diagram", id="research-framework", name="统一建模框架"
    )
    model = ET.SubElement(
        diagram,
        "mxGraphModel",
        dx="1400",
        dy="405",
        grid="0",
        gridSize="10",
        guides="1",
        tooltips="1",
        connect="1",
        arrows="1",
        fold="1",
        page="1",
        pageScale="1",
        pageWidth="1400",
        pageHeight="405",
        math="0",
        shadow="0",
    )
    root = ET.SubElement(model, "root")
    ET.SubElement(root, "mxCell", id="0")
    ET.SubElement(root, "mxCell", id="1", parent="0")

    input_centers = [290, 600, 910, 1220]
    for i, x in enumerate(input_centers, 1):
        edge(
            root,
            f"input_to_bus_{i}",
            [(x, 58), (x, 73), (755, 73), (755, 86)],
            color=LINE,
            width=1.4,
            dashed=True,
            arrow=False,
        )
    edge(root, "input_bus_to_core", [(755, 73), (755, 86)], color=NAVY, width=2.3)

    question_centers = [290, 600, 910, 1220]
    edge(
        root,
        "core_to_q1",
        [(755, 155), (755, 178), (290, 178), (290, 194)],
        color=NAVY,
        width=1.9,
    )
    for i, x in enumerate(question_centers[1:], 2):
        edge(
            root,
            f"core_to_q{i}",
            [(755, 178), (x, 178), (x, 194)],
            color=NAVY,
            width=1.9,
        )
    edge(root, "q2_to_q3_state", [(741, 246), (769, 246)], color=GOLD, width=3.2)
    edge(root, "questions_to_validation", [(755, 304), (755, 330)], color=NAVY, width=2.1)

    section_label(root, "label_input", "输入", 16)
    section_label(root, "label_core", "共享内核", 103)
    section_label(root, "label_questions", "四问", 226)
    section_label(root, "label_validation", "验证", 350)

    inputs = [
        "环境　<font color=\"#697482\">T∞(t), C∞(t)</font>",
        "几何　<font color=\"#697482\">L, R(t)</font>",
        "物性　<font color=\"#697482\">b(C), k(C), D(C,T)</font>",
        "判据　<font color=\"#697482\">maxᵣ C &lt; 0.15</font>",
    ]
    for i, value in enumerate(inputs):
        vertex(
            root,
            f"input_{i+1}",
            value,
            150 + i * 310,
            11,
            280,
            45,
            box_style(PALE, LINE, 20, False, 14),
        )

    vertex(
        root,
        "core_group",
        "",
        150,
        86,
        1210,
        69,
        "rounded=1;arcSize=8;html=1;fillColor=#EEF3F8;strokeColor=none;",
    )
    core_blocks = [
        ("径向热湿守恒", NAVY, "#FFFFFF"),
        ("对称中心 · Robin 边界", "#FFFFFF", INK),
        ("真实表面 · 环形有限体积", "#FFFFFF", INK),
        ("BDF · 稀疏 Jacobian", "#FFFFFF", INK),
    ]
    widths = [250, 280, 350, 250]
    xs = [170, 440, 740, 1110]
    for i, ((value, fill, color), x, width) in enumerate(zip(core_blocks, xs, widths)):
        style = box_style(fill, NAVY if i == 0 else LINE, 20, i == 0, 8)
        style = style.replace(f"fontColor={INK}", f"fontColor={color}")
        vertex(root, f"core_{i+1}", value, x, 99, width, 43, style)

    questions = [
        (
            "<b>Q1　解耦基线</b><br>"
            "<font color=\"#697482\">固定物性 · 30 min 场</font>"
        ),
        (
            "<b>Q2　双向耦合</b><br>"
            "<font color=\"#697482\">变物性 · 3 h 场</font>"
        ),
        (
            "<b>Q3　全域事件</b><br>"
            "<font color=\"#697482\">max C=0.15 · </font>"
            "<font color=\"#B88A3B\"><b>57.4740 h</b></font>"
        ),
        (
            "<b>Q4　收缩移动域</b><br>"
            "<font color=\"#697482\">ξ=r/R(t) · </font>"
            "<font color=\"#B88A3B\"><b>51.0920 h</b></font>"
        ),
    ]
    for i, value in enumerate(questions):
        x = 150 + i * 310
        vertex(
            root,
            f"question_{i+1}",
            value,
            x,
            195,
            280,
            94,
            box_style("#FFFFFF", NAVY, 21, False, 8),
        )
        vertex(
            root,
            f"question_rule_{i+1}",
            "",
            x,
            195,
            280,
            4,
            "rounded=0;html=1;fillColor=#314E73;strokeColor=none;",
        )

    state_notes = ["原始初态", "原始初态", "延续 Q2 状态", "原始初态"]
    for i, note in enumerate(state_notes):
        color = GOLD if i == 2 else MID
        vertex(
            root,
            f"state_note_{i+1}",
            note,
            150 + i * 310,
            293,
            280,
            20,
            text_style(18, color, i == 2),
        )

    vertex(
        root,
        "validation_group",
        "",
        150,
        331,
        1210,
        69,
        "rounded=1;arcSize=8;html=1;fillColor=#F6F7F9;strokeColor=none;",
    )
    validation = [
        ("数值一致性", "解析解 · 独立 FEM · 收敛 · 守恒"),
        ("输入与参数", "留出检验 · ±10% 灵敏度"),
        ("结构假设", "长期边界 · 端面情景 · 2×2 析因"),
    ]
    for i, (head, body) in enumerate(validation):
        x = 175 + i * 395
        vertex(
            root,
            f"validation_head_{i+1}",
            head,
            x,
            341,
            360,
            23,
            text_style(20, NAVY_DARK, True),
        )
        vertex(
            root,
            f"validation_body_{i+1}",
            body,
            x,
            365,
            360,
            22,
            text_style(18, MID),
        )
        if i < 2:
            vertex(
                root,
                f"validation_sep_{i+1}",
                "",
                x + 382,
                342,
                1,
                42,
                "rounded=0;html=1;fillColor=#D3D8DE;strokeColor=none;",
            )

    return mxfile


def main():
    tree = ET.ElementTree(make_graph())
    ET.indent(tree, space="  ")
    tree.write(OUT, encoding="utf-8", xml_declaration=False)
    metadata = {
        "title": "四问统一建模与验证框架",
        "style": "restrained conference figure",
        "canvas": [1400, 405],
        "palette": {
            "primary": NAVY,
            "event_accent": GOLD,
            "neutral": [INK, MID, LINE, PALE],
        },
        "layers": ["输入", "共享内核", "四问", "验证"],
        "state_semantics": {
            "independent": ["Q1", "Q2", "Q4"],
            "continued": "Q2 → Q3",
        },
        "key_outputs": {"Q3_h": 57.4740, "Q4_h": 51.0920},
        "editable": True,
        "embedded_raster": False,
    }
    META.write_text(
        json.dumps(metadata, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    print(json.dumps({"drawio": str(OUT), "metadata": str(META)}, ensure_ascii=False))


if __name__ == "__main__":
    main()
