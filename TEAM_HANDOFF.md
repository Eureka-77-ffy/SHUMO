# Workspace 队内交接说明

更新日期：2026-09-12  
项目：2026 年全国大学生数学建模竞赛 A 题——药材烘干

这份文件供接手论文的队友使用。建议第一次打开项目时，先读第 1—4 节；准备改模型或重新计算时，再读第 8—10 节。

## 1. 当前可用版本

论文主入口是 `document.tex`，当前正式 PDF 是：

- `document.pdf`：当前正式稿，58 页 A4；
- `document_final_reviewed.pdf`：与 `document.pdf` 完全相同的留档副本；
- `document_table_revision.pdf`：较早的表格调整版，已经过期，不要继续使用或提交；
- `AI工具使用详情.pdf`：按竞赛要求单独准备的 AI 使用详情文件。

当前正文、AI 使用声明和参考文献在前 30 页内结束，附录从第 30 页下半页开始，核心代码清单位于第 33—58 页。论文现有 11 张图、14 张表，编译日志中没有溢出、未定义引用或文献警告。

正式结果文件有两套相同副本：

- `results/final/result1.xlsx`—`result4.xlsx`：验收通过的正式结果；
- `results/submission/result1.xlsx`—`result4.xlsx`：提交用副本。

不要把 `inputs/official/result1.xlsx`—`result4.xlsx` 当成答案。那四个文件是赛题提供的空白模板，只用于规定工作表、表头和格式。

## 2. 已冻结的核心结论

除非模型、边界条件或参数发生实质变化，下列结果不要只因排版或措辞调整而改动。

### 第三问

- 全域判据：对整个径向区间搜索最大干基含水率，而不是只检查中心、表面或截面平均值；
- 连续临界时刻：`206906.4911824653 s = 57.47402532846259 h`；
- 论文题定结果：`57.4740 h`；
- 首个严格达标整数秒：`206907 s`；
- Q3 从 Q2 的未舍入全精度状态继续，不能从四位小数 Excel 结果续算；
- 1280→2560 网格变化约 0.076 s，文中采用约 0.1 s 作为保守数值尺度；
- 独立有限元与主解约 0.0021 s 的差异只说明数值一致性，不代表实际过程具有毫秒级精度。

第三问是全文最需要保持口径一致的部分。对应正文、图表和验证记录分别位于：

- `texfile/5MakeModel.tex`；
- `tables/q3_tables.tex`、`tables/q3_cross_validation.tex`；
- `figures/fig04_q3_endpoint.pdf`；
- `reports/q3_acceptance.md`；
- `reports/q3_cross_validation.md`；
- `reports/q3_event_checks.json`；
- `reports/q3_independent_verification.json`；
- `reports/q3_convergence.json`。

### 第四问

- 连续临界时刻：`183931.2887197449 s = 51.09202464437358 h`；
- 论文题定结果：`51.0920 h`；
- 首个严格达标整数秒：`183932 s`；
- Q4 按原始初态独立重算，不继承 Q3 末态；
- 半径采用附件 2 的分段线性观测轨迹，主结果不需要超出 72 h 数据范围；
- 水分方程按守恒干固体质量定义，干固体密度因子在归一化方程中约去；经验密度只作为有效热学系数使用；
- 四组因素对照是给定边界和参数下的条件比较，不应写成普遍因果结论，也不要把基于 A 组的差值称为标准析因“主效应”。

### 边界口径

- 附件 1 只观测到 4 h；主情景在 4 h 后保持 3—4 h 观测均值；
- 空气含湿量通过 `C_eq = beta * Y` 映射为等效平衡含水率，这是模型闭合假设，不是题面给出的实验定律；
- 文中已用更宽边界情景说明这一闭合对长期终点的影响；
- 主模型是径向一维、端面绝热且不透湿，端面暴露只作为结构情景审查，不能表述为实验误差上界。

## 3. 项目目录地图

| 路径 | 内容 | 修改建议 |
|---|---|---|
| `document.tex` | LaTeX 主入口、标题、宏包、各章节装配顺序 | 保留为唯一主入口 |
| `texfile/` | 摘要、问题分析、模型、误差、灵敏度、评价、AI 声明、参考文献、附录 | 正文修改主要发生在这里 |
| `tables/` | 论文表格的 `.tex` 源、数值 CSV 和 Markdown 对照 | 改数值前先查来源 |
| `figures/` | 正式 PDF/SVG/PNG/TIFF 图片及可编辑源文件 | 论文实际引用 PDF |
| `figures/scripts/` | 正式图和附图的 Python 生成程序 | 数据图修改从这里开始 |
| `figures/source_data/` | 附图使用的显式数据表 | 由脚本生成，不要手填 |
| `figures/archive/` | 旧版图片备份 | 仅供回退，不参与编译 |
| `figures/layout_previews/` | 早期版式样稿 | 不参与正式论文 |
| `inputs/official/` | 赛题 PDF、两个附件、四个空白结果模板的交接副本 | 原始输入，只读 |
| `inputs/reference/` | 前期调研报告副本 | 仅作参考 |
| `config/` | 模型总配置和 Q1/Q23/Q4 正式数值配置 | 修改会影响缓存指纹 |
| `src/` | 有限体积、边界、物性、耦合、事件和移动域核心代码 | 属于模型内核 |
| `scripts/` | 四问求解、导出和机制分析程序 | 裸跑默认参数不是正式配置 |
| `validation/` | 收敛、守恒、独立有限元、端面、留出和总验收 | 改模型后按依赖重跑 |
| `results/final/` | 正式全精度 NPZ 与四份正式工作簿 | 论文数据的权威来源 |
| `results/submission/` | 与正式结果哈希一致的提交副本 | 最终上传使用 |
| `results/cache/` | 主解、灵敏度、网格加密、独立解等缓存，约 954 MB | 不要随意删除或改名 |
| `results/alternatives/` | Q2 前三小时等备选产物 | 不是正式提交结果 |
| `reports/` | 各阶段验收、误差、图表和论文终检记录 | 查依据时先看这里 |
| `.mathmodel/paper/config.json` | 论文模板入口配置 | 已设为 CUMCM/A/中文 |
| `.venv/` | 当前机器的 Python 虚拟环境，约 196 MB | 换电脑后建议重建 |
| `requirements.txt` | 固定的 Python 依赖版本 | 新环境按此安装 |
| `book.bib` | 真实参考文献 BibTeX | 增删引用时同步编译 BibTeX |
| `cumcmthesis.cls` | 当前国赛 LaTeX 模板 | 非必要不要修改 |
| `review.md` | 评委视角问题和已完成修改记录 | 用于了解稿件薄弱处 |
| `README.md` | 数值复现命令和模型口径 | 需要重算时详细阅读 |

整个文件夹当前约 1.3 GB，主要空间来自 `results/cache/`、`.venv/` 和多格式图片。

## 4. 数据、模型和论文的关系

项目的主要数据流如下：

```text
inputs/official/
        │
        ├── config/model_config.json + src/ + scripts/
        │                         │
        │                         ├── Q1 独立计算
        │                         ├── Q2 原始初态计算 ──全精度状态──> Q3 事件续算
        │                         └── Q4 原始初态独立计算
        │
        └── results/cache/ ──> results/final/*.npz
                                    │
                                    ├── results/final/result*.xlsx
                                    ├── tables/*.csv、*.tex
                                    └── figures/*.pdf
                                               │
texfile/*.tex + tables/*.tex + figures/*.pdf + book.bib
                                               │
                                         document.pdf
```

需要牢记三个关系：

1. Q1、Q2、Q4 都从原始初态开始；只有 Q2→Q3 是同一条未舍入轨迹的延续。
2. Excel 文件是提交格式，不是后续求解初值。计算和作图优先读取 `results/final/*_full_precision.npz`。
3. `results/cache/` 的元数据绑定了输入、配置和源代码哈希。来源不一致时程序会拒绝复用，这是有意设计的保护机制。

## 5. 论文正文怎么改

`document.tex` 只负责装配。正文按以下顺序引入：

| 文件 | 内容 |
|---|---|
| `texfile/1abstract.tex` | 摘要和关键词 |
| `texfile/2ProblemRestatement.tex` | 问题重述和研究框架图 |
| `texfile/3ProblemAnalysis.tex` | 总体分析和输入尺度 |
| `texfile/4AssumptionAndSign.tex` | 假设、符号和几何图 |
| `texfile/5MakeModel.tex` | 四问模型、求解、表格与主要结果图 |
| `texfile/6ErrorAnalysis.tex` | 数值误差与验证 |
| `texfile/7SensitivityAnalysis.tex` | 参数、边界和结构敏感性 |
| `texfile/7ModelEvaluation.tex` | 优缺点、推广和结论 |
| `texfile/8AIUsageStatement.tex` | 正文内的 AI 使用声明 |
| `texfile/8Reference.tex` | 参考文献样式与 `book.bib` 引入 |
| `texfile/9Appendix.tex` | 附图和核心程序清单 |

常见修改位置：

- 改摘要：`texfile/1abstract.tex`；
- 改第三问或第四问论证：`texfile/5MakeModel.tex`；
- 改误差、交叉验证表述：`texfile/6ErrorAnalysis.tex`；
- 改边界情景和灵敏度：`texfile/7SensitivityAnalysis.tex`；
- 改结论与模型评价：`texfile/7ModelEvaluation.tex`；
- 改参考文献条目：`book.bib`；
- 改 AI 声明：正文改 `texfile/8AIUsageStatement.tex`，详情文件改 `AI工具使用详情.tex`。

改动正文时请保持这些写作口径：

- 用“本文”“本模型”“数值结果表明”，避免宣传式语言和空泛评价；
- 先给判据或依据，再给数字结论；
- 不把数值误差写成物理误差；
- 不把情景差写成置信区间或普遍因果效应；
- 图题只说明图中对象，解释放在正文；
- 避免连续使用“首先、其次、此外、综上”等模板化连接词；
- 不新增无法由题面、公式、计算或真实文献支持的结论。

## 6. 表格怎么改

论文中的正式表格由 `tables/*.tex` 引入。当前主要对应关系是：

- Q1：`tables/q1_tables.tex`；
- Q2：`tables/q2_tables.tex`，温度与含水率已经合并为一个双分面表；
- Q3：`tables/q3_tables.tex`；
- Q3 交叉验证：`tables/q3_cross_validation.tex`；
- Q4：`tables/q4_tables.tex`。

同目录下的 CSV 是数值对照和制图数据。只改列名、表题、宽度或排版时，可以直接改 `.tex`；如果改了数值，必须先确认它来自哪个正式 NPZ、报告或导出器，再同步正文、表格、图片和验收记录。

不要在 LaTeX 表格中手工修饰某个末位来“对齐”结论。当前四位小数是输出格式，不表示每个末位都有相同物理精度。

## 7. 图片怎么改

论文实际插入的是 `figures/*.pdf`。PNG 用于预览，SVG/draw.io/TikZ 用于继续编辑，TIFF 用于高分辨率留档。

| 图片 | 主题 | 主要源文件 |
|---|---|---|
| `fig00a_research_roadmap` | 研究框架图 | `figures/fig00a_research_roadmap.drawio`、`figures/scripts/make_research_framework.py` |
| `fig00b_model_geometry` | 几何、边界和移动坐标 | `figures/fig00b_model_geometry.tex` |
| `fig01_inputs_scales` | 输入与时间尺度 | `figures/scripts/make_formal_enriched_figures.py` |
| `fig02_q1_fields` | Q1 温湿场 | 同上 |
| `fig03_q2_coupling` | Q2 耦合机制 | 同上 |
| `fig04_q3_endpoint` | Q3 全域终点 | 同上 |
| `fig05_q4_shrinkage` | Q4 收缩效应 | 同上 |
| `fig06_numerical_validation` | 数值验证 | 同上 |
| `fig07_sensitivity_risk` | 灵敏度与风险量级 | 同上 |
| `figS01_input_holdout` | 输入留出检验 | `figures/scripts/make_formal_enriched_supplements.py` |
| `figS02_endface_assessment` | 端面情景审查 | 同上 |

当前统一配色为白底、深蓝灰主色，金色只强调临界事件。不要重新加入多色渐变、阴影、大段图内文字或无解释作用的装饰。

在原电脑路径仍有效时，数据图可用以下命令重绘：

```sh
.venv/bin/python figures/scripts/make_formal_enriched_figures.py
.venv/bin/python figures/scripts/make_formal_enriched_supplements.py
.venv/bin/python figures/scripts/qa_figures.py
```

框架图的 Python 脚本只生成可编辑 `.drawio` 和元数据，PDF/SVG/PNG/TIFF 需用 diagrams.net/draw.io 导出。几何图用 XeLaTeX 单独编译，再将单页 PDF 保存在同名位置。

如果只调整图片样式而不改数值，不能顺手覆盖工作簿。重绘后至少检查：图中文字、色标、图例、Q4 域外白色遮罩、PDF 单页性、论文实际缩放后的可读性。

## 8. 编译论文

### LaTeX 环境

建议安装完整 TeX Live 或 MacTeX，并使用 XeLaTeX。项目依赖中文排版和 `cumcmthesis.cls`。主论文的完整编译顺序是：

```sh
xelatex -interaction=nonstopmode -halt-on-error document.tex
bibtex document
xelatex -interaction=nonstopmode -halt-on-error document.tex
xelatex -interaction=nonstopmode -halt-on-error document.tex
```

修改参考文献后必须执行 BibTeX；只改普通正文时也建议按完整链编译，避免交叉引用沿用旧缓存。

AI 使用详情单独编译：

```sh
xelatex -interaction=nonstopmode -halt-on-error "AI工具使用详情.tex"
xelatex -interaction=nonstopmode -halt-on-error "AI工具使用详情.tex"
```

编译后检查：

```sh
rg -n "Warning|Overfull|Underfull|undefined" document.log
```

理想输出为空。还要人工浏览摘要页、框架图页、Q2 表格页、Q3/Q4 结果页、参考文献末页、正文与附录分界页。不要只看编译成功。

`clean.bat` 是模板自带的 Windows 清理脚本；清理辅助文件不会恢复正文，也不能替代重新编译。

## 9. Python 环境与原始路径

当前数值环境记录在 `reports/environment.json`，Python 为 3.9.6，固定依赖见 `requirements.txt`。`.venv/` 含当前机器路径，不保证换电脑后可用。建议在新电脑重新建立：

```sh
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
```

Windows 对应解释器通常是：

```text
.venv\Scripts\python.exe
```

所有 Python 命令都应从 Workspace 根目录运行。

当前 `config/model_config.json` 仍保存原计算机的绝对路径。这是为了让现有缓存、审计记录和最终结果保持原验证状态。换电脑后：

- 只改论文、表格排版或直接编辑现有图片：不需要改配置；
- 需要重绘会调用 `verify_sources()` 的数据图：必须先处理输入路径；
- 需要重跑模型：将 `source_directory` 改为 `inputs/official`，将 `research_report` 改为 `inputs/reference/research_report.md`，可参照 `config/model_config_portable.example.json`。

路径改变后，即使输入文件内容相同，配置哈希也会变化，旧缓存按设计会失配。正确做法是重新运行阶段 0 审计并重算相应链条。不要直接编辑 JSON 中的哈希，不要关闭来源校验，也不要把旧缓存重新贴标签冒充新结果。

另外，`config/model_config.json` 的 `status` 和部分早期数值字段保留了阶段 0 记录，不代表正式计算参数。正式参数以这三个文件为准：

- `config/q1_final.json`；
- `config/q23_final.json`；
- `config/q4_final.json`。

`scripts/question1.py`—`question4.py` 的命令行默认值只是试算参数，也不是最终验收参数。

## 10. 什么时候需要重算

### 不需要重算

以下修改一般只需重新编译论文：

- 改错别字、句式、章节衔接；
- 调整图表位置、标题和字号；
- 修改不涉及数值的摘要、评价、假设表述；
- 更新 AI 使用说明；
- 调整参考文献格式。

### 需要局部重绘或重新导出

- 改图的颜色、布局、标注：重绘对应图片并执行图片 QA；
- 改 Excel 格式但不改数值：使用对应导出器重新生成，再核对哈希和单元格；
- 改表格抽样时刻：从全精度 NPZ 重新抽取，不能从 PDF 或四位小数 Excel 抄回。

### 必须重算

以下修改会改变科学结果：

- 改物性公式、边界系数、平衡含水率映射；
- 改 4 h 后环境延拓方式；
- 改网格、时间步、求解器容差或全域事件定义；
- 改 Q4 半径插值、材料坐标或守恒方程；
- 改核心 `src/`、问题驱动脚本或正式配置。

重算依赖顺序：

```text
阶段0
├── Q1 → Q1 导出与验证
├── Q2 → Q2 验证 → Q3 → Q3 导出与验证
└── Q4 → Q4 因素对照、导出与验证
                       ↓
            最终四工作簿总验收
                       ↓
                 表格、图片、论文
```

完整命令已经按问题写在 `README.md`。不要在不了解依赖时一次性删除缓存后全跑：独立有限元、2560 网格和二维端面审查耗时较长，占用空间也较大。

最末的总验收命令是：

```sh
.venv/bin/python -W error validation/input_holdout.py
.venv/bin/python -W error validation/final_result_validation.py
.venv/bin/python -W error validation/final_validation_audit.py
```

其中 `final_result_validation.py` 会回读约 890 万个 Excel 数值格，运行时间明显长于普通文本检查。

## 11. 重要报告怎么查

遇到问题时，建议按下表查证，不要只凭正文判断。

| 需要确认的内容 | 报告 |
|---|---|
| 总体模型口径 | `reports/model_contract.md` |
| 物理机制与量级 | `reports/physics_and_scales.md` |
| Q1 验收 | `reports/q1_acceptance.md` |
| Q2 验收 | `reports/q2_acceptance.md` |
| Q3 终点与误差 | `reports/q3_acceptance.md` |
| Q3 独立交叉验证 | `reports/q3_cross_validation.md` |
| Q4 验收与因素对照 | `reports/q4_acceptance.md` |
| 边界宽情景 | `reports/boundary_closure_scenarios.md` |
| 端面假设审查 | `reports/endface_assessment.md` |
| 统一误差预算 | `reports/error_budget.md` |
| 四份工作簿逐单元格验收 | `reports/final_result_validation.md` |
| 最终验证包审计 | `reports/final_validation_audit.json` |
| 图片规范检查 | `reports/figure_qa.md` |
| 论文页数、引用和实页终检 | `reports/paper_qa.md` |
| 评委视角修改记录 | `review.md` |

JSON 是机器记录，Markdown 是便于阅读的解释。若两者有冲突，应先停下来核对生成脚本和时间戳，不能只选择更有利的数字。

## 12. 当前正式文件哈希

文件传输后可用 SHA-256 判断是否损坏或被意外覆盖。

```text
a62dc57d7c78e1abf791fd280749a4d57aff2f0b40613d4492c5081b434dd01e  document.pdf
a62dc57d7c78e1abf791fd280749a4d57aff2f0b40613d4492c5081b434dd01e  document_final_reviewed.pdf
bcad0ca11d996ae9022b9eed1ed526fc29590cdfd4144c439be7f6c6b56f5ac0  AI工具使用详情.pdf
11f29d96fdbbbf982aa5589beecf50126b52a81e279e841bf69ff8d5a30765c6  results/final/result1.xlsx
ead4006426658b090d564129cd705dbed57dcfb728ef4fe3dc1ae667ea6b2687  results/final/result2.xlsx
7f1efee8dc8af9a7a0cb3a5ea503b306b88d13965cf5c8efeb96ec03b57df4aa  results/final/result3.xlsx
8086731d59ba174376121e8ec103e585766034f40a86e8ad5366faf7f1029326  results/final/result4.xlsx
```

`results/submission/` 下四个同名文件应分别与上述正式工作簿哈希一致。

macOS/Linux 可执行：

```sh
shasum -a 256 document.pdf "AI工具使用详情.pdf" results/final/result*.xlsx results/submission/result*.xlsx
```

项目已提供关键文件校验清单；传输整个文件夹后，可从 Workspace 根目录一次检查：

```sh
shasum -a 256 -c reports/handoff_integrity.sha256
```

Windows PowerShell 可用 `Get-FileHash -Algorithm SHA256`。

## 13. 提交前检查

最后提交前至少完成以下检查：

1. 打开 `document.pdf`，确认摘要第一页、页码连续、无封面和目录；
2. 确认正文、AI 声明、参考文献仍在前 30 页，附录没有把正文内容挤出限制；
3. 检查标题、正文、PDF 属性和图片中没有学校、队员姓名等身份信息；
4. `document.tex` 中的 `MCxxxxxxx` 目前因 `withoutpreface` 不显示；若改回带封面的模板，必须正确处理报名号；
5. 确认 Q3 仍为 `57.4740 h / 206907 s`，Q4 仍为 `51.0920 h / 183932 s`；
6. 确认 `results/final/` 与 `results/submission/` 四份工作簿逐一同哈希；
7. 确认 `AI工具使用详情.pdf` 与正文 AI 声明口径一致；
8. 再运行一次日志检查并逐页浏览 PDF；
9. 检查竞赛系统当年的文件命名、上传入口和大小限制。当前 `result2.xlsx` 约 27 MB，若平台对压缩包或单文件有限制，要按当年官方说明处理，不能为压缩而删减结果行。

正式提交应以根目录的 `document.pdf`、`AI工具使用详情.pdf` 和 `results/submission/` 为准。旧 PDF、预览图片、缓存和报告都不应混入上传包，除非赛事明确要求。

## 14. 队内协作建议

这个目录目前不是 Git 仓库。多人继续修改时，至少采用以下做法之一：

- 每次修改前复制一份带日期的正文源文件或整个轻量论文目录；
- 建立 Git 仓库，但不要提交 `.venv/`、LaTeX 辅助文件和无必要的巨型缓存；
- 一个人负责正文，一个人负责数值，避免同时改同一 `.tex` 或同一生成脚本；
- 每次交接注明“改了什么、是否重算、正式 PDF 哈希、四份工作簿是否变化”。

推荐每轮修改都遵循：

```text
修改源文件 → 编译/生成 → 自动检查 → 实页浏览 → 记录变化 → 再交给下一位
```

如果只是继续润色论文，最安全的起点是 `document.tex`、`texfile/`、`tables/` 和 `figures/`；不要碰 `src/`、`config/`、`results/cache/`。如果确实要改模型，先保存当前正式 PDF、四份工作簿和本说明中的哈希，再按 `README.md` 的依赖顺序重新计算。
