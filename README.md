# 2026 国赛 A 题：药材烘干

队友接手本目录时，请先阅读 [Workspace 队内交接说明](TEAM_HANDOFF.md)。原题、附件和空白结果模板的交接副本位于 [inputs/](inputs/README.md)；换电脑后的路径迁移、论文编译、数值依赖、正式文件与提交检查均已在交接说明中列明。

当前已完成M0—M4数值总验收、11张论文图和整篇论文的终稿修改。Q1非线性独立抽检、端面条件审查、输入留出检验及result1—result4逐单元格核验均已通过。Q3/Q4主情景临界时间分别为57.4740/51.0920 h，首个严格达标整数秒分别为206907/183932 s。论文入口为`document.tex`，当前编译得到58页A4 PDF；结论、AI声明和参考文献位于前30页，附录从第30页下半页开始，完整核心程序位于第33--58页。模型章节已补齐移动域干固体密度约去条件，灵敏度章节新增平衡边界宽情景，第三问以约0.1 s作为保守数值尺度。终检记录见[论文终检报告](reports/paper_qa.md)，评委视角修改记录见[论文评审](review.md)，AI使用详情已按本次实际使用记录编译为PDF。

最终验证结论见 [四份结果工作簿最终验证报告](reports/final_result_validation.md)、[机器可读验证记录](reports/final_result_validation.json)、[输入留出检验](reports/input_holdout_validation.json) 和 [验证包回读审计](reports/final_validation_audit.json)。共逐项核对8,882,410个结果数值格与23,165个Q4域外空白格；原始同名工作簿只作版式模板，未当作正确答案。此前端面与误差证据仍见 [端面假设审查](reports/endface_assessment.md)、[统一误差预算](reports/error_budget.md) 和 [270条误差证据CSV](tables/unified_error_budget.csv)。

最新结果见 [Q3验收、误差与机制报告](reports/q3_acceptance.md)、[完整result2.xlsx](results/final/result2.xlsx)、[result3.xlsx](results/final/result3.xlsx) 和 [论文表5](tables/q3_tables.md)。前三小时分析与表3—4仍见 [Q2阶段报告](reports/q2_acceptance.md)。Q1结果见 [问题一报告](reports/q1_acceptance.md)、[result1.xlsx](results/final/result1.xlsx) 和 [论文表1—2](tables/q1_tables.md)。模型依据见 [阶段0验收](reports/stage0_acceptance.md)、[物理机制与尺度](reports/physics_and_scales.md) 和 [模型契约](reports/model_contract.md)。完整安排见 [执行计划](执行计划.md)。

Q4最新结果见 [收缩与因素分离验收报告](reports/q4_acceptance.md)、[result4.xlsx](results/final/result4.xlsx)、[论文表](tables/q4_tables.md) 和 [四组条件对照](tables/q4_factorial.csv)。Q4原始初态复算逐元素一致；半径无需超出附件72 h范围。固定半径与附录四物性的情景在129.8484 h达标，该结果用于给定物性与外生半径轨迹下的条件比较。

## 复现阶段0

项目已创建`.venv`独立环境，当前Python为3.9.6；实际版本记录于`reports/environment.json`，依赖固定在`requirements.txt`。重新建立环境时可执行：

```sh
uv venv --python python3 .venv
uv pip install --python .venv/bin/python -r requirements.txt
```

在本目录运行：

```sh
.venv/bin/python scripts/stage0_audit.py
```

脚本只读原始文件，重新生成`reports/`内的JSON/CSV数值审计，不修改人工撰写的模型说明或原始Excel。输入路径在`config/model_config.json`中；移动原件后须更新路径并重新记录校验值。

检查结果包括12组物性、9项代数关系和4个错误时间线配置。原始输入与配置/脚本SHA-256记录在`reports/input_manifest.json`，便于后续绑定计算缓存。

## 当前模型口径

Q1为解耦的热传导与水分扩散；Q2/Q3用附录三；Q4用附录四与指定R(t)。Q2和Q4均从原始初态重算。Q4的经验密度只用于有效热学系数，输运按守恒干质量定义；严格总质量密度解释与收缩数据兼容性不足的证明已写入模型契约。

问题一实际数值参数单列于`config/q1_final.json`：表面加密1280格、rtol=1e-10、atol=1e-12、最大步长1 s。Q23配置见`config/q23_final.json`：1280格、rtol=1e-10、温度/水分atol分别1e-11和1e-13，前三小时最大步长2 s，之后30 s。Q4配置见`config/q4_final.json`：1280格、相同温湿容差、最大步长30 s、分段线性半径且禁止未声明外推。各正式配置均另以2560格复核。原始`model_config.json`保留阶段0初拟值以维持输入版本可追溯。所有代表状态的R²/D仅为局部数量级，不是最终干燥时长。

## 复现问题一

已有有效缓存时，执行验收、工作簿导出和机制诊断：

```sh
.venv/bin/python -W error validation/kernel_checks.py
.venv/bin/python -W error scripts/export_question1.py
.venv/bin/python -W error scripts/q1_mechanisms.py
```

从头重算网格、时间和敏感性证据时，先依次运行：

```sh
.venv/bin/python -W error validation/q1_convergence.py
.venv/bin/python -W error validation/q1_time_and_sensitivity.py
```

然后执行上述三条验收/导出命令。脚本生成项目内缓存和结果，不覆盖用户提供的原始结果模板。若修改Q1使用的顶层`src/*.py`、主配置或Q1驱动，导出器会拒绝失配缓存，需重算相关证据。新增的`src/coupled/`单独纳入Q2版本，不使未改动的Q1失效。直接运行`question1.py`的默认参数只是单个试算，并非最终验收配置。

验收包括独立Bessel热场、常D水分退化基准、网格加密至2560格、容差和实际时间步比较、Radau与Kirchhoff对照、累计平衡及8个单因素敏感性情景。Excel共75,600个结果值已逐项回读。M4已补齐Q1原题非线性水分的独立弱形式抽检，2560单元细步解与主解的水分稠密剖面最大差约1.32×10⁻⁶，见`reports/q1_nonlinear_independent.json`。四位小数仍是显示要求，尚未证明所有末位稳定或完成实验有效性校准。

## 复现问题二前三小时

在Q1有效缓存存在的基础上，运行：

```sh
.venv/bin/python -W error validation/q2_kernel_checks.py
.venv/bin/python -W error validation/q2_convergence.py
.venv/bin/python -W error validation/q2_sensitivity.py
.venv/bin/python -W error validation/q2_independent_verification.py
.venv/bin/python -W error scripts/export_question2.py
```

独立校验脚本会复用自身源文件/配置/内容哈希全部匹配的有限元缓存；不匹配则重算。其他脚本重新计算对应场景。已有匹配证据时，只需最后一条导出命令即可重新审计、制表与导出。

此阶段的工作簿仅覆盖1—10800 s，故明确命名`result2_first3h.xlsx`并放在`results/alternatives/`。现已另行生成完整`results/final/result2.xlsx`，不能把三小时当成干燥终点。`q2_first3h_full_precision.npz`包含原始计时的每60 s内部状态和10800 s检查点，可用于同模型续算；禁止从四位小数工作簿续算。

此处Q2阶段报告仅评价前三小时；Q3全过程和Q4各有下面的独立验证链，不能互相替代。`reports/q2_acceptance.md`明确记录一个论文表温度值跨越舍入边界的情况，不作全部末位稳定的承诺。

## 复现Q3与完整Q2工作簿

先确保上节的Q2基线、各级网格和参数情景缓存存在且版本匹配，再运行：

```sh
.venv/bin/python -W error validation/q3_convergence.py
.venv/bin/python -W error validation/q3_event_checks.py
.venv/bin/python -W error validation/q3_independent_verification.py
.venv/bin/python -W error validation/q3_sensitivity.py
.venv/bin/python -W error validation/q3_error_budget.py
.venv/bin/python -W error scripts/export_question3.py
```

比较脚本复用来源、参数和文件内容校验值全部匹配的缓存，失配则重算。独立有限元只继承它自身的Q2节点轨迹，并按需生成2560单元检查点，不继承主程序状态。新增代码在`src/events/`，独立纳入Q3版本而不改变Q1/Q2核心。

完整Q2每张表206907行，Q3有3449行，均追加精确临界时刻；已逐项回读8,762,523个结果值。完整Q2工作簿约26 MB，读写耗时明显高于单个主求解场景。`q23_full_precision.npz`是两问共同的未舍入输出；内部检查点仍在`results/cache/`。题面完整Q2输出时域有解释空间，本项目沿已冻结计划取至Q3临界时刻，前三小时版本保留为备选。

57.4740 h依赖“4 h后保持3—4 h均值”和等效含水率映射等假设。主网格渐近终点误差尺度约0.10 s，不是严格上界；全部输出末位稳定未获证明，真实边界/几何/物性误差不能由数值收敛消除。

## 复现Q4与四组因素分离

按下列顺序运行；验证与情景程序复用来源、参数、内容哈希均匹配的缓存，失配则重算。全新复算使用单独标签，不覆盖已验收的主解缓存：

```sh
.venv/bin/python -W error validation/q4_kernel_checks.py
.venv/bin/python -W error validation/q4_convergence.py
.venv/bin/python -W error validation/q4_independent_verification.py
.venv/bin/python -W error scripts/factorial_q4.py
.venv/bin/python -W error validation/q4_sensitivity.py
.venv/bin/python -W error validation/q4_factorial_verification.py
.venv/bin/python -W error validation/q4_diagnostics.py
.venv/bin/python -W error scripts/question4.py --label q4_reproduction --n 1280 --max-step 30
.venv/bin/python -W error scripts/export_question4.py
```

已有全部匹配证据和复算缓存时，最后一条即可再次验收导出。`src/moving/`单独纳入Q4版本，没有修改Q1/Q23内核。独立移动域FEM从初态重算，不读取主程序状态；额外对B/C因素情景分别做空间、时间和独立方法比较。

result4共3066个数据行、44287个数值与23165个域外空白，均已逐项回读；21个固定物理半径列之外另有真实动态表面。论文表6为每6 h及临界时刻，共9行。末行是临界等号，不是严格达标；真实表面不能用最后一个有效固定半径列代替。

Q4主网格1280→2560终点变化−0.05797 s，二阶渐近误差尺度约0.0773 s；最细独立求解终点互差约0.00067 s，但可能含空间误差抵消，不声称毫秒级物理精度。归一化水分累计残差2.48×10⁻¹⁰，修正有效热方程收支相对残差5.44×10⁻¹⁰；后者不代表完整多相总能量守恒。

端面/一维几何适用范围、Q1非线性水分独立抽检、统一误差预算与四份结果总验收已补齐，复现方法见下节。下一步是生成科研图表，再将适用条件与证据写入论文并补充论断证据表，之后完成真实文献、LaTeX正文/PDF及八项风险终审。原始附件与模板未改动；`results/final/`是所选模型的数值交付目录，不代表论文已经达到最终提交状态。

## 复现端面审查与统一误差预算

保留既有四问验收产物，在项目目录执行：

```sh
.venv/bin/python -W error validation/endface_suite.py
.venv/bin/python -W error validation/m4_scope_audit.py
```

第一条核验/生成29个端面场景缓存（17组同网格配对）、Q1非线性独立抽检、二维内核/有限圆柱模态对照、重建与几何分辨率诊断，以及270行JSON/CSV误差台账。匹配缓存会复用；缺失则计算。来源或内容不匹配会明确拒绝，不能重标已有缓存绕过验证。全量二维从头计算比一维主解耗时明显更长。

第二条回查JSON/CSV逐行一致性、新报告来源哈希、原始8个文件和此前54项验收产物未变，生成`reports/m4_validation_audit.json`。人工报告不会由复现脚本覆盖；计算结果或口径变化后，应同步审阅报告而非只重新生成摘要。

端面审查属于几何情景分析，不是新的独立物性/边界校准。当前保留一维主模型及其限定；若把暴露端面升级为正式主模型，需要另行对二维绝对解达到所需精度并重新导出，不能直接用此次粗网格对照替换。

## 复现最终结果总验收

保持既有缓存、工作簿和M4证据不变，在项目目录运行：

```sh
.venv/bin/python -W error validation/input_holdout.py
.venv/bin/python -W error validation/final_result_validation.py
.venv/bin/python -W error validation/final_validation_audit.py
```

第二条会重新读取约890万个Excel数值格，因此明显慢于纯JSON审计。验收同时检查量纲与边界、物理范围、全域事件、守恒残差、独立FEM/解析基线、网格和时间误差证据、输入留出、单因素敏感性、边界/插值情景及极端反例。附件没有内部温湿观测，所以留出检验只验证输入插值；报告明确不将其写成PDE输出的实验交叉验证。进入制图后应读取`results/final/*_full_precision.npz`和`tables/`，不要从四位小数Excel反推曲线。
