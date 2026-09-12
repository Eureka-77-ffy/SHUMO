# 图表插入与论证指南

图注只保留短图题；面板解释与结论放在正文。当前项目没有论文入口 `document.tex`，因此本文件给出建议插入章节和可直接复用的 LaTeX 片段，不自动改写正文。

| 论文编号与文件 | 建议章节 | 短图题 | 正文承担的论证 |
|---|---|---|---|
| 图1 `fig00a_research_roadmap` | 问题分析 | 整体技术路线 | 四问的复杂度递进、独立初态口径及统一验证链 |
| 图2 `fig00b_model_geometry` | 模型假设与符号 | 几何边界与收缩坐标 | 一维中截面解释、端面条件、侧面边界与Q4材料坐标 |
| 图3 `fig01_inputs_scales` | 数据分析与尺度分析 | 环境驱动与特征尺度 | 环境在数小时内趋稳，但水分扩散和收缩的时间尺度更长，支持时变边界与分尺度建模 |
| 图4 `fig02_q1_fields` | 问题一 | 问题一温湿场演化 | 30 min 内热影响已进入内部，而明显失水仍集中在表层 |
| 图5 `fig03_q2_coupling` | 问题二 | 问题二耦合响应 | 温度接近平台不等于干燥完成；温升和低含水对局部扩散率作用方向相反 |
| 图6 `fig04_q3_endpoint` | 问题三 | 问题三全域干燥终点 | 后期由全域最大含水率控制，57.4740 h 是等号临界点，206907 s 才是首个核验严格达标整数秒 |
| 图7 `fig05_q4_shrinkage` | 问题四 | 问题四收缩域效应 | 收缩加快几何输运，但Q3—Q4差异还包含物性变化与交互作用 |
| 图8 `fig06_numerical_validation` | 模型检验 | 数值收敛与独立验证 | 数值误差、独立实现差异和守恒残差均低于主要情景差异，但不等于物理实验验证 |
| 图9 `fig07_sensitivity_risk` | 灵敏度与模型评价 | 参数敏感性与模型风险 | 终点对D最敏感；边界和结构情景的影响远大于纯数值误差，情景范围不是置信区间 |
| 图S1 `figS01_input_holdout` | 附录：输入检验 | 输入插值留出检验 | 线性插值优于持续值基线；PCHIP未表现出稳定优势；检验对象仅为输入结点 |
| 图S2 `figS02_endface_assessment` | 附录：端面审查 | 端面假设条件对照 | 端面影响在中截面、整根平均量和临界时间上的尺度不同 |

## LaTeX 片段

```latex
\begin{figure}[htbp]
  \centering
  \includegraphics[width=\textwidth]{figures/fig00a_research_roadmap.pdf}
  \caption{整体技术路线}
  \label{fig:research-roadmap}
\end{figure}
```

其余各图只需替换文件名、图题和标签：

| 文件名 | 图题 | 标签 |
|---|---|---|
| `fig00b_model_geometry.pdf` | 几何边界与收缩坐标 | `fig:model-geometry` |
| `fig01_inputs_scales.pdf` | 环境驱动与特征尺度 | `fig:inputs-scales` |
| `fig02_q1_fields.pdf` | 问题一温湿场演化 | `fig:q1-fields` |
| `fig03_q2_coupling.pdf` | 问题二耦合响应 | `fig:q2-coupling` |
| `fig04_q3_endpoint.pdf` | 问题三全域干燥终点 | `fig:q3-endpoint` |
| `fig05_q4_shrinkage.pdf` | 问题四收缩域效应 | `fig:q4-shrinkage` |
| `fig06_numerical_validation.pdf` | 数值收敛与独立验证 | `fig:numerical-validation` |
| `fig07_sensitivity_risk.pdf` | 参数敏感性与模型风险 | `fig:sensitivity-risk` |
| `figS01_input_holdout.pdf` | 输入插值留出检验 | `fig:input-holdout` |
| `figS02_endface_assessment.pdf` | 端面假设条件对照 | `fig:endface-assessment` |

## 使用边界

- 热图未平滑，也未生成虚构中间状态；为减小文件体积只做声明过的显示采样。
- Q4物理坐标图的白色部分表示该位置已经位于药材之外。
- 数值收敛线给出相邻离散之间的观测差，不是严格误差上界。
- 独立有限元对照验证计算实现，不替代内部温湿实验数据。
- 图S1只验证观测结点间插值，不能用于证明4 h后的环境平台延拓。
- 图S2是端面与侧面采用相同系数的条件对照，不是未知真实端面条件的误差上界。
