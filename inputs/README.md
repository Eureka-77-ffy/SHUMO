# 输入文件说明

本目录用于队内交接，保存模型复算所需的原始材料副本。原件内容未修改。

## `official/`

- `A题.pdf`：赛题原文。
- `附件1.xlsx`：环境温度与空气含湿量观测。
- `附件2.xlsx`：半径观测。
- `result1.xlsx`—`result4.xlsx`：赛题提供的空白结果模板。它们不是本队答案，不要填写或覆盖。

正式结果位于项目根目录下的 `results/final/`，待提交副本位于 `results/submission/`。

## `reference/`

- `research_report.md`：前期调研材料，仅作为模型依据整理的参考，不是论文正文。

## 换电脑后重跑模型

当前已经验收的 `config/model_config.json` 仍保留原计算机绝对路径，以维持既有缓存与审计记录的可追溯性。换电脑后若只改论文，不必改它。

若需要重新计算，请在项目根目录运行，并将 `config/model_config.json` 中的两个路径改为：

```json
"source_directory": "inputs/official",
"research_report": "inputs/reference/research_report.md"
```

可参照 `config/model_config_portable.example.json`。更改配置后应重新执行阶段0审计并重算相应问题；不要手工改缓存中的哈希或来源记录来绕过检查。

## 原始副本 SHA-256

```text
052d8014bff5727c019b72e44fdffaf5c145ce04050dd938baaf3527db331736  A题.pdf
7ef32870abeef420b89560b2530ff60dfe4255917805151d89988d0311af9dd7  附件1.xlsx
5563acbfa4b4afb10cc6c03e2207e5369bf39da27576672aff14cc5c32e704af  附件2.xlsx
23b261b295c1b787d000eebbca6521c37075107b6fcf78724f8d395ce1798ff4  result1.xlsx
23b261b295c1b787d000eebbca6521c37075107b6fcf78724f8d395ce1798ff4  result2.xlsx
07e4793d620a7f899804c0298d49a16a197960440fd47f8bb780c57ec27e2859  result3.xlsx
86e9300ffa3d30c43de895ea6723da943e85b8740b137bcae5af7107f076eeac  result4.xlsx
1b9bc78c09a31a03d93e29d914c050ecf19e4f245db5c6bd940777561de135d8  reference/research_report.md
```
