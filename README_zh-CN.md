# 研究分析代码

本仓库用于整理本研究数据分析和图件绘制所使用的 R 脚本。

论文题目：*Plastisphere formation reshapes the ecological toxicity of microplastics in Corbicula fluminea.*

## 已整理的脚本

| 图件 | 分析内容 | R 脚本 |
|---|---|---|
| Fig. 1 | 河蚬生理与生化指标 | `scripts/main_figures/fig1_physiology_biochemistry.R` |
| Fig. 5 | 肠道 16S α 多样性 | `scripts/main_figures/fig5_gut_16s_alpha_diversity.R` |

其他分析和绘图脚本将在整理完成后陆续补充。

## 输入数据

输入数据不随代码仓库上传。运行脚本前，请将对应数据文件放入 `data/input/`。

英文 README 中列出了后续重新整理时拟采用的标准化数据名称。目前 R 脚本仍引用原始工作文件名；待输入数据重新整理后，再统一调整脚本中的文件路径。实际保存数据时应保留相应格式扩展名，例如 `.csv` 或 `.xlsx`。

## 运行方法

1. 安装 R 和 RStudio。
2. 使用 RStudio 打开 `Plastisphere_Microplastic_Toxicity.Rproj`。
3. 将对应输入数据放入 `data/input/`。
4. 运行 `scripts/install_packages.R` 安装所需 R 包。
5. 从项目根目录运行相应脚本。生成的图件及其他结果保存在 `results/`。

## 仓库结构

```text
scripts/
  install_packages.R
  main_figures/
    fig1_physiology_biochemistry.R
    fig5_gut_16s_alpha_diversity.R
data/
  input/
results/
```

GitHub: https://github.com/Xiaokang-He/plastisphere_toxicity
