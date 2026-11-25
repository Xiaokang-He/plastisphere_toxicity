# 微塑料际形成重塑微塑料对河蚬的生态毒性：代码仓库说明

本仓库用于整理论文分析和绘图所使用的 R 脚本，目前仍在持续整理和补充中。

## 目前已整理的内容

目前仓库已收录以下分析脚本：

- 主图 Fig. 1：生理和生化指标分析。
- 主图 Fig. 5：肠道 16S α 多样性分析。
- RT-qPCR 验证分析。
- EPS 定量及 PARAFAC 荧光组分分析。
- 预实验虹吸效率时间序列分析。
- 肠道微生物组成分析，包括 *Acinetobacter* 和 *Burkholderia*。
- VFDB 相关丰度汇总及热图分析。
- RDA 拟合、envfit 分析及数据导出。
- Mantel 关联覆盖率图绘制。
- R 分析环境及依赖包安装脚本。

相关 R 脚本和说明文件仍在整理中，后续将继续补充。

## 生物信息学分析

转录组、肠道 16S rRNA 基因及宏基因组数据的生物信息学分析主要通过美吉云平台完成。所使用的软件工具及版本信息见论文 Methods 和 Supporting Information。

## 如何运行

1. 安装 R。
2. 在 RStudio 中打开 `Plastisphere_Microplastic_Toxicity.Rproj`。
3. 根据英文 README 中列出的文件名称，将相应输入数据放入 `data/input/`。
4. 运行 `scripts/install_packages.R` 安装依赖包。
5. 根据需要运行相应分析脚本，结果输出至 `results/`。

## 仓库目录

```text
scripts/
  install_packages.R
  main_figures/
  validation/
  ordination/
  supplementary/
data/
  input/
results/
```

论文题目：*Plastisphere formation reshapes the ecological toxicity of microplastics in Corbicula fluminea.*

GitHub: https://github.com/Xiaokang-He/plastisphere_toxicity
