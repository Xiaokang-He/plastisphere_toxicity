# 微塑料际形成重塑微塑料对河蚬的生态毒性：代码仓库说明

本仓库整理了论文部分统计分析、排序分析、Source Data 导出和绘图所用的 R 脚本。论文题目为：

*Plastisphere formation reshapes the ecological toxicity of microplastics in Corbicula fluminea.*

## 仓库范围

本仓库目前包含部分下游分析脚本，**尚不能视为覆盖论文全部自编分析流程的完整复现仓库**。转录组、16S 和宏基因组测序数据的主要生物信息学处理由美吉云平台完成，平台内部流程不属于本仓库脚本。

目前仍未找到生成微生物—宿主跨组学 Pearson Mantel 检验及 BH 校正结果的完整 R 脚本。仓库中的 Supplementary Fig. 23 脚本只根据已经报告的关联数量绘制覆盖率图，并不重新计算 Mantel 检验或 FDR。因此，在补齐原始 Mantel 检验代码之前，不能声称仓库完整复现了 Supplementary Fig. 22 及相关关联分析。

已找到并纳入 RDA 拟合及 envfit/Source Data 导出脚本；该脚本包含 Hellinger 转换、`vegan::rda()` 和 999 次置换的 `vegan::envfit()`，但尚未找到单独的原始 RDA 绘图脚本。

## 如何运行

1. 安装 R，并记录最终验证所用的准确版本。
2. 在 RStudio 中打开 `Plastisphere_Microplastic_Toxicity.Rproj`。
3. 将对应输入文件放入 `data/input/`。
4. 运行 `scripts/install_packages.R` 安装依赖包。
5. 根据英文 README 的输入文件表，逐个运行对应脚本。结果输出至 `results/`。

## 当前限制

- 目前只完成静态路径和内容检查，尚未在完整 R 环境及所需输入数据上逐个运行验证。
- 未包含美吉云平台的测序数据处理流程。
- 完整 Mantel 检验及 BH 校正代码仍缺失；网络拓扑和覆盖率绘图脚本不能替代原始检验代码。
- RDA 拟合和 envfit/Source Data 导出脚本已纳入，但原始 RDA 绘图脚本尚未找到。
- 粒径/羰基指数六面板脚本暂不纳入，因为需要先与最终羰基指数数据版本核对。
- KEGG Level 3 PCA 脚本暂不纳入，因为保存的版本还会计算当前 SI 未描述的 PERMANOVA、PERMDISP 和 ANOSIM，需要先核对。
- 未提供原始输入数据。README 中列出的是原始分析工作表文件名；论文 Source Data 工作簿不一定能直接替代这些输入文件。
- 正式公开前还需要记录 R 版本、包版本和 `sessionInfo()`，并确认每个脚本均对应最终稿件中的图表。
- 尚未选择代码复用许可，不要未经作者确认添加 LICENSE。

**注意：本仓库仍为待核查初稿。补齐 Mantel 代码并完成实际运行验证后，才适合作为最终公开代码仓库。**