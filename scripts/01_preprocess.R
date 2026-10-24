## =============================================================================
## MTBLS17 非靶向 LC-MS 代谢组学分析 —— 上游数据预处理
##
## 流程：读入原始数据 → 峰检测 → 保留时间校正 → 跨样本对齐
##       → 缺失值填补 → 检出率过滤 → 技术重复合并 → QC 指标输出
## =============================================================================

library(xcms)
library(tidyverse)
library(MsExperiment)
library(BiocParallel)

## ===== 路径配置(按需修改) =====
data_dir <- "C:/mtbls17/data"
out_dir  <- "C:/mtbls17"
dir.create(out_dir, showWarnings = FALSE)

## ===== 整理样本清单，解析实验分组 =====
# 列出所有文件
files <- list.files(data_dir, pattern = "^Exp1_.*POS\\.[Cc][Dd][Ff]$",
                    full.names = TRUE)
# 建立样本元数据表 pd
pd <- data.frame(sample_name = sub("\\.[Cc][Dd][Ff]$", "", basename(files)),
  stringsAsFactors = FALSE)
# 从文件名提取分组
pd$group <- str_match(pd$sample_name, "^Exp1_(HCC|CRR|QC)")[, 2]
# 验证分组解析是否有NA
stopifnot(!any(is.na(pd$group))) 
print(table(pd$group))  
# 行名设为样本名
rownames(pd) <- pd$sample_name

## ===== 读入质谱数据 =====
raw <- readMsExperiment(spectraFiles = files, sampleData = pd)

## ===== 峰检测 =====
xdata <- findChromPeaks(raw, param = CentWaveParam(
  ppm = 25, peakwidth = c(5, 60), prefilter = c(3, 100),
  snthresh = 5, mzdiff = 0.01))
# 查看每样本检出峰数
tab <- table(chromPeaks(xdata)[, "sample"])
print(summary(as.numeric(tab)))

saveRDS(xdata, file.path(out_dir, "xdata_after_findChromPeaks.rds"))

## ===== 保留时间校正 =====
xdata <- adjustRtime(xdata, param = ObiwarpParam(binSize = 0.1))
saveRDS(xdata,file.path(out_dir, "xdata_after_adjustRtime.rds"))

## ===== 跨样本峰对齐=====
xdata <- groupChromPeaks(xdata, param = PeakDensityParam(
  sampleGroups = sampleData(xdata)$group, minFraction = 0.5, bw = 5, binSize = 0.01))
saveRDS(xdata, file.path(out_dir, "xdata_after_groupChromPeaks.rds"))

## ===== 缺失值填补 =====
xdata <- fillChromPeaks(xdata,param = ChromPeakAreaParam())

## ===== 导出特征矩阵 =====
ft <- featureValues(xdata, value = "into")
ft_def <- featureDefinitions(xdata)
saveRDS(list(xdata = xdata, ft = ft, ft_def = ft_def), file.path(out_dir, "xdata_final.rds"))

# 导出 CSV 
ft_def_df <- as.data.frame(ft_def)
ft_def_df <- ft_def_df[, !vapply(ft_def_df, is.list, logical(1))]
write.csv(cbind(ft_def_df, ft), file.path(out_dir, "feature_table.csv"),
          row.names = FALSE)

## ===== 检出率过滤 =====
grp <- pd$group[match(colnames(ft), pd$sample_name)]
pres <- apply(ft, 1, function(x)
  max(mean(!is.na(x[grp == "HCC"])), mean(!is.na(x[grp == "CRR"]))))
ft_f <- ft[pres > 0.8, ]

## ===== 缺失值统计填补 =====
feat_min <- apply(ft_f, 1, min, na.rm = TRUE)
ft_imp <- ft_f
na_idx <- is.na(ft_imp)
ft_imp[na_idx] <- rep(feat_min, ncol(ft_imp))[na_idx] / 2

saveRDS(list(ft = ft_imp, ft_def = ft_def[pres > 0.8, ]),
        file.path(out_dir, "feature_table_filtered.rds"))

## ===== 合并技术重复(每样本中位数) =====
cols <- colnames(ft_imp)
cols_clean <- sub("\\.[Cc][Dd][Ff]$", "", cols)
sample_id <- sub("[ab]_POS$", "", cols_clean)

stopifnot(all(table(sample_id) == 2))

df <- as.data.frame(t(ft_imp))
df$.id <- sample_id
agg <- aggregate(. ~ .id, data = df, FUN = median)
ft_sample <- t(agg[, -1])
colnames(ft_sample) <- agg$.id
rownames(ft_sample) <- rownames(ft_imp)
saveRDS(ft_sample, file.path(out_dir, "feature_table_per_sample.rds"))

## ===== 计算技术重复相关系数 =====
# 计算每个样本两次进样间的 Spearman 相关系数，
a_cols <- cols[grepl("a_POS", cols)]
b_cols <- sub("a_POS", "b_POS", a_cols)
cor_ab <- mapply(function(x, y) cor(ft_imp[, x], ft_imp[, y],
                                    use = "pairwise.complete.obs",
                                    method = "spearman"),
                 a_cols, b_cols)
saveRDS(cor_ab, file.path(out_dir, "cor_ab_technical_replicates.rds"))

## ===== 记录分析环境 =====
writeLines(capture.output(sessionInfo()),
           file.path(out_dir, "sessionInfo_01.txt"))
