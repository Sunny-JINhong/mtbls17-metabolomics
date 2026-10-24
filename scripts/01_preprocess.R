library(xcms)
library(tidyverse)
library(MsExperiment)
library(BiocParallel)

## ===== 路径配置 =====
data_dir <- "C:/mtbls17/data"
out_dir  <- "C:/mtbls17"
dir.create(out_dir, showWarnings = FALSE)

files <- list.files(data_dir, pattern = "^Exp1_.*POS\\.[Cc][Dd][Ff]$",
                    full.names = TRUE)

pd <- data.frame(sample_name = sub("\\.[Cc][Dd][Ff]$", "", basename(files)),
  stringsAsFactors = FALSE)

pd$group <- str_match(pd$sample_name, "^Exp1_(HCC|CRR|QC)")[, 2]

stopifnot(!any(is.na(pd$group))) 
print(table(pd$group))  
rownames(pd) <- pd$sample_name

#读入质谱数据
raw <- readMsExperiment(spectraFiles = files, sampleData = pd)
#峰检测
xdata <- findChromPeaks(raw, param = CentWaveParam(
  ppm = 25, peakwidth = c(5, 60), prefilter = c(3, 100),
  snthresh = 5, mzdiff = 0.01))

tab <- table(chromPeaks(xdata)[, "sample"])
print(summary(as.numeric(tab)))
saveRDS(xdata, file.path(out_dir, "xdata_after_findChromPeaks.rds"))

#校正保留时间
xdata <- adjustRtime(xdata, param = ObiwarpParam(binSize = 0.1))
saveRDS(xdata,file.path(out_dir, "xdata_after_adjustRtime.rds"))
#峰对齐
xdata <- groupChromPeaks(xdata, param = PeakDensityParam(
  sampleGroups = sampleData(xdata)$group, minFraction = 0.5, bw = 5, binSize = 0.01))
saveRDS(xdata, file.path(out_dir, "xdata_after_groupChromPeaks.rds"))
#填补缺失值
xdata <- fillChromPeaks(xdata,param = ChromPeakAreaParam())
#导出特征表
ft <- featureValues(xdata, value = "into")
ft_def <- featureDefinitions(xdata)
saveRDS(list(xdata = xdata, ft = ft, ft_def = ft_def), file.path(out_dir, "xdata_final.rds"))
ft_def_df <- as.data.frame(ft_def)
ft_def_df <- ft_def_df[, !vapply(ft_def_df, is.list, logical(1))]
write.csv(cbind(ft_def_df, ft), file.path(out_dir, "feature_table.csv"),
          row.names = FALSE)

grp <- pd$group[match(colnames(ft), pd$sample_name)]
pres <- apply(ft, 1, function(x)
  max(mean(!is.na(x[grp == "HCC"])), mean(!is.na(x[grp == "CRR"]))))
ft_f <- ft[pres > 0.8, ]

feat_min <- apply(ft_f, 1, min, na.rm = TRUE)
ft_imp <- ft_f
na_idx <- is.na(ft_imp)
ft_imp[na_idx] <- rep(feat_min, ncol(ft_imp))[na_idx] / 2

saveRDS(list(ft = ft_imp, ft_def = ft_def[pres > 0.8, ]),
        file.path(out_dir, "feature_table_filtered.rds"))
writeLines(capture.output(sessionInfo()),
            file.path(out_dir, "sessionInfo_01.txt"))
