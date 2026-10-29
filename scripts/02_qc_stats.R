## =============================================================================
## 质量控制与统计分析
## 输入：feature_table_per_sample.rds / cor_ab_technical_replicates.rds
## =============================================================================

library(pheatmap)
library(ropls)
library(tidyverse)
library(ggrepel)

out_dir <- "F:/study/meta"
fig_dir <- "F:/study/meta/figures"
dir.create(fig_dir, showWarnings = TRUE, recursive = TRUE)

# == 读入数据 ==
ft_sample_all <- readRDS(file.path(out_dir, "feature_table_per_sample.rds"))
cor_ab <- readRDS(file.path(out_dir, "cor_ab_technical_replicates.rds"))
#print(summary(cor_ab))

# == 剔除技术重复相关性低的样本 ==
threshold   <- 0.8
low_samples <- names(cor_ab)[cor_ab < threshold]  
low_ids     <- sub("[ab]_POS.[Cc][Dd][Ff]$", "", low_samples)
#cat("低于阈值", threshold, "的样本:", paste(low_ids, collapse = ", "), "\n")

ft_sample <- ft_sample_all[, !colnames(ft_sample_all) %in% low_ids, drop = FALSE]
saveRDS(ft_sample, file.path(out_dir, "feature_table_per_sample_clean.rds"))

## ===== 图 1：技术重复相关性分布 =====
cor_df <- data.frame(correlation = cor_ab)
p <- ggplot(cor_df, aes(x = correlation)) +
  geom_histogram(binwidth = 0.01, fill = "#4C72B0", color = "white", alpha = 0.85) +coord_cartesian(xlim = c(0.80, 1)) +
  labs(title = "技术重复Spearman相关系数分布（a vs b）",
       x = "Spearman 相关系数",y = "样本数" ) +
  theme_bw(base_size = 13)
ggsave(file.path(fig_dir, "technical_replicates_Spearman.pdf"),
       p, width = 8.5, height = 8, dpi = 300)

## ===== 图 2：样本间相关性热图 =====
# 计算样本间 Spearman 相关矩阵
cor_mat <- cor(ft_sample, method = "spearman")

# 样本分组注释
anno <- data.frame(group = str_match(colnames(ft_sample), "_(HCC|CRR)_")[, 2])
rownames(anno) <- colnames(ft_sample)

anno_colors <- list(group = c(HCC = "#E64B35", CRR = "#4DBBD5"))

# 画热图 
off <- cor_mat[upper.tri(cor_mat)]
lo <- quantile(off, 0.02)
p_heatmap <- pheatmap(cor_mat,
                      annotation_col  = anno,
                      annotation_row  = anno,
                      annotation_colors = anno_colors,
                      clustering_method = "complete",
                      show_colnames     = TRUE,
                      show_rownames     = TRUE,
                      fontsize_row  = 6,
                      fontsize_col  = 6,
                      angle_col     = 90,
                      treeheight_row    = 30,
                      treeheight_col    = 30,
                      color             = colorRampPalette(c("#3B4CC0", "#F7F7F7", "#B40426")
                                                           )(100),
                      breaks            = seq(lo, 1, length.out = 101), 
#                      main              = "样本间相关性热图（Spearman）",
)
ggsave(file.path(fig_dir, "sample_Spearman_heatmap.pdf"),
       p_heatmap, width = 20, height = 20, dpi = 300)

## ===== 图 3：PCA =====
# 计算log10值
mat <- t(log10(ft_sample))

# Pareto 缩放：每个特征减去均值、除以标准差的平方根
center <- colMeans(mat)
sc     <- sqrt(apply(mat, 2, sd))
mat_sc <- scale(mat, center = center, scale = sc)

pca <- prcomp(mat_sc)
ve  <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
print(ve[1:10])

scores <- as.data.frame(pca$x)
scores$sample <- rownames(scores)
scores$group  <- ifelse(grepl("_HCC_", scores$sample), "HCC", "CRR")

p_pca <- ggplot(scores, aes(PC1, PC2, color = group)) +
  stat_ellipse(level = 0.95, linewidth = 0.8, show.legend = FALSE) +
  geom_point(size = 2.5, alpha = 0.8) +
  scale_color_manual(values = c(HCC = "#E64B35", CRR = "#4DBBD5")) +
  labs(x = paste0("PC1 (", ve[1], "%)"),
       y = paste0("PC2 (", ve[2], "%)"),
       color = "分组") +
  theme_bw(base_size = 13)

ggsave(file.path(fig_dir, "fig_PCA.pdf"),
       p_pca, width = 10, height = 10, dpi = 300)

## ===== 图 4：OPLS-DA（含复核）=====
grp <- ifelse(grepl("_HCC_", rownames(mat_sc)), "HCC", "CRR")  #构造分组向量 (Y向量)
set.seed(123)  #固定随机数种子

pdf(NULL) 
opls_res <- opls(mat_sc, grp, predI = 1, orthoI = 1)  #建 OPLS-DA 模型
# 置换检验（200 次）
perm_res <- opls(mat_sc, grp, predI = 1, orthoI = 1, permI = 200)  #做 200 次置换检验
# 复核
opls_res2 <- opls(mat, grp, predI = 1, orthoI = 1, scaleC = "pareto")
dev.off()

print(opls_res2)

scores_opls <- data.frame(
  t1    = opls_res@scoreMN[, 1],
  to1   = opls_res@orthoScoreMN[, 1],
  sample = rownames(opls_res@scoreMN)
)
scores_opls$group <- ifelse(grepl("_HCC_", scores_opls$sample), "HCC", "CRR")

p_opls <- ggplot(scores_opls, aes(t1, to1, color = group)) +
  stat_ellipse(level = 0.95, linewidth = 0.8, show.legend = FALSE) +
  geom_point(size = 2.5, alpha = 0.8) +
  scale_color_manual(values = c(HCC = "#E64B35", CRR = "#4DBBD5")) +
  scale_x_continuous(limits = c(-30, 30)) +
  scale_y_continuous(limits = c(-100, 100)) +
  labs(x = "t_pred",
       y = "t_ortho",
       color = "Group") +
  theme_bw(base_size = 13)

ggsave(file.path(fig_dir, "fig_oplsda_scores_ggplot.pdf"), p_opls,
       width = 10, height = 10, dpi = 300)
       

## ===== 图 5：火山图 =====
ft_def    <- readRDS(file.path(out_dir, "feature_table_filtered.rds"))$ft_def
stopifnot(nrow(ft_def) == nrow(ft_sample))

fc_dir_cut <- log2(1.5)
p_cut <- 0.05

grp <- ifelse(grepl("_HCC_", colnames(ft_sample)), "HCC", "CRR")
hcc_idx <- which(grp == "HCC")
crr_idx <- which(grp == "CRR")

pval <- apply(ft_sample, 1, function(x)
  t.test(x[hcc_idx], x[crr_idx])$p.value)

fc <- rowMeans(ft_sample[, hcc_idx]) /
  rowMeans(ft_sample[, crr_idx])

ft_def_df <- as.data.frame(ft_def)

res <- data.frame(
  feature  = rownames(ft_sample),
  mz       = ft_def_df$mzmed,
  rt       = ft_def_df$rtmed,
  log2FC   = log2(fc),
  pval     = pval,
  padj     = p.adjust(pval, method = "BH"),
  row.names = NULL
)

write.table(res %>% select(m.z = mz, p.value = pval),
            file.path(out_dir, "mummichog_input_full.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

#上传mummichog_input_full.txt至metaboanalyst网站，选 "Functional Analysis [LC-MS1]"模块，p-value cutoff： 0.05
#下载分析结果到本地

deg <- res %>%
  filter(pval < 0.05 & abs(log2FC) > fc_dir_cut) %>%
  arrange(pval) %>%
  mutate(evidence = ifelse(padj < 0.05, "strict", "loose"))

write.csv(deg,file.path(out_dir, "differential_features_38.csv"),
          row.names = FALSE)
# 合并注释
mum <- read.csv("F:/study/meta/Functional_Analysis_[LC-MS1]/mummichog_matched_compound_all.csv", check.names = FALSE)
# 每个特征取Mass.Diff最小的候选
mum_top <- mum %>%
  arrange(`Mass.Diff`) %>%
  group_by(`Query.Mass`) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  select(mz = `Query.Mass`, compound_id = `Matched.Compound`,
         form = `Matched.Form`, mass_diff = `Mass.Diff`)

deg$mz_r <- round(deg$mz, 4)
mum_top$mz_r <- round(mum_top$mz, 4)

mum_top <- mum_top %>%
  arrange(mass_diff) %>%
  group_by(mz_r) %>%
  slice_head(n = 1) %>%
  ungroup()
deg_annot <- deg %>%
  left_join(mum_top, by = "mz_r") %>%
  select(-mz_r)

# 将化合物编号转化为化合物名称
# compound_id 是 KEGG 编号（如 C05466），用 KEGG API 查名称
ids <- unique(deg_annot$compound_id)
ids <- ids[!is.na(ids)]

get_name <- function(id) {
  name <-  tryCatch({
    txt <- readLines(paste0("https://rest.kegg.jp/get/", id), warn = FALSE)
    name_line <- grep("^NAME", txt, value = TRUE)
    if (length(name_line) == 0) return (NA_character_)
    gsub("^NAME\\s+", "", name_line[1])
  }, error = function(e) NA_character_)
  as.character(name)
}

id2name <- data.frame(
  compound_id   = ids,
  compound_name = vapply(ids, get_name, character(1)),
  stringsAsFactors = FALSE
)

deg_annot <- deg_annot %>%
  left_join(id2name, by = "compound_id")

write.csv(deg_annot,
          file.path(out_dir, "differential_features_38_annotated.csv"),
          row.names = FALSE)

res$regulation <- ifelse(res$pval < p_cut & res$log2FC >  fc_dir_cut, "上调",
                         ifelse(res$pval < p_cut & res$log2FC < -fc_dir_cut, "下调",
                                "无显著变化"))
res$regulation <- factor(res$regulation,
                         levels = c("上调", "下调", "无显著变化"))
label_df <- res %>%
  left_join(deg_annot %>% select(feature, compound_name), by = "feature") %>%
  filter(!is.na(compound_name)) %>%
  arrange(pval) %>%
  mutate(label = paste0(gsub(";", "", compound_name), " (", feature, ")"))

p_vol <- ggplot(res, aes(log2FC, -log10(pval), color = regulation)) +
  geom_point(alpha = 0.85, size = 2) +
  geom_text_repel(data = label_df,
                  aes(label = label),
                  size = 3.2, max.overlaps = 30,
                  box.padding = 0.5, segment.alpha = 0.5,
                  show.legend = FALSE) +
  scale_color_manual(values = c("上调" = "#E64B35",
                                "下调" = "#4DBBD5",
                                "无显著变化" = "grey70")) +
  geom_vline(xintercept = c(-fc_dir_cut, fc_dir_cut), linetype = 3, color = "grey40") +
  geom_hline(yintercept = -log10(p_cut), linetype = 3, color = "grey40") +
  labs(x = "log2 Fold Change (HCC/CRR)",
       y = "-log10 pvalue",color = NULL) +
  theme_bw(base_size = 13)
ggsave(file.path(fig_dir, "fig_volcano.pdf"), p_vol,
       width = 7, height = 6, dpi = 300)
  
## ===== 记录环境 =====
writeLines(capture.output(sessionInfo()),
           file.path(out_dir, "sessionInfo_02.txt"))
