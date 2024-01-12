#------------查看DEG和DSG在表达量上是否存在差异----------------
#从volcano_data中找到MeanTPM_Queen和MeanTPM_Worker表达量最小的DSG，设置为两个阈值
DSG_threshold_queen <- min(volcano_data[volcano_data$GeneType == "DSG", "Mean_TPM_Queen"])
DSG_threshold_worker <- min(volcano_data[volcano_data$GeneType == "DSG", "Mean_TPM_Worker"])
#根据两个阈值来过滤volcano_data
volcano_data <- volcano_data[volcano_data$Mean_TPM_Queen >= DSG_threshold_queen & volcano_data$Mean_TPM_Worker >= DSG_threshold_worker, ]
#从volcanodata读取DEG和DSG的基因列表
DEG_DSG_list <- volcano_data[volcano_data$GeneType == "DEG&DSG", "gene"]
DEG_list <- c(volcano_data[volcano_data$GeneType == "DEG", "gene"], DEG_DSG_list)
DSG_list <- c(volcano_data[volcano_data$GeneType == "DSG", "gene"], DEG_DSG_list)
Not_significant_list <- volcano_data[volcano_data$GeneType == "Not significant", "gene"]
#从volcanodata读取queen上调的基因列表
DEG_up_list <- volcano_data[volcano_data$GeneType == "DEG" & volcano_data$log2FoldChange > 0, "gene"]
#如果DEG&DSGDSG上调，加入到DEG_up_list
DEG_up_list <- c(DEG_up_list, volcano_data[volcano_data$GeneType == "DEG&DSG" & volcano_data$log2FoldChange > 0, "gene"])
#worker上调的基因列表
DEG_down_list <- volcano_data[volcano_data$GeneType == "DEG" & volcano_data$log2FoldChange < 0, "gene"]
#如果DEG&DSGDSG在worker上调，加入到DEG_down_list
DEG_down_list <- c(DEG_down_list, volcano_data[volcano_data$GeneType == "DEG&DSG" & volcano_data$log2FoldChange < 0, "gene"])
#从volcanodata读取queen上调和worker的DPSI上调的基因列表
DSG_up_list <- volcano_data[volcano_data$GeneType == "DSG" & volcano_data$DPSI > 0, "gene"]
DSG_down_list <- volcano_data[volcano_data$GeneType == "DSG" & volcano_data$DPSI < 0, "gene"]
#如果DEG&DSGDSG在queen上调，加入到DSG_up_list
DSG_up_list <- c(DSG_up_list, volcano_data[volcano_data$GeneType == "DEG&DSG" & volcano_data$DPSI > 0, "gene"])
#如果DEG&DSGDSG在worker上调，加入到DSG_down_list
DSG_down_list <- c(DSG_down_list, volcano_data[volcano_data$GeneType == "DEG&DSG" & volcano_data$DPSI < 0, "gene"])

#从gene_tpm中提取DEG和DSG的基因表达量
DEG_tpm <- gene_tpm[DEG_list, colnames(gene_tpm) %in% rownames(sample_info)]
DEG_up_tpm <- gene_tpm[DEG_up_list, colnames(gene_tpm) %in% rownames(sample_info)]
DEG_down_tpm <- gene_tpm[DEG_down_list, colnames(gene_tpm) %in% rownames(sample_info)]
DSG_tpm <- gene_tpm[DSG_list, colnames(gene_tpm) %in% rownames(sample_info)]
DSG_up_tpm <- gene_tpm[DSG_up_list, colnames(gene_tpm) %in% rownames(sample_info)]
DSG_down_tpm <- gene_tpm[DSG_down_list, colnames(gene_tpm) %in% rownames(sample_info)]
Not_significant_tpm <- gene_tpm[Not_significant_list, colnames(gene_tpm) %in% rownames(sample_info)]

combined_df <- bind_rows(
  mutate(DEG_up_tpm, gene_type="worker-DEG"),
  mutate(DEG_down_tpm, gene_type="queen-DEG"),
  mutate(DSG_up_tpm, gene_type="queen-DSG"),
  mutate(DSG_down_tpm, gene_type="worker-DSG"),
  mutate(Not_significant_tpm, gene_type="Not_significant")
)
#combineddf增加一列gene
combined_df$gene <- rownames(combined_df)
long_df <- gather(combined_df, "sample", "tpm", -gene, -gene_type)
#longdf根据sample_info加一列sample_type
long_df$sample_type <- sample_info[long_df$sample, "sample_type"]
#对Longdf的gene_type进行因子化
long_df$gene_type <- factor(long_df$gene_type, levels=c("worker-DEG", "queen-DEG","worker-DSG", "queen-DSG", "Not_significant"))
#对longdf按照gene_type排序
long_df <- arrange(long_df, gene_type)
#绘制箱型图
# 对比1: DEG_up和not significant在worker样本中的表达量对比
worker_df <- subset(long_df, sample_type == "worker" & gene_type %in% c("worker-DEG", "queen-DEG","worker-DSG", "queen-DSG", "Not_significant"))
# 对比2: DEG_down和not significant在queen样本中的表达量对比
queen_df <- subset(long_df, sample_type == "queen" & gene_type %in% c("worker-DEG", "queen-DEG","worker-DSG", "queen-DSG", "Not_significant"))

# 绘制箱型图
#设置颜色
box_colors <- c("worker-DEG" = "darkblue", "queen-DEG" = "#f5c242", "worker-DSG" = "#56B4E9", "queen-DSG" = "#F0E442", "Not_significant" = "grey")
custom_title_part <- "Expression in Queen Samples"
title <- get_species_title(new_abbreviation, custom_title_part)
p <- ggplot(queen_df, aes(x = gene_type, y = log1p(tpm), fill = gene_type)) +
  geom_boxplot() +
  scale_fill_manual(values=box_colors) +
  #stat_summary(fun.y = mean, geom = "point", shape = 20, size = 3, color = "darkred") +
  theme(plot.title = element_markdown()) +
  guides(fill=FALSE) +
  labs(title = title,
       x = "Gene Type",
       y = "log(TPM + 1)")
# 使用stat_compare_means进行两两比较并添加到图上
p <- p + stat_compare_means(comparisons = list( c("queen-DEG", "Not_significant"), 
                                                c("queen-DEG", "queen-DSG"), 
                                                c("queen-DEG", "worker-DEG"), 
                                                c("queen-DEG", "worker-DSG"), 
                                                c("queen-DSG", "Not_significant"), 
                                                c("queen-DSG", "worker-DEG"), 
                                                c("queen-DSG", "worker-DSG"), 
                                                c("worker-DEG", "Not_significant"), 
                                                c("worker-DEG", "worker-DSG"), 
                                                c("worker-DSG", "Not_significant")), 
                            label = "p.signif", method = "wilcox.test", hide.ns = TRUE)
p

custom_title_part <- "Expression in Worker Samples"
title <- get_species_title(new_abbreviation, custom_title_part)
p <- ggplot(worker_df, aes(x = gene_type, y = log1p(tpm), fill = gene_type)) +
  geom_boxplot() +
  scale_fill_manual(values=box_colors) +
  theme(plot.title = element_markdown()) +
  guides(fill=FALSE) +
  labs(title = title,
       x = "Gene Type",
       y = "log(TPM + 1)")
p <- p + stat_compare_means(comparisons = list( c("queen-DEG", "Not_significant"), 
                                                c("queen-DEG", "queen-DSG"), 
                                                c("queen-DEG", "worker-DEG"), 
                                                c("queen-DEG", "worker-DSG"), 
                                                c("queen-DSG", "Not_significant"), 
                                                c("queen-DSG", "worker-DEG"), 
                                                c("queen-DSG", "worker-DSG"), 
                                                c("worker-DEG", "Not_significant"), 
                                                c("worker-DEG", "worker-DSG"), 
                                                c("worker-DSG", "Not_significant")), 
                            label = "p.signif", method = "wilcox.test", hide.ns = TRUE)
p
