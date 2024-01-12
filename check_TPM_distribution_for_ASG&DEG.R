###Followed by overlap_checking scripts

# ----------查看DSG和非DSG在queen和worekr的绝对表达量----------------
#读取基因绝对表达量矩阵
#gene_tpm <- read.csv(gene_tpm_filename, row.names = 1)
#根据volcano_data的基因筛选基因绝对表达量矩阵
#gene_tpm <- gene_tpm[rownames(gene_tpm) %in% volcano_data$gene, ]
#过滤在queen的组平均tpm小于1的基因，但若有一组平均tpm大于1，则保留
#gene_tpm <- gene_tpm[rowMeans(gene_tpm[, dominant_cols]) >= 1, ]
#与worker组平均tpm小于1的基因
#gene_tpm <- gene_tpm[rowMeans(gene_tpm[, subordinate_cols]) >= 1, ]
#过滤其它类型的样本
gene_tpm <- gene_tpm[, colnames(gene_tpm) %in% rownames(sample_info)]
#计算queen和worker的平均表达量
mean_tpm_queen <- rowMeans(gene_tpm[, sample_info$sample_type == "queen"], na.rm = TRUE)
mean_tpm_worker <- rowMeans(gene_tpm[, sample_info$sample_type == "worker"], na.rm = TRUE)
#volcano data新增平均表达量列
volcano_data <- cbind(volcano_data, Mean_TPM_Queen = mean_tpm_queen[match(volcano_data$gene, names(mean_tpm_queen))], Mean_TPM_Worker = mean_tpm_worker[match(volcano_data$gene, names(mean_tpm_worker))])

# 绘图，横坐标为queen平均表达量，纵坐标为worker平均表达量
custom_title_part <- "Gene Expression Profile"
title <- get_species_title(new_abbreviation, custom_title_part)
p <- ggplot(volcano_data, aes(x=Mean_TPM_Queen, y=Mean_TPM_Worker, color=GeneType)) +
  xlim(0, 2000) +
  ylim(0, 2000) +
  facet_zoom(xlim = c(1, 50), ylim = c(1, 50), horizontal = FALSE) +
  #geom_point(data=subset(volcano_data, GeneType=="Not significant")) +
  geom_point() +
  #geom_point(data=subset(volcano_data, GeneType=="DSG")) +
  #geom_point(data=subset(volcano_data, GeneType=="DEG&DSG")) +
  scale_color_manual(values=color_blind_friendly) +
  guides(color=guide_legend(title="Gene Type")) +
  labs(title=title,
       x="Mean TPM (Queen)",
       y="Mean TPM (Worker)") +
  theme(plot.title = element_markdown())
print(p)

custom_title_part <- "Expression distribution in queen and worker"
title <- get_species_title(new_abbreviation, custom_title_part)
p <- ggplot(volcano_data, aes(x=Mean_TPM_Queen, y=Mean_TPM_Worker, color=GeneType)) +
  xlim(0, 50) +
  ylim(0, 50) +
  #facet_zoom(xlim = c(0, 200), ylim = c(0, 200), horizontal = FALSE) +
  #geom_point() +
  #geom_point(data=subset(volcano_data, GeneType=="Not significant")) +
  geom_point(data=subset(volcano_data, GeneType!="Not significant")) +
  #geom_point(data=subset(volcano_data, GeneType=="DSG")) +
  #geom_point(data=subset(volcano_data, GeneType=="DEG&DSG")) +
  scale_color_manual(values=color_blind_friendly) +
  guides(color=guide_legend(title="Gene Type")) +
  labs(title=title,
       x="Mean TPM (Queen)",
       y="Mean TPM (Worker)") +
  theme(legend.position="bottom",plot.title = element_markdown())
print(p)

# 使用ggMarginal添加边缘直方图
ggMarginal(p, type = "histogram", alpha = 0.4, fill = "GeneType", margins = "both", bins = 20,
           groupColour = TRUE, groupFill = TRUE)
