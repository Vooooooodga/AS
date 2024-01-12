library(dplyr)
library(tidyr)
library(ggplot2)
library(ggtext)
library(ggforce)
library(khroma)
library(ggExtra)
library(ggridges)
library(ggpubr)
#color_blind_friendly <- color("muted")
change_filename <- function(original_name, new_abbreviation) {
  return(gsub("xx", new_abbreviation, original_name))
}

#设置缩写，物种名与social levels对应关系
species_info <- data.frame(
  abbreviation = c("Am", "Bt", "Ed", "Cc", "Cj", "Mg"),
  full_name = c("Apis mellifera", "Bombus terrestris", "Euglossa dilemma", 
                "Ceratina calcarata", "Ceratina japonica", "Megaloptera genalis"),
  social_level = c("eusocial", "eusocial", "sub-social", "sub-social", "sub-social", "sub-social")
)
get_species_title <- function(abbreviation, custom_title_part) {
  species_row <- species_info[species_info$abbreviation == abbreviation, ]
  color <- ifelse(species_row$social_level == "eusocial", "#f5c242", "#4f71be") # Use hex color codes for HTML
  full_name_colored <- sprintf("<span style='color: %s;'><i>%s</i></span>", color, species_row$full_name)
  title <- paste(custom_title_part, "in", full_name_colored)
  return(title)
}

# 加载文件名
new_abbreviation <- "Am"  # 你想要的新缩写
gene_count_filename <- change_filename("gene_count_matrix_xx.csv", new_abbreviation)
gene_tpm_filename <- change_filename("gene_tpm_matrix_xx.csv", new_abbreviation)
transcript_tpm_filename <- change_filename("tid2gid_tpm_matrix_xx.csv", new_abbreviation)
transcript_count_filename <- change_filename("transcript_count_matrix_xx.csv", new_abbreviation)
dominant_cols_filename <- change_filename("xx.queen.txt", new_abbreviation)
subordinate_cols_filename <- change_filename("xx.worker.txt", new_abbreviation)
duplicated_genes_filename <- change_filename("duplicates.xx.txt", new_abbreviation)
MXE_filename <- change_filename("MXE.MATS.JC.xx.txt", new_abbreviation)
SE_filename <- change_filename("SE.MATS.JC.xx.txt", new_abbreviation)
# 读取 MXE 和 SE 文件
mxe_data <- read.table(MXE_filename, header = TRUE, sep = "\t", stringsAsFactors = FALSE) 
se_data <- read.table(SE_filename, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# 添加一个新列来表示事件类型
mxe_data <- mxe_data %>% mutate(Event_Type = "MXE")
se_data <- se_data %>% mutate(Event_Type = "SE")

# 确保两个数据框具有相同的列（如果没有，则可能需要选择和/或重命名列）
# 你可能需要根据你的实际数据进行调整
common_columns <- intersect(names(mxe_data), names(se_data))
mxe_data <- mxe_data[, common_columns]
se_data <- se_data[, common_columns]
# 合并数据框
data <- rbind(mxe_data, se_data)
#找到不重复的GeneID
#unique_genes <- unique(data$GeneID)
#data <- data %>% filter(FDR<0.05)
#unique_genes_sig <- unique(data$GeneID)
# 保存合并的数据到文件
#write.table(combined_data, "combined_AS_data.Cj.txt", sep = "\t", row.names = FALSE)

# 加载数据
#data <- read.table("combined_AS_data.Cj.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# 保留FDR < 0.05，且在个基因中DPSI绝对值最大的事件
data <- data %>%
  filter(FDR < 0.05) %>%  # Filter events with FDR < 0.05
  mutate(AbsDPSI = abs(as.numeric(IncLevelDifference))) %>%  # Calculate absolute DPSI and convert factors to numeric if necessary
  group_by(GeneID) %>%  # Group by gene
  filter(AbsDPSI == max(AbsDPSI)) %>%  # Filter for the row with the max absolute DPSI for each gene
  ungroup()

# 函数将字符串转换为数值向量
#str_to_num_vector <- function(str) {
#  as.numeric(unlist(strsplit(str, ",")))
#}
#过滤出显著的 ΔPSI 事件，例如 FDR < 0.05
#significant_events <- data %>% filter(FDR < 0.05)
# 保存结果到文件
#write.table(data, "DPSI_gene.Bt.txt", sep = "\t", row.names = FALSE)

# 可选：绘制 ΔPSI 分布的直方图
ggplot(data, aes(x = IncLevelDifference)) +
  geom_histogram(binwidth = 0.01, fill = "skyblue", color = "black") +
  labs(title = "Distribution of ΔPSI", x = "ΔPSI", y = "Frequency") +
  theme_minimal()

#--------表达差异分析--------
library(DESeq2)


# 加载数据
gene_count <- read.csv(gene_count_filename, row.names = 1)
dominant_cols <- readLines(dominant_cols_filename)
subordinate_cols <- readLines(subordinate_cols_filename)


# 创建样本信息DataFrame
sample_info <- data.frame(sample_type = rep(NA, ncol(gene_count)))
rownames(sample_info) <- colnames(gene_count)
sample_info$sample_type[rownames(sample_info) %in% dominant_cols] <- "queen"
sample_info$sample_type[rownames(sample_info) %in% subordinate_cols] <- "worker"
# 删除sample——type为NA的行，保持sample_info的格式不变
sample_info <- sample_info[!is.na(sample_info$sample_type), , drop = FALSE]

#读取基因绝对表达量矩阵
gene_tpm <- read.csv(gene_tpm_filename, row.names = 1)
#过滤在queen的组平均tpm小于1的基因，但若有一组平均tpm大于1，则保留
gene_tpm <- gene_tpm[rowMeans(gene_tpm[, dominant_cols]) >= 1, ]

#与worker组平均tpm小于1的基因
gene_tpm <- gene_tpm[rowMeans(gene_tpm[, subordinate_cols]) >= 1, ]


###### 计算每个基因在多少样本中表达量大于等于10
#expressed_in_samples <- rowSums(gene_count >= 20) 

# 计算样本总数的10%
#threshold <- ncol(gene_count) * 0.10 

# 过滤低表达基因
#filtered_gene_count <- gene_count[expressed_in_samples >= threshold, ]

#保存gene_tpm中的基因
filtered_gene_count <- gene_count[rownames(gene_count) %in% rownames(gene_tpm), ]
# 除去sampletype之外的列
filtered_gene_count <- filtered_gene_count[, colnames(filtered_gene_count) %in% rownames(sample_info)]
# 更新DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = filtered_gene_count,
                              colData = sample_info,
                              design = ~ sample_type)

# 运行DESeq2分析
dds <- DESeq(dds)
res <- results(dds)
# 检查 'res' 结果中的 'padj' 是否存在并非空
if("padj" %in% names(res)) {
  # 如果存在NA值，将它们替换为1（表示非显著）
  res$padj[is.na(res$padj)] <- 1
  # 计算significance
  res$significance <- -log10(res$padj)
} else {
  stop("Adjusted p-values (padj) are missing from the DESeq2 results.")
}
# 为了避免后续错误，我们不直接替换res的列，而是创建新的变量
log2FoldChange <- res$log2FoldChange
padj <- res$padj

# 使用默认值替换NA值
log2FoldChange[is.na(log2FoldChange)] <- 0
padj[is.na(padj)] <- 1
significance <- -log10(padj)

# 添加gene列
res$gene <- rownames(res)



#生成volcano_data，合并res与significant_events的DPSI
volcano_data <- cbind(as.data.frame(res), data$IncLevelDifference[match(rownames(res), data$GeneID)],data$FDR[match(rownames(res), data$GeneID)])
# volcano_data的最后两列名字改为DPSI和FDR-DPSI
colnames(volcano_data)[ncol(volcano_data)-1] <- "DPSI"
colnames(volcano_data)[ncol(volcano_data)] <- "FDR-DPSI"
# DPSI的NA值替换为0
volcano_data$DPSI[is.na(volcano_data$DPSI)] <- 0
#FDR的NA为1
volcano_data$`FDR-DPSI`[is.na(volcano_data$`FDR-DPSI`)] <- 1
#新列为显著的 ΔPSI 事件，例如 |ΔPSI| >= 0.1，FDR < 0.05
volcano_data$DSG <- ifelse(abs(volcano_data$DPSI) >= 0.00001 & volcano_data$`FDR-DPSI` < 0.05, TRUE, FALSE)  

# 设置DEG显著性阈值
log2fc_threshold <- 1
pvalue_threshold <- 0.05
significance_threshold <- -log10(pvalue_threshold)

# 新列为显著的 DEG 事件，例如 LFC >= 1，FDR < 0.05
volcano_data$DEG <- ifelse(volcano_data$significance > significance_threshold & abs(volcano_data$log2FoldChange) >= log2fc_threshold, TRUE, FALSE)
# merge geen type into one column, 根据significant和DSG将基因分为四类
volcano_data$GeneType <- ifelse(volcano_data$DEG == TRUE & volcano_data$DSG == TRUE, "DEG&DSG", 
                                ifelse(volcano_data$DEG == TRUE & volcano_data$DSG == FALSE, "DEG", 
                                       ifelse(volcano_data$DEG == FALSE & volcano_data$DSG == TRUE, "DSG", "Not significant")))
volcano_data$GeneType <- factor(volcano_data$GeneType, levels = c("Not significant", "DEG", "DSG", "DEG&DSG"))
##统计DEG的数量，在queen中显著上调的基因数量
DEG_in_queen <- sum(volcano_data$DEG[volcano_data$log2FoldChange>0] == TRUE)
#在worker中显著上调的基因数量
DEG_in_worker <- sum(volcano_data$DEG[volcano_data$log2FoldChange<0] == TRUE)
#DSG
DSG_in_queen <- sum(volcano_data$DSG[volcano_data$DPSI>0] == TRUE)
DSG_in_worker <- sum(volcano_data$DSG[volcano_data$DPSI<0] == TRUE)
#按照GeneType对volcano_data进行排序
volcano_data <- volcano_data[order(volcano_data$GeneType), ]
#保存volcano_data
#write.csv(volcano_data, file = "volcano_data.csv", row.names = FALSE)
custom_title_part <- "Differential Gene Expression in"
title <- get_species_title(new_abbreviation, custom_title_part)
#color_blind_friendly <- color("muted")
color_blind_friendly <- c("Not significant"="grey", "DEG"="#332288", "DSG"="#DDCC77", "DEG&DSG"="#AA4499")
p <- ggplot(volcano_data, aes(x=log2FoldChange, y=significance, color=GeneType)) +
  geom_point() +
  scale_color_manual(values=color_blind_friendly) +
  guides(color=guide_legend(title="GeneType")) +
  geom_hline(yintercept = significance_threshold, linetype="dashed", color = "black") +  # 添加 p 值阈值线
  geom_vline(xintercept=c(-log2fc_threshold, log2fc_threshold), linetype="dashed", color = "black") +  # 添加 log2FC 阈值线
  theme(plot.title = element_markdown()) +
  labs(title=title,
       x="Log2 Fold Change",
       y="-Log10 P-value")

print(p)

















