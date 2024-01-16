#!/usr/bin/env Rscript
library(DESeq2)
library(ggplot2)
library(plotly)
library(optparse)

## 测试
# Rscript Deseq2.R -i 'input_file/expression_matrix_deseq2.csv' -o 'output_file/deseq2.tsv' -n 3

# 定义命令行参数
option_list <- list(
  make_option(c("-i", "--input"), type = "character", default = "", help = "输入文件路径"),
  make_option(c("-o", "--output"), type = "character", default = "", help = "输出文件路径"),
  make_option(c("-n", "--repetition"), type = "integer", default = 3, help = "样本重复数") 
)

# 解析命令行参数
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# 检查输入文件是否存在
if (!file.exists(opt$input)) {
  stop("输入文件不存在: ", opt$input)
}

# 读取基因表达谱
count_df <- read.csv(file = opt$input, header = TRUE, row.names = 1)
rownames(count_df) <- gsub("-", ".", rownames(count_df))

# 过滤掉垃圾基因
count_df <- count_df[rowSums(count_df) > 10 & apply(count_df,1,function(x){ all(x > 0) }),]

# 根据重复数分组
n <- opt$repetition
colnames_split <- strsplit(colnames(count_df), "_")
group_names <- unique(sapply(colnames_split, `[`, 1))

if (length(group_names) != 2) {
  stop("输入的count表必须只包含两组数据")
}

# 创建样本分组信息
sample_info <- data.frame(row.names = colnames(count_df))
sample_info$group <- rep(group_names, each = n)

# 创建dds对象用于差异表达分析
dds <- DESeqDataSetFromMatrix(countData = count_df,
                              colData = DataFrame(sample_info),
                              design = ~ group)

# 执行DESeq分析
dds <- DESeq(dds)

# 获取差异表达结果
res <- results(dds)

# 保存结果为TSV文件
# write.table(as.data.frame(res), file = opt$output, sep = "\t", quote = FALSE, row.names = TRUE)
# write.csv(as.data.frame(res), file = opt$output, row.names = TRUE)

# 将差异表达分析的结果转换为data.frame
res_df <- as.data.frame(res)

# 将行名（基因名）作为新列“Gene”添加到data.frame
res_df$Gene <- rownames(res_df)

# 重新排列列，使得"Gene"列成为第一列
res_df <- res_df[, c("Gene", setdiff(names(res_df), "Gene"))]

# 保存结果为TSV文件，第一列列名为"Gene"
write.table(res_df, file = opt$output, sep = "\t", quote = FALSE, row.names = FALSE)

