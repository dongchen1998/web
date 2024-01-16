#!/usr/bin/env Rscript
library(optparse)
library(clusterProfiler)
library(dplyr)


## 测试
#  Rscript kegg_enrich.R -p 0.05 -i input_file/gene_list_nc.txt -o output_file/kegg.tsv -k input_file/kegg_pathway_nc.txt -f input_file/kegg_gene_nc.txt


# 设置命令行选项
option_list <- list(
  make_option(c("-p", "--pvalue"), type = "numeric", default = 0.05, help = "pvalue阈值", metavar = "pvalue"),
  make_option(c("-i", "--input"), type = "character", default = NULL, help = "输入文件路径", metavar = "file"),
  make_option(c("-o", "--output"), type = "character", default = NULL, help = "输出文件路径", metavar = "file"),
  make_option(c("-k", "--kegg_file1"), type = "character", default = NULL, help = "kegg通路文件", metavar = "file"),
  make_option(c("-f", "--kegg_file2"), type = "character", default = NULL, help = "kegg通路基因文件", metavar = "file")

)

# 解析命令行参数
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# 检查输入文件是否存在
if (!file.exists(opt$input)) {
  stop("输入文件不存在: ", opt$input)
}

# 背景文件，需要师哥你修改一下路径

ko2name <- read.delim(file = opt$kegg_file1,stringsAsFactors = FALSE)
ko2gene <- read.delim(file = opt$kegg_file2,stringsAsFactors = FALSE)

# 读取gene list
genelist <- read.csv(file = opt$input, row.names = 1)
genelist <- as.character(rownames(genelist))

# 如果genelist中有重复的基因，只保留一个
genelist <- unique(genelist)

# kegg
enrich_kegg <- enricher(genelist,
  TERM2GENE = ko2gene,
  TERM2NAME = ko2name,
  pAdjustMethod = "BH", # 使用FDR进行校正
  pvalueCutoff = opt$pvalue,
  qvalueCutoff = 1
)

# 保存结果到指定的输出文件
write.table(enrich_kegg, file = opt$output, sep = "\t", quote = FALSE, row.names = FALSE)
