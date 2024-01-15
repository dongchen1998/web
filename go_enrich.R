#!/usr/bin/env Rscript
library(optparse)
library(clusterProfiler)
library(dplyr)


## 测试
#  Rscript go_enrich.R -p 0.05 -i 'input_file/gene_list_nc.txt' -o 'output_file/go.tsv' -g 'enrich_background_file/go_gene_nc.txt'

# 设置命令行选项
option_list <- list(
  make_option(c("-p", "--pvalue"), type = "numeric", default = 0.05, help = "pvalue阈值", metavar = "pvalue"),
  make_option(c("-i", "--input"), type = "character", default = NULL, help = "输入文件路径", metavar = "file"),
  make_option(c("-g", "--go_file"), type = "character", default = NULL, help = "GO文件", metavar = "file"),
  make_option(c("-o", "--output"), type = "character", default = NULL, help = "输出文件路径", metavar = "file")
)

# 解析命令行参数
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# 检查输入文件是否存在
if (!file.exists(opt$input)) {
  stop("输入文件不存在: ", opt$input)
}

# 背景文件，需要师哥你修改一下路径
go2name <- read.delim('/Users/dongjiacheng/Desktop/coder/web/enrich_background_file/go2name.txt', stringsAsFactors=FALSE)
go2gene_df <- read.delim(opt$go_file, stringsAsFactors=FALSE)


# 读取gene list
genelist <- read.csv(file = opt$input, row.names = 1)
genelist <- as.character(rownames(genelist))

# 如果genelist中有重复的基因，只保留一个
genelist <- unique(genelist)

go2gene_list <- split(go2gene_df, f = go2gene_df$CLASS)

enrich_MF = enricher(genelist, TERM2GENE=go2gene_list[['MF']][c(1,2)], TERM2NAME=go2name, pvalueCutoff = 1, qvalueCutoff = 1)

enrich_BP = enricher(genelist, TERM2GENE=go2gene_list[['BP']][c(1,2)], TERM2NAME=go2name, pvalueCutoff = 1, qvalueCutoff = 1)

enrich_CC = enricher(genelist, TERM2GENE=go2gene_list[['CC']][c(1,2)], TERM2NAME=go2name, pvalueCutoff = 1, qvalueCutoff = 1)

# 合并三个表
enrich_MF = as.data.frame(enrich_MF)
enrich_BP = as.data.frame(enrich_BP)
enrich_CC = as.data.frame(enrich_CC)
enrich_all <- bind_rows(
  mutate(enrich_MF, category = "MF"),
  mutate(enrich_BP, category = "BP"),
  mutate(enrich_CC, category = "CC")
)

# 保存结果到指定的输出文件
write.table(enrich_all, file = opt$output, sep = "\t", quote = FALSE, row.names = FALSE)