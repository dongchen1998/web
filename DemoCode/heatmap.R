library(optparse)
library(pheatmap)
library(dplyr)
library(ggplot2)

# Rscript heatmap.R --input input_file/expression_matrix_heatmap.csv --output output_file/heatmap.png --color_up "default" --color_down "default" --color_mid "default" --show_border TRUE --scale "z-score" --cluster_rows TRUE --cluster_cols FALSE --cellwidth 20 --cellheight 20 --fontsize 10


# 定义命令行选项
option_list <- list(
    make_option(c("--input"), type = "character", default = NULL, help = "输入文件"),
    make_option(c("--output"), type = "character", default = NULL, help = "输出文件"),
    make_option(c("--color_up"), type = "character", default = "default", help = "高表达量颜色"),
    make_option(c("--color_down"), type = "character", default = "default", help = "低表达量颜色"),
    make_option(c("--color_mid"), type = "character", default = "default", help = "中等表达量颜色"),
    make_option(c("--show_border"), type = "logical", default = FALSE, help = "是否显示边框"),
    make_option(c("--cluster_rows"), type = "logical", default = TRUE, help = "是否对基因进行聚类"),
    make_option(c("--cluster_cols"), type = "logical", default = FALSE, help = "是否对样本进行聚类"),
    make_option(c("--scale"), type = "character", default = "log2", help = "标准化方式"),
    make_option(c("--cellwidth"), type = "numeric", default = 20, help = "单元格宽度"),
    make_option(c("--cellheight"), type = "numeric", default = 20, help = "单元格高度"),
    make_option(c("--fontsize"), type = "numeric", default = 10, help = "字体大小")
)

# 解析命令行选项
args <- parse_args(OptionParser(option_list = option_list))

row_zscores <- function(x) {
  row_means <- rowMeans(x)
  row_sds <- apply(x, 1, sd)
  scale(x, center = row_means, scale = row_sds)
}

draw_heatmap <- function(args) {
    data <- read.csv(args$input) # 读取csv文件
    data <- as.data.frame(data) # 转换为data.frame
    row.names(data) <- make.unique(data[,1]) # 使用make.unique确保行名的唯一性
    data <- data[,-1]
    # data <- row_zscores(data) # z-score标准化
    
    # 如果上下调与中等表达量参数为default，则运行
    if(args$color_up == "default" & args$color_down == "default" & args$color_mid == "default") {
      p <- pheatmap(data,
          # border_color = ifelse(args$show_border, args$color_border, NA),
          border = args$show_border,
          breaks = seq(-2, 2, length.out = 100),
          scale = "row",
          cluster_rows = args$cluster_rows,
          cluster_cols = args$cluster_cols,
          cellwidth = args$cellwidth, 
          cellheight = args$cellheight,
          fontsize = args$fontsize,
          legend = TRUE,
          silent = TRUE
    )
    }
    else {
      my_colors <- colorRampPalette(c(args$color_down, args$color_mid, args$color_up))(n = 256)
      p <- pheatmap(data,
        color = my_colors,
        # border_color = ifelse(args$show_border, args$color_border, NA),
        border = args$show_border,
        breaks = seq(-2, 2, length.out = 256),
        scale = "row",
        cluster_rows = args$cluster_rows,
        cluster_cols = args$cluster_cols,
        cellwidth = args$cellwidth, 
        cellheight = args$cellheight,
        fontsize = args$fontsize,
        legend = TRUE,
        silent = TRUE
    )
    }
    return(p)
}

draw_heatmap_log <- function(args) {
    data <- read.csv(args$input) # 读取csv文件
    data <- as.data.frame(data) # 转换为data.frame
    row.names(data) <- make.unique(data[,1]) # 使用make.unique确保行名的唯一性
    data <- data[,-1]
    data <- log2(data+1) # log2化

    # 如果上下调与中等表达量参数为default，则运行
    if(args$color_up == "default" & args$color_down == "default" & args$color_mid == "default") {
      p <- pheatmap(data,
          # border_color = ifelse(args$show_border, args$color_border, NA),
          border = args$show_border,
          cluster_rows = args$cluster_rows,
          cluster_cols = args$cluster_cols,
          cellwidth = args$cellwidth, 
          cellheight = args$cellheight,
          fontsize = args$fontsize,
          legend = TRUE,
          silent = TRUE
    )
    }
    else {
      my_colors <- colorRampPalette(c(args$color_down, args$color_mid, args$color_up))(n = 256)
      p <- pheatmap(data,
          color = my_colors,
          # border_color = ifelse(args$show_border, args$color_border, NA),
          border = args$show_border,
          cluster_rows = args$cluster_rows,
          cluster_cols = args$cluster_cols,
          cellwidth = args$cellwidth, 
          cellheight = args$cellheight,
          fontsize = args$fontsize,
          legend = TRUE,
          silent = TRUE
      )
    return(p)
}
}

# 如果标准化方式为z-score，则运行
if(args$scale == "z-score")
    p <- draw_heatmap(args)

# 如果标准化方式为log2，则运行
if(args$scale == "log2")
    p <- draw_heatmap_log(args)

# 保存图片
ggsave(args$output, plot = p, width = 10, height = 10, dpi = 600)