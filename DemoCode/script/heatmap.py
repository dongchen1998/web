#-*- coding:utf-8 -*-

import subprocess
import os


def run_heatmap(workdir, input_path, output_path, color_up="default", color_down="default", color_mid="default", show_border=False, scale='row', cluster_rows=True, cluster_cols=False, cellwidth="20", cellheight="20", fontsize="10"):
    """根据输入的表达矩阵，生成热图

    Args:
        workdir (str): 工作目录
        input_path (str): 输入的表达矩阵文件路径
        output_path (str): 输出的热图文件路径
        color_up (str): 颜色上限. Defaults to "default".
        color_down (str): 颜色下限. Defaults to "default".
        color_mid (str): 颜色中间值. Defaults to "default".
        show_border (str): 是否显示边框. Defaults to 'FALSE'.
        scale (str): 标准化方式，可选值为"row"、"log2" Defaults to 'row'.
        cluster_rows (str): 是否对行进行聚类. Defaults to 'TRUE'.
        cluster_cols (str): 是否对列进行聚类. Defaults to 'FALSE'.
        cellwidth (str): 每个单元格的宽度. Defaults to "20".
        cellheight (str): 每个单元格的高度. Defaults to "20".
        fontsize (str): 字体大小. Defaults to "10".
    """

    # R脚本的路径
    script_path =os.path.join(workdir,'heatmap.R')

    # Rscript heatmap.R --input input_file/expression_matrix_heatmap.csv --output output_file/heatmap.png --color_up "default" --color_down "default" --color_mid "default" --show_border TRUE --scale "z-score" --cluster_rows TRUE --cluster_cols FALSE --cellwidth 20 --cellheight 20 --fontsize 10
    cmd = [
        'Rscript', script_path,
        '--input', input_path,
        '--output', output_path,
        '--color_up', color_up,
        '--color_down', color_down,
        '--color_mid', color_mid,
        '--show_border', str(show_border).upper(),
        '--scale', scale,
        '--cluster_rows', str(cluster_rows).upper(),
        '--cluster_cols', str(cluster_cols).upper(),
        '--cellwidth', str(cellwidth),  # 转换为字符串
        '--cellheight', str(cellheight),  # 转换为字符串
        '--fontsize', str(fontsize)  # 转换为字符串
    ]

    # 执行R脚本并捕获输出
    try:
        result = subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        return result.stdout
    except subprocess.CalledProcessError as e:
        return e.stderr


if __name__ == "__main__":

    # 示例调用
    run_heatmap(
        '/Users/dongjiacheng/Desktop/Github/omic_analysis/heatmap',
        '/Users/dongjiacheng/Desktop/Github/omic_analysis/heatmap/input_file/expression_matrix_heatmap.csv',
        '/Users/dongjiacheng/Desktop/Github/omic_analysis/heatmap/output_file/heatmap.png',
        color_up="#F76809",
        color_down="#0766AD",
        color_mid="#FFFFFF",
        show_border=True,
        scale='log2',
        cluster_rows=True,
        cluster_cols=False,
        cellwidth="20",
        cellheight="20",
        fontsize="10"
    )



