#-*- coding:utf-8 -*-

import subprocess
import os

def run_pca(workdir, input_count_path, input_sample_path, output_png_path,output_hmtl_path):
    """
    根据输入的表达量矩阵和样本信息表，运行R脚本，生成PCA分析的结果

    Args:
        input_count_path: 输入的表达量矩阵路径
        input_sample_path: 输入的样本信息表路径
        output_png_path: 输出的PCA分析结果图路径
    """

    # R脚本的路径，需要师哥你改路径
    script_path = os.path.join(workdir, 'PCA.R')

    # Rscript PCA.R -c 'input_file/expression_matrix.csv' -s 'input_file/sample_info.csv' -o 'output_file/pca.png' -t 'output_file/pca.html'
    cmd = [
        'Rscript', script_path,
        '--input_count', input_count_path,
        '--input_sample', input_sample_path,
        '--output_png', output_png_path,
        '--output_html', output_hmtl_path,
    ]

    # 执行R脚本并捕获输出
    try:
        result = subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        return result.stdout
    except subprocess.CalledProcessError as e:
        return e.stderr


if __name__ == "__main__":

    # 示例调用
    run_pca(
        '/Users/dongjiacheng/Desktop/Github/omic_analysis/pca_analysis',
        '/Users/dongjiacheng/Desktop/Github/omic_analysis/pca_analysis/input_file/expression_matrix.csv',
        '/Users/dongjiacheng/Desktop/Github/omic_analysis/pca_analysis/input_file/sample_info.csv',
        '/Users/dongjiacheng/Desktop/Github/omic_analysis/pca_analysis/output_file/pca.png',
        '/Users/dongjiacheng/Desktop/Github/omic_analysis/pca_analysis/output_file/pca.html',
    )