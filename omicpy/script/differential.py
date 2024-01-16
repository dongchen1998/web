#-*- coding:utf-8 -*-

import pandas as pd
import numpy as np
import subprocess
import plotly.graph_objects as go
import os


def run_deseq(workdir, input_path, output_path, repetition):
    """
    运行R脚本，将输入文件进行差异分析，生成差异分析结果表格

    Args:
        workdir (str): 工作目录
        input_path (str): 输入文件路径
        output_path (str): 输出文件路径
        repetition (int): 重复次数
    """

    # R脚本的路径，需要师哥你改路径
    script_path = os.path.join(workdir,'Deseq2.R')

    # Rscript Deseq2.R -i 'input_file/expression_matrix_deseq2.csv' -o 'output_file/deseq2.tsv' -n 3

    cmd = [
        'Rscript', script_path,
        '--input', input_path,
        '--output', output_path,
        '--repetition', str(repetition),
    ]
    
    # 执行R脚本并捕获输出
    try:
        result = subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        return result.stdout
    except subprocess.CalledProcessError as e:
        return e.stderrs
    

def plot_volcano(workdir, deseq2_relust, genelist=None, width=1200, height=900, p_threshold=0.05, logFC_threshold=1, color_schemes=1, bubble_size=5, opacity=0.8, up_donw_info=True, x_fix=False, y_fix=False):
    """根据差异分析结果绘制火山图

    Args:
        workdir (str): 工作目录
        deseq2_relust (pd.DataFrame): 差异分析结果表格
        genelist (list): 需要标记的基因列表. Defaults to None.
        width (int): 图像宽度. Defaults to 1200.
        height (int): 图像高度. Defaults to 900.
        p_threshold (float): 显著性阈值. Defaults to 0.05.
        logFC_threshold (int): logFC阈值. Defaults to 1.
        color_schemes (int): 颜色方案. Defaults to 1.
        bubble_size (int): 气泡大小. Defaults to 5.
        opacity (float): 气泡透明度. Defaults to 0.8.
        up_donw_info (bool): 是否显示上调和下调基因的数量. Defaults to True.
        x_fix (bool): 是否固定x轴范围. Defaults to False.
        y_fix (bool): 是否固定y轴范围. Defaults to False.
    """

    # 读取数据
    deseq2_relust = pd.read_csv(deseq2_relust, sep='\t')

    # 预处理
    deseq2_relust = deseq2_relust[['Gene', 'log2FoldChange', 'padj']].copy()
    deseq2_relust.columns = ['Symbol', 'logFC', 'P.adjust']
    deseq2_relust.dropna(inplace=True)

    # 对logFC范围进行修正
    if x_fix == True: 
        deseq2_relust['logFC'] = np.clip(deseq2_relust['logFC'], -6, 6)
        
    deseq2_relust['P.adjust'] = deseq2_relust['P.adjust'].replace(0, 1e-300)
    pvalue = -np.log10(deseq2_relust['P.adjust'])

    # 对y轴范围进行修正
    if y_fix == True:
        pvalue= np.clip(pvalue, 0, 100)
    
    # 提取数据
    gene_names = deseq2_relust['Symbol'].values
    logFC = deseq2_relust['logFC'].values

    # 根据阈值筛选差异显著的基因，给差异最显著的基因添加标签
    significant = (np.abs(logFC) > logFC_threshold) & (deseq2_relust['P.adjust'] < p_threshold)
    upregulated = significant & (logFC > 0)
    downregulated = significant & (logFC < 0)
    nonsignificant = ~significant

    # 颜色设置
    colors = {
        1: ('#f08d1a', '#7fa4ca'),
        2: ('#f26c6a', '#54a857'),
        3: ('#c42121', '#15609b'),
        4: ('#df7415', '#3d9241')
    }
    up_color, down_color = colors.get(color_schemes, ('#f08d1a', '#7fa4ca'))

    # 绘图
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=logFC[upregulated], y=pvalue[upregulated], mode='markers',
                            marker=dict(color=up_color, size=bubble_size, sizemode='area',symbol='circle',opacity=opacity,
                                        # line_dict=1
                                        line=dict(color='black',width=0.4)), name='Up',
                            text=gene_names[upregulated]))
    fig.add_trace(go.Scatter(x=logFC[downregulated], y=pvalue[downregulated], mode='markers',
                            marker=dict(color=down_color, size=bubble_size, sizemode='area',symbol='circle',opacity=opacity,
                                        line=dict(color='black',width=0.4)), name='Down',
                            text=gene_names[downregulated]))
    fig.add_trace(go.Scatter(x=logFC[nonsignificant], y=pvalue[nonsignificant], mode='markers',
                            marker=dict(color='#A9A9A9', size=bubble_size, sizemode='area',symbol='circle',opacity=opacity,
                                        line=dict(color='black',width=0.4)), name='Nonsignificant',
                            text=gene_names[nonsignificant]))

    if up_donw_info == True:
    # 计算上调和下调基因的数量
        upregulated_num = np.sum(upregulated)
        downregulated_num = np.sum(downregulated)
    
        # 在图像中右上角添加上调基因的数量注释
        fig.add_annotation(xref="paper", yref="paper",
                        x=1, y=1, showarrow=False,
                        xanchor='right', yanchor='top',
                        text="Up: {} genes".format(upregulated_num),
                        font=dict(size=16),
                        align="right",
                        bgcolor="white",
                        borderpad=4)
        
        # 在图像中左上角添加下调基因的数量注释
        fig.add_annotation(xref="paper", yref="paper",
                        x=0, y=1, showarrow=False,
                        xanchor='left', yanchor='top',
                        text="Down: {} genes".format(downregulated_num),
                        font=dict(size=16),
                        align="left",
                        bgcolor="white",
                        borderpad=4)
        
    # 添加火山图的阈值线
    x_min = np.min(logFC)
    x_max = np.max(logFC)
    fig.update_layout(shapes=[
        dict(type="line", x0=x_min, x1=x_max, y0=-np.log10(p_threshold), y1=-np.log10(p_threshold), line=dict(color="Black", width=1, dash="dash")),
        dict(type="line", x0=logFC_threshold, x1=logFC_threshold, y0=0, y1=max(pvalue)+10, line=dict(color="Black", width=1, dash="dash")),
        dict(type="line", x0=-logFC_threshold, x1=-logFC_threshold, y0=0, y1=max(pvalue)+10, line=dict(color="Black", width=1, dash="dash"))
    ])

    # 添加基因名的注释
    if genelist is not None:
        highlighted_genes = np.isin(gene_names, genelist)
        fig.add_trace(go.Scatter(
            x=logFC[highlighted_genes], 
            y=pvalue[highlighted_genes], 
            mode='markers+text',
            marker=dict(color='red', size=bubble_size*1.2, symbol='circle', opacity=1),
            # 字体大小
            textfont=dict(size=12),
            name='Marker Genes',
            text=gene_names[highlighted_genes],
            textposition="top center"
        ))
        
    # 设置图像布局，并限制 y 轴的范围
    fig.update_layout(
        xaxis_title='log2 Fold Change',
        yaxis_title='-log10(p-value)',
        title='DE Analysis Volcano Plot',
        template="plotly_white",
        height=height,
        width=width,
    )

    if y_fix == True:
        fig.update_layout(yaxis_range=[0, 110])

    # 保存上调和下调基因列表
    upregulated_genes_path = os.path.join(workdir, 'output_file', 'upregulated_genes.txt')
    downregulated_genes_path = os.path.join(workdir, 'output_file', 'downregulated_genes.txt')
    deseq2_relust[upregulated]['Symbol'].to_csv(upregulated_genes_path, index=False, header=False)
    deseq2_relust[downregulated]['Symbol'].to_csv(downregulated_genes_path, index=False, header=False)

    # 保存为png，scale设置为4
    output_image_path = os.path.join(workdir, "output_file", "volcano.png")
    fig.write_image(output_image_path, scale=4)
    # 保存为html
    output_html_path = os.path.join(workdir, "output_file", "volcano.html")
    fig.write_html(output_html_path)

    # 返回路径
    return upregulated_genes_path, downregulated_genes_path, output_image_path, output_html_path



if __name__ == "__main__":

    # 运行差异分析
    run_deseq(
        "/Users/dongjiacheng/Desktop/Github/omic_analysis/differential_analysis/",
        "/Users/dongjiacheng/Desktop/Github/omic_analysis/differential_analysis/input_file/expression_matrix_deseq2.csv",
        "/Users/dongjiacheng/Desktop/Github/omic_analysis/differential_analysis/output_file/deseq2.tsv",
        3,
        )
    
    # 绘制火山图
    plot_volcano(
        "/Users/dongjiacheng/Desktop/Github/omic_analysis/differential_analysis/",
        "/Users/dongjiacheng/Desktop/Github/omic_analysis/differential_analysis/output_file/deseq2.tsv",
        genelist=['MYCTH_98562', 'MYCTH_2121938'],
        width=1200,
        height=900,
        p_threshold=0.05,
        logFC_threshold=1,
        color_schemes=1,
        bubble_size=8,
        opacity=0.8,
        up_donw_info=True,
        x_fix=False,
        y_fix=False,
    )

        