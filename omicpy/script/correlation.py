#-*- coding:utf-8 -*-

import pandas as pd
import plotly.express as px
import plotly.graph_objects as corr_network
import networkx as nx
import os



def plot_corr_sample(workdir, gene_exp_input_path, width=900, height=600, color='geyser', method='pearson'):
    """
    根据输入文件路径中的标准化基因表达水平矩阵，绘制样本相关性热图
    
    Args:
        gene_exp_input_path (str): 标准化基因表达水平矩阵的文件路径
        workdir (str): 输出文件的工作目录
        width (int): 图片宽度
        height (int): 图片高度
        color (str): 颜色方案
        method (str): 计算相关性的方法，pearson或spearman
    """
    # 读取数据
    df_corr = pd.read_csv(gene_exp_input_path)
    df_corr = df_corr.iloc[:, 1:]
    df_corr = df_corr.round(2)
    # 计算df_corr的样本数量
    sample_num = df_corr.shape[1]

    # 选择计算相关性的方法
    if method == 'pearson':
        # 计算样本之间的相关性（使用皮尔森系数）
        correlation_matrix = df_corr.corr(method='pearson')
    else:
        correlation_matrix = df_corr.corr(method='spearman')

    # 创建注释文本，用于标记相关性值
    annotations = []
    font_size_control = int(10 - sample_num/10)
    for i, row in enumerate(correlation_matrix.values):
        for j, value in enumerate(row):
            if abs(value) > 0.95:  # 用星星标记大于0.8的相关性值
                text = "★★"
                font_size = font_size_control
            elif abs(value) > 0.8:  # 显示大于0.6的相关性值
                text = "★"
                font_size = font_size_control
            else:
                continue  # 不显示小于或等于0.7的值

            annotations.append(dict(
                x=j, y=i, text=text, xref='x1', yref='y1',
                font=dict(color='white' if abs(value) > 0.5 else 'black', size=font_size),
                showarrow=False))
            
    # 绘制热图
    fig = px.imshow(correlation_matrix,
                    color_continuous_scale=color,
                    zmin=-1, zmax=1,
                    labels=dict(color='Correlation'))
    
    # 更新布局以将图例放置在图的右侧
    fig.update_layout(
        width=width,
        height=height,
        annotations=annotations,
    )
    # 添加解释性标注到图外的空白处
    fig.add_annotation(
        xref="paper", yref="paper",
        x=1.05, y=1.07,  # 在图外的空白处定位
        text="★★:Correlation > 0.95<br>★: Correlation > 0.8",
        showarrow=False,
        align="left"
    )
    # 保存为高清图片
    output_image_path = os.path.join(workdir, "output_file", "corr_sample.png")
    fig.write_image(output_image_path, scale=4)

    # 保存为html文件
    output_html_path = os.path.join(workdir, "output_file", "corr_sample.html")
    fig.write_html(output_html_path)

    # 返回路径
    return output_image_path, output_html_path


def plot_corr_gene(workdir, gene_exp_input_path, width=900, height=600, color='geyser', method='pearson'):
    """
    根据输入文件路径中的标准化基因表达水平矩阵，绘制基因相关性热图
    
    Args:
        gene_exp_input_path (str): 标准化基因表达水平矩阵的文件路径
        workdir (str): 输出文件的工作目录
        width (int): 图片宽度
        height (int): 图片高度
        color (str): 颜色方案
        method (str): 计算相关性的方法，pearson或spearman
    """
    # 读取数据
    df_corr = pd.read_csv(gene_exp_input_path)
    # 预处理,只保留前30个基因
    gene_ids = df_corr.iloc[:30, 0]
    df_corr = df_corr.iloc[:30, 1:]
    df_corr = df_corr.round(2)

    # 计算df_corr有多少行
    sample_num = df_corr.shape[0]
    
    # 选择计算相关性的方法
    if method == 'pearson':
        # 计算样本之间的相关性（使用皮尔森系数）
        correlation_matrix = df_corr.T.corr(method='pearson')
    else:
        correlation_matrix = df_corr.T.corr(method='spearman')
    # 创建注释文本，用于标记相关性值
    annotations = []
    font_size_control = 10 - abs(int(sample_num/5))
    for i, row in enumerate(correlation_matrix.values):
        for j, value in enumerate(row):
            if abs(value) > 0.9:  # 用星星标记大于0.8的相关性值
                text = "★★"
                font_size = font_size_control
            elif abs(value) > 0.7:  # 显示大于0.7的相关性值
                text = "★"
                font_size = font_size_control
            else:
                continue  # 不显示小于或等于0.7的值

            annotations.append(dict(
                x=j, y=i, text=text, xref='x1', yref='y1',
                font=dict(color='white' if abs(value) > 0.5 else 'black', size=font_size),
                showarrow=False))
                 
    # 绘制相关性矩阵的热图
    fig = px.imshow(correlation_matrix,
                    color_continuous_scale=color,
                    zmin=-1,  # 设置颜色比例尺的最小值
                    zmax=1,  # 设置颜色比例尺的最大值
                    x=gene_ids, # 设置x轴为基因id
                    y=gene_ids) # 设置y轴为基因id
    # 更新布局以将图例放置在图的右侧
    fig.update_layout(
        width=width,
        height=height,
        annotations=annotations,
    )
    # 添加解释性标注到图外的空白处
    fig.add_annotation(
        xref="paper", yref="paper",
        x=1.05, y=1.07,  # 在图外的空白处定位
        text="★★:Correlation > 0.8<br>★: Correlation > 0.6",
        showarrow=False,
        align="left"
    )
    # 保存为高清图片
    output_image_path = os.path.join(workdir, "output_file", "corr_gene.png")
    fig.write_image(output_image_path, scale=4)

    # 保存为html文件
    output_html_path = os.path.join(workdir, "output_file", "corr_gene.html")
    fig.write_html(output_html_path)

    # 返回路径
    return output_image_path, output_html_path


def create_gene_network(workdir, gene_exp_input_path, width=1200, height=900, bubble_size=4, threshold=0.6, k_value=0.5, iterations=10, color='gnbu', method='pearson'):
    """
    根据输入文件路径中的标准化基因表达水平矩阵，绘制基因相关性网络图
    
    Args:
        gene_exp_input_path (str): 标准化基因表达水平矩阵的文件路径
        workdir (str): 输出文件的工作目录
        ...其他参数...
    """
    # 读取数据
    df_corr = pd.read_csv(gene_exp_input_path)
    # 预处理
    df_corr = df_corr.rename(columns={df_corr.columns[0]: 'gene_id'})
    df_corr = df_corr.set_index('gene_id')
    df_corr = df_corr.round(2)

    # df_corr 只保留前100个基因
    df_corr = df_corr.iloc[:100, :]

    # 选择计算相关性的方法
    if method == 'pearson':
        # 计算样本之间的相关性（使用皮尔森系数）
        correlation_matrix = df_corr.transpose().corr(method='pearson')
    else:
        correlation_matrix = df_corr.transpose().corr(method='spearman')
        
    # 使用networkx创建一个网络图
    G = nx.Graph()
    for gene1 in df_corr.index:
        for gene2 in df_corr.index:
            if gene1 != gene2:
                G.add_edge(gene1, gene2, weight=correlation_matrix.loc[gene1, gene2])

    threshold = threshold # 相关性系数
    edges = [(u, v) for (u, v, d) in G.edges(data=True) if abs(d['weight']) > threshold]
    G = G.edge_subgraph(edges).copy()  # 使用edges创建一个新的图，并使用copy()避免状态问题
    # pos设置,k越小则点越紧,iterations越大则点越稳定
    pos = nx.spring_layout(G, k=k_value, iterations=iterations)
    # 将位置作为节点属性添加到G中
    for node in G.nodes():
        G.nodes[node]['pos'] = pos[node]

    # 使用plotly创建网络图
    edge_x = []
    edge_y = []
    for edge in G.edges():
        x0, y0 = G.nodes[edge[0]]['pos']
        x1, y1 = G.nodes[edge[1]]['pos']
        edge_x.append(x0)
        edge_x.append(x1)
        edge_x.append(None)
        edge_y.append(y0)
        edge_y.append(y1)
        edge_y.append(None)
    
    # 分别创建正相关和负相关的边
    edge_x_pos, edge_y_pos = [], []
    edge_x_neg, edge_y_neg = [], []

    # 根据权重将边分为正负两组
    for edge in G.edges(data=True):
        x0, y0 = pos[edge[0]]
        x1, y1 = pos[edge[1]]
        if edge[2]['weight'] > 0:
            edge_x_pos.extend([x0, x1, None])
            edge_y_pos.extend([y0, y1, None])
        else:
            edge_x_neg.extend([x0, x1, None])
            edge_y_neg.extend([y0, y1, None])

    edge_trace_pos = corr_network.Scatter(
        x=edge_x_pos, y=edge_y_pos,
        line=dict(width=0.3, color='red'),
        hoverinfo='none',
        mode='lines'
    )
    
    edge_trace_neg = corr_network.Scatter(
        x=edge_x_neg, y=edge_y_neg,
        line=dict(width=0.3, color='blue'),
        hoverinfo='none',
        mode='lines'
    )

    node_x = []
    node_y = []
    for node in G.nodes():
        x, y = G.nodes[node]['pos']
        node_x.append(x)
        node_y.append(y)

    node_trace = corr_network.Scatter(
        x=node_x, y=node_y,
        mode='markers',
        hoverinfo='text',
        marker=dict(
            showscale=True,
            colorscale=color,
            colorbar=dict(
                thickness=15,
                title='Number of Significantly Correlated Nodes',
                xanchor='left',
                titleside='right'
            ),
            line_width=2))

    node_adjacencies = []
    node_text = []
    node_sizes = []  # 添加一个列表来存储基于节点连接数的大小

    # 根据基因数量设置节点大小
    if df_corr.shape[0] > 50:
        # 计算每个节点的连接数并设置节点大小
        for node, adjacencies in enumerate(G.adjacency()):
            node_adjacencies.append(len(adjacencies[1]))
            node_text.append(adjacencies[0])
            node_degree = len(adjacencies[1])
            scaled_size = 15 + (node_degree * bubble_size/10)
            node_sizes.append(scaled_size)
    
    else:
        # 计算每个节点的连接数并设置节点大小
        for node, adjacencies in enumerate(G.adjacency()):
            node_adjacencies.append(len(adjacencies[1]))
            node_text.append(adjacencies[0])
            node_degree = len(adjacencies[1])
            scaled_size = 15 + (node_degree * bubble_size/2)
            node_sizes.append(scaled_size)

    # norm = plt.Normalize(vmin=min(node_adjacencies), vmax=max(node_adjacencies)) # 将连接数映射到0-1范围  
    node_trace.marker.color = node_adjacencies
    node_trace.text = node_text
    node_trace.marker.size = node_sizes  # 更新marker的大小

    # 创建图形
    fig = corr_network.Figure(data = [edge_trace_pos, edge_trace_neg, node_trace],
                    layout=corr_network.Layout(
                        title='',
                        showlegend=False,
                        hovermode='closest',
                        margin=dict(b=0, l=0, r=0, t=40),
                        annotations=[
                            dict(
                                text="",
                                showarrow=False,
                                xref="paper", yref="paper",
                                x=0.005, y=-0.002)
                        ],
                        xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                        yaxis=dict(showgrid=False, zeroline=False, showticklabels=False))
                    )
    # 定义fig的布局，设置宽度和高度
    fig.update_layout(
        autosize=False,
        width=width, 
        height=height,
        template="plotly_white"
    )
    # 保存为高清图片
    output_image_path = os.path.join(workdir, "output_file", "corr_gene_network.png")
    fig.write_image(output_image_path, scale=4)

    # 保存为html文件
    output_html_path = os.path.join(workdir, "output_file", "corr_gene_network.html")
    fig.write_html(output_html_path)

    # 返回路径
    return output_image_path, output_html_path



if __name__ == "__main__":

    # 样本相关性热图
    plot_corr_sample(
                '/Users/dongjiacheng/Desktop/Github/omic_analysis/correlation_analysis',
                '/Users/dongjiacheng/Desktop/Github/omic_analysis/correlation_analysis/input_file/expression_matrix_corr.csv',
                width=900,
                height=600,
                color='geyser',
                method='pearson')
    
    # 基因相关性热图
    plot_corr_gene(
                '/Users/dongjiacheng/Desktop/Github/omic_analysis/correlation_analysis',
                '/Users/dongjiacheng/Desktop/Github/omic_analysis/correlation_analysis/input_file/expression_matrix_corr.csv',
                width=900,
                height=600,
                color='geyser',
                method='pearson')
    
    # 基因相关性网络图
    create_gene_network(
                '/Users/dongjiacheng/Desktop/Github/omic_analysis/correlation_analysis',
                '/Users/dongjiacheng/Desktop/Github/omic_analysis/correlation_analysis/input_file/expression_matrix_corr.csv',
                width=900,
                height=600,
                bubble_size=4,
                threshold=0.6,
                k_value=0.5,
                iterations=10,
                color='geyser',
                method='pearson')