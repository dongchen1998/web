#-*- coding:utf-8 -*-

import base64
import tempfile
import streamlit as st
import pandas as pd
import os


from script.correlation import plot_corr_sample, plot_corr_gene, create_gene_network
from script.pca import run_pca
from script.differential import run_deseq, plot_volcano
from script.enrich import run_enrich, plot_kegg_chart, plot_go_chart
from script.heatmap import run_heatmap


# streamlit run app.py

# 设置工作目录
workdir = '/Users/dongjiacheng/Desktop/web'

# 设置页面标题
st.title("Bulk RNA-seq 分析工具")
st.markdown("###### 欢迎，这是一个基于Python的组学分析工具，可以生成高质量的图片与交互式的页面")

# 分割线
st.markdown("---")

# 设置侧边栏
st.sidebar.markdown("### 选择功能页面")
page = st.sidebar.radio("功能列表", ["相关性分析", "主成分分析", "差异分析", "富集分析", "热图绘制"])


# 文件下载功能
def get_binary_file_downloader_html(bin_file, file_label='File'):
    with open(bin_file, 'rb') as f:
        data = f.read()
    bin_str = base64.b64encode(data).decode()
    href = f'<a href="data:file/octet-stream;base64,{bin_str}" download="{os.path.basename(bin_file)}">{file_label}</a>'
    return href


@st.cache_resource
def get_corr_data_sample_csv():
    corr_data_path = os.path.join(workdir, 'input_file', 'expression_matrix_corr.csv')
    with open(corr_data_path, "r", encoding="utf-8") as f:  # Open as regular text file
        data = f.read()
    return data

def corr_page():
    # 下载示例CSV文件
    corr_file = get_corr_data_sample_csv()
    st.download_button(
        label="下载示例文件",
        data=corr_file,
        file_name='expression_matrix_corr.csv',
        mime="text/csv;charset=utf-8",
    )
    
    uploaded_file = st.file_uploader("上传标准化后的表达转录水平数据", type=["csv"])
    if uploaded_file is not None:
        data = pd.read_csv(uploaded_file)
        st.write(data)
            
        # 用户输入区域
        st.sidebar.markdown("### 选择你的参数")
        width = st.sidebar.number_input("设置图片宽度", 900, 2000, 900)
        height = st.sidebar.number_input("设置图片高度", 600, 2000, 600)
        color = st.sidebar.selectbox("选择颜色方案", ['gnbu', 'geyser', 'portland', 'Temps'])
        method = st.sidebar.selectbox("选择相关性计算方法", ['pearson', 'spearman'])

        if st.button("生成样本相关性热图"):
            with tempfile.NamedTemporaryFile(delete=False, suffix='.csv') as tmp_file:
                tmp_file.write(uploaded_file.getvalue())
                tmp_file_path = tmp_file.name
            # Assuming plot_corr_sample returns a Plotly figure object and not file paths
            output_image_path, output_html_path = plot_corr_sample(workdir, tmp_file_path, width, height, color, method)
            st.image(output_image_path)
            st.markdown(get_binary_file_downloader_html(output_image_path, 'Download Image'), unsafe_allow_html=True)
            st.markdown(get_binary_file_downloader_html(output_html_path, 'Download HTML'), unsafe_allow_html=True)
            os.remove(tmp_file_path)

        if st.button("生成基因相关性热图"):
            with tempfile.NamedTemporaryFile(delete=False, suffix='.csv') as tmp_file:
                tmp_file.write(uploaded_file.getvalue())
                tmp_file_path = tmp_file.name
            # Assuming plot_corr_gene returns a Plotly figure object and not file paths
            output_image_path, output_html_path = plot_corr_gene(workdir, tmp_file_path, width, height, color, method)
            st.image(output_image_path)
            st.markdown(get_binary_file_downloader_html(output_image_path, 'Download Image'), unsafe_allow_html=True)
            st.markdown(get_binary_file_downloader_html(output_html_path, 'Download HTML'), unsafe_allow_html=True)
            os.remove(tmp_file_path)

        if st.button("生成基因网络图"):
            with tempfile.NamedTemporaryFile(delete=False, suffix='.csv') as tmp_file:
                tmp_file.write(uploaded_file.getvalue())
                tmp_file_path = tmp_file.name
            # Assuming create_gene_network returns a Plotly figure object and not file paths
            output_image_path, output_html_path = create_gene_network(workdir, tmp_file_path, width, height)
            st.image(output_image_path)
            st.markdown(get_binary_file_downloader_html(output_image_path, 'Download Image'), unsafe_allow_html=True)
            st.markdown(get_binary_file_downloader_html(output_html_path, 'Download HTML'), unsafe_allow_html=True)
            os.remove(tmp_file_path)


@st.cache_resource
def get_pca_data_sample1_csv():
    pca_data_path = os.path.join(workdir, 'input_file', 'expression_matrix_pca.csv')
    with open(pca_data_path, "r", encoding="utf-8") as f:
        data = f.read()
    return data

@st.cache_resource  
def get_pca_data_sample2_csv():
    pca_data_path = os.path.join(workdir, 'input_file', 'sample_info_pca.csv')
    with open(pca_data_path, "r", encoding="utf-8") as f:
        data = f.read()
    return data


def pca_page():
    st.markdown("如果你不知道上传什么文件，可以点击下方的下载示例CSV文件按钮")
    col1, col2 = st.columns(2)

    with col1:
        st.download_button(
            label="下载基因表达水平数据(ReadCount)",
            data=get_pca_data_sample1_csv(),
            file_name="expression_matrix_pca.csv",
            mime="text/csv;charset=utf-8",
        )

    with col2:
        st.download_button(
            label="下载样本分组数据",
            data=get_pca_data_sample2_csv(),
            file_name="sample_info_pca.csv",
            mime="text/csv;charset=utf-8",
        )
    # 文件上传
    uploaded_file1 = st.file_uploader("上传表达转录水平数据", type=["csv"])
    uploaded_file2 = st.file_uploader("上传样本分组数据", type=["csv"])
    
    if uploaded_file1 and uploaded_file2:
        data1 = pd.read_csv(uploaded_file1, index_col=0)
        data2 = pd.read_csv(uploaded_file2, index_col=0)
        st.write(data1)
        st.write(data2)
    
        if st.button("生成主成分分析图"):
            if uploaded_file1 and uploaded_file2:
                with tempfile.NamedTemporaryFile(delete=False, suffix='.csv') as tmp_file1, tempfile.NamedTemporaryFile(delete=False, suffix='.csv') as tmp_file2:
                    tmp_file1.write(uploaded_file1.getvalue())
                    tmp_file2.write(uploaded_file2.getvalue())
                    tmp_file_path1 = tmp_file1.name
                    tmp_file_path2 = tmp_file2.name
                # Assuming run_pca returns the output of the R script
                output_image_path = os.path.join(workdir,'output_file' ,'pca.png')
                output_html_path = os.path.join(workdir,'output_file' ,'pca.html')
                output = run_pca(workdir, tmp_file_path1, tmp_file_path2, output_image_path, output_html_path)
                st.image(output_image_path)
                st.markdown(get_binary_file_downloader_html(output_image_path, 'Download Image'), unsafe_allow_html=True)
                st.markdown(get_binary_file_downloader_html(output_html_path, 'Download HTML'), unsafe_allow_html=True)
                os.remove(tmp_file_path1)
                os.remove(tmp_file_path2)


@st.cache_resource
def get_diff_data_sample_csv():
    diff_data_path = os.path.join(workdir, 'input_file', 'expression_matrix_deseq2.csv')
    with open(diff_data_path, "r", encoding="utf-8") as f:
        data = f.read()
    return data
        
def diff_page():

    # 下载示例CSV文件
    diff_file = get_diff_data_sample_csv()
    st.download_button(
        label="下载示例文件",
        data=diff_file,
        file_name='expression_matrix_deseq2.csv',
        mime="text/csv;charset=utf-8",
    )

    uploaded_file = st.file_uploader("上传标准化后的表达转录水平数据", type=["csv"])
    if uploaded_file is not None:
        data = pd.read_csv(uploaded_file)
        st.write(data)
        
        # 用户输入区域
        st.sidebar.markdown("### 选择你的参数")
        width = int(st.sidebar.number_input("设置图片宽度", 900, 2000, 1200))
        height = st.sidebar.number_input("设置图片高度", 600, 2000, 1000)
        repetition = st.sidebar.number_input("设置重复次数", 3)
        p_threshold = st.sidebar.number_input("设置p值阈值", 0.05)
        logFC_threshold = st.sidebar.number_input("设置log2fc阈值", 1)
        color_schemes = st.sidebar.selectbox("选择颜色方案", ['1', '2', '3', '4'])
        bubble_size = st.sidebar.number_input("设置气泡大小", 1, 20, 8)
        opacity = st.sidebar.number_input("设置透明度", 0.1, 1.0, 0.8)
        # 是否显示上下调基因信息
        up_donw_info = st.sidebar.checkbox("是否显示上下调基因信息", value=True)
        # 让用户输入基因列表
        genelist = st.sidebar.text_input("输入基因列表", value='')


        if st.button("进行差异基因分析"):
            with tempfile.NamedTemporaryFile(delete=False, suffix='.csv') as tmp_file:
                tmp_file.write(uploaded_file.getvalue())
                tmp_file_path = tmp_file.name
            # Assuming plot_volcano returns a Plotly figure object and not file paths
            output_deseq2_result_path = os.path.join(workdir, 'output_file', 'deseq2_result.tsv')
            run_deseq(workdir, tmp_file_path,output_deseq2_result_path, repetition)
            output_image_path = os.path.join(workdir, 'output_file', 'volcano.png')
            output_html_path = os.path.join(workdir, 'output_file', 'volcano.html')
            output_up_genelist_path = os.path.join(workdir, 'output_file', 'upregulated_genes.txt')
            output_down_genelist_path = os.path.join(workdir, 'output_file', 'downregulated_genes.txt')
            plot_volcano(workdir, output_deseq2_result_path, genelist, width, height, p_threshold, logFC_threshold, color_schemes, bubble_size, opacity, up_donw_info)
            st.write(pd.read_csv(output_deseq2_result_path, sep='\t'))
            st.image(output_image_path)
            st.markdown(get_binary_file_downloader_html(output_image_path, 'Download Image'), unsafe_allow_html=True)
            st.markdown(get_binary_file_downloader_html(output_html_path, 'Download HTML'), unsafe_allow_html=True)
            st.markdown(get_binary_file_downloader_html(output_up_genelist_path, 'Download Up Genelist'), unsafe_allow_html=True)
            st.markdown(get_binary_file_downloader_html(output_down_genelist_path, 'Download Down Genelist'), unsafe_allow_html=True)
            os.remove(tmp_file_path)



@st.cache_resource
def get_enrich_data_sample_genelist_csv():
    enrich_data_path = os.path.join(workdir, 'input_file', 'gene_list_nc.txt')
    with open(enrich_data_path, "r", encoding="utf-8") as f:
        data = f.read()
    return data

@st.cache_resource
def get_enrich_data_sample_go_csv():
    enrich_data_path = os.path.join(workdir, 'input_file', 'go_gene_nc.txt')
    with open(enrich_data_path, "r", encoding="utf-8") as f:
        data = f.read()
    return data

@st.cache_resource
def get_enrich_data_sample_kegg1_csv():
    enrich_data_path = os.path.join(workdir, 'input_file', 'kegg_pathway_nc.txt')
    with open(enrich_data_path, "r", encoding="utf-8") as f:
        data = f.read()
    return data

@st.cache_resource
def get_enrich_data_sample_kegg2_csv():
    enrich_data_path = os.path.join(workdir, 'input_file', 'kegg_gene_nc.txt')
    with open(enrich_data_path, "r", encoding="utf-8") as f:
        data = f.read()
    return data

def enrich_page():
    # 下载示例CSV文件
    enrich_file1 = get_enrich_data_sample_genelist_csv()
    enrich_file2 = get_enrich_data_sample_go_csv()
    enrich_file3 = get_enrich_data_sample_kegg1_csv()
    enrich_file4 = get_enrich_data_sample_kegg2_csv()

    col1, col2, col3, col4= st.columns(4)
    with col1:
        st.download_button(
            label="下载基因列表文件",
            data=enrich_file1,
            file_name='gene_list_nc.txt',
            mime="text/txt;charset=utf-8",
        )
    with col2:
        st.download_button(
            label="下载GO基因注释",
            data=enrich_file2,
            file_name='go_gene_nc.txt',
            mime="text/txt;charset=utf-8",
        )
    with col3:
        st.download_button(
            label="下载KEGG通路注释",
            data=enrich_file3,
            file_name='kegg_pathway_nc.txt',
            mime="text/txt;charset=utf-8",
        )
    with col4:
        st.download_button(
            label="下载KEGG基因注释",
            data=enrich_file4,
            file_name='kegg_gene_nc.txt',
            mime="text/txt;charset=utf-8",
        )
    # 文件上传
    uploaded_file1 = st.file_uploader("上传基因列表文件", type=["txt"])
    uploaded_file2 = st.file_uploader("上传GO注释文件", type=["txt"])
    uploaded_file3 = st.file_uploader("上传KEGG通路注释文件", type=["txt"])
    uploaded_file4 = st.file_uploader("上传KEGG通路基因注释文件", type=["txt"])

    if uploaded_file1 and uploaded_file2 and uploaded_file3:
        data1 = pd.read_csv(uploaded_file1,sep='\t')
        data2 = pd.read_csv(uploaded_file2,sep='\t')
        data3 = pd.read_csv(uploaded_file3,sep='\t')
        data4 = pd.read_csv(uploaded_file4,sep='\t')
        st.write(data1)
        st.write(data2)
        st.write(data3)
        st.write(data4)

        # 用户输入区域
        st.sidebar.markdown("### 选择你的参数")
        width = int(st.sidebar.number_input("设置图片宽度", 900, 2000, 1200))
        height = st.sidebar.number_input("设置图片高度", 600, 2000, 1000)
        enrich_type = st.sidebar.selectbox("选择富集分析类型", ['GO', 'KEGG'])
        p_adjust = st.sidebar.number_input("设置p值阈值", 0.01, 1.0, 0.05)
        font_size = st.sidebar.number_input("设置字体大小", 15)
        chart_num = st.sidebar.number_input("设置最大显示数量", 30)
        chart_size = st.sidebar.number_input("设置图表显示大小", 30)
        pic_type = st.sidebar.selectbox("选择图表类型", ['bubble', 'bar'])
        color = st.sidebar.selectbox("选择颜色方案", ['gnbu', 'geyser', 'portland', 'Temps'])
        function_type = st.sidebar.selectbox("选择功能类型", ['All', 'CC', 'MF', 'BP'])

        if st.button("进行基因富集分析"):
            if enrich_type == 'GO':
                with tempfile.NamedTemporaryFile(delete=False, suffix='.txt') as tmp_file1, tempfile.NamedTemporaryFile(delete=False, suffix='.txt') as tmp_file2, tempfile.NamedTemporaryFile(delete=False, suffix='.txt') as tmp_file3:
                    tmp_file1.write(uploaded_file1.getvalue())
                    tmp_file2.write(uploaded_file2.getvalue())
                    tmp_file3.write(uploaded_file3.getvalue())
                    tmp_file_path1 = tmp_file1.name
                    tmp_file_path2 = tmp_file2.name
                    tmp_file_path3 = tmp_file3.name
                output_go_enrich_result_path = os.path.join(workdir, 'output_file', 'go.tsv')
                output_go_image_path = os.path.join(workdir, 'output_file', 'go.png')
                output_go_html_path = os.path.join(workdir, 'output_file', 'go.html')

                run_enrich(workdir, tmp_file_path1, output_go_enrich_result_path, p_adjust, enrich_type, tmp_file_path2, tmp_file_path3)
                plot_go_chart(workdir, output_go_enrich_result_path, width, height, p_adjust,font_size, chart_num, chart_size, pic_type, color, function_type)

                st.write(pd.read_csv(output_go_enrich_result_path, sep='\t'))
                st.image(output_go_image_path)
                st.markdown(get_binary_file_downloader_html(output_go_image_path, 'Download Image'), unsafe_allow_html=True)
                st.markdown(get_binary_file_downloader_html(output_go_html_path, 'Download HTML'), unsafe_allow_html=True)
                os.remove(tmp_file_path1)
                os.remove(tmp_file_path2)
                os.remove(tmp_file_path3)


            if enrich_type == 'KEGG':
                with tempfile.NamedTemporaryFile(delete=False, suffix='.txt') as tmp_file1, tempfile.NamedTemporaryFile(delete=False, suffix='.txt') as tmp_file2, tempfile.NamedTemporaryFile(delete=False, suffix='.txt') as tmp_file3, tempfile.NamedTemporaryFile(delete=False, suffix='.txt') as tmp_file4:
                    tmp_file1.write(uploaded_file1.getvalue())
                    tmp_file2.write(uploaded_file2.getvalue())
                    tmp_file3.write(uploaded_file3.getvalue())
                    tmp_file4.write(uploaded_file4.getvalue())
                    tmp_file_path1 = tmp_file1.name
                    tmp_file_path2 = tmp_file2.name
                    tmp_file_path3 = tmp_file3.name
                    tmp_file_path4 = tmp_file4.name
                output_kegg_enrich_result_path = os.path.join(workdir, 'output_file', 'kegg.tsv')
                output_kegg_image_path = os.path.join(workdir, 'output_file', 'kegg.png')
                output_kegg_html_path = os.path.join(workdir, 'output_file', 'kegg.html')

                run_enrich(workdir, tmp_file_path1, output_kegg_enrich_result_path, p_adjust, enrich_type, tmp_file_path3, tmp_file_path4)
                plot_kegg_chart(workdir, output_kegg_enrich_result_path, width, height, p_adjust, font_size, chart_num, chart_size, color, pic_type, function_type)

                st.write(pd.read_csv(output_kegg_enrich_result_path, sep='\t'))
                st.image(output_kegg_image_path)
                st.markdown(get_binary_file_downloader_html(output_kegg_image_path, 'Download Image'), unsafe_allow_html=True)
                st.markdown(get_binary_file_downloader_html(output_kegg_html_path, 'Download HTML'), unsafe_allow_html=True)
                os.remove(tmp_file_path1)
                os.remove(tmp_file_path2)
                os.remove(tmp_file_path3)
                os.remove(tmp_file_path4)


@st.cache_resource
def heatmap_data_sample_csv():
    heatmap_data_path = os.path.join(workdir, 'input_file', 'expression_matrix_heatmap.csv')
    with open(heatmap_data_path, "r", encoding="utf-8") as f:
        data = f.read()
    return data

def heatmap_page():
    
    # 下载示例CSV文件
    heatmap_file = heatmap_data_sample_csv()
    st.download_button(
        label="下载示例文件",
        data=heatmap_file,
        file_name='expression_matrix_heatmap.csv',
        mime="text/csv;charset=utf-8",
    )

    uploaded_file = st.file_uploader("上传标准化后的表达转录水平数据", type=["csv"])
    if uploaded_file is not None:
        data = pd.read_csv(uploaded_file)
        st.write(data)
        
        # 用户输入区域
        st.sidebar.markdown("### 选择你的参数")
        color_down = st.sidebar.text_input("设置下调颜色", value='default')
        color_up = st.sidebar.text_input("设置上调颜色", value='default')
        color_mid = st.sidebar.text_input("设置中性颜色", value='default')
        show_border = st.sidebar.checkbox("是否显示边框", value=True)
        scale = st.sidebar.selectbox("选择标准化方式", ['log2', 'z-score'])
        cluster_row = st.sidebar.checkbox("是否行聚类分析", value=True)
        cluster_cols = st.sidebar.checkbox("是否列聚类分析", value=False)
        cellwidth = st.sidebar.number_input("设置单元格宽度", 20)
        cellheight = st.sidebar.number_input("设置单元格高度", 20)
        fontsize = st.sidebar.number_input("设置字体大小", 10)

        if st.button("生成热图"):
            with tempfile.NamedTemporaryFile(delete=False, suffix='.csv') as tmp_file:
                tmp_file.write(uploaded_file.getvalue())
                tmp_file_path = tmp_file.name
            
            output_image_path = os.path.join(workdir, 'output_file', 'heatmap.png')
            run_heatmap(workdir, tmp_file_path, output_image_path, color_up, color_down, color_mid, show_border, scale, cluster_row, cluster_cols, cellwidth, cellheight, fontsize)
            st.image(output_image_path)
            st.markdown(get_binary_file_downloader_html(output_image_path, 'Download Image'), unsafe_allow_html=True)
            os.remove(tmp_file_path)


if page == "相关性分析":
    corr_page()
elif page == "主成分分析":
    pca_page() 
elif page == "差异分析":
    diff_page()
elif page == "富集分析":
    enrich_page()
elif page == "热图绘制":
    heatmap_page()
else:
    corr_page()  # 默认加载相关性分析界面



