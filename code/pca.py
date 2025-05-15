import pandas as pd
import numpy as np
from sklearn.impute import SimpleImputer
from sklearn.decomposition import PCA
#import matplotlib.pyplot as plt

"""
py=3.5.3
"""

# ----------------------
# 数据预处理与特征工程
# ----------------------
def preprocess_genotype_matrix(input_file, missing_strategy='most_frequent'):
    """
    处理多等位基因型矩阵，包含：
    1. 转置矩阵（样本×标记）
    2. 处理缺失值
    3. 多等位基因独热编码
    
    参数：
    input_file: 输入TSV文件路径
    missing_strategy: 缺失值填充策略（默认众数填充）
    
    返回：
    处理后的特征矩阵（DataFrame）
    """
    # 读取原始数据
    raw_df = pd.read_csv(input_file, sep='\t', index_col='Marker')

    # 定义缺失值
    raw_df.replace(["N", "-"], np.nan, inplace=True)
    
    # 缺失值处理
    imputer = SimpleImputer(strategy=missing_strategy)
    imputed_df = pd.DataFrame(
        imputer.fit_transform(raw_df),
        columns=raw_df.columns,
        index=raw_df.index
    )
    
    # 多等位基因独热编码
    encoded_dfs = []
    for marker in imputed_df.columns:
        # 为每个标记生成独热编码
        dummies = pd.get_dummies(
            imputed_df[marker], 
            prefix=marker,
            dummy_na=False  # 不单独生成缺失值列
        )
        encoded_dfs.append(dummies)
    
    encoded_dfs=pd.concat(encoded_dfs, axis=1)

    # 输出处理后的矩阵
    encoded_dfs.to_csv(processed_file, index=True, sep="\t")

    return encoded_dfs

# ----------------------
# PCA分析与可视化
# ----------------------
def perform_pca(feature_matrix, n_components=5):
    """
    执行PCA分析并可视化结果
    
    参数：
    feature_matrix: 预处理后的特征矩阵
    n_components: 保留的主成分数
    """
    # 执行PCA
    pca = PCA(n_components=n_components)
    pca_result = pca.fit_transform(feature_matrix)
    
    # 构建结果DataFrame
    result_df = pd.DataFrame(
        pca_result,
        columns=['PC{}'.format(i+1) for i in range(n_components)],
        index=feature_matrix.index
    )
    
    # 可视化
    #plt.figure(figsize=(10, 8))
    #plt.scatter(result_df['PC1'], result_df['PC2'], alpha=0.7)
    
    # 添加方差解释率
    #var_ratio = pca.explained_variance_ratio_
    #plt.xlabel(f'PC1 ({var_ratio[0]:.1%})')
    #plt.ylabel(f'PC2 ({var_ratio[1]:.1%})')
    
    # 添加样本标签
    #for sample in result_df.index:
    #    plt.annotate(sample, (result_df.loc[sample, 'PC1'], result_df.loc[sample, 'PC2']),
    #                 textcoords="offset points", xytext=(0,5), ha='center')
    
    #plt.title('PCA of Multi-allelic Genotypes')
    #plt.grid(True)
    #plt.show()

    # 输出特征值矩阵
    result_df.to_csv(score_file, index=True, sep="\t")

# ----------------------
# 主程序执行
# ----------------------
if __name__ == "__main__":
    # 输入文件路径（示例）
    input_tsv = "../../rawdata/data/kasp_matrix.t.tsv"
    # 处理后的基因型矩阵
    processed_file = "kasp_matrix.t.process.tsv"
    # 特征矩阵
    score_file = "scores.tsv"

    # 执行预处理
    feature_matrix = preprocess_genotype_matrix(input_tsv)
    print("生成特征矩阵维度: {}".format(feature_matrix.shape))
    
    # 执行PCA分析
    perform_pca(feature_matrix)
