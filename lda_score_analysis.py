#!/usr/bin/env python3
"""
宏基因组俩俩之间的LDA score计算与可视化
计算所有宏基因组特征之间的LDA score，筛选前30个进行柱状图可视化
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.preprocessing import StandardScaler
from typing import Optional, Tuple, Dict, Any, List
from dataset import BioSmokeDataset
from split_group import prepare_crc_smoke_groups
import os
from itertools import combinations

def get_metagenome_columns(df: pd.DataFrame, prefix: str = 'tax_') -> List[str]:
    """
    获取宏基因组特征列
    
    Parameters:
    -----------
    df : pd.DataFrame
        包含宏基因组数据的DataFrame
    prefix : str
        宏基因组列的前缀
    
    Returns:
    --------
    List[str]
        宏基因组特征列名列表
    """
    return [col for col in df.columns if col.startswith(prefix)]


def calculate_pairwise_lda_scores(df: pd.DataFrame, 
                                 metagenome_cols: List[str],
                                 group_col: str = 'group',
                                 min_samples_per_group: int = 5) -> pd.DataFrame:
    """
    计算宏基因组俩俩之间的LDA score
    
    Parameters:
    -----------
    df : pd.DataFrame
        包含分组和宏基因组数据的DataFrame
    metagenome_cols : List[str]
        宏基因组特征列列表
    group_col : str
        分组列名
    min_samples_per_group : int
        每组最少样本数要求
    
    Returns:
    --------
    pd.DataFrame
        包含LDA score结果的DataFrame
    """
    
    # 筛选有足够样本的分组
    group_counts = df[group_col].value_counts()
    valid_groups = group_counts[group_counts >= min_samples_per_group].index.tolist()
    
    if len(valid_groups) < 2:
        raise ValueError(f"需要至少2个分组，每个分组至少有{min_samples_per_group}个样本")
    
    print(f"有效分组: {valid_groups}")
    print(f"每个分组的样本数: {group_counts[valid_groups].to_dict()}")
    
    # 筛选有效分组的数据
    valid_df = df[df[group_col].isin(valid_groups)].copy()
    
    # 标准化宏基因组数据
    scaler = StandardScaler()
    metagenome_data = scaler.fit_transform(valid_df[metagenome_cols])
    
    # 准备结果存储
    lda_results = []
    
    # 计算所有宏基因组特征对之间的LDA score
    print(f"计算{len(metagenome_cols)}个宏基因组特征之间的LDA score...")
    
    for i, j in combinations(range(len(metagenome_cols)), 2):
        feature1 = metagenome_cols[i]
        feature2 = metagenome_cols[j]
        
        # 提取两个特征的数据
        X = metagenome_data[:, [i, j]]
        y = valid_df[group_col].values
        
        # 检查是否有足够的数据点
        if np.sum(~np.isnan(X).any(axis=1)) < min_samples_per_group * len(valid_groups):
            continue
        
        # 移除NaN值
        valid_mask = ~np.isnan(X).any(axis=1)
        X_clean = X[valid_mask]
        y_clean = y[valid_mask]
        
        if len(np.unique(y_clean)) < 2:
            continue
        
        try:
            # 执行LDA分析
            lda = LinearDiscriminantAnalysis()
            lda.fit(X_clean, y_clean)
            
            # 计算LDA score（使用判别函数的系数绝对值之和）
            lda_score = np.sum(np.abs(lda.coef_[0]))
            
            # 计算解释方差比例
            explained_variance_ratio = lda.explained_variance_ratio_[0] if len(lda.explained_variance_ratio_) > 0 else 0
            
            # 存储结果
            lda_results.append({
                'feature1': feature1,
                'feature2': feature2,
                'lda_score': lda_score,
                'explained_variance': explained_variance_ratio,
                'n_samples': len(X_clean),
                'n_groups': len(np.unique(y_clean))
            })
            
        except Exception as e:
            # 如果LDA计算失败，跳过这个特征对
            continue
    
    # 转换为DataFrame并排序
    if lda_results:
        lda_df = pd.DataFrame(lda_results)
        lda_df = lda_df.sort_values('lda_score', ascending=False)
        print(f"成功计算了{len(lda_df)}对宏基因组特征的LDA score")
        return lda_df
    else:
        raise ValueError("未能计算任何LDA score，请检查数据质量")


def filter_top_lda_scores(lda_df: pd.DataFrame, top_n: int = 30) -> pd.DataFrame:
    """
    筛选前N个LDA score最高的特征对
    
    Parameters:
    -----------
    lda_df : pd.DataFrame
        包含LDA score结果的DataFrame
    top_n : int
        要筛选的前N个特征对
    
    Returns:
    --------
    pd.DataFrame
        前N个LDA score最高的特征对
    """
    return lda_df.head(top_n).copy()


def plot_lda_scores_bar(lda_df: pd.DataFrame,
                       title: str = 'Top LDA Scores Between Metagenomic Features',
                       figsize: Tuple[int, int] = (16, 12),
                       colors: Optional[List[str]] = None,
                       save_path: Optional[str] = None,
                       dpi: int = 300) -> plt.Figure:
    """
    绘制LDA score柱状图
    
    Parameters:
    -----------
    lda_df : pd.DataFrame
        包含LDA score结果的DataFrame
    title : str
        图表标题
    figsize : Tuple[int, int]
        图表尺寸
    colors : Optional[List[str]]
        颜色列表
    save_path : Optional[str]
        保存路径
    dpi : int
        保存图片的DPI
    
    Returns:
    --------
    plt.Figure
        生成的图表对象
    """
    
    # 创建特征对标签
    feature_pairs = [f"{row['feature1']}\nvs\n{row['feature2']}" 
                     for _, row in lda_df.iterrows()]
    
    # 默认颜色
    if colors is None:
        colors = ['#FF6B6B', '#4ECDC4', '#45B7D1', '#96CEB4', '#FECA57', 
                  '#FF9FF3', '#54A0FF', '#5F27CD', '#00D2D3', '#FF9F43']
    
    # 美化样式
    sns.set(style='whitegrid', context='notebook', font_scale=1.1)
    fig, ax = plt.subplots(figsize=figsize)
    
    # 创建柱状图
    bars = ax.bar(range(len(lda_df)), lda_df['lda_score'], 
                  color=colors * (len(lda_df) // len(colors) + 1)[:len(lda_df)],
                  alpha=0.8, edgecolor='black', linewidth=0.5)
    
    # 设置图表标题和标签
    ax.set_title(title, fontsize=18, fontweight='bold', pad=20)
    ax.set_xlabel('Metagenomic Feature Pairs', fontsize=14)
    ax.set_ylabel('LDA Score', fontsize=14)
    
    # 设置x轴标签
    ax.set_xticks(range(len(lda_df)))
    ax.set_xticklabels(feature_pairs, rotation=45, ha='right', fontsize=10)
    
    # 添加数值标签
    for i, bar in enumerate(bars):
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2., height + 0.01,
                f'{height:.3f}', ha='center', va='bottom', fontsize=9, fontweight='bold')
    
    # 添加网格线
    ax.grid(True, alpha=0.3, axis='y')
    
    # 调整布局
    plt.tight_layout()
    
    # 保存图片
    if save_path:
        os.makedirs(os.path.dirname(save_path), exist_ok=True)
        plt.savefig(save_path, dpi=dpi, bbox_inches='tight')
        print(f"图表已保存至: {save_path}")
    
    return fig


def analyze_metagenome_lda_scores(df: pd.DataFrame,
                                 metagenome_prefix: str = 'tax_',
                                 group_col: str = 'group',
                                 top_n: int = 30,
                                 min_samples_per_group: int = 5,
                                 output_dir: Optional[str] = None) -> Dict[str, Any]:
    """
    宏基因组LDA score分析流水线
    
    Parameters:
    -----------
    df : pd.DataFrame
        原始数据
    metagenome_prefix : str
        宏基因组列前缀
    group_col : str
        分组列名
    top_n : int
        要显示的前N个LDA score
    min_samples_per_group : int
        每组最少样本数要求
    output_dir : Optional[str]
        输出目录
    
    Returns:
    --------
    Dict[str, Any]
        分析结果
    """
    
    results = {}
    
    # 1. 准备分组数据
    print("步骤1: 准备CRC吸烟与非吸烟分组...")
    df_grouped = prepare_crc_smoke_groups(df)
    results['group_counts'] = df_grouped[group_col].value_counts().to_dict()
    print(f"分组统计: {results['group_counts']}")
    
    # 2. 获取宏基因组特征
    print("步骤2: 获取宏基因组数据...")
    metagenome_cols = get_metagenome_columns(df_grouped, metagenome_prefix)
    results['n_metagenomes'] = len(metagenome_cols)
    print(f"宏基因组特征数量: {len(metagenome_cols)}")
    
    # 3. 计算LDA score
    print("步骤3: 计算宏基因组俩俩之间的LDA score...")
    lda_scores_df = calculate_pairwise_lda_scores(
        df_grouped, metagenome_cols, group_col, min_samples_per_group
    )
    results['lda_scores'] = lda_scores_df
    print(f"成功计算了{len(lda_scores_df)}对特征的LDA score")
    
    # 4. 筛选前N个LDA score
    print(f"步骤4: 筛选前{top_n}个LDA score...")
    top_lda_df = filter_top_lda_scores(lda_scores_df, top_n)
    results['top_lda_scores'] = top_lda_df
    print(f"筛选出前{len(top_lda_df)}个LDA score最高的特征对")
    
    # 5. 显示前10个结果
    print("\n前10个LDA score最高的特征对:")
    print(top_lda_df.head(10)[['feature1', 'feature2', 'lda_score', 'explained_variance']])
    
    # 6. 创建输出目录
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
        
        # 保存完整的LDA score结果
        lda_output_path = os.path.join(output_dir, 'complete_lda_scores.csv')
        lda_scores_df.to_csv(lda_output_path, index=False)
        print(f"完整LDA score结果已保存至: {lda_output_path}")
        
        # 保存前N个LDA score结果
        top_lda_output_path = os.path.join(output_dir, f'top_{top_n}_lda_scores.csv')
        top_lda_df.to_csv(top_lda_output_path, index=False)
        print(f"前{top_n}个LDA score结果已保存至: {top_lda_output_path}")
    
    return results


def main():
    """主函数"""
    print("=" * 80)
    print("宏基因组俩俩之间的LDA score计算与可视化")
    print("计算所有宏基因组特征之间的LDA score，筛选前30个进行柱状图可视化")
    print("=" * 80)
    
    # 加载数据
    ds = BioSmokeDataset()
    df = ds.merge_to_dataframe()
    
    print(f"总样本数: {df.shape[0]}")
    print(f"总特征数: {df.shape[1]}")
    
    # 创建输出目录
    output_dir = './results/lda_score_analysis'
    os.makedirs(output_dir, exist_ok=True)
    
    # 执行分析
    results = analyze_metagenome_lda_scores(
        df,
        metagenome_prefix='tax_',
        group_col='group',
        top_n=30,
        min_samples_per_group=5,
        output_dir=output_dir
    )
    
    # 绘制柱状图
    print("\n步骤5: 绘制LDA score柱状图...")
    plot_save_path = os.path.join(output_dir, 'top_30_lda_scores_barplot.png')
    
    fig = plot_lda_scores_bar(
        results['top_lda_scores'],
        title='Top 30 LDA Scores Between Metagenomic Features\n(CRC Smoking vs Non-smoking vs Control)',
        figsize=(18, 10),
        save_path=plot_save_path
    )
    
    print("\n分析结果汇总:")
    print(f"分组统计: {results['group_counts']}")
    print(f"分析的宏基因组特征数量: {results['n_metagenomes']}")
    print(f"计算的LDA score特征对数量: {len(results['lda_scores'])}")
    print(f"显示的前N个LDA score数量: {len(results['top_lda_scores'])}")
    print(f"柱状图已保存至: {plot_save_path}")


if __name__ == '__main__':
    main()