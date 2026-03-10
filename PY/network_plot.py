#!/usr/bin/env python3
"""
CRC早期与晚期吸烟与非吸烟的宏基因组相关性网络图分析
分组：早期（tnm0, tnm1）和晚期（tnm2, tnm3, tnm4）
对比：CRC早期吸烟、CRC早期非吸烟、CRC晚期吸烟、CRC晚期非吸烟四组
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import networkx as nx
from typing import Optional, Tuple, Dict, Any, List
from dataset import BioSmokeDataset
import os
from split_group import prepare_early_late_groups


def get_metagenome_columns(df: pd.DataFrame) -> List[str]:
    """获取宏基因组列"""
    return [col for col in df.columns if col.startswith('tax_')]


def calculate_mean_expression(df: pd.DataFrame, metagenome_cols: List[str]) -> pd.Series:
    """计算每个宏基因组的平均表达量"""
    return df[metagenome_cols].mean()


def calculate_expression_difference(mean_early_smoking: pd.Series, 
                                   mean_early_nonsmoking: pd.Series, 
                                   mean_late_smoking: pd.Series, 
                                   mean_late_nonsmoking: pd.Series) -> pd.DataFrame:
    """计算四组之间的表达量差异"""
    # 计算所有可能的组间差异
    differences = {
        'early_smoking_vs_early_nonsmoking': abs(mean_early_smoking - mean_early_nonsmoking),
        'early_smoking_vs_late_smoking': abs(mean_early_smoking - mean_late_smoking),
        'early_smoking_vs_late_nonsmoking': abs(mean_early_smoking - mean_late_nonsmoking),
        'early_nonsmoking_vs_late_smoking': abs(mean_early_nonsmoking - mean_late_smoking),
        'early_nonsmoking_vs_late_nonsmoking': abs(mean_early_nonsmoking - mean_late_nonsmoking),
        'late_smoking_vs_late_nonsmoking': abs(mean_late_smoking - mean_late_nonsmoking)
    }
    
    # 创建差异DataFrame
    diff_df = pd.DataFrame(differences)
    
    # 计算每个宏基因组的最大差异
    diff_df['max_difference'] = diff_df.max(axis=1)
    
    # 按最大差异排序
    diff_df = diff_df.sort_values('max_difference', ascending=False)
    
    return diff_df


def select_top_metagenomes(diff_df: pd.DataFrame, n: int = 12) -> List[str]:
    """选择差异最大的前n个宏基因组"""
    return diff_df.head(n).index.tolist()


def calculate_correlation_matrix(df: pd.DataFrame, metagenome_cols: List[str]) -> pd.DataFrame:
    """计算宏基因组之间的相关性矩阵"""
    return df[metagenome_cols].corr()


def calculate_significant_correlations(corr_matrix: pd.DataFrame, n_samples: int, alpha: float = 0.005) -> pd.DataFrame:
    """计算显著的相关性（调整后的p值 < 0.005）"""
    import scipy.stats as stats
    
    # 创建空的p值矩阵
    p_matrix = pd.DataFrame(index=corr_matrix.index, columns=corr_matrix.columns)
    
    # 计算每个相关系数的p值
    for i, metagenome1 in enumerate(corr_matrix.index):
        for j, metagenome2 in enumerate(corr_matrix.columns):
            if i < j:
                correlation = corr_matrix.iloc[i, j]
                # 计算p值
                t_stat = correlation * np.sqrt((n_samples - 2) / (1 - correlation**2))
                p_value = 2 * (1 - stats.t.cdf(abs(t_stat), df=n_samples-2))
                p_matrix.iloc[i, j] = p_value
                p_matrix.iloc[j, i] = p_value
    
    # 对角线设为1（自相关）
    np.fill_diagonal(p_matrix.values, 1)
    
    # 多重检验校正（Bonferroni校正）
    num_tests = len(corr_matrix.index) * (len(corr_matrix.index) - 1) / 2
    adjusted_alpha = alpha / num_tests
    
    # 创建显著相关性矩阵
    significant_corr = corr_matrix.copy()
    significant_corr[p_matrix > adjusted_alpha] = 0
    
    return significant_corr


def create_correlation_network(corr_matrix: pd.DataFrame, n_samples: int, alpha: float = 0.005) -> nx.Graph:
    """创建相关性网络图（仅显示显著相关性）"""
    G = nx.Graph()
    
    # 计算显著相关性
    significant_corr = calculate_significant_correlations(corr_matrix, n_samples, alpha)
    
    # 添加节点
    for metagenome in corr_matrix.index:
        G.add_node(metagenome)
    
    # 添加边（仅添加显著相关的）
    for i, metagenome1 in enumerate(corr_matrix.index):
        for j, metagenome2 in enumerate(corr_matrix.columns):
            if i < j:  # 避免重复添加边
                correlation = significant_corr.iloc[i, j]
                if correlation != 0:  # 只添加显著相关的
                    G.add_edge(metagenome1, metagenome2, weight=correlation)
    
    return G


def plot_correlation_network(G: nx.Graph, 
                             group_name: str, 
                             title: str, 
                             group_data: pd.DataFrame, 
                             top_metagenomes: List[str],
                             figsize: Tuple[int, int] = (16, 14),
                             colors: Optional[Dict[str, str]] = None,
                             save_path: Optional[str] = None,
                             dpi: int = 300) -> plt.Figure:
    """
    绘制相关性网络图，模仿violin_plot的风格
    
    Parameters:
    -----------
    G : nx.Graph
        相关性网络
    group_name : str
        分组名称
    title : str
        图表标题
    group_data : pd.DataFrame
        分组数据，用于计算相对丰度
    top_metagenomes : List[str]
        前12个差异最大的宏基因组
    figsize : Tuple[int, int]
        图表尺寸
    colors : Optional[Dict[str, str]]
        颜色映射
    save_path : Optional[str]
        保存路径
    dpi : int
        保存图片的DPI
    
    Returns:
    --------
    plt.Figure
        生成的图表对象
    """
    # 默认颜色映射（与violin_plot保持一致）
    if colors is None:
        colors = {
            'early_smoking': '#FF6B6B',      # 红色 - 早期吸烟
            'early_nonsmoking': '#4ECDC4',   # 青色 - 早期非吸烟
            'late_smoking': '#FFD166',       # 黄色 - 晚期吸烟
            'late_nonsmoking': '#06D6A0'     # 绿色 - 晚期非吸烟
        }
    
    # 美化样式并创建图形（去掉背景网格）
    sns.set(style='white', context='notebook', font_scale=1.2)
    fig, ax = plt.subplots(figsize=figsize)
    
    # 使用圆形布局，形成12边形
    pos = nx.circular_layout(G)
    
    # 计算log10相对丰度作为节点大小
    node_sizes = []
    log_abundances = []
    for metagenome in G.nodes():
        if metagenome in group_data.columns:
            # 计算平均丰度
            mean_abundance = group_data[metagenome].mean()
            # 计算log10相对丰度（加1e-10避免log(0)）
            log_abundance = np.log10(mean_abundance + 1e-10)
            log_abundances.append(log_abundance)
            # 缩放节点大小
            node_size = max(300, min(1800, (log_abundance + 10) * 120))
            node_sizes.append(node_size)
        else:
            log_abundances.append(0)
            node_sizes.append(600)
    
    # 绘制节点（所有组使用相同颜色）
    node_color = '#6C757D'  # 灰色 - 所有组使用相同颜色
    nx.draw_networkx_nodes(G, pos, node_size=node_sizes, node_color=node_color, alpha=0.8, ax=ax, edgecolors='black', linewidths=1.5)
    
    # 绘制边（根据相关性值设置颜色和宽度）
    edges = G.edges(data=True)
    edge_weights = [d['weight'] for u, v, d in edges]
    edge_colors = []
    edge_widths = []
    for u, v, d in edges:
        weight = d['weight']
        if weight > 0:
            edge_colors.append('#DC3545')  # 正相关 - 红色
        else:
            edge_colors.append('#007BFF')  # 负相关 - 蓝色
        edge_widths.append(8)
    
    nx.draw_networkx_edges(G, pos, edgelist=edges, edge_color=edge_colors, width=edge_widths, alpha=0.8, ax=ax)
    
    # 绘制节点标签（只显示简化的名称，显示在节点下方）
    node_labels = {}
    label_pos = {}
    for node in G.nodes():
        # 简化节点标签，只显示物种名称
        if 'tax_' in node:
            # 提取物种名称（去掉'tax_'前缀）
            species_name = node.replace('tax_', '')
            # 全部显示，不截断
            node_labels[node] = species_name
        else:
            node_labels[node] = node
        # 调整标签位置到节点下方
        x, y = pos[node]
        label_pos[node] = (x, y - 0.1)  # 向下偏移0.1
    
    nx.draw_networkx_labels(G, label_pos, labels=node_labels, font_size=10, font_weight='bold', ax=ax, verticalalignment='top')
    
    # 添加相关性颜色条
    from matplotlib import cm
    from matplotlib.colors import LinearSegmentedColormap
    
    # 创建自定义颜色映射
    colors = ['#007BFF', '#FFFFFF', '#DC3545']  # 蓝色 -> 白色 -> 红色
    cmap = LinearSegmentedColormap.from_list('correlation', colors, N=256)
    
    # 创建颜色条
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=-1, vmax=1))
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax, orientation='vertical', pad=0.02, aspect=40)
    cbar.set_label('Correlation Coefficient', fontsize=12, fontweight='bold')
    cbar.ax.tick_params(labelsize=10)
    
    # 添加节点大小图例
    from matplotlib.lines import Line2D
    
    # 计算节点大小示例
    min_size = min(node_sizes)
    max_size = max(node_sizes)
    mid_size = (min_size + max_size) / 2
    
    # 计算对应的log10丰度
    min_log = min(log_abundances)
    max_log = max(log_abundances)
    mid_log = (min_log + max_log) / 2
    
    # 创建节点大小图例
    node_legend_elements = [
        Line2D([0], [0], marker='o', color='w', markerfacecolor=node_color, markersize=np.sqrt(min_size)/2, label=f'Low Abundance\n(log10: {min_log:.2f})'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor=node_color, markersize=np.sqrt(mid_size)/2, label=f'Medium Abundance\n(log10: {mid_log:.2f})'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor=node_color, markersize=np.sqrt(max_size)/2, label=f'High Abundance\n(log10: {max_log:.2f})')
    ]
    
    # 添加图例
    node_legend = ax.legend(handles=node_legend_elements, loc='upper left', fontsize=10, title='Node Size (Abundance)')
    ax.add_artist(node_legend)
    
    # 添加网络统计信息
    num_nodes = len(G.nodes())
    num_edges = len(G.edges())
    avg_degree = 2 * num_edges / num_nodes if num_nodes > 0 else 0
    
    stats_text = f'Nodes: {num_nodes}\nEdges: {num_edges}\nAverage Degree: {avg_degree:.2f}\nSignificance: p < 0.005 (adjusted)'
    ax.text(0.02, 0.98, stats_text, transform=ax.transAxes,
            fontsize=12, fontweight='bold',
            verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8, edgecolor='gray', linewidth=1))
    
    # 设置图表标题
    ax.set_title(title, fontsize=20, fontweight='bold', pad=30)
    
    # 去除图片外围方框
    ax.axis('off')
    
    # 调整布局，确保所有内容都能显示
    plt.tight_layout()
    plt.subplots_adjust(left=0.05, right=0.95, top=0.9, bottom=0.1)
    
    # 保存图片
    if save_path:
        plt.savefig(save_path, dpi=dpi, bbox_inches='tight')
        print(f"图片已保存至: {save_path}")
    
    return fig

def analyze_specific_metagenome_correlation(df: pd.DataFrame, 
                                   specific_metagenomes: List[str],
                                   output_dir: Optional[str] = None,
                                   alpha: float = 0.005) -> Dict[str, Any]:
    """
    特定宏基因组相关性分析流水线
    
    Parameters:
    -----------
    df : pd.DataFrame
        原始数据
    specific_metagenomes : List[str]
        要分析的特定宏基因组列表
    output_dir : Optional[str]
        输出目录
    alpha : float
        显著性水平阈值
    
    Returns:
    --------
    Dict[str, Any]
        分析结果
    """
    results = {}
    
    # 1. 准备分组（早期与晚期）
    print("步骤1: 准备CRC早期与晚期吸烟与非吸烟分组...")
    df_grouped = prepare_early_late_groups(df)
    results['group_counts'] = df_grouped['group'].value_counts().to_dict()
    print(f"分组统计: {results['group_counts']}")
    
    # 2. 过滤出存在的特定宏基因组
    print("步骤2: 过滤特定宏基因组...")
    available_metagenomes = [mg for mg in specific_metagenomes if mg in df_grouped.columns]
    results['specific_metagenomes'] = available_metagenomes
    print(f"可用的特定宏基因组数量: {len(available_metagenomes)}")
    for mg in available_metagenomes:
        print(f"  - {mg}")
    
    if len(available_metagenomes) < 2:
        print("错误: 可用的特定宏基因组数量不足，无法进行相关性分析")
        return results
    
    # 3. 计算每组的平均表达量
    print("步骤3: 计算每组的平均表达量...")
    group_order = ['CRC_early_smoking', 'CRC_early_nonsmoking', 'CRC_late_smoking', 'CRC_late_nonsmoking']
    mean_expressions = {}
    
    for group in group_order:
        if group in df_grouped['group'].values:
            group_data = df_grouped[df_grouped['group'] == group]
            mean_expressions[group] = calculate_mean_expression(group_data, available_metagenomes)
    
    # 4. 计算表达量差异
    print("步骤4: 计算特定宏基因组的表达量差异...")
    if all(group in mean_expressions for group in group_order):
        diff_df = calculate_expression_difference(
            mean_expressions['CRC_early_smoking'],
            mean_expressions['CRC_early_nonsmoking'],
            mean_expressions['CRC_late_smoking'],
            mean_expressions['CRC_late_nonsmoking']
        )
        
        results['expression_differences'] = diff_df
        print("特定宏基因组表达量差异计算完成")
    else:
        print("警告: 某些分组数据缺失，无法计算表达量差异")
    
    # 5. 为每组创建相关性网络图
    print("步骤5: 为每组创建特定宏基因组相关性网络图...")
    colors = {
        'early_smoking': '#FF6B6B',      # 红色 - 早期吸烟
        'early_nonsmoking': '#4ECDC4',   # 青色 - 早期非吸烟
        'late_smoking': '#FFD166',       # 黄色 - 晚期吸烟
        'late_nonsmoking': '#06D6A0'     # 绿色 - 晚期非吸烟
    }
    
    network_results = {}
    
    for group in group_order:
        if group in df_grouped['group'].values:
            print(f"\n处理分组: {group}")
            
            # 获取该组数据
            group_data = df_grouped[df_grouped['group'] == group]
            
            # 计算相关性矩阵
            corr_matrix = calculate_correlation_matrix(group_data, available_metagenomes)
            
            # 创建相关性网络（仅显示显著相关性）
            n_samples = len(group_data)
            G = create_correlation_network(corr_matrix, n_samples=n_samples, alpha=alpha)
            
            # 生成图表标题
            title_map = {
                'CRC_early_smoking': 'Early Stage Smoking',
                'CRC_early_nonsmoking': 'Early Stage Non-Smoking',
                'CRC_late_smoking': 'Late Stage Smoking',
                'CRC_late_nonsmoking': 'Late Stage Non-Smoking'
            }
            title = f'Specific Metagenome Correlation Network: {title_map.get(group, group)}'
            
            # 保存路径
            if output_dir:
                save_path = f"{output_dir}/specific_metagenome_correlation_{group}.png"
            else:
                save_path = f"specific_metagenome_correlation_{group}.png"
            
            # 绘制网络图
            fig = plot_correlation_network(
                G, 
                group, 
                title, 
                group_data=group_data,
                top_metagenomes=available_metagenomes,
                save_path=save_path,
                colors=colors
            )
            
            network_results[group] = {
                'figure': fig,
                'network': G,
                'correlation_matrix': corr_matrix,
                'save_path': save_path,
                'n_samples': n_samples,
                'num_edges': len(G.edges()),
                'num_nodes': len(G.nodes())
            }
            
            print(f"分组 {group}: {len(G.nodes())} 个节点, {len(G.edges())} 条边")
        else:
            print(f"警告: 分组 {group} 数据缺失")
    
    results['network_results'] = network_results
    
    print("\n特定宏基因组相关性分析完成!")
    return results
    


def analyze_metagenome_correlation(df: pd.DataFrame, 
                                   output_dir: Optional[str] = None, 
                                   top_n: int = 16) -> Dict[str, Any]:
    """
    宏基因组相关性网络分析流水线
    
    Parameters:
    -----------
    df : pd.DataFrame
        原始数据
    output_dir : Optional[str]
        输出目录
    top_n : int
        选择差异最大的前n个宏基因组
    
    Returns:
    --------
    Dict[str, Any]
        分析结果
    """
    results = {}
    
    # 1. 准备分组（早期与晚期）
    print("步骤1: 准备CRC早期与晚期吸烟与非吸烟分组...")
    df_grouped = prepare_early_late_groups(df)
    results['group_counts'] = df_grouped['group'].value_counts().to_dict()
    print(f"分组统计: {results['group_counts']}")
    
    # 2. 获取宏基因组列
    print("步骤2: 获取宏基因组数据...")
    metagenome_cols = get_metagenome_columns(df_grouped)
    print(f"宏基因组数量: {len(metagenome_cols)}")
    
    # 3. 计算每组的平均表达量
    print("步骤3: 计算每组的平均表达量...")
    group_order = ['CRC_early_smoking', 'CRC_early_nonsmoking', 'CRC_late_smoking', 'CRC_late_nonsmoking']
    mean_expressions = {}
    
    for group in group_order:
        if group in df_grouped['group'].values:
            group_data = df_grouped[df_grouped['group'] == group]
            mean_expressions[group] = calculate_mean_expression(group_data, metagenome_cols)
    
    # 4. 计算表达量差异并选择前16个差异最大的宏基因组
    print("步骤4: 计算表达量差异并选择前16个差异最大的宏基因组...")
    if all(group in mean_expressions for group in group_order):
        diff_df = calculate_expression_difference(
            mean_expressions['CRC_early_smoking'],
            mean_expressions['CRC_early_nonsmoking'],
            mean_expressions['CRC_late_smoking'],
            mean_expressions['CRC_late_nonsmoking']
        )
        
        top_metagenomes = select_top_metagenomes(diff_df, top_n)
        print(f"选择的前{top_n}个差异最大的宏基因组:")
        for mg in top_metagenomes:
            print(f"  - {mg}")
        
        results['top_metagenomes'] = top_metagenomes
        results['expression_differences'] = diff_df
    else:
        print("警告: 某些分组数据缺失，无法计算表达量差异")
        # 如果数据缺失，使用所有宏基因组
        top_metagenomes = metagenome_cols[:top_n]
        results['top_metagenomes'] = top_metagenomes
    
    # 5. 为每组创建相关性网络图
    print("步骤5: 为每组创建相关性网络图...")
    colors = {
        'early_smoking': '#FF6B6B',      # 红色 - 早期吸烟
        'early_nonsmoking': '#4ECDC4',   # 青色 - 早期非吸烟
        'late_smoking': '#FFD166',       # 黄色 - 晚期吸烟
        'late_nonsmoking': '#06D6A0'     # 绿色 - 晚期非吸烟
    }
    
    network_results = {}
    
    for group in group_order:
        if group in df_grouped['group'].values:
            print(f"\n处理分组: {group}")
            
            # 获取该组数据
            group_data = df_grouped[df_grouped['group'] == group]
            
            # 计算相关性矩阵
            corr_matrix = calculate_correlation_matrix(group_data, top_metagenomes)
            
            # 创建相关性网络（仅显示显著相关性）
            n_samples = len(group_data)
            G = create_correlation_network(corr_matrix, n_samples=n_samples, alpha=0.8)
            
            # 生成图表标题
            title_map = {
                'early_smoking': 'Early Stage Smoking',
                'early_nonsmoking': 'Early Stage Non-Smoking',
                'late_smoking': 'Late Stage Smoking',
                'late_nonsmoking': 'Late Stage Non-Smoking'
            }
            title = f'Metagenome Correlation Network: {title_map.get(group, group)}'
            
            # 保存路径
            if output_dir:
                save_path = f"{output_dir}/metagenome_correlation_{group}.png"
            else:
                save_path = f"metagenome_correlation_{group}.png"
            
            # 绘制网络图
            fig = plot_correlation_network(
                G, 
                group, 
                title, 
                group_data=group_data,
                top_metagenomes=top_metagenomes,
                save_path=save_path,
                colors=colors
            )
            
            network_results[group] = {
                'figure': fig,
                'network': G,
                'correlation_matrix': corr_matrix,
                'save_path': save_path
            }
        else:
            print(f"警告: 分组 {group} 数据缺失")
    
    results['network_results'] = network_results
    
    # 6. 保存差异分析结果
    if output_dir and 'expression_differences' in results:
        output_file = f"{output_dir}/metagenome_expression_differences.csv"
        results['expression_differences'].to_csv(output_file)
        results['differences_file'] = output_file
        print(f"表达量差异结果已保存至: {output_file}")
    
    print("\n宏基因组相关性网络分析完成!")
    return results


def main():
    """主函数"""
    print("=" * 80)
    print("CRC早期与晚期吸烟与非吸烟的宏基因组相关性网络分析")
    print("分组: 早期吸烟、早期非吸烟、晚期吸烟、晚期非吸烟")
    print("分析: 每组内宏基因组之间的相关性网络图")
    print("选择: 表达量差值最大的前16个宏基因组")
    print("=" * 80)
    
    # 加载数据
    ds = BioSmokeDataset()
    df = ds.merge_to_dataframe()
    
    print(f"总样本数: {df.shape[0]}")
    print(f"总特征数: {df.shape[1]}")
    
    # 检查TNM分期数据
    print(f"\nTNM分期分布:")
    print(df['tnm_stage'].value_counts())
    
    # 检查吸烟状态
    print(f"\n吸烟状态分布:")
    print(df['smoking_label'].value_counts())

    # 创建输出目录
    output_dir = './results/network_plots'
    os.makedirs(output_dir, exist_ok=True)
    
    # 执行分析
    results = analyze_metagenome_correlation(
        df,
        output_dir=output_dir,
        top_n=12
    )

    specific_metagenomes = [
            # 有益菌种
            'tax_s__Faecalibacterium_prausnitzii',  # 抗炎菌，CRC中减少
            'tax_s__Bacteroides_fragilis',          # 机会致病菌
            'tax_s__Bacteroides_vulgatus',          # 常见肠道菌
            'tax_s__Bifidobacterium_bifidum',       # 益生菌
            'tax_s__Akkermansia_muciniphila',       # 黏液降解菌
            
            # 致病相关菌种
            'tax_s__Escherichia_coli',              # 机会致病菌
            'tax_s__Fusobacterium_nucleatum',       # CRC相关菌
            'tax_s__Peptostreptococcus_stomatis',   # CRC相关菌
            'tax_s__Parvimonas_micra',              # CRC相关菌
            'tax_s__Streptococcus_gallolyticus',    # CRC相关菌
            
            # 其他重要菌种
            'tax_s__Prevotella_copri',              # 代谢相关
            'tax_s__Roseburia_intestinalis',        # 丁酸产生菌
            'tax_s__Ruminococcus_bromii',           # 淀粉降解菌
            'tax_s__Alistipes_putredinis',          # 肠道菌
            'tax_s__Coprococcus_comes',             # 丁酸产生菌
            
            # 从差异分析中发现的重要菌种
            'tax_s__Bacteroides_cellulosilyticus',
            'tax_s__Bacteroides_finegoldii',
            'tax_s__Bacteroides_plebeius',
            'tax_s__Megamonas_funiformis',
            'tax_s__Parabacteroides_merdae',
            'tax_s__Ruminococcaceae_bacterium',
            'tax_s__Clostridium_sp._AM42-36',
            'tax_s__Bacteroides_stercoris',
            'tax_s__Coprobacter_fastidiosus',
            'tax_s__Oscillibacter_sp.',
            'tax_s__Ruminococcus_callidus',
            'tax_s__Ruminococcaceae_bacterium_AF10-16',
            'tax_s__Shigella_sonnei',
            'tax_s__Shigella_flexneri',
            'tax_s__Prevotellamassilia_timonensis',
            'tax_s__Oscillibacter_sp._ER4',
            'tax_s__Clostridiales_bacterium_42_27',
            'tax_s__Ruminococcus_lactaris',
            'tax_s__uncultured_Faecalibacterium_sp.',
            'tax_s__Prevotella_sp._Marseille-P4119'
        ]

    analyze_specific_metagenome_correlation(
        df,
        specific_metagenomes,
        output_dir=output_dir
    )
    
    print("\n分析结果汇总:")
    print(f"分组统计: {results['group_counts']}")
    print(f"选择的宏基因组数量: {len(results['top_metagenomes'])}")
    print(f"生成的网络图数量: {len(results['network_results'])}")


if __name__ == '__main__':
    main()
