import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from typing import Optional, Tuple, Dict, Any
from dataset import BioSmokeDataset
import os
from split_group import prepare_crc_diff_groups

def calculate_alpha_diversity(df: pd.DataFrame, 
                            bacteria_prefix: str = 'tax_s__',
                            method: str = 'shannon') -> pd.DataFrame:
    """
    计算Alpha多样性指数
    
    Parameters:
    -----------
    df : pd.DataFrame
        包含细菌数据的DataFrame
    bacteria_prefix : str
        细菌列的前缀
    method : str
        多样性指数计算方法：'shannon', 'simpson', 'richness'
    
    Returns:
    --------
    pd.DataFrame
        包含Alpha多样性指数的DataFrame
    """
    # 识别细菌数据列
    bacteria_columns = [col for col in df.columns if col.startswith(bacteria_prefix)]
    
    if method == 'shannon':
        def calculate_shannon(row):
            abundances = row[bacteria_columns].values.astype(float)
            abundances = abundances[abundances > 0]
            
            if len(abundances) == 0:
                return 0
            
            total = np.sum(abundances)
            relative_abundances = abundances / total
            shannon_index = -np.sum(relative_abundances * np.log(relative_abundances))
            return shannon_index
        
        df['alpha_diversity'] = df.apply(calculate_shannon, axis=1)
        df['diversity_type'] = 'Shannon Index'
        
    elif method == 'simpson':
        def calculate_simpson(row):
            abundances = row[bacteria_columns].values.astype(float)
            abundances = abundances[abundances > 0]
            
            if len(abundances) == 0:
                return 0
            
            total = np.sum(abundances)
            relative_abundances = abundances / total
            simpson_index = 1 - np.sum(relative_abundances ** 2)
            return simpson_index
        
        df['alpha_diversity'] = df.apply(calculate_simpson, axis=1)
        df['diversity_type'] = 'Simpson Index'
        
    elif method == 'richness':
        def calculate_richness(row):
            abundances = row[bacteria_columns].values.astype(float)
            # 计算观测到的物种数（非零值）
            richness = np.sum(abundances > 0)
            return richness
        
        df['alpha_diversity'] = df.apply(calculate_richness, axis=1)
        df['diversity_type'] = 'Species Richness'
    
    else:
        raise ValueError(f"不支持的多样性指数计算方法: {method}。请选择 'shannon', 'simpson' 或 'richness'")
    
    return df


def calculate_early_late_pvalues(df: pd.DataFrame, group_col: str = 'group', value_col: str = 'alpha_diversity') -> Dict[str, float]:
    """
    计算CTRL, CRC_welldiff, CRC_poordiff三组之间的p值
    
    Parameters:
    -----------
    df : pd.DataFrame
        包含分组和Alpha多样性数据的DataFrame
    group_col : str
        分组列名
    value_col : str
        多样性值列名
    
    Returns:
    --------
    Dict[str, float]
        p值结果字典
    """
    pvalues = {}
    
    groups = ['CTRL', 'CRC_welldiff', 'CRC_poordiff']
    
    group_data = {}
    for group in groups:
        if group in df[group_col].values:
            group_data[group] = df[df[group_col] == group][value_col].dropna().values
        else:
            group_data[group] = np.array([])
    
    print("\n计算CTRL, CRC_welldiff, CRC_poordiff三组之间的p值...")
    
    print("\n1. CRC_welldiff vs CTRL:")
    if len(group_data['CRC_welldiff']) > 0 and len(group_data['CTRL']) > 0:
        stat, p = stats.mannwhitneyu(group_data['CRC_welldiff'], group_data['CTRL'], alternative='two-sided')
        pvalues['CRC_welldiff_vs_CTRL'] = p
        print(f"  CRC_welldiff vs CTRL: p = {p:.6f}")
    
    print("\n2. CRC_poordiff vs CTRL:")
    if len(group_data['CRC_poordiff']) > 0 and len(group_data['CTRL']) > 0:
        stat, p = stats.mannwhitneyu(group_data['CRC_poordiff'], group_data['CTRL'], alternative='two-sided')
        pvalues['CRC_poordiff_vs_CTRL'] = p
        print(f"  CRC_poordiff vs CTRL: p = {p:.6f}")
    
    print("\n3. CRC_poordiff vs CRC_welldiff:")
    if len(group_data['CRC_poordiff']) > 0 and len(group_data['CRC_welldiff']) > 0:
        stat, p = stats.mannwhitneyu(group_data['CRC_poordiff'], group_data['CRC_welldiff'], alternative='two-sided')
        pvalues['CRC_poordiff_vs_welldiff'] = p
        print(f"  CRC_poordiff vs CRC_welldiff: p = {p:.6f}")
    
    return pvalues


def plot_early_late_alpha_diversity_violin(df: pd.DataFrame, 
                                          group_col: str = 'group',
                                          value_col: str = 'alpha_diversity',
                                          title: str = 'Alpha Diversity',
                                          figsize: Tuple[int, int] = (14, 10),
                                          colors: Optional[Dict[str, str]] = None,
                                          save_path: Optional[str] = None,
                                          dpi: int = 300) -> plt.Figure:
    """
    绘制CRC分化程度的Alpha多样性小提琴图
    
    Parameters:
    -----------
    df : pd.DataFrame
        包含分组和Alpha多样性数据的DataFrame
    group_col : str
        分组列名
    value_col : str
        多样性值列名
    title : str
        图表标题
    figsize : Tuple[int, int]
        图表尺寸
    colors : Optional[Dict[str, str]]
        自定义颜色映射
    save_path : Optional[str]
        保存路径
    dpi : int
        保存图片的DPI
    
    Returns:
    --------
    plt.Figure
        生成的图表对象
    """
    # 定义3个分组的顺序
    group_order = ['CTRL', 'CRC_welldiff', 'CRC_poordiff']

    plot_group_name = {
        'CTRL': 'CTRL',
        'CRC_welldiff': 'CRC-Well',
        'CRC_poordiff': 'CRC-Poor'
    }
    
    # 默认颜色映射
    if colors is None:
        colors = {
            'CRC_poordiff': '#D7263D',
            'CRC_welldiff': '#F18F01',
            'CTRL': '#2E86AB'
        }
    
    # 筛选出3个分组的数据
    plot_data = df[df[group_col].isin(group_order)].copy()
    
    if plot_data.empty:
        raise ValueError("没有找到有效的分组数据，请检查数据中的分组标签")
    
    # 计算p值
    pvalues = calculate_early_late_pvalues(plot_data, group_col, value_col)
    
    # 美化样式并创建图形（去掉背景网格）
    sns.set(style='white', context='notebook', font_scale=1.1)
    fig, ax = plt.subplots(figsize=figsize)

    # 小提琴图：去掉内部线条并移除边缘线，按宽度缩放
    sns.violinplot(x=group_col, y=value_col, data=plot_data,
                   order=group_order, palette=colors, inner=None,
                   cut=0, scale='width', linewidth=0, ax=ax)

    # 叠加箱线图，细窄以便显示在小提琴中央
    sns.boxplot(x=group_col, y=value_col, data=plot_data,
                order=group_order, width=0.14,showcaps=True, showfliers=False,
                boxprops=dict(facecolor='none', edgecolor='black', linewidth=1.2),
                whiskerprops=dict(color='black', linewidth=1.2),
                capprops=dict(color='black', linewidth=1.2),
                medianprops=dict(color='black', linewidth=1.6),
                zorder=2, ax=ax)

    # 添加散点分布（抖动），放在小提琴和箱线之上
    sns.stripplot(x=group_col, y=value_col, data=plot_data,
                  order=group_order, color='gray', alpha=0.4,
                  size=4, jitter=True, zorder=3, ax=ax)

    # 设置图表标题和标签
    ax.set_title(title, fontsize=18, fontweight='bold', pad=25)
    ax.set_xlabel('Group', fontsize=14)
    ax.set_ylabel('Alpha Diversity Index (Shannon)', fontsize=14)

    # 获取y轴范围
    y_min, y_max = ax.get_ylim()
    y_range = y_max - y_min

    # 添加每个组的样本数量标签
    for i, group in enumerate(group_order):
        if group in plot_data[group_col].values:
            group_data = plot_data[plot_data[group_col] == group]
            count = len(group_data)
            mean_value = group_data[value_col].mean()
            median_value = group_data[value_col].median()

            # 在图表顶部添加样本数量
            ax.text(i, y_max - y_range * 0.03, f'n={count}',
                    ha='center', va='bottom', fontsize=10, fontweight='bold')

            # 在图表底部添加统计信息
            ax.text(i, y_min + y_range * 0.06,
                    f'Mean: {mean_value:.2f}\nMedian: {median_value:.2f}',
                    ha='center', va='top', fontsize=9,
                    bbox=dict(facecolor='white', alpha=0.7, edgecolor='none'))

    same_stage_y = y_max + y_range * 0.1
    
    if 'CRC_poordiff_vs_welldiff' in pvalues:
        p = pvalues['CRC_poordiff_vs_welldiff']
        x_pos_mid = 1.5
        p_text = f'p = {p:.3f}' if p >= 0.001 else 'p < 0.001'
        print(f"  CRC_poordiff vs CRC_welldiff: p = {p:.6f}")
        ax.text(x_pos_mid, same_stage_y+0.0, p_text,
            ha='center', va='bottom', fontsize=10, fontweight='bold')
        y_line = same_stage_y - y_range * 0.04
        bracket_h = y_range * 0.02
        ax.plot([1, 1], [y_line, y_line + bracket_h], color='black', linewidth=1.5)
        ax.plot([1, 2], [y_line + bracket_h, y_line + bracket_h], color='black', linewidth=1.5)
        ax.plot([2, 2], [y_line, y_line + bracket_h], color='black', linewidth=1.5)
    
    if 'CRC_poordiff_vs_CTRL' in pvalues:
        p = pvalues['CRC_poordiff_vs_CTRL']
        x_pos_mid = 1.0
        p_text = f'p = {p:.3f}' if p >= 0.001 else 'p < 0.001'
        print(f"  CRC_poordiff vs CTRL: p = {p:.6f}")
        ax.text(x_pos_mid, same_stage_y + y_range * 0.1, p_text,
            ha='center', va='bottom', fontsize=10, fontweight='bold')
        y_line = same_stage_y + y_range * 0.06
        bracket_h = y_range * 0.02
        ax.plot([0, 0], [y_line, y_line + bracket_h], color='black', linewidth=1.5)
        ax.plot([0, 2], [y_line + bracket_h, y_line + bracket_h], color='black', linewidth=1.5)
        ax.plot([2, 2], [y_line, y_line + bracket_h], color='black', linewidth=1.5)
    
    if 'CRC_welldiff_vs_CTRL' in pvalues:
        p = pvalues['CRC_welldiff_vs_CTRL']
        x_pos_mid = 0.5
        p_text = f'p = {p:.3f}' if p >= 0.001 else 'p < 0.001'
        print(f"  CRC_welldiff vs CTRL: p = {p:.6f}")
        ax.text(x_pos_mid, same_stage_y+0.0, p_text,
            ha='center', va='bottom', fontsize=10, fontweight='bold')
        y_line = same_stage_y - y_range * 0.04
        bracket_h = y_range * 0.02
        ax.plot([0, 0], [y_line, y_line + bracket_h], color='black', linewidth=1.5)
        ax.plot([0, 1], [y_line + bracket_h, y_line + bracket_h], color='black', linewidth=1.5)
        ax.plot([1, 1], [y_line, y_line + bracket_h], color='black', linewidth=1.5)

    # 调整y轴范围以容纳所有p值标注
    ax.set_ylim(y_min, y_max + y_range * 0.25)
    
    # 旋转x轴标签
    plt.xticks(ticks=range(len(group_order)), labels=[plot_group_name[g] for g in group_order], rotation=45, ha='right', fontsize=12)
    plt.yticks(fontsize=12)
    plt.tight_layout()
    
    # 保存图片
    if save_path:
        plt.savefig(save_path, dpi=dpi, bbox_inches='tight')
        print(f"图片已保存至: {save_path}")
    
    return fig, pvalues


def analyze_early_late_alpha_diversity(df: pd.DataFrame,
                                      bacteria_prefix: str = 'tax_s__',
                                      diversity_method: str = 'shannon',
                                      output_dir: Optional[str] = None) -> Dict[str, Any]:
    """
    CRC分化程度的Alpha多样性分析流水线
    
    Parameters:
    -----------
    df : pd.DataFrame
        原始数据
    bacteria_prefix : str
        细菌列前缀
    diversity_method : str
        多样性计算方法
    output_dir : Optional[str]
        输出目录
    
    Returns:
    --------
    Dict[str, Any]
        分析结果
    """
    results = {}
    
    # 1. 准备分组（分化程度）
    print("步骤1: 准备CRC分化程度分组...")
    df_grouped = prepare_crc_diff_groups(df)
    results['group_counts'] = df_grouped['group'].value_counts().to_dict()
    print(f"分组统计: {results['group_counts']}")
    
    # 2. 计算Alpha多样性
    print(f"步骤2: 计算{diversity_method}多样性指数...")
    df_diversity = calculate_alpha_diversity(df_grouped, bacteria_prefix, diversity_method)
    results['diversity_method'] = diversity_method
    
    # 3. 描述性统计
    print("步骤3: 计算描述性统计...")
    descriptive_stats = {}
    group_order = ['CTRL', 'CRC_welldiff', 'CRC_poordiff']
    
    for group in group_order:
        if group in df_diversity['group'].values:
            group_data = df_diversity[df_diversity['group'] == group]['alpha_diversity']
            descriptive_stats[group] = {
                'n': len(group_data),
                'mean': float(group_data.mean()),
                'std': float(group_data.std()),
                'median': float(group_data.median()),
                'min': float(group_data.min()),
                'max': float(group_data.max()),
                'q1': float(group_data.quantile(0.25)),
                'q3': float(group_data.quantile(0.75))
            }
    
    results['descriptive_stats'] = descriptive_stats
    
    # 4. 绘制小提琴图并计算p值
    print("步骤4: 绘制CRC分化程度小提琴图...")
    title_map = {
        'shannon': 'Shannon Index',
        'simpson': 'Simpson Index',
        'richness': 'Species Richness'
    }
    
    if output_dir:
        save_path = f"{output_dir}/alpha_diversity_{diversity_method}_violin.png"
    else:
        save_path = f"alpha_diversity_{diversity_method}_violin.png"
    
    fig, pvalues = plot_early_late_alpha_diversity_violin(
        df_diversity,
        save_path=save_path
    )
    
    results['figure'] = fig
    results['pvalues'] = pvalues
    results['save_path'] = save_path
    
    # 5. 保存数据
    if output_dir:
        output_file = f"{output_dir}/diff_alpha_diversity_results_{diversity_method}.csv"
        df_diversity[['group', 'alpha_diversity', 'diversity_type']].to_csv(output_file, index=True)
        results['data_file'] = output_file
        print(f"数据已保存至: {output_file}")        
    
    print("CRC分化程度分析完成!")
    return results


def main():
    """主函数"""
    print("=" * 80)
    print("CRC分化程度的Alpha多样性分析")
    print("分组: CTRL, CRC_welldiff (中-高分化), CRC_poordiff (低分化)")
    print("对比: CTRL vs CRC_welldiff, CTRL vs CRC_poordiff, CRC_welldiff vs CRC_poordiff")
    print("=" * 80)
    
    ds = BioSmokeDataset()
    ds.preprocess_taxonomy_data()
    ds.preprocess_metabolomics_data()
    df = ds.merge_to_dataframe()
    
    print(f"总样本数: {df.shape[0]}")
    print(f"总特征数: {df.shape[1]}")
    
    # 检查分化程度数据
    print(f"\n分化程度分布:")
    print(df['differentiation'].value_counts())

    output_dir = './results/violin_plots'
    os.makedirs(output_dir, exist_ok=True)
    
    # 执行分析
    results = analyze_early_late_alpha_diversity(
        df=df,
        bacteria_prefix='tax_s__',
        diversity_method='shannon',
        output_dir=output_dir
    )
    
    print("\n" + "=" * 80)
    print("分析完成!")
    print("=" * 80)


if __name__ == "__main__":
    main()
