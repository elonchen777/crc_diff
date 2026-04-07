import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from typing import Optional, Tuple, Dict, Any
from dataset import BioSmokeDataset
from split_group import prepare_crc_diff_groups
import os


GROUP_COLORS = {
    'CRC_poordiff': '#D7263D',
    'CRC_welldiff': '#F18F01',
    'CTRL': '#2E86AB'
}

GROUP_NAMES = {
    'CTRL': 'CTRL',
    'CRC_welldiff': 'CRC-Well',
    'CRC_poordiff': 'CRC-Poor'
}

GROUP_ORDER = ['CTRL', 'CRC_welldiff', 'CRC_poordiff']

QUADRANT_COLORS = {
    'FP_high_PS_low': '#2166AC',
    'FP_low_PS_high': '#B2182B',
    'FP_high_PS_high': '#67A9CF',
    'FP_low_PS_low': '#D6604D'
}

QUADRANT_NAMES = {
    'FP_high_PS_low': 'FP↑ PS↓',
    'FP_low_PS_high': 'FP↓ PS↑',
    'FP_high_PS_high': 'FP↑ PS↑',
    'FP_low_PS_low': 'FP↓ PS↓'
}


def calculate_fp_ps_ratio(df: pd.DataFrame, 
                          fp_col: str = 'tax_s__Faecalibacterium_prausnitzii',
                          ps_col: str = 'tax_s__Peptostreptococcus_stomatis',
                          epsilon: float = 1e-6) -> pd.DataFrame:
    """
    计算F. prausnitzii / P. stomatis比值
    
    Parameters:
    -----------
    df : pd.DataFrame
        包含物种丰度数据的DataFrame
    fp_col : str
        F. prausnitzii列名
    ps_col : str
        P. stomatis列名
    epsilon : float
        避免除零的小常数
    
    Returns:
    --------
    pd.DataFrame
        添加了FP丰度、PS丰度和FP/PS比值的DataFrame
    """
    df_result = df.copy()
    
    if fp_col not in df.columns:
        raise ValueError(f"未找到F. prausnitzii列: {fp_col}")
    if ps_col not in df.columns:
        raise ValueError(f"未找到P. stomatis列: {ps_col}")
    
    df_result['FP_abundance'] = pd.to_numeric(df_result[fp_col], errors='coerce').fillna(0)
    df_result['PS_abundance'] = pd.to_numeric(df_result[ps_col], errors='coerce').fillna(0)
    
    df_result['FP_PS_ratio'] = df_result['FP_abundance'] / (df_result['PS_abundance'] + epsilon)
    
    df_result['FP_PS_log_ratio'] = np.log10(df_result['FP_PS_ratio'] + epsilon)
    
    print(f"\nFP/PS比值计算完成:")
    print(f"  FP丰度范围: {df_result['FP_abundance'].min():.6f} - {df_result['FP_abundance'].max():.6f}")
    print(f"  PS丰度范围: {df_result['PS_abundance'].min():.6f} - {df_result['PS_abundance'].max():.6f}")
    print(f"  FP/PS比值范围: {df_result['FP_PS_ratio'].min():.6f} - {df_result['FP_PS_ratio'].max():.6f}")
    
    return df_result


def perform_quadrant_analysis(df: pd.DataFrame,
                               fp_col: str = 'FP_abundance',
                               ps_col: str = 'PS_abundance',
                               method: str = 'median') -> pd.DataFrame:
    """
    进行四象限分析
    
    Parameters:
    -----------
    df : pd.DataFrame
        包含FP和PS丰度的DataFrame
    fp_col : str
        FP丰度列名
    ps_col : str
        PS丰度列名
    method : str
        分组方法：'median' 或 'tertile'
    
    Returns:
    --------
    pd.DataFrame
        添加了四象限分类的DataFrame
    """
    df_result = df.copy()
    
    if method == 'median':
        fp_threshold = df_result[fp_col].median()
        ps_threshold = df_result[ps_col].median()
        print(f"\n使用中位数分组:")
        print(f"  FP阈值: {fp_threshold:.6f}")
        print(f"  PS阈值: {ps_threshold:.6f}")
        
        df_result['FP_group'] = (df_result[fp_col] > fp_threshold).astype(int)
        df_result['PS_group'] = (df_result[ps_col] > ps_threshold).astype(int)
        
    elif method == 'tertile':
        fp_q1, fp_q2 = df_result[fp_col].quantile([0.33, 0.67])
        ps_q1, ps_q2 = df_result[ps_col].quantile([0.33, 0.67])
        print(f"\n使用三分位数分组:")
        print(f"  FP阈值: {fp_q1:.6f}, {fp_q2:.6f}")
        print(f"  PS阈值: {ps_q1:.6f}, {ps_q2:.6f}")
        
        df_result['FP_group'] = pd.cut(df_result[fp_col], 
                                        bins=[-np.inf, fp_q1, fp_q2, np.inf],
                                        labels=[0, 1, 2])
        df_result['PS_group'] = pd.cut(df_result[ps_col], 
                                        bins=[-np.inf, ps_q1, ps_q2, np.inf],
                                        labels=[0, 1, 2])
    else:
        raise ValueError(f"不支持的分组方法: {method}")
    
    conditions = [
        (df_result['FP_group'] == 1) & (df_result['PS_group'] == 0),
        (df_result['FP_group'] == 0) & (df_result['PS_group'] == 1),
        (df_result['FP_group'] == 1) & (df_result['PS_group'] == 1),
        (df_result['FP_group'] == 0) & (df_result['PS_group'] == 0)
    ]
    
    choices = ['FP_high_PS_low', 'FP_low_PS_high', 'FP_high_PS_high', 'FP_low_PS_low']
    df_result['quadrant'] = np.select(conditions, choices, default='unknown')
    
    quadrant_counts = df_result['quadrant'].value_counts()
    print(f"\n四象限分布:")
    for q in choices:
        if q in quadrant_counts.index:
            print(f"  {QUADRANT_NAMES[q]}: {quadrant_counts[q]} ({quadrant_counts[q]/len(df_result)*100:.1f}%)")
    
    return df_result


def plot_fp_ps_scatter_quadrant(df: pd.DataFrame,
                                  figsize: Tuple[int, int] = (10, 8),
                                  save_path: Optional[str] = None,
                                  dpi: int = 300) -> plt.Figure:
    """
    绘制FP vs PS散点图（四象限分析）
    
    Parameters:
    -----------
    df : pd.DataFrame
        包含FP、PS丰度和四象限分类的DataFrame
    figsize : Tuple[int, int]
        图表尺寸
    save_path : Optional[str]
        保存路径
    dpi : int
        DPI
    
    Returns:
    --------
    plt.Figure
        图表对象
    """
    sns.set(style='white', context='notebook', font_scale=1.1)
    fig, ax = plt.subplots(figsize=figsize)
    
    fp_median = df['FP_abundance'].median()
    ps_median = df['PS_abundance'].median()
    
    for quadrant in ['FP_high_PS_low', 'FP_low_PS_high', 'FP_high_PS_high', 'FP_low_PS_low']:
        subset = df[df['quadrant'] == quadrant]
        if len(subset) > 0:
            ax.scatter(subset['FP_abundance'], subset['PS_abundance'],
                      c=QUADRANT_COLORS[quadrant], 
                      label=f"{QUADRANT_NAMES[quadrant]} (n={len(subset)})",
                      alpha=0.6, s=50, edgecolors='white', linewidth=0.5)
    
    ax.axvline(x=fp_median, color='gray', linestyle='--', linewidth=1.5, alpha=0.7)
    ax.axhline(y=ps_median, color='gray', linestyle='--', linewidth=1.5, alpha=0.7)
    
    ax.set_xlabel('F. prausnitzii Abundance', fontsize=14, fontweight='bold')
    ax.set_ylabel('P. stomatis Abundance', fontsize=14, fontweight='bold')
    ax.set_title('FP vs PS Quadrant Analysis', fontsize=16, fontweight='bold', pad=20)
    
    ax.legend(loc='upper right', fontsize=10, frameon=True, fancybox=True, shadow=True)
    
    ax.grid(True, alpha=0.3, linestyle='-', linewidth=0.5)
    ax.set_axisbelow(True)
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=dpi, bbox_inches='tight')
        print(f"散点图已保存至: {save_path}")
    
    return fig


def plot_fp_ps_ratio_by_group(df: pd.DataFrame,
                               figsize: Tuple[int, int] = (12, 8),
                               save_path: Optional[str] = None,
                               dpi: int = 300) -> plt.Figure:
    """
    绘制FP/PS比值在三个分组中的分布（箱线图+小提琴图）
    
    Parameters:
    -----------
    df : pd.DataFrame
        包含FP/PS比值和分组信息的DataFrame
    figsize : Tuple[int, int]
        图表尺寸
    save_path : Optional[str]
        保存路径
    dpi : int
        DPI
    
    Returns:
    --------
    plt.Figure
        图表对象
    """
    sns.set(style='white', context='notebook', font_scale=1.1)
    fig, ax = plt.subplots(figsize=figsize)
    
    plot_data = df[df['group'].isin(GROUP_ORDER)].copy()
    
    sns.violinplot(x='group', y='FP_PS_log_ratio', data=plot_data,
                   order=GROUP_ORDER, palette=GROUP_COLORS, inner=None,
                   cut=0, scale='width', linewidth=0, ax=ax)
    
    sns.boxplot(x='group', y='FP_PS_log_ratio', data=plot_data,
                order=GROUP_ORDER, width=0.14, showcaps=True, showfliers=False,
                boxprops=dict(facecolor='none', edgecolor='black', linewidth=1.2),
                whiskerprops=dict(color='black', linewidth=1.2),
                capprops=dict(color='black', linewidth=1.2),
                medianprops=dict(color='black', linewidth=1.6),
                zorder=2, ax=ax)
    
    sns.stripplot(x='group', y='FP_PS_log_ratio', data=plot_data,
                  order=GROUP_ORDER, color='gray', alpha=0.4,
                  size=4, jitter=True, zorder=3, ax=ax)
    
    ax.set_xlabel('Group', fontsize=14, fontweight='bold')
    ax.set_ylabel('log10(FP/PS Ratio)', fontsize=14, fontweight='bold')
    ax.set_title('FP/PS Ratio Across Differentiation Groups', fontsize=16, fontweight='bold', pad=20)
    
    y_min, y_max = ax.get_ylim()
    y_range = y_max - y_min
    
    for i, group in enumerate(GROUP_ORDER):
        if group in plot_data['group'].values:
            group_data = plot_data[plot_data['group'] == group]
            count = len(group_data)
            mean_value = group_data['FP_PS_log_ratio'].mean()
            median_value = group_data['FP_PS_log_ratio'].median()
            
            ax.text(i, y_max - y_range * 0.03, f'n={count}',
                    ha='center', va='bottom', fontsize=10, fontweight='bold')
            
            ax.text(i, y_min + y_range * 0.06,
                    f'Mean: {mean_value:.2f}\nMedian: {median_value:.2f}',
                    ha='center', va='top', fontsize=9,
                    bbox=dict(facecolor='white', alpha=0.7, edgecolor='none'))
    
    pvalues = {}
    group_data_dict = {}
    for group in GROUP_ORDER:
        if group in plot_data['group'].values:
            group_data_dict[group] = plot_data[plot_data['group'] == group]['FP_PS_log_ratio'].dropna().values
    
    comparisons = [
        ('CRC_welldiff', 'CTRL', 0.5),
        ('CRC_poordiff', 'CTRL', 1.0),
        ('CRC_poordiff', 'CRC_welldiff', 1.5)
    ]
    
    same_stage_y = y_max + y_range * 0.1
    
    for i, (g1, g2, x_pos) in enumerate(comparisons):
        if g1 in group_data_dict and g2 in group_data_dict:
            stat, p = stats.mannwhitneyu(group_data_dict[g1], group_data_dict[g2], alternative='two-sided')
            pvalues[f'{g1}_vs_{g2}'] = p
            
            p_text = f'p = {p:.3f}' if p >= 0.001 else 'p < 0.001'
            
            if i == 0:
                y_offset = 0
            elif i == 1:
                y_offset = y_range * 0.1
            else:
                y_offset = 0
            
            ax.text(x_pos, same_stage_y + y_offset, p_text,
                    ha='center', va='bottom', fontsize=10, fontweight='bold')
            
            y_line = same_stage_y + y_offset - y_range * 0.04
            bracket_h = y_range * 0.02
            
            idx1 = GROUP_ORDER.index(g1)
            idx2 = GROUP_ORDER.index(g2)
            
            ax.plot([idx1, idx1], [y_line, y_line + bracket_h], color='black', linewidth=1.5)
            ax.plot([idx1, idx2], [y_line + bracket_h, y_line + bracket_h], color='black', linewidth=1.5)
            ax.plot([idx2, idx2], [y_line, y_line + bracket_h], color='black', linewidth=1.5)
    
    ax.set_ylim(y_min, y_max + y_range * 0.25)
    
    plt.xticks(ticks=range(len(GROUP_ORDER)), 
               labels=[GROUP_NAMES[g] for g in GROUP_ORDER], 
               rotation=45, ha='right', fontsize=12)
    plt.yticks(fontsize=12)
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=dpi, bbox_inches='tight')
        print(f"比值分布图已保存至: {save_path}")
    
    return fig, pvalues


def plot_quadrant_distribution_by_group(df: pd.DataFrame,
                                          figsize: Tuple[int, int] = (14, 6),
                                          save_path: Optional[str] = None,
                                          dpi: int = 300) -> plt.Figure:
    """
    绘制四象限在三个分组中的分布（堆叠条形图）
    
    Parameters:
    -----------
    df : pd.DataFrame
        包含四象限分类和分组信息的DataFrame
    figsize : Tuple[int, int]
        图表尺寸
    save_path : Optional[str]
        保存路径
    dpi : int
        DPI
    
    Returns:
    --------
    plt.Figure
        图表对象
    """
    sns.set(style='white', context='notebook', font_scale=1.1)
    fig, axes = plt.subplots(1, 2, figsize=figsize)
    
    plot_data = df[df['group'].isin(GROUP_ORDER)].copy()
    
    crosstab = pd.crosstab(plot_data['group'], plot_data['quadrant'], normalize='index') * 100
    crosstab = crosstab.reindex(GROUP_ORDER)
    
    quadrant_order = ['FP_high_PS_low', 'FP_low_PS_high', 'FP_high_PS_high', 'FP_low_PS_low']
    crosstab = crosstab[quadrant_order]
    
    crosstab.plot(kind='bar', stacked=True, ax=axes[0],
                  color=[QUADRANT_COLORS[q] for q in quadrant_order],
                  edgecolor='white', linewidth=1)
    
    axes[0].set_xlabel('Group', fontsize=12, fontweight='bold')
    axes[0].set_ylabel('Percentage (%)', fontsize=12, fontweight='bold')
    axes[0].set_title('Quadrant Distribution by Group', fontsize=14, fontweight='bold', pad=15)
    axes[0].set_xticklabels([GROUP_NAMES[g] for g in GROUP_ORDER], rotation=45, ha='right')
    axes[0].legend([QUADRANT_NAMES[q] for q in quadrant_order], 
                   loc='upper right', fontsize=9, frameon=True)
    axes[0].grid(axis='y', alpha=0.3, linestyle='-', linewidth=0.5)
    
    crosstab_counts = pd.crosstab(plot_data['group'], plot_data['quadrant'])
    crosstab_counts = crosstab_counts.reindex(GROUP_ORDER)
    crosstab_counts = crosstab_counts[quadrant_order]
    
    crosstab_counts.plot(kind='bar', stacked=True, ax=axes[1],
                         color=[QUADRANT_COLORS[q] for q in quadrant_order],
                         edgecolor='white', linewidth=1)
    
    axes[1].set_xlabel('Group', fontsize=12, fontweight='bold')
    axes[1].set_ylabel('Count', fontsize=12, fontweight='bold')
    axes[1].set_title('Quadrant Counts by Group', fontsize=14, fontweight='bold', pad=15)
    axes[1].set_xticklabels([GROUP_NAMES[g] for g in GROUP_ORDER], rotation=45, ha='right')
    axes[1].legend([QUADRANT_NAMES[q] for q in quadrant_order], 
                   loc='upper right', fontsize=9, frameon=True)
    axes[1].grid(axis='y', alpha=0.3, linestyle='-', linewidth=0.5)
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=dpi, bbox_inches='tight')
        print(f"四象限分布图已保存至: {save_path}")
    
    return fig


def analyze_fp_ps_ratio(df: pd.DataFrame,
                        output_dir: Optional[str] = None) -> Dict[str, Any]:
    """
    FP/PS比值分析流水线
    
    Parameters:
    -----------
    df : pd.DataFrame
        原始数据
    output_dir : Optional[str]
        输出目录
    
    Returns:
    --------
    Dict[str, Any]
        分析结果
    """
    results = {}
    
    print("=" * 80)
    print("FP/PS比值分析与四象限分析")
    print("=" * 80)
    
    print("\n步骤1: 准备CRC分化程度分组...")
    df_grouped = prepare_crc_diff_groups(df)
    results['group_counts'] = df_grouped['group'].value_counts().to_dict()
    print(f"分组统计: {results['group_counts']}")
    
    print("\n步骤2: 计算FP/PS比值...")
    df_ratio = calculate_fp_ps_ratio(df_grouped)
    
    print("\n步骤3: 进行四象限分析...")
    df_quadrant = perform_quadrant_analysis(df_ratio, method='median')
    
    print("\n步骤4: 统计分析...")
    descriptive_stats = {}
    for group in GROUP_ORDER:
        if group in df_quadrant['group'].values:
            group_data = df_quadrant[df_quadrant['group'] == group]['FP_PS_log_ratio']
            descriptive_stats[group] = {
                'n': len(group_data),
                'mean': float(group_data.mean()),
                'std': float(group_data.std()),
                'median': float(group_data.median()),
                'min': float(group_data.min()),
                'max': float(group_data.max())
            }
    results['descriptive_stats'] = descriptive_stats
    
    print("\n步骤5: 绘制图表...")
    
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
        
        scatter_path = f"{output_dir}/fp_ps_scatter_quadrant.png"
        fig_scatter = plot_fp_ps_scatter_quadrant(df_quadrant, save_path=scatter_path)
        results['scatter_path'] = scatter_path
        
        ratio_path = f"{output_dir}/fp_ps_ratio_by_group.png"
        fig_ratio, pvalues = plot_fp_ps_ratio_by_group(df_quadrant, save_path=ratio_path)
        results['ratio_path'] = ratio_path
        results['pvalues'] = pvalues
        
        quadrant_path = f"{output_dir}/quadrant_distribution_by_group.png"
        fig_quadrant = plot_quadrant_distribution_by_group(df_quadrant, save_path=quadrant_path)
        results['quadrant_path'] = quadrant_path
        
        output_file = f"{output_dir}/fp_ps_ratio_analysis.csv"
        df_quadrant[['group', 'FP_abundance', 'PS_abundance', 'FP_PS_ratio', 
                     'FP_PS_log_ratio', 'quadrant']].to_csv(output_file, index=True)
        results['data_file'] = output_file
        print(f"\n数据已保存至: {output_file}")
    
    print("\n" + "=" * 80)
    print("FP/PS比值分析完成!")
    print("=" * 80)
    
    return results


def main():
    """主函数"""
    print("=" * 80)
    print("F. prausnitzii / P. stomatis 比值分析")
    print("分析FP和PS的生态竞争关系及其与CRC分化的关联")
    print("=" * 80)
    
    ds = BioSmokeDataset()
    ds.preprocess_taxonomy_data()
    ds.preprocess_metabolomics_data()
    df = ds.merge_to_dataframe()
    
    print(f"\n总样本数: {df.shape[0]}")
    print(f"总特征数: {df.shape[1]}")
    
    output_dir = './results/fp_ps_ratio_analysis'
    
    results = analyze_fp_ps_ratio(df=df, output_dir=output_dir)
    
    print("\n关键发现:")
    print("-" * 80)
    
    if 'descriptive_stats' in results:
        print("\nFP/PS比值（log10转换）统计:")
        for group, stats in results['descriptive_stats'].items():
            print(f"  {GROUP_NAMES[group]}: Mean={stats['mean']:.2f}, "
                  f"Median={stats['median']:.2f}, n={stats['n']}")
    
    if 'pvalues' in results:
        print("\n组间比较（Mann-Whitney U检验）:")
        for comparison, p in results['pvalues'].items():
            sig = '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else 'ns'
            print(f"  {comparison}: p={p:.4f} {sig}")


if __name__ == "__main__":
    main()
