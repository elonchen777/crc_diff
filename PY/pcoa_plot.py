import os
import pandas as pd
import numpy as np
from dataset import BioSmokeDataset
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import pairwise_distances
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import Ellipse
try:
    from scipy.stats import chi2
except Exception:
    chi2 = None
import warnings
warnings.filterwarnings('ignore')
from split_group import prepare_crc_diff_groups


def pcoa_and_plot(df: pd.DataFrame, prefix: str, out_png: str, distance_metric: str = 'braycurtis'):
    """
    执行PCoA分析并绘图
    
    参数:
        df: 包含特征和分组的数据框
        prefix: 特征前缀 ('tax_' for species, 'met_' for metabolites)
        out_png: 输出图片路径
        distance_metric: 距离度量方法 ('braycurtis', 'euclidean', 'jaccard', etc.)
    """
    # 提取特征
    features = [c for c in df.columns if c.startswith(prefix)]
    if len(features) == 0:
        print(f'没有找到以 "{prefix}" 开头的特征列')
        return
    
    print(f"分析 {prefix} 数据: {len(features)} 个特征, {df.shape[0]} 个样本")
    
    X = df[features].values.astype(float)
    
    # 标准化数据
    X_scaled = StandardScaler().fit_transform(X)
    
    # 计算距离矩阵
    print(f"计算距离矩阵 (metric: {distance_metric})...")
    try:
        D = pairwise_distances(X_scaled, metric=distance_metric)
    except Exception as e:
        print(f"距离计算失败: {e}, 使用欧氏距离")
        D = pairwise_distances(X_scaled, metric='euclidean')
    
    # PCoA (经典多维尺度分析)
    print("执行PCoA分析...")
    
    # 中心化距离矩阵
    n = D.shape[0]
    H = np.eye(n) - np.ones((n, n)) / n
    B = -0.5 * H @ (D ** 2) @ H
    
    # 特征值分解
    try:
        eigvals, eigvecs = np.linalg.eigh(B)
    except Exception as e:
        print(f"特征值分解失败: {e}")
        return
    
    # 按特征值降序排列
    idx = eigvals.argsort()[::-1]
    eigvals = eigvals[idx]
    eigvecs = eigvecs[:, idx]
    
    # 取前两个主坐标
    pcoa_coords = eigvecs[:, :2] * np.sqrt(np.maximum(eigvals[:2], 0))
    
    # 创建结果数据框
    pcoa_df = pd.DataFrame(pcoa_coords, index=df.index, columns=['PCo1', 'PCo2'])
    pcoa_df['group'] = df['group'].values
    
    # 计算解释方差比例
    total_variance = np.sum(np.maximum(eigvals, 0))
    explained_variance = np.maximum(eigvals[:2], 0) / total_variance if total_variance > 0 else [0, 0]
    
    print(f"PCo1解释方差: {explained_variance[0]:.4f}")
    print(f"PCo2解释方差: {explained_variance[1]:.4f}")
    
    # 绘图
    plt.figure(figsize=(10, 8))
    ax = plt.gca()
    
    # 定义三组的颜色
    groups = ['CTRL', 'CRC_welldiff', 'CRC_poordiff']
    colors = {
        'CRC_poordiff': '#FF6B6B',
        'CRC_welldiff': '#FFD166',
        'CTRL': '#06D6A0'
    }
    
    # 绘制散点图
    for group in groups:
        if group in pcoa_df['group'].values:
            group_data = pcoa_df[pcoa_df['group'] == group]
            ax.scatter(group_data['PCo1'], group_data['PCo2'], 
                      color=colors[group], s=60, label=group, alpha=0.7, linewidth=0)
    
    # 置信椭圆：95% 的 chi2 临界值（2维）约为 5.991
    chi2_val = 5.991
    if chi2 is not None:
        try:
            chi2_val = float(chi2.ppf(0.95, 2))
        except Exception:
            chi2_val = 5.991
    
    # 为每组添加置信椭圆和均值标记
    for group in groups:
        if group in pcoa_df['group'].values:
            group_data = pcoa_df[pcoa_df['group'] == group][['PCo1', 'PCo2']]
            if group_data.shape[0] >= 2:  # 至少需要2个点才能计算椭圆
                mean = group_data.mean().values
                
                # 计算协方差并绘制椭圆
                cov = np.cov(group_data.values.T)
                try:
                    vals, vecs = np.linalg.eigh(cov)
                except Exception:
                    vals = np.array([0.0, 0.0])
                    vecs = np.eye(2)
                
                # 将特征值按降序排列
                order = vals.argsort()[::-1]
                vals = vals[order]
                vecs = vecs[:, order]
                
                angle = np.degrees(np.arctan2(vecs[1, 0], vecs[0, 0]))
                width, height = 2 * np.sqrt(np.clip(vals * chi2_val, 0, None))
                
                ellipse = Ellipse(xy=mean, width=width, height=height, angle=angle,
                                  edgecolor='none', facecolor=colors[group], 
                                  lw=0, alpha=0.3)
                ax.add_patch(ellipse)
    
    # 设置图表属性
    if prefix == 'tax_':
        title = 'PCoA of Species (Bray-Curtis Distance)'
        xlabel = 'PCo1 - Species'
        ylabel = 'PCo2 - Species'
        filename_suffix = 'species'
    else:
        title = 'PCoA of Metabolites (Bray-Curtis Distance)'
        xlabel = 'PCo1 - Metabolites'
        ylabel = 'PCo2 - Metabolites'
        filename_suffix = 'metabolites'
    
    ax.set_xlabel(f'{xlabel} ({explained_variance[0]:.1%})', fontsize=12)
    ax.set_ylabel(f'{ylabel} ({explained_variance[1]:.1%})', fontsize=12)
    ax.set_title(f'{title}\nCRC Patients by Differentiation', fontsize=14, fontweight='bold')
    
    # 添加图例
    ax.legend(loc='upper right', fontsize=10, framealpha=0.9)
    
    # 添加网格
    ax.grid(True, alpha=0.3, linestyle='--')
    
    plt.tight_layout()
    
    # 保存图片
    plt.savefig(out_png, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f'已保存: {out_png}')
    
    # 返回PCoA结果
    return pcoa_df, explained_variance


def calculate_permanova(df: pd.DataFrame, prefix: str, group_col: str = 'group'):
    """
    计算PERMANOVA统计检验（组间差异显著性）
    
    参数:
        df: 包含特征和分组的数据框
        prefix: 特征前缀
        group_col: 分组列名
    
    返回:
        F统计量和p值
    """
    try:
        from skbio.stats.distance import permanova
        from skbio import DistanceMatrix
        
        # 提取特征
        features = [c for c in df.columns if c.startswith(prefix)]
        if len(features) == 0:
            return None, None
        
        X = df[features].values.astype(float)
        X_scaled = StandardScaler().fit_transform(X)
        
        # 计算距离矩阵
        D = pairwise_distances(X_scaled, metric='braycurtis')
        
        # 创建距离矩阵对象
        dm = DistanceMatrix(D, ids=[str(i) for i in range(len(D))])
        
        # 执行PERMANOVA
        groups = df[group_col].values
        result = permanova(dm, groups, permutations=6000)

        print(result)
        
        return result['test statistic'], result['p-value']
    
    except ImportError:
        print("警告: 未安装skbio库，无法计算PERMANOVA")
        return None, None
    except Exception as e:
        print(f"PERMANOVA计算失败: {e}")
        return None, None


def pcoa_pairwise_comparison(df: pd.DataFrame, prefix: str, group1: str, group2: str, 
                            output_dir: str = 'figures', distance_metric: str = 'braycurtis'):
    """
    执行俩俩对比的PCoA分析并绘图
    
    参数:
        df: 包含特征和分组的数据框
        prefix: 特征前缀 ('tax_' for species, 'met_' for metabolites)
        group1: 第一组名称
        group2: 第二组名称
        output_dir: 输出目录
        distance_metric: 距离度量方法
    """
    # 筛选两组数据
    df_pair = df[df['group'].isin([group1, group2])].copy()
    
    if df_pair.empty:
        print(f"警告: 没有找到 {group1} 或 {group2} 的数据")
        return None, None, None
    
    print(f"\n分析 {group1} vs {group2}: {len(df_pair)} 个样本")
    
    # 提取特征
    features = [c for c in df_pair.columns if c.startswith(prefix)]
    if len(features) == 0:
        print(f'没有找到以 "{prefix}" 开头的特征列')
        return None, None, None
    
    X = df_pair[features].values.astype(float)
    
    # 标准化数据
    X_scaled = StandardScaler().fit_transform(X)
    
    # 计算距离矩阵
    try:
        D = pairwise_distances(X_scaled, metric=distance_metric)
    except Exception as e:
        print(f"距离计算失败: {e}, 使用欧氏距离")
        D = pairwise_distances(X_scaled, metric='euclidean')
    
    # PCoA (经典多维尺度分析)
    n = D.shape[0]
    H = np.eye(n) - np.ones((n, n)) / n
    B = -0.5 * H @ (D ** 2) @ H
    
    # 特征值分解
    try:
        eigvals, eigvecs = np.linalg.eigh(B)
    except Exception as e:
        print(f"特征值分解失败: {e}")
        return None, None, None
    
    # 按特征值降序排列
    idx = eigvals.argsort()[::-1]
    eigvals = eigvals[idx]
    eigvecs = eigvecs[:, idx]
    
    # 取前两个主坐标
    pcoa_coords = eigvecs[:, :2] * np.sqrt(np.maximum(eigvals[:2], 0))
    
    # 创建结果数据框
    pcoa_df = pd.DataFrame(pcoa_coords, index=df_pair.index, columns=['PCo1', 'PCo2'])
    pcoa_df['group'] = df_pair['group'].values
    
    # 计算解释方差比例
    total_variance = np.sum(np.maximum(eigvals, 0))
    explained_variance = np.maximum(eigvals[:2], 0) / total_variance if total_variance > 0 else [0, 0]
    
    # 绘图
    plt.figure(figsize=(10, 8))
    ax = plt.gca()
    
    # 定义颜色
    colors = {
        'CRC_poordiff': '#FF6B6B',
        'CRC_welldiff': '#FFD166',
        'CTRL': '#06D6A0'
    }
    
    # 绘制散点图
    for group in [group1, group2]:
        if group in pcoa_df['group'].values:
            group_data = pcoa_df[pcoa_df['group'] == group]
            ax.scatter(group_data['PCo1'], group_data['PCo2'], 
                      color=colors.get(group, '#808080'), s=80, label=group, 
                      alpha=0.8, linewidth=0)
    
    # 置信椭圆：95% 的 chi2 临界值（2维）约为 5.991
    chi2_val = 5.991
    if chi2 is not None:
        try:
            chi2_val = float(chi2.ppf(0.95, 2))
        except Exception:
            chi2_val = 5.991
    
    # 为每组添加置信椭圆和均值标记
    for group in [group1, group2]:
        if group in pcoa_df['group'].values:
            group_data = pcoa_df[pcoa_df['group'] == group][['PCo1', 'PCo2']]
            if group_data.shape[0] >= 2:  # 至少需要2个点才能计算椭圆
                mean = group_data.mean().values
                
                # 计算协方差并绘制椭圆
                cov = np.cov(group_data.values.T)
                try:
                    vals, vecs = np.linalg.eigh(cov)
                except Exception:
                    vals = np.array([0.0, 0.0])
                    vecs = np.eye(2)
                
                # 将特征值按降序排列
                order = vals.argsort()[::-1]
                vals = vals[order]
                vecs = vecs[:, order]
                
                angle = np.degrees(np.arctan2(vecs[1, 0], vecs[0, 0]))
                width, height = 2 * np.sqrt(np.clip(vals * chi2_val, 0, None))
                
                ellipse = Ellipse(xy=mean, width=width, height=height, angle=angle,
                                  edgecolor='none', facecolor=colors.get(group, '#808080'), 
                                  lw=0, alpha=0.3)
                ax.add_patch(ellipse)
    
    # 设置图表属性
    if prefix == 'tax_':
        title = f'PCoA of Species: {group1} vs {group2}'
        xlabel = 'PCo1 - Species'
        ylabel = 'PCo2 - Species'
        filename_suffix = 'species'
    else:
        title = f'PCoA of Metabolites: {group1} vs {group2}'
        xlabel = 'PCo1 - Metabolites'
        ylabel = 'PCo2 - Metabolites'
        filename_suffix = 'metabolites'
    
    ax.set_xlabel(f'{xlabel} ({explained_variance[0]:.1%})', fontsize=12)
    ax.set_ylabel(f'{ylabel} ({explained_variance[1]:.1%})', fontsize=12)
    ax.set_title(f'{title}\n(Bray-Curtis Distance)', fontsize=14, fontweight='bold')
    
    # 添加图例
    ax.legend(loc='upper right', fontsize=11, framealpha=0.9)
    
    # 添加网格
    # ax.grid(True, alpha=0.3, linestyle='--')
    
    # 计算PERMANOVA p值
    F_stat, p_value = calculate_permanova(df_pair, prefix, 'group')
    
    # 添加p值标注
    if F_stat is not None and p_value is not None:
        significance = "***" if p_value < 0.001 else "**" if p_value < 0.01 else "*" if p_value < 0.05 else "ns"
        p_text = f'PERMANOVA: F = {F_stat:.2f}, p = {p_value:.4f} ({significance})'
        ax.text(0.02, 0.98, p_text, transform=ax.transAxes,
                fontsize=11, fontweight='bold', va='top', ha='left',
                bbox=dict(facecolor='white', alpha=0.8, edgecolor='black'))
    
    plt.tight_layout()
    
    # 生成文件名
    safe_group1 = group1.replace('CRC_', '').replace('_', '-')
    safe_group2 = group2.replace('CRC_', '').replace('_', '-')
    out_png = os.path.join(output_dir, f'pcoa_{filename_suffix}_{safe_group1}_vs_{safe_group2}.png')
    
    # 保存图片
    plt.savefig(out_png, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f'  已保存: {out_png}')
    print(f'  解释方差: PCo1 = {explained_variance[0]:.3f}, PCo2 = {explained_variance[1]:.3f}')
    if F_stat is not None and p_value is not None:
        print(f'  PERMANOVA: F = {F_stat:.4f}, p = {p_value:.6f}')
    
    return pcoa_df, explained_variance, (F_stat, p_value)


def main():    
    # 加载数据
    ds = BioSmokeDataset()
    ds.preprocess_taxonomy_data()
    ds.preprocess_metabolomics_data()
    merged = ds.merge_to_dataframe()
    
    print(f"总样本数: {merged.shape[0]}")
    print(f"总特征数: {merged.shape[1]}")
    
    # 准备分组
    merged = prepare_crc_diff_groups(merged)
    print(f"\n过滤后样本数: {merged.shape[0]}")
    
    # 打印各组样本数
    print("\n各组样本数:")
    groups = ['CTRL', 'CRC_welldiff', 'CRC_poordiff']
    group_counts = merged['group'].value_counts()
    for group in groups:
        if group in group_counts:
            print(f"  {group}: {group_counts[group]} 样本")
        else:
            print(f"  {group}: 0 样本")
    
    # 创建输出目录
    output_dir = 'results/pcoa_plots'
    os.makedirs(output_dir, exist_ok=True)  
    
    # 执行PCoA分析
    print("\n" + "=" * 80)
    print("1. 物种PCoA分析 (Species)")
    print("=" * 80)
    pcoa_species, var_species = pcoa_and_plot(
        merged, 
        'tax_', 
        os.path.join(output_dir, 'pcoa_species.png'),
        distance_metric='braycurtis'
    )
    
    print("\n" + "=" * 80)
    print("2. 代谢物PCoA分析 (Metabolites)")
    print("=" * 80)
    pcoa_metabolites, var_metabolites = pcoa_and_plot(
        merged, 
        'met_', 
        os.path.join(output_dir, 'pcoa_metabolites.png'),
        distance_metric='braycurtis'
    )
    
    # 计算PERMANOVA统计检验
    print("\n" + "=" * 80)
    print("3. PERMANOVA统计检验 (组间差异显著性)")
    print("=" * 80)
    
    # 保存PCoA坐标数据
    if pcoa_species is not None:
        pcoa_species.to_csv(os.path.join(output_dir, 'pcoa_species_coordinates.csv'))
        print(f"\n物种PCoA坐标已保存至: {output_dir}/pcoa_species_coordinates.csv")
    
    if pcoa_metabolites is not None:
        pcoa_metabolites.to_csv(os.path.join(output_dir, 'pcoa_metabolites_coordinates.csv'))
        print(f"代谢物PCoA坐标已保存至: {output_dir}/pcoa_metabolites_coordinates.csv")
    
    # 定义所有可能的俩俩对比组合
    pairwise_comparisons = [
        ('CRC_poordiff', 'CRC_welldiff', '低分化CRC vs 中高分化CRC'),
        ('CRC_poordiff', 'CTRL', '低分化CRC vs 健康对照'),
        ('CRC_welldiff', 'CTRL', '中高分化CRC vs 健康对照')
    ]
    
    # 存储俩俩对比结果
    pairwise_results = {
        'species': {},
        'metabolites': {}
    }
    
    # 执行物种俩俩对比分析
    print("\n4.1 物种俩俩对比PCoA分析:")
    for group1, group2, description in pairwise_comparisons:
        print(f"\n  {description}:")
        result = pcoa_pairwise_comparison(
            merged, 'tax_', group1, group2, output_dir, 'braycurtis'
        )
        if result[0] is not None:
            pairwise_results['species'][f"{group1}_vs_{group2}"] = {
                'description': description,
                'explained_variance': result[1],
                'F_stat': result[2][0] if result[2][0] is not None else None,
                'p_value': result[2][1] if result[2][1] is not None else None
            }
    
    # 执行代谢物俩俩对比分析
    print("\n4.2 代谢物俩俩对比PCoA分析:")
    for group1, group2, description in pairwise_comparisons:
        print(f"\n  {description}:")
        result = pcoa_pairwise_comparison(
            merged, 'met_', group1, group2, output_dir, 'braycurtis'
        )
        if result[0] is not None:
            pairwise_results['metabolites'][f"{group1}_vs_{group2}"] = {
                'description': description,
                'explained_variance': result[1],
                'F_stat': result[2][0] if result[2][0] is not None else None,
                'p_value': result[2][1] if result[2][1] is not None else None
            }

if __name__ == '__main__':
    main()
