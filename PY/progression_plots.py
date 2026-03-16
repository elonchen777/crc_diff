import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from scipy import stats
from dataset import BioSmokeDataset
from split_group import prepare_crc_diff_groups

def plot_microbial_progression(df, target_species, output_dir):
    """
    绘制指定物种的 Heatmap 和 Boxplot 展现疾病演变
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # 1. 准备分组
    if 'group' not in df.columns:
        df = prepare_crc_diff_groups(df)
    
    # 过滤掉不在列表中的样本（如果有的话）并设置顺序
    group_order = ['CTRL', 'CRC_welldiff', 'CRC_poordiff']
    df = df[df['group'].isin(group_order)].copy()
    df['group'] = pd.Categorical(df['group'], categories=group_order, ordered=True)

    # 映射分组名称
    group_mapping = {
        "CTRL": "CTRL",
        "CRC_welldiff": "CRC-Well",
        "CRC_poordiff": "CRC-Poor"
    }
    df['display_group'] = df['group'].map(group_mapping)
    display_group_order = [group_mapping[g] for g in group_order]
    df['display_group'] = pd.Categorical(df['display_group'], categories=display_group_order, ordered=True)
    
    # 2. 匹配物种列名
    species_cols = [col for col in df.columns if col.startswith('tax_s__')]
    matched_cols = []
    for s in target_species:
        matches = [col for col in species_cols if s.lower() in col.lower()]
        if matches:
            matched_cols.append(matches[0])
        else:
            print(f"Warning: Species {s} not found in dataset.")
            
    if not matched_cols:
        print("No target species found. Exiting.")
        return

    # --- 绘图 1: Heatmap (Z-score normalized) ---
    plt.figure(figsize=(10, 6))
    
    # 按分组排序样本
    plot_df = df.sort_values('display_group')
    heatmap_data = plot_df[matched_cols].T
    
    # Z-score 标准化
    heatmap_data_z = heatmap_data.apply(lambda x: (x - x.mean()) / x.std(), axis=1)
    
    # 绘制热图，不显示样本名称
    sns.heatmap(heatmap_data_z, cmap='RdBu_r', center=0, 
                xticklabels=False, yticklabels=[s.replace('tax_s__', '') for s in matched_cols])
    
    # 在顶部添加分组背景色/标签
    ax = plt.gca()
    n_samples = plot_df.shape[0]
    group_counts = plot_df['display_group'].value_counts().reindex(display_group_order)
    
    start = 0
    colors = ['#2E86AB', '#F18F01', '#D7263D'] # CTRL, Well, Poor
    for i, count in enumerate(group_counts):
        if count > 0:
            ax.axvspan(start, start + count, color=colors[i], alpha=0.3, label=display_group_order[i])
            plt.text(start + count/2, -0.5, display_group_order[i], ha='center', fontweight='bold')
            start += count

    plt.title("Microbial Progression Heatmap (Z-score)", fontsize=14)
    plt.savefig(output_dir / "microbial_progression_heatmap.png", dpi=300, bbox_inches='tight')
    plt.close()

    # --- 绘图 2: Boxplots (Grid) ---
    n_taxa = len(matched_cols)
    cols = 5
    rows = (n_taxa + cols - 1) // cols
    
    fig, axes = plt.subplots(rows, cols, figsize=(20, 4 * rows))
    axes = axes.flatten()
    
    for i, col in enumerate(matched_cols):
        ax = axes[i]
        # Remove outliers by setting showfliers=False
        sns.boxplot(data=df, x='display_group', y=col, order=display_group_order, 
                    palette=colors, ax=ax, showfliers=False)
        # Stripplot with small size and alpha for individual distributions
        sns.stripplot(data=df, x='display_group', y=col, order=display_group_order, 
                      color='black', size=1.5, alpha=0.2, ax=ax)
        
        # Clean up labels and title
        clean_name = col.replace('tax_s__', '').replace('_', ' ')
        ax.set_title(clean_name, fontsize=11, style='italic', fontweight='bold')
        ax.set_xlabel('')
        ax.set_ylabel('log1p(Rel. Abund.)', fontsize=9)
        ax.tick_params(axis='both', labelsize=8)
        
        # Adjust y-limit to focus on box part
        group_stats = df.groupby('display_group')[col].describe()
        upper_whisker = (group_stats['75%'] + 1.5 * (group_stats['75%'] - group_stats['25%'])).max()
        y_max = upper_whisker * 1.2
        ax.set_ylim(bottom=df[col].min() * 0.9, top=y_max)

        # --- Add P-values (Pairwise Wilcoxon/Mann-Whitney) ---
        # Define comparisons: (group1, group2, x1, x2, y_offset_factor)
        comparisons = [
            (group_mapping['CTRL'], group_mapping['CRC_welldiff'], 0, 1, 0.05),
            (group_mapping['CRC_welldiff'], group_mapping['CRC_poordiff'], 1, 2, 0.05),
            (group_mapping['CTRL'], group_mapping['CRC_poordiff'], 0, 2, 0.18)
        ]
        
        y_range = y_max - df[col].min()
        base_y = upper_whisker
        
        for g1, g2, x1, x2, offset in comparisons:
            data1 = df[df['display_group'] == g1][col]
            data2 = df[df['display_group'] == g2][col]
            
            if len(data1) > 0 and len(data2) > 0:
                stat, p = stats.mannwhitneyu(data1, data2, alternative='two-sided')
                
                # Get star symbol
                if p < 0.001: star = '***'
                elif p < 0.01: star = '**'
                elif p < 0.05: star = '*'
                else: star = None  # Skip if not significant
                
                # Only draw bracket and text if significant
                if star:
                    y_line = base_y + y_range * offset
                    h = y_range * 0.02
                    ax.plot([x1, x1, x2, x2], [y_line, y_line + h, y_line + h, y_line], 
                            lw=1, c='black')
                    ax.text((x1 + x2) / 2, y_line + h, star, 
                            ha='center', va='bottom', fontsize=9, 
                            color='red', fontweight='bold')

        # Update y-limit to accommodate brackets
        ax.set_ylim(top=base_y + y_range * 0.35)
    # 移除多余的子图
    for j in range(i + 1, len(axes)):
        fig.delaxes(axes[j])

    # 在图底部添加 P 值显著性图例
    fig.text(
        0.5, - 0.02, 
        'Significance: *** p < 0.001, ** p < 0.01, * p < 0.05', 
        ha='center', fontsize=10, style='italic', fontweight='bold'
    )
        
    plt.tight_layout()
    plt.savefig(output_dir / "taxa_progression_boxplot.png", dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Plots saved to {output_dir}")

if __name__ == "__main__":
    # 加载数据
    ds = BioSmokeDataset()
    ds.preprocess_taxonomy_data(remove_low_expression=False, transform=True)
    merged = ds.intersection_to_dataframe()
    
    target_species = [
        'Peptostreptococcus_stomatis',
        'Porphyromonas_gingivalis',
        'Prevotella_intermedia',
        'Fusobacterium_periodonticum',
        'Campylobacter_rectus',
        'Faecalibacterium_prausnitzii',
        'Roseburia_intestinalis',
        'Eubacterium_rectale',
        'Coprococcus_comes',
        'Ruminococcus_lactaris'
    ]
    
    plot_microbial_progression(merged, target_species, "results/progression_plots")
