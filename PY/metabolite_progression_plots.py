import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import re
from pathlib import Path
from scipy import stats
from dataset import BioSmokeDataset
from split_group import prepare_crc_diff_groups


def _normalize_metabolite_name(name):
    return re.sub(r'[^a-z0-9]+', '', str(name).lower())


def plot_metabolite_progression(df, target_metabolites, output_dir):
    """
    绘制指定代谢物的 Heatmap 和 Boxplot 展现疾病演变
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

    # 2. 匹配代谢物列名
    metabolite_cols = [col for col in df.columns if col.startswith('met_')]
    matched_cols = []
    metabolite_lookup = {}
    for col in metabolite_cols:
        raw_name = col[4:].strip().strip('"')
        metabolite_lookup.setdefault(_normalize_metabolite_name(raw_name), []).append(col)

    for m in target_metabolites:
        normalized_name = _normalize_metabolite_name(m)

        exact_matches = [
            col for col in metabolite_cols
            if col[4:].strip().strip('"').lower() == str(m).strip().strip('"').lower()
        ]
        if len(exact_matches) == 1:
            matched_cols.append(exact_matches[0])
            continue

        normalized_matches = metabolite_lookup.get(normalized_name, [])
        if len(normalized_matches) == 1:
            matched_cols.append(normalized_matches[0])
            continue

        fallback_matches = [
            col for col in metabolite_cols
            if normalized_name and normalized_name in _normalize_metabolite_name(col[4:])
        ]
        if len(fallback_matches) == 1:
            matched_cols.append(fallback_matches[0])
        elif len(fallback_matches) > 1:
            print(f"Warning: Metabolite {m} matched multiple columns: {fallback_matches}; skipped.")
        else:
            print(f"Warning: Metabolite {m} not found in dataset.")

    if not matched_cols:
        print("No target metabolites found. Exiting.")
        return

    # --- 绘图 1: Heatmap (Z-score normalized) ---
    plt.figure(figsize=(10, 6))

    # 按分组排序样本
    plot_df = df.sort_values('display_group')
    heatmap_data = plot_df[matched_cols].T

    # Z-score 标准化
    heatmap_data_z = heatmap_data.apply(lambda x: (x - x.mean()) / x.std(), axis=1)

    # 绘制热图，不显示样本名称
    sns.heatmap(
        heatmap_data_z,
        cmap='RdBu_r',
        center=0,
        xticklabels=False,
        yticklabels=[m.replace('met_', '') for m in matched_cols],
    )

    # 在顶部添加分组背景色/标签
    ax = plt.gca()
    group_counts = plot_df['display_group'].value_counts().reindex(display_group_order)

    start = 0
    colors = ['#2E86AB', '#F18F01', '#D7263D']  # CTRL, Well, Poor
    for i, count in enumerate(group_counts):
        if count > 0:
            ax.axvspan(start, start + count, color=colors[i], alpha=0.3, label=display_group_order[i])
            plt.text(start + count / 2, -0.5, display_group_order[i], ha='center', fontweight='bold')
            start += count

    plt.title("Metabolite Progression Heatmap (Z-score)", fontsize=14)
    plt.savefig(output_dir / "metabolite_progression_heatmap.png", dpi=300, bbox_inches='tight')
    plt.close()

    # --- 绘图 2: Barplots (Grid) ---
    n_metabs = len(matched_cols)
    cols = 5
    rows = (n_metabs + cols - 1) // cols

    fig, axes = plt.subplots(rows, cols, figsize=(20, 4 * rows))
    axes = axes.flatten()

    for i, col in enumerate(matched_cols):
        ax = axes[i]
        sns.barplot(
            data=df,
            x='display_group',
            y=col,
            order=display_group_order,
            palette=colors,
            ax=ax,
            estimator=np.mean,
            errorbar='sd', # 增加标准差误差棒 (Standard Deviation)
            errwidth= 1,    # 增加误差棒宽度
            capsize=0.1,   # 增加误差棒顶端横线
        )

        # Clean up labels and title
        import textwrap
        clean_name = col.replace('met_', '').replace('_', ' ')
        wrapped_name = "\n".join(textwrap.wrap(clean_name, width=25)) # 换行长度设为25
        ax.set_title(wrapped_name, fontsize=11, fontweight='bold', pad=25) # 增加 pad 防止重叠
        ax.set_xlabel('')
        ax.set_ylabel('log1p(Abund.)', fontsize=9)
        ax.tick_params(axis='both', labelsize=8)

        # Use data range for annotation placement
        group_stats = df.groupby('display_group', observed=False)[col].agg(['mean', 'std'])
        y_min = df[col].min()
        y_max = (group_stats['mean'] + group_stats['std']).max()
        y_range = y_max - y_min
        if y_range == 0:
            y_range = 1.0

        # --- Add P-values (Pairwise Wilcoxon/Mann-Whitney) ---
        comparisons = [
            (group_mapping['CTRL'], group_mapping['CRC_welldiff'], 0, 1, 0.1),  # 增加 offset
            (group_mapping['CRC_welldiff'], group_mapping['CRC_poordiff'], 1, 2, 0.1),
            (group_mapping['CTRL'], group_mapping['CRC_poordiff'], 0, 2, 0.3),  # 顶层连线间距拉大
        ]

        # 统一基准高度，使用数据最大值加上一定的 padding
        base_y = y_max + y_range * 0.05 

        for g1, g2, x1, x2, offset in comparisons:
            data1 = df[df['display_group'] == g1][col]
            data2 = df[df['display_group'] == g1][col] # 注意：之前代码这里可能有误，应为 g2，但我先按逻辑修复

            # 重新获取正确的数据用于统计
            data1 = df[df['display_group'] == g1][col]
            data2 = df[df['display_group'] == g2][col]

            if len(data1) > 0 and len(data2) > 0:
                stat, p = stats.mannwhitneyu(data1, data2, alternative='two-sided')

                if p < 0.001:
                    star = '***'
                elif p < 0.01:
                    star = '**'
                elif p < 0.05:
                    star = '*'
                else:
                    star = None

                if star:
                    y_line = base_y + y_range * offset
                    h = y_range * 0.03 # 增加连线转角的垂直长度
                    ax.plot([x1, x1, x2, x2], [y_line, y_line + h, y_line + h, y_line], lw=1, c='black')
                    ax.text(
                        (x1 + x2) / 2,
                        y_line + h + y_range * 0.01, # 文字再往上移一点点，避免压线
                        star,
                        ha='center',
                        va='bottom',
                        fontsize=10,
                        color='red',
                        fontweight='bold',
                    )
        
        # 动态调整 y 轴范围，确保 P 值标注不被切掉
        ax.set_ylim(y_min, base_y + y_range * 0.5)

        # Let Matplotlib autoscale based on drawn annotations

    # 移除多余的子图
    for j in range(i + 1, len(axes)):
        fig.delaxes(axes[j])

    # 在图底部添加 P 值显著性图例
    fig.text(
        0.5, 0.02, 
        'Significance: *** p < 0.001, ** p < 0.01, * p < 0.05', 
        ha='center', fontsize=10, style='italic', fontweight='bold'
    )

    # 调整子图间距
    plt.subplots_adjust(hspace=0.6, wspace=0.4, bottom=0.08) # 留出底部空间给图例
    # plt.tight_layout() # tight_layout 有时会和 subplots_adjust 冲突，先注释掉看看效果
    plt.savefig(output_dir / "metabolite_progression_barplot.png", dpi=300, bbox_inches='tight')
    plt.close()

    print(f"Plots saved to {output_dir}")


if __name__ == "__main__":
    # 加载数据
    ds = BioSmokeDataset()
    ds.preprocess_metabolomics_data(transform=True)
    merged = ds.intersection_to_dataframe()

    target_metabolites = [
        "SQDG 26:2; SQDG(13:1/13:1)",
        "Cytosine",
        "Perfluorooctanesulfonic acid",
        "Methyl dihydrojasmonate",
        "Pyrogallol-2-O-sulphate",
        "5'-(3',4'-Dihydroxyphenyl)-gamma-valerolactone sulfate",
        "2-Hydroxy-4,7-dimethoxy-2H-1,4-benzoxazin-3(4H)-one",
        "trans-3,5-Dimethoxy-4-hydroxycinnamaldehyde",
        "(R)-3-Hydroxy-5-phenylpentanoic acid",
        "N-Methyl-D-glucamine",
        "Chenodeoxycholic acid sulfate",
        "Lucidenic acid F",
        "Demissidine",
        "Alpha-Hydroxyisobutyric acid",
        "Pyrocatechol",
        "Gentisic acid",
        "D-Galacturonic acid",
        "1,3-Dimethyluric acid",
        "4-Hydroxy-5-(phenyl)-valeric acid-O-sulphate",
        "Cholesterol"
    ]

    plot_metabolite_progression(merged, target_metabolites, "results/progression_plots")
