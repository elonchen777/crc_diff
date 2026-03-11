import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from typing import Dict, List, Optional, Tuple, Set
import warnings
warnings.filterwarnings('ignore')
from dataset import BioSmokeDataset
from pathlib import Path
from split_group import prepare_crc_diff_groups

def calculate_unified_group_correlations(
    df: pd.DataFrame, 
    top_n_species: int = 10, 
    top_n_metabolites: int = 20,
    method: str = 'spearman',
    fixed_species: Optional[List[str]] = None,
    fixed_metabolites: Optional[List[str]] = None
) -> Tuple[Dict[str, pd.DataFrame], Dict[str, pd.DataFrame], List[str], List[str]]:
    """
    Calculate correlations between species and metabolites within each group using common species and metabolites
    
    Args:
        df: DataFrame containing group information and features
        top_n_species: Number of top species to select for visualization
        top_n_metabolites: Number of top metabolites to select for visualization
        method: Correlation calculation method ('spearman' or 'pearson')
        fixed_species: Optional list of fixed species names (without tax_ prefix) to use
        fixed_metabolites: Optional list of fixed metabolite names (without met_ prefix) to use
        
    Returns:
        Tuple: (correlations_by_group, pvalues_by_group, common_species, common_metabolites)
        - correlations_by_group: Dictionary with group names as keys and correlation matrices as values
        - pvalues_by_group: Dictionary with group names as keys and p-value matrices as values
        - common_species: List of commonly selected species
        - common_metabolites: List of commonly selected metabolites
    """
    # Separate species and metabolite features
    species_cols = [col for col in df.columns if col.startswith('tax_')]
    
    metabolite_cols = [col for col in df.columns if col.startswith('met_')]
    
    if not species_cols:
        raise ValueError("Species feature columns (tax_*) not found")
    if not metabolite_cols:
        raise ValueError("Metabolite feature columns (met_*) not found")
    
    # Get group information
    if 'group' not in df.columns:
        df = prepare_crc_diff_groups(df)
    
    groups = df['group'].unique()
    print(f"Found {len(groups)} groups: {groups}")
    
    # Step 1: Calculate mean values for species and metabolites in each group
    group_species_means = {}
    group_metabolite_means = {}
    
    for group in groups:
        group_data = df[df['group'] == group]
        
        if len(group_data) < 3:
            print(f"  Group {group} has insufficient samples ({len(group_data)}), skipping")
            continue
        
        # Extract species and metabolite data
        species_data = group_data[species_cols]
        metabolite_data = group_data[metabolite_cols]
        
        # Calculate mean for each feature
        group_species_means[group] = species_data.mean(axis=0)
        group_metabolite_means[group] = metabolite_data.mean(axis=0)
    
    # Step 2: Select common species and metabolites
    # Use fixed species if provided, otherwise select by ranking
    if fixed_species:
        # Match fixed species names with dataset columns (add tax_ prefix)
        common_species = []
        missing_species = []
        for species_name in fixed_species:
            matched_cols = [col for col in species_cols if species_name.lower() in col.lower()]
            if matched_cols:
                common_species.append(matched_cols[0])
            else:
                missing_species.append(species_name)
        
        if missing_species:
            raise ValueError(f"Species not found in dataset: {missing_species}")
        
        print(f"\nUsing fixed species: {[s.replace('tax_', '') for s in common_species]}")
    else:
        # Calculate average ranking of species across all groups
        all_species_ranks = pd.DataFrame()
        for group, means in group_species_means.items():
            ranked = means.rank(ascending=False, method='first')
            all_species_ranks[group] = ranked
        
        avg_species_ranks = all_species_ranks.mean(axis=1)
        common_species = avg_species_ranks.nsmallest(top_n_species).index.tolist()
    
    # Same method for metabolites
    # Use fixed metabolites if provided, otherwise select by ranking
    if fixed_metabolites:
        # Match fixed metabolite names with dataset columns (add met_ prefix)
        common_metabolites = []
        missing_metabolites = []
        for metabolite_name in fixed_metabolites:
            # Try exact match first, then partial match
            matched_cols = [col for col in metabolite_cols if metabolite_name.lower() in col.lower()]
            if matched_cols:
                common_metabolites.append(matched_cols[0])
            else:
                missing_metabolites.append(metabolite_name)
        
        if missing_metabolites:
            raise ValueError(f"Metabolites not found in dataset: {missing_metabolites}")
        
        print(f"\nUsing fixed metabolites: {[m.replace('met_', '') for m in common_metabolites]}")
    else:
        all_metabolite_ranks = pd.DataFrame()
        for group, means in group_metabolite_means.items():
            ranked = means.rank(ascending=False, method='first')
            all_metabolite_ranks[group] = ranked
        
        avg_metabolite_ranks = all_metabolite_ranks.mean(axis=1)
        common_metabolites = avg_metabolite_ranks.nsmallest(top_n_metabolites).index.tolist()
    
    print(f"\nCommonly selected species (top {top_n_species}): {[s.replace('tax_', '')[:20] for s in common_species]}")
    print(f"Commonly selected metabolites (top {top_n_metabolites}): {[m.replace('met_', '')[:20] for m in common_metabolites]}")
    
    # Step 3: Calculate correlation matrix and p-value matrix for each group
    correlations_by_group = {}
    pvalues_by_group = {}
    
    for group in groups:
        print(f"\nProcessing group: {group}")
        group_data = df[df['group'] == group]
        
        if len(group_data) < 3:
            print(f"  Group {group} has insufficient samples ({len(group_data)}), skipping correlation calculation")
            continue
        
        # Extract selected feature data
        selected_species = group_data[common_species]
        selected_metabolites = group_data[common_metabolites]
        
        # Create correlation matrix and p-value matrix
        correlation_matrix = pd.DataFrame(
            index=[s.replace('tax_', '')[:30] for s in common_species],
            columns=[m.replace('met_', '')[:30] for m in common_metabolites]
        )
        
        pvalue_matrix = pd.DataFrame(
            index=[s.replace('tax_', '')[:30] for s in common_species],
            columns=[m.replace('met_', '')[:30] for m in common_metabolites]
        )
        
        # Calculate Spearman or Pearson correlation
        for i, species_col in enumerate(common_species):
            species_values = selected_species[species_col].values
            # Ensure species values are numeric
            if species_values.dtype == object:
                try:
                    species_values = pd.to_numeric(species_values, errors='coerce')
                except:
                    print(f"  Warning: Cannot convert species {species_col} to numeric, skipping")
                    continue
                    
            for j, metabolite_col in enumerate(common_metabolites):
                metabolite_values = selected_metabolites[metabolite_col].values
                # Ensure metabolite values are numeric
                if metabolite_values.dtype == object:
                    try:
                        metabolite_values = pd.to_numeric(metabolite_values, errors='coerce')
                    except:
                        print(f"  Warning: Cannot convert metabolite {metabolite_col} to numeric, skipping")
                        correlation_matrix.iloc[i, j] = np.nan
                        pvalue_matrix.iloc[i, j] = np.nan
                        continue
                
                try:
                    # Check for NaN values
                    species_nan_count = np.isnan(species_values).sum()
                    metabolite_nan_count = np.isnan(metabolite_values).sum()
                    
                    if species_nan_count > 0 or metabolite_nan_count > 0:
                        raise ValueError(f"NaN values found: species '{species_col}' has {species_nan_count} NaNs, metabolite '{metabolite_col}' has {metabolite_nan_count} NaNs")
                    
                    # Check for constant values (no variance)
                    if np.std(species_values) == 0:
                        raise ValueError(f"Species '{species_col}' has constant values (no variance)")
                    if np.std(metabolite_values) == 0:
                        raise ValueError(f"Metabolite '{metabolite_col}' has constant values (no variance)")
                    
                    if method == 'spearman':
                        corr, p_value = stats.spearmanr(species_values, metabolite_values, nan_policy='omit')
                    else:  # pearson
                        corr, p_value = stats.pearsonr(species_values, metabolite_values)
                    
                    # Check if result is NaN
                    if np.isnan(corr) or np.isnan(p_value):
                        raise ValueError(f"Calculation returned NaN: corr={corr}, p_value={p_value}")
                    
                    # Store correlation coefficient and p-value
                    correlation_matrix.iloc[i, j] = corr
                    pvalue_matrix.iloc[i, j] = p_value
                except Exception as e:
                    raise ValueError(f"Error calculating correlation between '{species_col}' and '{metabolite_col}': {e}")
        
        # Convert to numeric type
        correlation_matrix = correlation_matrix.astype(float)
        pvalue_matrix = pvalue_matrix.astype(float)
        
        # Store results
        correlations_by_group[group] = correlation_matrix
        pvalues_by_group[group] = pvalue_matrix
        
        print(f"  Group {group} correlation matrix shape: {correlation_matrix.shape}")
        print(f"  Correlation range: [{correlation_matrix.min().min():.3f}, {correlation_matrix.max().max():.3f}]")
        print(f"  Significant correlations (p < 0.05): {(pvalue_matrix < 0.05).sum().sum()}")
    
    return correlations_by_group, pvalues_by_group, common_species, common_metabolites


def visualize_unified_correlations_heatmap(
    correlations_by_group: Dict[str, pd.DataFrame],
    pvalues_by_group: Dict[str, pd.DataFrame],
    group_order: List[str] = None,
    figsize: Tuple[int, int] = (14, 16),
    cmap: str = 'RdBu_r',
    vmin: float = -1,
    vmax: float = 1,
    title: str = "Species-Metabolite Correlation Heatmap",
    output_file: str = None
):
    """
    Visualize three-group correlation matrix heatmap using common species and metabolites with asterisk annotations
    
    Args:
        correlations_by_group: Dictionary containing correlation matrices for each group
        pvalues_by_group: Dictionary containing p-value matrices for each group
        group_order: Order of group display (if not provided, uses dictionary order)
        figsize: Overall figure size
        cmap: Color map
        vmin: Color bar minimum value
        vmax: Color bar maximum value
        title: Overall title
        output_file: Output file path (saves image if provided)
    """
    if not correlations_by_group or not pvalues_by_group:
        print("No correlation data available")
        return
    
    # Determine group order - 1x3 layout for CTRL, CRC_welldiff, CRC_poordiff
    if group_order is None:
        # Order for 1x3 layout: CTRL, CRC_welldiff (中-高分化), CRC_poordiff (低分化)
        group_order = ['CTRL', 'CRC_welldiff', 'CRC_poordiff']
    
    # Ensure all groups have matrices with same dimensions
    first_group = group_order[0]
    first_corr = correlations_by_group[first_group]
    n_species = len(first_corr.index)
    n_metabolites = len(first_corr.columns)
    
    # Verify all matrix dimensions are consistent
    for group in group_order:
        corr_matrix = correlations_by_group[group]
        pvalue_matrix = pvalues_by_group[group]
        if len(corr_matrix.index) != n_species or len(corr_matrix.columns) != n_metabolites:
            raise ValueError(f"Matrix dimensions for group {group} are inconsistent with other groups")
    
    # Create figure with tighter layout
    fig = plt.figure(figsize=figsize)
    
    # Create 3x1 grid of subplots with reduced spacing (vertical layout)
    gs = plt.GridSpec(3, 1, figure=fig, wspace=0.02, hspace=0.08,
                      left=0.18, right=0.78, bottom=0.08, top=0.92)
    
    # Create axes for each subplot
    axes = []
    for i in range(3):
        ax = fig.add_subplot(gs[i, 0])
        axes.append(ax)
    axes = np.array(axes)
    
    # Define asterisk annotation function
    def get_star_symbol(p_value):
        if p_value < 0.001:
            return '***'
        elif p_value < 0.01:
            return '**'
        elif p_value < 0.05:
            return '*'
        else:
            return ''
    
    # Create list of metabolite names for x-axis
    metabolite_names = correlations_by_group[first_group].columns.tolist()
    species_names = correlations_by_group[first_group].index.tolist()
    species_names = [s.replace('tax_s__', '')[:25] for s in species_names]
    
    # Define colors for labels
    label_colors = {
        'CTRL': '#2E86AB',  # Green
        'CRC_welldiff': '#A23B72',  # Yellow
        'CRC_poordiff': '#F18F01'   # Red
    }
    
    # Plot heatmap for each group
    for idx, group_name in enumerate(group_order):
        row = idx
        ax = axes[row]
        corr_matrix = correlations_by_group[group_name]
        pvalue_matrix = pvalues_by_group[group_name]
        
        # Plot heatmap
        im = ax.imshow(corr_matrix.values, cmap=cmap, aspect='auto', vmin=vmin, vmax=vmax)
        
        # 设置刻度但不显示标签（内部热图）
        ax.set_xticks(range(n_metabolites))
        ax.set_yticks(range(n_species))
        
        # 隐藏所有内部热图的坐标轴标签
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        
        # 添加刻度线但不显示刻度值
        ax.tick_params(axis='both', which='both', length=0)

        ax.axis('off')
        
        # Add grid lines
        ax.set_xticks(np.arange(-0.5, n_metabolites, 1), minor=True)
        ax.set_yticks(np.arange(-0.5, n_species, 1), minor=True)
        ax.grid(which="minor", color="gray", linestyle='-', linewidth=0.1)
        ax.tick_params(which="minor", size=0)
        
        # Add asterisk annotations for significant correlations
        for i in range(n_species):
            for j in range(n_metabolites):
                value = corr_matrix.iloc[i, j]
                p_value = pvalue_matrix.iloc[i, j]
                
                # Skip NaN values
                if np.isnan(value) or np.isnan(p_value):
                    continue
                
                # Get asterisk symbol
                star_symbol = get_star_symbol(p_value)
                
                # Add asterisk if significant
                if star_symbol:
                    # Determine text color
                    color = 'white' if abs(value) > 0.6 else 'black'
                    
                    # Add asterisk annotation
                    ax.text(j, i, star_symbol, 
                           ha='center', va='center', 
                           color=color, fontsize=6, fontweight='bold')
    
    # Add shared colorbar - 调整位置
    cbar_ax = fig.add_axes([0.82, 0.25, 0.015, 0.5])  # [left, bottom, width, height]
    cbar = fig.colorbar(im, cax=cbar_ax)
    cbar.set_label('Correlation Coefficient', rotation=270, labelpad=15, fontsize=11)
    
    # Add overall title
    fig.suptitle(title, fontsize=14, fontweight='bold', y=0.96)
    
    # ================== 添加右侧纵向标签 ==================
    # 获取子图的边界框
    ax0_bbox = axes[0].get_position()
    ax1_bbox = axes[1].get_position()
    ax2_bbox = axes[2].get_position()
    
    # CTRL 标签 - 右侧竖向
    ctrl_y_center = (ax0_bbox.y0 + ax0_bbox.y1) / 2
    fig.text(ax2_bbox.x1 + 0.01, ctrl_y_center, 'CTRL', 
             ha='left', va='center', rotation='vertical', fontsize=11, fontweight='bold')
    
    # 添加背景色块
    rect_ctrl = plt.Rectangle(
        (ax2_bbox.x1 + 0.005, ax0_bbox.y0), 0.015, ax0_bbox.height,
        facecolor=label_colors['CTRL'], edgecolor='black', linewidth=0,
        transform=fig.transFigure, clip_on=False, zorder=0
    )
    fig.add_artist(rect_ctrl)
    
    # CRC_welldiff (中-高分化) 标签 - 右侧竖向
    welldiff_y_center = (ax1_bbox.y0 + ax1_bbox.y1) / 2
    fig.text(ax2_bbox.x1 + 0.01, welldiff_y_center, 'CRC WELLDIFF', 
             ha='left', va='center', rotation='vertical', fontsize=11, fontweight='bold')
    
    # 添加背景色块
    rect_welldiff = plt.Rectangle(
        (ax2_bbox.x1 + 0.005, ax1_bbox.y0), 0.015, ax1_bbox.height,
        facecolor=label_colors['CRC_welldiff'], edgecolor='black', linewidth=0,
        transform=fig.transFigure, clip_on=False, zorder=0
    )
    fig.add_artist(rect_welldiff)
    
    # CRC_poordiff (低分化) 标签 - 右侧竖向
    poordiff_y_center = (ax2_bbox.y0 + ax2_bbox.y1) / 2
    fig.text(ax2_bbox.x1 + 0.01, poordiff_y_center, 'CRC POORDIFF', 
             ha='left', va='center', rotation='vertical', fontsize=11, fontweight='bold')
    
    # 添加背景色块
    rect_poordiff = plt.Rectangle(
        (ax2_bbox.x1 + 0.005, ax2_bbox.y0), 0.015, ax2_bbox.height,
        facecolor=label_colors['CRC_poordiff'], edgecolor='black', linewidth=0,
        transform=fig.transFigure, clip_on=False, zorder=0
    )
    fig.add_artist(rect_poordiff)
    
    # ================== 添加左侧纵向坐标轴标签（物种名称）==================
    # 在最左侧的热图显示y轴标签（顶部热图）
    for i, species_name in enumerate(species_names):
        y_pos = ax0_bbox.y0 + (ax0_bbox.height * (i + 0.4) / n_species)
        fig.text(ax0_bbox.x0 - 0.01, y_pos, species_name, 
                 ha='right', va='center', fontsize=9, rotation=0)
    
    # ================== 添加底部横向坐标轴标签（代谢物名称）==================
    # 在底部热图显示x轴标签
    for j, metabolite_name in enumerate(metabolite_names):
        x_pos = ax2_bbox.x0 + (ax2_bbox.width * (j + 0.5) / n_metabolites)
        fig.text(x_pos, ax2_bbox.y0 - 0.01, metabolite_name, 
                 ha='center', va='top', fontsize=9, rotation=90)
    
    # ================== 添加图例说明 ==================
    legend_text = "Significance levels: *p < 0.05, **p < 0.01, ***p < 0.001"
    fig.text(0.5, 0.03, legend_text, ha='center', fontsize=10, style='italic', transform=fig.transFigure)
    
    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"Three-group comparison heatmap saved to: {output_file}")
    
    # plt.show()
    
    # Return statistics
    stats_info = {}
    for group_name in group_order:
        corr_matrix = correlations_by_group[group_name]
        pvalue_matrix = pvalues_by_group[group_name]
        
        stats_info[group_name] = {
            'matrix_shape': corr_matrix.shape,
            'mean_correlation': corr_matrix.values[~np.isnan(corr_matrix.values)].mean(),
            'median_correlation': np.median(corr_matrix.values[~np.isnan(corr_matrix.values)]),
            'min_correlation': corr_matrix.values[~np.isnan(corr_matrix.values)].min(),
            'max_correlation': corr_matrix.values[~np.isnan(corr_matrix.values)].max(),
            'strong_pos_corr': ((corr_matrix.values > 0.5) & (~np.isnan(corr_matrix.values))).sum(),
            'strong_neg_corr': ((corr_matrix.values < -0.5) & (~np.isnan(corr_matrix.values))).sum(),
            'sig_corr_p_05': (pvalue_matrix < 0.05).sum().sum(),
            'sig_corr_p_01': (pvalue_matrix < 0.01).sum().sum(),
            'sig_corr_p_001': (pvalue_matrix < 0.001).sum().sum(),
        }
    
    return stats_info

def analyze_unified_group_correlations(
    df: pd.DataFrame,
    top_n_species: int = 10,
    top_n_metabolites: int = 20,
    method: str = 'spearman',
    output_prefix: str = 'unified_group_correlations',
    fixed_species: Optional[List[str]] = None,
    fixed_metabolites: Optional[List[str]] = None
):
    """
    Complete analysis of group correlations with visualization (using common species and metabolites)
    
    Args:
        df: Input DataFrame
        top_n_species: Number of top species to select
        top_n_metabolites: Number of top metabolites to select
        method: Correlation calculation method
        output_prefix: Output file prefix
        fixed_species: Optional list of fixed species names to use
        fixed_metabolites: Optional list of fixed metabolite names to use
    """
    print("="*80)
    print("Starting unified group correlation analysis (using common species and metabolites)")
    print("="*80)
    
    # 1. Calculate unified correlations (including p-values)
    correlations, pvalues, common_species, common_metabolites = calculate_unified_group_correlations(
        df, 
        top_n_species=top_n_species, 
        top_n_metabolites=top_n_metabolites,
        method=method,
        fixed_species=fixed_species,
        fixed_metabolites=fixed_metabolites
    )
    
    if not correlations:
        print("No valid correlation matrices calculated")
        return None, None, None, None
    
    # 2. Visualize three-group comparison heatmap (with asterisk annotations)
    print("\n" + "="*80)
    print("Generating three-group comparison heatmap (with asterisk annotations)")
    print("="*80)
    
    # Define group order: CTRL, CRC_welldiff, CRC_poordiff
    group_order = ['CTRL', 'CRC_welldiff', 'CRC_poordiff']
    
    # Only keep existing groups
    group_order = [g for g in group_order if g in correlations]
    
    stats_info = visualize_unified_correlations_heatmap(
        correlations,
        pvalues,
        group_order=group_order,
        figsize=(14, 16),
        cmap='RdBu_r',
        title=f"{method.capitalize()} Correlation - CRC Differentiation Comparison",
        output_file=f"{output_prefix}_three_groups_heatmap.png"
    )
    
    # 3. Output statistical information
    print("\n" + "="*80)
    print("Correlation Statistics")
    print("="*80)
    
    for group_name in group_order:
        if group_name in stats_info:
            stats = stats_info[group_name]
            print(f"\nGroup: {group_name}")
            print(f"  Matrix shape: {stats['matrix_shape']}")
            print(f"  Mean correlation: {stats['mean_correlation']:.3f}")
            print(f"  Median correlation: {stats['median_correlation']:.3f}")
            print(f"  Min correlation: {stats['min_correlation']:.3f}")
            print(f"  Max correlation: {stats['max_correlation']:.3f}")
            print(f"  Strong positive correlations (>0.5): {stats['strong_pos_corr']}")
            print(f"  Strong negative correlations (<-0.5): {stats['strong_neg_corr']}")
            print(f"  Significant correlations (p < 0.05): {stats['sig_corr_p_05']}")
            print(f"  Significant correlations (p < 0.01): {stats['sig_corr_p_01']}")
            print(f"  Significant correlations (p < 0.001): {stats['sig_corr_p_001']}")
    
    # 4. Save correlation matrices and p-value matrices to CSV files
    print("\n" + "="*80)
    print("Saving correlation matrices and p-value matrices to CSV files")
    print("="*80)
    
    for group_name in group_order:
        if group_name in correlations:
            # Save correlation matrix
            corr_filename = f"{output_prefix}_{group_name}_correlations.csv"
            correlations[group_name].to_csv(corr_filename)
            print(f"  Group {group_name} correlation matrix saved to: {corr_filename}")
            
            # Save p-value matrix
            pval_filename = f"{output_prefix}_{group_name}_pvalues.csv"
            pvalues[group_name].to_csv(pval_filename)
            print(f"  Group {group_name} p-value matrix saved to: {pval_filename}")
    
    # 5. Identify strongest correlation pairs (considering significance)
    print("\n" + "="*80)
    print("Identifying Strongest Significant Correlation Pairs (Top 10)")
    print("="*80)
    
    for group_name in group_order:
        if group_name in correlations and group_name in pvalues:
            print(f"\nGroup: {group_name}")
            
            # Flatten correlation matrix and p-value matrix
            corr_pairs = []
            corr_matrix = correlations[group_name]
            pvalue_matrix = pvalues[group_name]
            
            for i, species in enumerate(corr_matrix.index):
                for j, metabolite in enumerate(corr_matrix.columns):
                    corr_value = corr_matrix.iloc[i, j]
                    p_value = pvalue_matrix.iloc[i, j]
                    
                    # Skip NaN values
                    if np.isnan(corr_value) or np.isnan(p_value):
                        continue
                    
                    corr_pairs.append({
                        'species': species,
                        'metabolite': metabolite,
                        'correlation': corr_value,
                        'abs_correlation': abs(corr_value),
                        'p_value': p_value,
                        'significant': p_value < 0.05
                    })
            
            # Sort by absolute correlation value, considering only significant correlations
            corr_pairs_df = pd.DataFrame(corr_pairs)
            significant_pairs = corr_pairs_df[corr_pairs_df['significant']]
            
            if len(significant_pairs) > 0:
                top_pos = significant_pairs.sort_values('correlation', ascending=False).head(5)
                top_neg = significant_pairs.sort_values('correlation', ascending=True).head(5)
                
                print("  Strongest significant positive correlations:")
                for _, row in top_pos.iterrows():
                    star_symbol = '***' if row['p_value'] < 0.001 else ('**' if row['p_value'] < 0.01 else '*')
                    print(f"    {row['species']} ↔ {row['metabolite']}: {row['correlation']:.3f} {star_symbol}")
                
                print("  Strongest significant negative correlations:")
                for _, row in top_neg.iterrows():
                    star_symbol = '***' if row['p_value'] < 0.001 else ('**' if row['p_value'] < 0.01 else '*')
                    print(f"    {row['species']} ↔ {row['metabolite']}: {row['correlation']:.3f} {star_symbol}")
            else:
                print("  No significant correlations found (p < 0.05)")
    
    # 6. Save common species and metabolites lists
    print("\n" + "="*80)
    print("Saving Common Species and Metabolites Lists")
    print("="*80)
    
    # Save species list
    species_df = pd.DataFrame({
        'original_name': common_species,
        'display_name': [s.replace('tax_', '')[:30] for s in common_species]
    })
    species_filename = f"{output_prefix}_common_species.csv"
    species_df.to_csv(species_filename, index=False)
    print(f"  Common species list saved to: {species_filename}")
    
    # Save metabolites list
    metabolites_df = pd.DataFrame({
        'original_name': common_metabolites,
        'display_name': [m.replace('met_', '')[:30] for m in common_metabolites]
    })
    metabolites_filename = f"{output_prefix}_common_metabolites.csv"
    metabolites_df.to_csv(metabolites_filename, index=False)
    print(f"  Common metabolites list saved to: {metabolites_filename}")
    
    print("\n" + "="*80)
    print("Analysis Complete!")
    print("="*80)
    
    return correlations, pvalues, common_species, common_metabolites


# Usage example
if __name__ == "__main__":
    ds = BioSmokeDataset()
    ds.preprocess_taxonomy_data(remove_low_expression = False)
    ds.preprocess_metabolomics_data(relative_abund=False)
    merged = ds.intersection_to_dataframe()

    out_dir = Path('results/heatmap_plots')
    out_dir.mkdir(parents=True, exist_ok=True)
    
    # Execute unified correlation analysis
    fixed_species_list = [
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
    
    fixed_metabolites_list = [
        'SQDG 26:2; SQDG(13:1/13:1)',
        'Cytosine',
        'Perfluorooctanesulfonic acid',
        'Methyl dihydrojasmonate',
        'Pyrogallol-2-O-sulphate',
        "5'-(3',4'-Dihydroxyphenyl)-gamma-valerolactone sulfate",
        '2-Hydroxy-4,7-dimethoxy-2H-1,4-benzoxazin-3(4H)-one',
        'trans-3,5-Dimethoxy-4-hydroxycinnamaldehyde',
        '(R)-3-Hydroxy-5-phenylpentanoic acid',
        'N-Methyl-D-glucamine',
        'Chenodeoxycholic acid sulfate',
        'Creatinine',
        'Lucidenic acid F',
        'Demissidine',
        'Alpha-Hydroxyisobutyric acid',
        'Pyrocatechol',
        'Gentisic acid',
        'D-Galacturonic acid',
        '1,3-Dimethyluric acid',
        '4-Hydroxy-5-(phenyl)-valeric acid-O-sulphate'
    ]
    
    correlations, pvalues, common_species, common_metabolites = analyze_unified_group_correlations(
        merged,
        top_n_species=10,      # Select top 10 species
        top_n_metabolites=20,  # Select top 20 metabolites
        method='spearman',     # Use Spearman correlation
        output_prefix=str(out_dir / 'unified_group_correlations'),
        fixed_species=fixed_species_list,
        fixed_metabolites=fixed_metabolites_list
    )