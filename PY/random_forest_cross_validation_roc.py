from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_curve, auc, confusion_matrix
from sklearn.model_selection import train_test_split

from dataset import BioSmokeDataset
from split_group import prepare_crc_diff_groups

GROUP_MAP = {
    'CTRL': 'CTRL',
    'CRC_welldiff': 'CRC-Well',
    'CRC_poordiff': 'CRC-Poor',
}

FIXED_SPECIES_LIST = [
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

FIXED_METABOLITES_LIST = [
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

# Optional switch: if True, only fixed species/metabolites are used as biomarkers.
USE_FIXED_BIOMARKERS = True


def match_fixed_features(columns: list[str], target_list: list[str], prefix: str) -> list[str]:
    matched = []
    missing = []

    for target in target_list:
        # Match by substring to tolerate source column naming details.
        matches = [column for column in columns if column.startswith(prefix) and target.lower() in column.lower()]
        if matches:
            matched.append(matches[0])
        else:
            missing.append(target)

    if missing:
        raise ValueError(f'以下 biomarker 在数据中未找到: {missing}')

    return matched

def load_and_prepare_data(use_fixed_biomarkers: bool = USE_FIXED_BIOMARKERS) -> pd.DataFrame:
    ds = BioSmokeDataset()
    ds.preprocess_taxonomy_data(transform=True)
    ds.preprocess_metabolomics_data(relative_abund=False, transform=True)
    merged = ds.intersection_to_dataframe()

    grouped = prepare_crc_diff_groups(merged)
    grouped['group'] = grouped['group'].map(GROUP_MAP)
    grouped = grouped[grouped['group'].isin(['CTRL', 'CRC-Well', 'CRC-Poor'])].copy()
    
    # Feature extraction
    feature_cols = [c for c in grouped.columns if c.startswith('tax_') or c.startswith('met_')]

    if use_fixed_biomarkers:
        matched_species = match_fixed_features(feature_cols, FIXED_SPECIES_LIST, 'tax_')
        matched_metabolites = match_fixed_features(feature_cols, FIXED_METABOLITES_LIST, 'met_')
        feature_cols = matched_species + matched_metabolites

    x = grouped[feature_cols].apply(pd.to_numeric, errors='coerce').fillna(0.0)
    x = x.loc[:, x.var(axis=0) > 0].copy()
    feature_names = x.columns.tolist()
    grouped = pd.concat([grouped[['group']], x], axis=1)
    if use_fixed_biomarkers:
        print(f'Using fixed biomarkers only: {len(feature_names)} features')
    return grouped, feature_names

def get_subset(df: pd.DataFrame, groups: list, positive_class: str = 'CRC') -> tuple[pd.DataFrame, pd.Series]:
    """Retrieve features and binary labels (CTRL=0, other=1) for selected groups."""
    subset = df[df['group'].isin(groups)].copy()
    y = (subset['group'] != 'CTRL').astype(int)
    x = subset.drop(columns=['group'])
    return x, y

def plot_cross_val_confusion_matrices(model, model_name, tasks, test_df, out_path):
    # Set global plotting standards for this figure
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans']
    
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))
    fig.suptitle(f'Validation Confusion Matrices (Model: {model_name})', fontsize=18, fontweight='bold', y=1.05)
    
    for ax, (val_name, groups) in zip(axes, tasks.items()):
        x_test, y_test = get_subset(test_df, groups)
        y_pred = model.predict(x_test)
        cm = confusion_matrix(y_test, y_pred)
        
        # Plot with 0,0 at top-left and 1,1 at bottom-right (default behavior)
        sns.heatmap(cm, annot=True, fmt='d', cmap='Blues', 
                    xticklabels=['CTRL', 'CRC'], yticklabels=['CTRL', 'CRC'],
                    annot_kws={'size': 16, 'fontweight': 'bold'}, cbar=False, 
                    linewidths=1, linecolor='white', ax=ax)
        
        ax.set_title(f'Test Set: {val_name}', fontsize=15, fontweight='bold', pad=12)
        ax.set_ylabel('True Category', fontsize=13, fontweight='bold')
        ax.set_xlabel('Predicted Category', fontsize=13, fontweight='bold')
        
        # Remove spines/outlines accurately
        for _, spine in ax.spines.items():
            spine.set_visible(False)
        ax.grid(False) # Ensure no grid lines
            
    plt.tight_layout()
    plt.savefig(out_path, dpi=300, bbox_inches='tight', transparent=False)
    plt.close()

def plot_importance_heatmap(models, feature_names, out_path, top_n=30):
    plt.rcParams['font.family'] = 'sans-serif'
    
    imp_dict = {}
    for task_name, model in models.items():
        imp_dict[task_name] = model.feature_importances_
        
    imp_df = pd.DataFrame(imp_dict, index=feature_names)
    
    imp_df['mean_imp'] = imp_df.mean(axis=1)
    imp_df = imp_df.sort_values(by='mean_imp', ascending=False)
    
    top_imp_df = imp_df.head(top_n).drop(columns=['mean_imp'])
    
    biomarker_types = ['species' if feature.startswith('tax_') else 'metabolite' for feature in top_imp_df.index]
    clean_labels = [
        f.replace('tax_', '').replace('met_', '').replace('s__', '')[:40]
        for f in top_imp_df.index
    ]
    top_imp_df.index = clean_labels

    # Horizontal layout: biomarkers on x-axis, models on y-axis.
    top_imp_df = top_imp_df.T
    fig, ax = plt.subplots(figsize=(max(10, top_n * 0.45), 4.2))
    
    # Use a professional, linear sequential colormap e.g., 'rocket_r' or 'mako_r' from Seaborn
    sns.heatmap(top_imp_df, cmap='flare', linewidths=0.5, linecolor='white',
                cbar_kws={'shrink': 0.8}, ax=ax)
    
    ax.set_title('Top Biomarker Importances Across Models', fontsize=16, fontweight='bold', pad=18)
    
    plt.yticks(rotation=0, fontsize=11, fontweight='bold')
    plt.xticks(rotation=45, ha='right', fontsize=10)

    # Species names stay italic; metabolites and model labels remain normal.
    for label, biomarker_type in zip(ax.get_xticklabels(), biomarker_types):
        label.set_fontstyle('italic' if biomarker_type == 'species' else 'normal')
    for label in ax.get_yticklabels():
        label.set_fontstyle('normal')
    
    # Remove spines for the heatmap and ensure grid is off
    sns.despine(ax=ax, left=True, bottom=True, top=True, right=True)
    ax.grid(False)
        
    plt.tight_layout()
    plt.savefig(out_path, dpi=300, bbox_inches='tight')
    plt.close()

def main():
    sns.set_theme(style='whitegrid', context='talk')
    out_dir = Path(__file__).resolve().parents[1] / 'results' / 'random_forest_cross_val_roc'
    out_dir.mkdir(parents=True, exist_ok=True)

    df, feature_names = load_and_prepare_data()

    # Split the original dataset to fix the Train and Test pools for all combinations
    # This prevents data leakage (same sample being in A's train and B's test)
    train_df, test_df = train_test_split(df, test_size=0.3, random_state=42, stratify=df['group'])
    
    tasks = {
        'CTRL_vs_CRC-Well': ['CTRL', 'CRC-Well'],
        'CTRL_vs_CRC-Poor': ['CTRL', 'CRC-Poor'],
        'CTRL_vs_CRC_All': ['CTRL', 'CRC-Well', 'CRC-Poor']
    }
    
    # Prepare Train Models
    models = {}
    results_summary = []

    for task_name, groups in tasks.items():
        x_train, y_train = get_subset(train_df, groups)
        
        forest = RandomForestClassifier(
            n_estimators=1000,
            random_state=42,
            class_weight='balanced_subsample',
            min_samples_leaf=2,
            n_jobs=-1
        )
        forest.fit(x_train, y_train)
        models[task_name] = forest
        print(f"Trained {task_name} (Train size: {len(x_train)})")

    # Final Summary for ChatGPT analysis
    summary_data = []

    # ==== Define Top-Tier Journal Colors ====
    journal_colors = ['#E64B35', '#4DBBD5', '#00A087']  # NPG (Nature Publishing Group) style

    # Plot ROC for each model using all three validation sets
    for model_name, model in models.items():
        fig, ax = plt.subplots(figsize=(7, 6))
        
        # Save model feature importances
        importance = model.feature_importances_
        top_indices = np.argsort(importance)[::-1][:50] # Top 50 features
        
        for idx, (val_name, groups) in enumerate(tasks.items()):
            x_test, y_test = get_subset(test_df, groups)
            
            # Predict
            y_proba = model.predict_proba(x_test)[:, 1]
            fpr, tpr, _ = roc_curve(y_test, y_proba)
            roc_auc = auc(fpr, tpr)
            
            # Predict labels for accuracy
            y_pred = model.predict(x_test)
            accuracy = (y_pred == y_test).mean()

            # Record metrics for summary CSV
            for rank, feat_idx in enumerate(top_indices):
                summary_data.append({
                    'Trained_Model': model_name,
                    'Validation_Set': val_name,
                    'AUC': roc_auc,
                    'Accuracy': accuracy,
                    'Feature_Rank': rank + 1,
                    'Feature_Name': feature_names[feat_idx],
                    'Importance_Score': importance[feat_idx]
                })

            ax.plot(fpr, tpr, lw=2.5, color=journal_colors[idx % len(journal_colors)], 
                    label=f'{val_name} (AUC = {roc_auc:.3f})')
            
        ax.plot([0, 1], [0, 1], color='#7E7E7E', lw=2, linestyle=(0, (5, 5)))
        ax.set_xlim([-0.02, 1.02])
        ax.set_ylim([-0.02, 1.05])
        ax.set_xlabel('False Positive Rate', fontsize=14, fontweight='bold')
        ax.set_ylabel('True Positive Rate', fontsize=14, fontweight='bold')
        ax.set_title(f'ROC - Trained on {model_name}', fontsize=16, fontweight='bold', pad=15)
        ax.tick_params(axis='both', which='major', labelsize=12)
        
        # Clean spines (Nature/Cell style)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_linewidth(1.5)
        ax.spines['bottom'].set_linewidth(1.5)
        
        # Grid lines
        ax.grid(True, linestyle='--', alpha=0.3, color='#B0B0B0')
        
        plt.legend(loc="lower right", fontsize=11, frameon=False)
        
        output_path = out_dir / f'ROC_{model_name}.png'
        plt.savefig(output_path, dpi=300, bbox_inches='tight', transparent=False)
        plt.close()
        print(f"Saved plot: {output_path}")

        # ==== New: Plot Confusion Matrix Grid ====
        cm_path = out_dir / f'CM_{model_name}_across_all_vals.png'
        plot_cross_val_confusion_matrices(model, model_name, tasks, test_df, cm_path)
        print(f"Saved CM Grid for {model_name}")

    # ==== New: Plot Unified Biomarker Importance Heatmap ====
    # The feature names are the same across all models, so we can extract it from the last processed subset
    unified_heatmap_path = out_dir / 'Biomarker_Importance_Unified_Heatmap.png'
    plot_importance_heatmap(models, feature_names, unified_heatmap_path, top_n=30)
    print(f"Saved Unified Biomarker Importance Heatmap: {unified_heatmap_path}")

    # ==== New: Save CSV Summary for AI Analysis ====
    summary_df = pd.DataFrame(summary_data)
    summary_csv_path = out_dir / 'RF_Cross_Validation_Summary.csv'
    summary_df.to_csv(summary_csv_path, index=False)
    print(f"Saved AI-Ready Summary CSV: {summary_csv_path}")

if __name__ == '__main__':
    main()
