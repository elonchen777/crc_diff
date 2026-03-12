from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.metrics import auc, classification_report, confusion_matrix, roc_curve
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.preprocessing import label_binarize

from dataset import BioSmokeDataset
from split_group import prepare_crc_diff_groups


GROUP_MAP = {
    'CTRL': 'CTRL',
    'CRC_welldiff': 'WELL',
    'CRC_poordiff': 'POOR',
}

LABEL_MAP = {
    'CTRL': 0,
    'WELL': 1,
    'POOR': 2,
}

CLASS_ORDER = ['CTRL', 'WELL', 'POOR']
CLASS_COLORS = {
    'CTRL': '#2E86AB',
    'WELL': '#F4A261',
    'POOR': '#D1495B',
}


def clean_feature_name(feature_name: str) -> str:
    if feature_name.startswith('tax_'):
        return feature_name.replace('tax_', '', 1)
    if feature_name.startswith('met_'):
        return feature_name.replace('met_', '', 1)
    return feature_name


def feature_type(feature_name: str) -> str:
    if feature_name.startswith('tax_'):
        return 'Microbe'
    if feature_name.startswith('met_'):
        return 'Metabolite'
    return 'Other'


def load_analysis_data() -> pd.DataFrame:
    ds = BioSmokeDataset()
    ds.preprocess_taxonomy_data(transform=True)
    ds.preprocess_metabolomics_data(relative_abund=False, transform=True)
    merged = ds.intersection_to_dataframe()

    grouped = prepare_crc_diff_groups(merged)
    grouped['group'] = grouped['group'].map(GROUP_MAP)
    grouped = grouped[grouped['group'].isin(CLASS_ORDER)].copy()

    return grouped


def prepare_features(df: pd.DataFrame) -> tuple[pd.DataFrame, pd.Series]:
    feature_cols = [
        column for column in df.columns
        if column.startswith('tax_') or column.startswith('met_')
    ]
    if not feature_cols:
        raise ValueError('未找到 tax_ 或 met_ 特征列。')

    x = df[feature_cols].apply(pd.to_numeric, errors='coerce').fillna(0.0)
    x = x.loc[:, x.var(axis=0) > 0].copy()
    y = df['group'].map(LABEL_MAP)

    if y.isna().any():
        raise ValueError('存在无法映射到标签编码的 group。')

    return x, y.astype(int)


def build_cv(y: pd.Series) -> StratifiedKFold:
    min_class_count = int(y.value_counts().min())
    if min_class_count < 2:
        raise ValueError('至少每组需要 2 个样本才能进行交叉验证。')

    n_splits = min(5, min_class_count)
    if n_splits < 3:
        n_splits = 2

    return StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=42)


def select_feature_count(feature_count: int) -> int:
    return max(20, min(30, feature_count)) if feature_count >= 20 else feature_count


def fit_models(x: pd.DataFrame, y: pd.Series):
    selected_feature_count = select_feature_count(x.shape[1])
    selector = SelectKBest(score_func=f_classif, k=selected_feature_count)
    forest = RandomForestClassifier(
        n_estimators=1000,
        random_state=42,
        class_weight='balanced_subsample',
        min_samples_leaf=2,
        n_jobs=-1,
    )

    cv = build_cv(y)
    x_selected = selector.fit_transform(x, y)

    oof_pred = cross_val_predict(forest, x_selected, y, cv=cv, method='predict', n_jobs=1)
    oof_proba = cross_val_predict(forest, x_selected, y, cv=cv, method='predict_proba', n_jobs=1)

    forest.fit(x_selected, y)

    selected_columns = x.columns[selector.get_support()].tolist()
    importance_df = pd.DataFrame(
        {
            'feature': selected_columns,
            'importance': forest.feature_importances_,
        }
    ).sort_values('importance', ascending=False)
    importance_df['feature_type'] = importance_df['feature'].map(feature_type)
    importance_df['display_name'] = importance_df['feature'].map(clean_feature_name)

    return {
        'selector': selector,
        'forest': forest,
        'cv': cv,
        'oof_pred': oof_pred,
        'oof_proba': oof_proba,
        'importance_df': importance_df,
        'selected_columns': selected_columns,
    }


def plot_feature_importance(importance_df: pd.DataFrame, output_dir: Path) -> pd.DataFrame:
    top20 = importance_df.head(20).copy()
    top20 = top20.iloc[::-1]

    color_map = {'Microbe': '#287271', 'Metabolite': '#E76F51', 'Other': '#666666'}
    colors = top20['feature_type'].map(color_map)

    plt.figure(figsize=(11, 8))
    plt.barh(top20['display_name'], top20['importance'], color=colors)
    plt.xlabel('Random Forest Importance')
    plt.ylabel('Feature')
    plt.title('Top 20 Feature Importance')
    plt.tight_layout()
    plt.savefig(output_dir / 'feature_importance_top20.png', dpi=300, bbox_inches='tight')
    plt.close()

    return top20.iloc[::-1]


def plot_roc_curves(y: pd.Series, oof_proba: np.ndarray, output_dir: Path) -> pd.DataFrame:
    y_binary = label_binarize(y, classes=[0, 1, 2])
    roc_rows = []

    plt.figure(figsize=(8, 7))
    for class_index, class_name in enumerate(CLASS_ORDER):
        fpr, tpr, _ = roc_curve(y_binary[:, class_index], oof_proba[:, class_index])
        roc_auc = auc(fpr, tpr)
        roc_rows.append({'class': class_name, 'auc': roc_auc})
        plt.plot(
            fpr,
            tpr,
            color=CLASS_COLORS[class_name],
            linewidth=2.2,
            label=f'{class_name} vs rest (AUC = {roc_auc:.3f})',
        )

    plt.plot([0, 1], [0, 1], linestyle='--', color='#888888', linewidth=1)
    plt.xlim(0, 1)
    plt.ylim(0, 1.02)
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('One-vs-Rest ROC Curves')
    plt.legend(loc='lower right')
    plt.tight_layout()
    plt.savefig(output_dir / 'roc_curves_ovr.png', dpi=300, bbox_inches='tight')
    plt.close()

    return pd.DataFrame(roc_rows)


def plot_confusion(y: pd.Series, oof_pred: np.ndarray, output_dir: Path) -> pd.DataFrame:
    cm = confusion_matrix(y, oof_pred, labels=[0, 1, 2])
    cm_df = pd.DataFrame(cm, index=CLASS_ORDER, columns=CLASS_ORDER)

    plt.figure(figsize=(7, 6))
    sns.heatmap(cm_df, annot=True, fmt='d', cmap='Blues', cbar=False)
    plt.xlabel('Predicted Label')
    plt.ylabel('True Label')
    plt.title('Confusion Matrix')
    plt.tight_layout()
    plt.savefig(output_dir / 'confusion_matrix.png', dpi=300, bbox_inches='tight')
    plt.close()

    return cm_df


def plot_biomarker_panel(biomarkers_df: pd.DataFrame, output_dir: Path) -> None:
    plot_df = biomarkers_df.iloc[::-1].copy()
    color_map = {'Microbe': '#287271', 'Metabolite': '#E76F51', 'Other': '#666666'}

    plt.figure(figsize=(10, 6.5))
    plt.barh(
        plot_df['display_name'],
        plot_df['importance'],
        color=plot_df['feature_type'].map(color_map),
    )
    plt.xlabel('Random Forest Importance')
    plt.ylabel('Biomarker')
    plt.title('Top 10 Microbe-Metabolite Biomarker Panel')
    plt.tight_layout()
    plt.savefig(output_dir / 'biomarker_panel_top10.png', dpi=300, bbox_inches='tight')
    plt.close()


def save_reports(
    output_dir: Path,
    df: pd.DataFrame,
    y: pd.Series,
    model_results: dict,
    roc_df: pd.DataFrame,
    confusion_df: pd.DataFrame,
) -> None:
    importance_df = model_results['importance_df']
    biomarkers_df = importance_df.head(10).copy()
    biomarkers_df['label_code'] = biomarkers_df['feature'].map(lambda _: '')

    top20_df = importance_df.head(20).copy()
    top20_df.to_csv(output_dir / 'feature_importance_top20.csv', index=False)
    importance_df.to_csv(output_dir / 'feature_importance_all_selected.csv', index=False)
    biomarkers_df.to_csv(output_dir / 'biomarker_panel_top10.csv', index=False)
    roc_df.to_csv(output_dir / 'roc_auc_summary.csv', index=False)
    confusion_df.to_csv(output_dir / 'confusion_matrix.csv')

    metadata_df = df[['group']].copy()
    metadata_df['label_code'] = y.values
    metadata_df['predicted_label_code'] = model_results['oof_pred']
    metadata_df['predicted_group'] = metadata_df['predicted_label_code'].map({v: k for k, v in LABEL_MAP.items()})
    proba_df = pd.DataFrame(
        model_results['oof_proba'],
        index=df.index,
        columns=[f'prob_{class_name}' for class_name in CLASS_ORDER],
    )
    pd.concat([metadata_df, proba_df], axis=1).to_csv(output_dir / 'cross_validated_predictions.csv')

    report_text = classification_report(
        y,
        model_results['oof_pred'],
        labels=[0, 1, 2],
        target_names=CLASS_ORDER,
        digits=4,
        zero_division=0,
    )
    accuracy = float((y.values == model_results['oof_pred']).mean())
    selected_feature_count = len(model_results['selected_columns'])

    summary_lines = [
        'Random Forest multi-class classification summary',
        '',
        f'Samples: {len(df)}',
        f'Selected features: {selected_feature_count}',
        f'Groups: {df["group"].value_counts().to_dict()}',
        f'Cross-validated accuracy: {accuracy:.4f}',
        '',
        'Label encoding:',
        'CTRL = 0',
        'WELL = 1',
        'POOR = 2',
        '',
        'Classification report:',
        report_text,
    ]
    (output_dir / 'random_forest_summary.txt').write_text('\n'.join(summary_lines), encoding='utf-8')


def main() -> None:
    sns.set_theme(style='whitegrid', context='talk')

    project_root = Path(__file__).resolve().parents[1]
    output_dir = project_root / 'results' / 'random_forest'
    output_dir.mkdir(parents=True, exist_ok=True)

    df = load_analysis_data()
    x, y = prepare_features(df)
    model_results = fit_models(x, y)

    plot_feature_importance(model_results['importance_df'], output_dir)
    roc_df = plot_roc_curves(y, model_results['oof_proba'], output_dir)
    confusion_df = plot_confusion(y, model_results['oof_pred'], output_dir)
    plot_biomarker_panel(model_results['importance_df'].head(10), output_dir)
    save_reports(output_dir, df, y, model_results, roc_df, confusion_df)

    print(f'随机森林结果已输出到: {output_dir}')
    print(f'使用特征数: {len(model_results["selected_columns"])}')
    print('标签编码: CTRL=0, WELL=1, POOR=2')


if __name__ == '__main__':
    main()