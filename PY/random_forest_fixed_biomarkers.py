from pathlib import Path

import pandas as pd

from random_forest_diff import (
    CLASS_ORDER,
    LABEL_MAP,
    clean_feature_name,
    load_analysis_data,
    plot_biomarker_panel,
    plot_confusion,
    plot_feature_importance,
    plot_roc_curves,
    save_reports,
)
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import cross_val_predict


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


def feature_type(feature_name: str) -> str:
    if feature_name.startswith('tax_'):
        return 'Microbe'
    if feature_name.startswith('met_'):
        return 'Metabolite'
    return 'Other'


def match_fixed_features(columns: list[str], target_list: list[str], prefix: str) -> list[str]:
    matched = []
    missing = []

    for target in target_list:
        matches = [column for column in columns if column.startswith(prefix) and target.lower() in column.lower()]
        if matches:
            matched.append(matches[0])
        else:
            missing.append(target)

    if missing:
        raise ValueError(f'以下 biomarker 在数据中未找到: {missing}')

    return matched


def prepare_fixed_biomarker_features(df: pd.DataFrame) -> tuple[pd.DataFrame, pd.Series, list[str]]:
    feature_columns = [column for column in df.columns if column.startswith('tax_') or column.startswith('met_')]
    matched_species = match_fixed_features(feature_columns, FIXED_SPECIES_LIST, 'tax_')
    matched_metabolites = match_fixed_features(feature_columns, FIXED_METABOLITES_LIST, 'met_')

    selected_columns = matched_species + matched_metabolites
    x = df[selected_columns].apply(pd.to_numeric, errors='coerce').fillna(0.0)
    y = df['group'].map(LABEL_MAP)

    if y.isna().any():
        raise ValueError('存在无法映射到标签编码的 group。')

    return x, y.astype(int), selected_columns


def build_cv(y: pd.Series):
    min_class_count = int(y.value_counts().min())
    if min_class_count < 2:
        raise ValueError('至少每组需要 2 个样本才能进行交叉验证。')

    n_splits = min(5, min_class_count)
    if n_splits < 3:
        n_splits = 2

    from sklearn.model_selection import StratifiedKFold

    return StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=42)


def fit_fixed_biomarker_model(x: pd.DataFrame, y: pd.Series) -> dict:
    forest = RandomForestClassifier(
        n_estimators=1000,
        random_state=42,
        class_weight='balanced_subsample',
        min_samples_leaf=2,
        n_jobs=-1,
    )

    cv = build_cv(y)
    oof_pred = cross_val_predict(forest, x, y, cv=cv, method='predict', n_jobs=1)
    oof_proba = cross_val_predict(forest, x, y, cv=cv, method='predict_proba', n_jobs=1)
    forest.fit(x, y)

    importance_df = pd.DataFrame(
        {
            'feature': x.columns,
            'importance': forest.feature_importances_,
        }
    ).sort_values('importance', ascending=False)
    importance_df['feature_type'] = importance_df['feature'].map(feature_type)
    importance_df['display_name'] = importance_df['feature'].map(clean_feature_name)

    return {
        'forest': forest,
        'cv': cv,
        'oof_pred': oof_pred,
        'oof_proba': oof_proba,
        'importance_df': importance_df,
        'selected_columns': list(x.columns),
    }


def save_selected_biomarker_list(output_dir: Path, selected_columns: list[str]) -> None:
    biomarker_df = pd.DataFrame(
        {
            'feature': selected_columns,
            'feature_type': [feature_type(column) for column in selected_columns],
            'display_name': [clean_feature_name(column) for column in selected_columns],
        }
    )
    biomarker_df.to_csv(output_dir / 'selected_fixed_biomarkers.csv', index=False)


def main() -> None:
    project_root = Path(__file__).resolve().parents[1]
    output_dir = project_root / 'results' / 'random_forest_fixed_biomarkers'
    output_dir.mkdir(parents=True, exist_ok=True)

    df = load_analysis_data()
    x, y, selected_columns = prepare_fixed_biomarker_features(df)
    model_results = fit_fixed_biomarker_model(x, y)

    plot_feature_importance(model_results['importance_df'], output_dir)
    roc_df = plot_roc_curves(y, model_results['oof_proba'], output_dir)
    confusion_df = plot_confusion(y, model_results['oof_pred'], output_dir)
    plot_biomarker_panel(model_results['importance_df'].head(10), output_dir)
    save_reports(output_dir, df, y, model_results, roc_df, confusion_df)
    save_selected_biomarker_list(output_dir, selected_columns)

    print(f'固定 biomarker 随机森林结果已输出到: {output_dir}')
    print(f'菌数: {len(FIXED_SPECIES_LIST)}, 代谢物数: {len(FIXED_METABOLITES_LIST)}')
    print(f'总特征数: {len(selected_columns)}')
    print(f'分组顺序: {CLASS_ORDER}')


if __name__ == '__main__':
    main()