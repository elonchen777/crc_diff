from itertools import combinations
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score, classification_report, roc_auc_score
from sklearn.model_selection import train_test_split

from random_forest_diff import CLASS_ORDER, clean_feature_name, load_analysis_data


def feature_type(feature_name: str) -> str:
    if feature_name.startswith('tax_'):
        return 'Microbe'
    if feature_name.startswith('met_'):
        return 'Metabolite'
    return 'Other'


def prepare_feature_matrix(df: pd.DataFrame) -> pd.DataFrame:
    feature_cols = [
        column for column in df.columns
        if column.startswith('tax_') or column.startswith('met_')
    ]
    if not feature_cols:
        raise ValueError('未找到 tax_ 或 met_ 特征列。')

    x = df[feature_cols].apply(pd.to_numeric, errors='coerce').fillna(0.0)
    x = x.loc[:, x.var(axis=0) > 0].copy()
    return x


def run_pairwise_random_forest(
    df: pd.DataFrame,
    pair: tuple[str, str],
    output_dir: Path,
) -> pd.DataFrame:
    pair_name = f'{pair[0]}_vs_{pair[1]}'
    pair_dir = output_dir / pair_name
    pair_dir.mkdir(parents=True, exist_ok=True)

    subset = df[df['group'].isin(pair)].copy()
    if subset.empty:
        raise ValueError(f'{pair_name} 没有样本。')

    x = prepare_feature_matrix(subset)
    y = subset['group'].copy()

    # 保证每一类在训练/验证中都能分到样本。
    min_count = int(y.value_counts().min())
    if min_count < 2:
        raise ValueError(f'{pair_name} 至少每组需要 2 个样本。')

    x_train, x_test, y_train, y_test = train_test_split(
        x,
        y,
        test_size=0.3,
        random_state=42,
        stratify=y,
    )

    forest = RandomForestClassifier(
        n_estimators=1000,
        random_state=42,
        class_weight='balanced_subsample',
        min_samples_leaf=2,
        n_jobs=-1,
    )
    forest.fit(x_train, y_train)

    y_pred = forest.predict(x_test)
    y_proba = forest.predict_proba(x_test)
    accuracy = accuracy_score(y_test, y_pred)
    auc = roc_auc_score(y_test, y_proba[:, 1])

    labels = list(forest.classes_)
    report = classification_report(y_test, y_pred, target_names=labels, digits=4, zero_division=0)
    summary_lines = [
        f'Pair: {pair_name}',
        f'Samples: {len(subset)}',
        f'Train samples: {len(x_train)}',
        f'Test samples: {len(x_test)}',
        f'Class counts: {y.value_counts().to_dict()}',
        f'Accuracy: {accuracy:.4f}',
        f'ROC AUC: {auc:.4f}',
        '',
        'Classification report:',
        report,
    ]
    (pair_dir / 'pairwise_summary.txt').write_text('\n'.join(summary_lines), encoding='utf-8')

    importance_df = pd.DataFrame(
        {
            'feature': x.columns,
            'importance': forest.feature_importances_,
        }
    ).sort_values('importance', ascending=False)
    importance_df['feature_type'] = importance_df['feature'].map(feature_type)
    importance_df['display_name'] = importance_df['feature'].map(clean_feature_name)
    importance_df.to_csv(pair_dir / 'feature_importance_all.csv', index=False)
    importance_df.head(30).to_csv(pair_dir / 'feature_importance_top30.csv', index=False)

    return importance_df[['feature', 'display_name', 'importance']].assign(pair=pair_name)


def save_importance_heatmap(pairwise_importance_df: pd.DataFrame, output_dir: Path, top_n: int = 30) -> None:
    matrix = pairwise_importance_df.pivot_table(
        index='display_name',
        columns='pair',
        values='importance',
        aggfunc='mean',
        fill_value=0.0,
    )

    # 取跨三组比较平均 importance 最高的特征，保证图可读。
    matrix['mean_importance'] = matrix.mean(axis=1)
    matrix = matrix.sort_values('mean_importance', ascending=False).head(top_n)
    matrix = matrix.drop(columns=['mean_importance'])

    matrix.to_csv(output_dir / 'pairwise_importance_heatmap_matrix_top30.csv')

    plt.figure(figsize=(10, max(7, 0.32 * len(matrix))))
    sns.heatmap(matrix, cmap='YlOrRd', linewidths=0.35, linecolor='white')
    plt.title('Pairwise Random Forest Feature Importance Heatmap (Top 30)')
    plt.xlabel('Pairwise Comparison')
    plt.ylabel('Feature')
    plt.tight_layout()
    plt.savefig(output_dir / 'pairwise_importance_heatmap_top30.png', dpi=300, bbox_inches='tight')
    plt.close()


def main() -> None:
    sns.set_theme(style='whitegrid', context='talk')

    project_root = Path(__file__).resolve().parents[1]
    output_dir = project_root / 'results' / 'random_forest_pairwise'
    output_dir.mkdir(parents=True, exist_ok=True)

    df = load_analysis_data()
    all_pairwise_importance = []

    for pair in combinations(CLASS_ORDER, 2):
        pair_df = run_pairwise_random_forest(df, pair, output_dir)
        all_pairwise_importance.append(pair_df)

    pairwise_importance_df = pd.concat(all_pairwise_importance, ignore_index=True)
    pairwise_importance_df.to_csv(output_dir / 'pairwise_importance_long_format.csv', index=False)
    save_importance_heatmap(pairwise_importance_df, output_dir, top_n=30)

    print(f'两两随机森林结果已输出到: {output_dir}')
    print('完成比较: CTRL_vs_CRC-Well, CTRL_vs_CRC-Poor, CRC-Well_vs_CRC-Poor')
    print('训练/验证比例: 0.7 / 0.3')


if __name__ == '__main__':
    main()
