"""
Convert `merged_dataset.csv` into two TSVs usable by MaAsLin2/MaAsLin3 (features and metadata).

参考 adonis_analysis.R 的分组逻辑，生成适用于 Maaslin 分析的输入文件。

Produces:
- dataset/maaslin_features.tsv  (samples x features)
- dataset/maaslin_ko_features.tsv  (samples x KO features)
- dataset/maaslin_metabolite_features.tsv  (samples x metabolite features)
- dataset/maaslin_metadata.tsv  (samples x covariates)

Usage:
    python merged2maaslin.py --input dataset/merged_dataset.csv --outdir dataset
"""

from pathlib import Path
import argparse
import sys

import pandas as pd

def prepare_crc_diff_groups(df: pd.DataFrame) -> pd.DataFrame:
    def _label(row):
        crc = int(row.get('crc_label', 0))
        diff = int(row.get('differentiation', 0))
        
        # Only consider CRC patients (crc_label == 1)
        if crc == 1:
            if diff == 1:
                return f'CRC_poordiff'
            else:
                return f'CRC_welldiff'
        else:
            return f'CTRL'
    
    df = df.copy()
    df['group'] = df.apply(_label, axis=1)
    df = df[df['group'].notna()].copy()
    return df

def _clean_names(cols):
    """Sanitize column names: remove tabs/newlines and make unique."""
    out = []
    seen = {}
    for c in cols:
        s = str(c)
        s = s.replace('\t', ' ').replace('\n', ' ').replace('\r', ' ').strip()
        # replace characters that R may treat specially (e.g. '#' is comment.char)
        s = s.replace('#', '_')
        if s == '':
            s = 'NA'
        if s in seen:
            seen[s] += 1
            s = f"{s}__{seen[s]}"
        else:
            seen[s] = 0
        out.append(s)
    return out

def main():
    parser = argparse.ArgumentParser(
        description='Convert merged_dataset.csv to MaAsLin input format'
    )
    parser.add_argument(
        '--input', 
        type=str, 
        default='dataset/merged_dataset_raw.csv',
        help='Path to merged_dataset.csv'
    )
    parser.add_argument(
        '--outdir', 
        type=str, 
        default='dataset/maaslin',
        help='Output directory for TSV files'
    )
    args = parser.parse_args()
    
    input_path = Path(args.input)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    
    if not input_path.exists():
        print(f"Input file not found: {input_path}", file=sys.stderr)
        sys.exit(1)
    
    # 读取合并数据
    print(f"读取数据: {input_path}")
    df = pd.read_csv(input_path, low_memory=False)
    print(f"原始数据形状: {df.shape}")
    
    # 检查必要的列
    required_cols = ['SAMPLE_ID', 'crc_label', 'tnm_stage', 'smoking_label', 'age', 'gender_label','differentiation']
    missing_cols = [c for c in required_cols if c not in df.columns]
    if missing_cols:
        print(f"警告: 缺少列 {missing_cols}")
        print(f"可用列: {df.columns.tolist()[:20]}...")
    
    # 只保留 CRC 患者数据 (crc_label = 1)
    # print("\n筛选 CRC 患者数据 (crc_label = 1)...")
    # df = df[df['crc_label'] == 1].copy()
    # print(f"筛选后数据形状: {df.shape}")
    
    # 提取分类学特征列 (tax_ 开头的列)
    taxonomy_cols = [c for c in df.columns if c.startswith('tax_')]
    print(f"\n找到 {len(taxonomy_cols)} 个分类学特征")

    # 提取KO特征列 (kegg_ 开头的列)
    ko_cols = [c for c in df.columns if c.startswith('kegg_')]
    print(f"找到 {len(ko_cols)} 个KO特征")
    
    # 提取代谢物特征列 (met_ 开头的列)
    met_cols = [c for c in df.columns if c.startswith('met_')]
    print(f"找到 {len(met_cols)} 个代谢物特征")
    
    if len(taxonomy_cols) == 0:
        print("错误: 未找到 tax_ 开头的分类学特征列", file=sys.stderr)
        sys.exit(1)
    
    # 准备特征矩阵 (samples x features)
    # 提取 SAMPLE_ID 和分类学特征
    features_df = df[['SAMPLE_ID'] + taxonomy_cols].copy()
    
    # 清理特征列名
    features_df.columns = ['SAMPLE_ID'] + _clean_names(taxonomy_cols)
    
    # 设置 SAMPLE_ID 为索引
    features_df = features_df.set_index('SAMPLE_ID')
    
    # 移除全零特征（物种）
    nonzero_features = features_df.columns[features_df.sum() > 0]
    features_df = features_df[nonzero_features]
    print(f"移除全零特征后: {features_df.shape[1]} 个特征")
    
    # 准备元数据
    # 选择用于 Maaslin 分析的协变量
    metadata_cols = ['SAMPLE_ID', 'group', 'age', 'gender_label', 'crc_label', 
                     'tnm_stage', 'smoking_label','differentiation','diff_stage']
    
    # 准备分组列
    df = prepare_crc_diff_groups(df)

    df['diff_stage'] = df['group'].map({'CRC_poordiff': '2', 'CRC_welldiff': '1', 'CTRL': '0'})
    
    # 只保留存在的列
    available_meta_cols = [c for c in metadata_cols if c in df.columns]
    metadata_df = df[available_meta_cols].copy()
    
    # 清理列名
    metadata_df.columns = _clean_names(metadata_df.columns)
    
    # 设置 SAMPLE_ID 为索引
    metadata_df = metadata_df.set_index('SAMPLE_ID')
    
    # 确保索引为字符串类型
    features_df.index = features_df.index.astype(str)
    metadata_df.index = metadata_df.index.astype(str)
    
    # 确保样本匹配
    common_samples = features_df.index.intersection(metadata_df.index)
    print(f"\n共同样本数: {len(common_samples)}")
    
    features_df = features_df.loc[common_samples]
    metadata_df = metadata_df.loc[common_samples]

    # 准备KO特征矩阵 (samples x KO features)
    ko_features_df = None
    if len(ko_cols) > 0:
        ko_features_df = df[['SAMPLE_ID'] + ko_cols].copy()
        ko_features_df.columns = ['SAMPLE_ID'] + _clean_names(ko_cols)
        ko_features_df = ko_features_df.set_index('SAMPLE_ID')

        # 转换为数值并移除全零特征
        ko_features_df = ko_features_df.apply(pd.to_numeric, errors='coerce').fillna(0)
        nonzero_ko_features = ko_features_df.columns[ko_features_df.sum() > 0]
        ko_features_df = ko_features_df[nonzero_ko_features]
        print(f"移除全零KO特征后: {ko_features_df.shape[1]} 个特征")

        ko_features_df.index = ko_features_df.index.astype(str)
        ko_features_df = ko_features_df.loc[ko_features_df.index.intersection(metadata_df.index)]
        ko_features_df = ko_features_df.loc[metadata_df.index.intersection(ko_features_df.index)]
        ko_features_df.index.name = 'SAMPLE_ID'
    
    # 准备代谢物特征矩阵 (samples x met features)
    met_features_df = None
    if len(met_cols) > 0:
        met_features_df = df[['SAMPLE_ID'] + met_cols].copy()
        met_features_df.columns = ['SAMPLE_ID'] + _clean_names(met_cols)
        met_features_df = met_features_df.set_index('SAMPLE_ID')

        # 转换为数值并移除全零特征
        met_features_df = met_features_df.apply(pd.to_numeric, errors='coerce').fillna(0)
        nonzero_met_features = met_features_df.columns[met_features_df.sum() > 0]
        met_features_df = met_features_df[nonzero_met_features]
        print(f"移除全零代谢物特征后: {met_features_df.shape[1]} 个特征")

        met_features_df.index = met_features_df.index.astype(str)
        met_features_df = met_features_df.loc[met_features_df.index.intersection(metadata_df.index)]
        met_features_df = met_features_df.loc[metadata_df.index.intersection(met_features_df.index)]
        met_features_df.index.name = 'SAMPLE_ID'
    
    # 写入输出文件
    feat_out = outdir / "maaslin_features.tsv"
    ko_feat_out = outdir / "maaslin_ko_features.tsv"
    met_feat_out = outdir / "maaslin_metabolite_features.tsv"
    meta_out = outdir / "maaslin_metadata.tsv"
    
    # 确保索引标签兼容 R 的 read.table(row.names=1)
    features_df.index.name = 'SAMPLE_ID'
    metadata_df.index.name = 'SAMPLE_ID'
    
    features_df.to_csv(feat_out, sep="\t", index=True, index_label='SAMPLE_ID', na_rep='')
    metadata_df.to_csv(meta_out, sep="\t", index=True, index_label='SAMPLE_ID', na_rep='')
    if ko_features_df is not None:
        ko_features_df.to_csv(ko_feat_out, sep="\t", index=True, index_label='SAMPLE_ID', na_rep='')
    if met_features_df is not None:
        met_features_df.to_csv(met_feat_out, sep="\t", index=True, index_label='SAMPLE_ID', na_rep='')
    
    print(f"\n输出文件:")
    print(f"  Features: {feat_out} (samples x features: {features_df.shape})")
    if ko_features_df is not None:
        print(f"  KO Features: {ko_feat_out} (samples x KO features: {ko_features_df.shape})")
    else:
        print(f"  KO Features: 未导出（未找到 kegg_ 开头特征列）")
    if met_features_df is not None:
        print(f"  Metabolite Features: {met_feat_out} (samples x met features: {met_features_df.shape})")
    else:
        print(f"  Metabolite Features: 未导出（未找到 met_ 开头特征列）")
    print(f"  Metadata: {meta_out} (samples x covariates: {metadata_df.shape})")
    
    # 打印分组信息供参考
    print("\n分组信息（用于 Maaslin 分析）:")
    if 'group' in metadata_df.columns:
        print(metadata_df['group'].value_counts())


if __name__ == "__main__":
    main()
