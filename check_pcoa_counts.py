import pandas as pd
import numpy as np
from dataset import BioSmokeDataset

def check_sample_counts():
    """检查met和tax特征矩阵的样本数量差异"""
    print("=" * 80)
    print("检查PCoA图中met和tax点数量不一致的原因")
    print("=" * 80)
    
    # 加载数据
    ds = BioSmokeDataset()
    merged = ds.merge_to_dataframe()
    
    print(f"总样本数: {merged.shape[0]}")
    print(f"总特征数: {merged.shape[1]}")
    
    # 提取tax和met特征
    tax_cols = [c for c in merged.columns if str(c).startswith('tax_')]
    met_cols = [c for c in merged.columns if str(c).startswith('met_')]
    
    print(f"\ntax特征数量: {len(tax_cols)}")
    print(f"met特征数量: {len(met_cols)}")
    
    # 检查每个样本的tax和met特征是否有值
    tax_sum = merged[tax_cols].sum(axis=1)
    met_sum = merged[met_cols].sum(axis=1)
    
    print(f"\ntax特征总和为0的样本数: {(tax_sum == 0).sum()}")
    print(f"met特征总和为0的样本数: {(met_sum == 0).sum()}")
    
    # 检查哪些样本在tax或met中为0
    tax_zero_samples = merged.index[tax_sum == 0].tolist()
    met_zero_samples = merged.index[met_sum == 0].tolist()
    
    print(f"\ntax特征为0的样本: {len(tax_zero_samples)}个")
    if tax_zero_samples:
        print(f"  示例: {tax_zero_samples[:5]}")
    
    print(f"\nmet特征为0的样本: {len(met_zero_samples)}个")
    if met_zero_samples:
        print(f"  示例: {met_zero_samples[:5]}")
    
    # 检查两个特征集都为0的样本
    both_zero = merged.index[(tax_sum == 0) & (met_sum == 0)].tolist()
    print(f"\ntax和met特征都为0的样本: {len(both_zero)}个")
    
    # 检查分组信息
    from pcoa_plot import prepare_early_late_groups
    grouped = prepare_early_late_groups(merged)
    
    print(f"\n分组后的样本数: {grouped.shape[0]}")
    print(f"分组信息:")
    print(grouped['group'].value_counts())
    
    # 检查分组后tax和met的特征
    grouped_tax_sum = grouped[tax_cols].sum(axis=1)
    grouped_met_sum = grouped[met_cols].sum(axis=1)
    
    print(f"\n分组后tax特征总和为0的样本数: {(grouped_tax_sum == 0).sum()}")
    print(f"分组后met特征总和为0的样本数: {(grouped_met_sum == 0).sum()}")
    
    # 检查每个分组的样本数量
    for group_name in ['CRC_early_smoking', 'CRC_early_nonsmoking', 'CRC_late_smoking', 'CRC_late_nonsmoking']:
        group_mask = grouped['group'] == group_name
        group_samples = grouped[group_mask]
        if len(group_samples) > 0:
            group_tax_zero = (group_samples[tax_cols].sum(axis=1) == 0).sum()
            group_met_zero = (group_samples[met_cols].sum(axis=1) == 0).sum()
            print(f"\n{group_name}:")
            print(f"  样本数: {len(group_samples)}")
            print(f"  tax为0的样本: {group_tax_zero}")
            print(f"  met为0的样本: {group_met_zero}")
    
    return merged, grouped, tax_zero_samples, met_zero_samples

if __name__ == '__main__':
    check_sample_counts()