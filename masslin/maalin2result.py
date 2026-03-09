import pandas as pd
import os

def simple_filter_smoke_species(all_results_file, species_abund_file, output_prefix="smoke_filtered"):
    """
    简化版本：筛选smoke相关菌种
    """
    
    # 1. 读取数据
    all_results = pd.read_csv(all_results_file, sep='\t')
    species_abund = pd.read_csv(species_abund_file, sep='\t')
    print(f"读取all_results.tsv，共 {len(all_results)} 行")
    print(f"读取taxonomy_Species_abund.txt，共 {len(species_abund)} 行")
    
    # 2. 找出smoke相关的特征
    smoke_features = all_results[all_results['metadata'] == 'smoke']
    print(f"smoke相关特征共 {len(smoke_features)} 个")

    # 2.1 筛选pval < 0.05
    smoke_features = smoke_features[smoke_features['qval'] < 0.05]
    print(f"pval < 0.05 共 {len(smoke_features)} 个")

    # 2.2 根据coef绝对值大小排序
    smoke_features = pd.DataFrame(smoke_features)
    smoke_features = smoke_features.sort_values(by='coef', key=abs, ascending=False)
    
    # 3. 获取物种名列名称
    species_col = species_abund.columns[0]
    
    # 4. 筛选物种丰度数据
    # 使用部分匹配（因为物种名可能不完全一致）
    mask = species_abund[species_col].apply(
        lambda x: any(feature in str(x) for feature in smoke_features['feature'])
    )
    
    filtered_data = species_abund[mask]
    
    # 5. 保存结果
    filtered_data.to_csv(f"dataset/masslin/{output_prefix}_abundance.txt", sep='\t', index=False)
    smoke_features.to_csv(f"dataset/masslin/{output_prefix}_species_list.csv", sep='\t', index=False)
    
    print(f"找到 {len(smoke_features)} 个smoke相关菌种")
    print(f"筛选出 {len(filtered_data)} 行丰度数据")
    print(f"结果已保存到 {output_prefix}_abundance.txt 和 {output_prefix}_species_list.txt")
    
    return filtered_data

# 使用示例
if __name__ == "__main__":
    # 直接指定文件路径
    all_results_path = "results/all_results.tsv"
    species_abund_path = "dataset/taxonomy_Species_abund.txt"
    
    if os.path.exists(all_results_path) and os.path.exists(species_abund_path):
        simple_filter_smoke_species(all_results_path, species_abund_path)
    else:
        print("请指定正确的文件路径")