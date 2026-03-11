import pandas as pd
import os

def simple_filter_diff_stage_species(all_results_file, species_abund_file, output_prefix="diff_stage_filtered"):
    """
    简化版本：筛选diff_stage相关菌种
    """
    
    # 1. 读取数据
    all_results = pd.read_csv(all_results_file, sep='\t')
    species_abund = pd.read_csv(species_abund_file, sep='\t')
    print(f"读取all_results.tsv，共 {len(all_results)} 行")
    print(f"读取taxonomy_Species_abund.txt，共 {len(species_abund)} 行")
    
    # 2. 找出diff_stage相关的特征
    diff_stage_features = all_results[all_results['metadata'] == 'diff_stage']
    print(f"diff_stage相关特征共 {len(diff_stage_features)} 个")

    # 2.1 筛选pval < 0.05
    # diff_stage_features = diff_stage_features[diff_stage_features['qval'] < 0.05]
    # print(f"pval < 0.05 共 {len(diff_stage_features)} 个")

    # 2.2 根据coef绝对值大小排序
    diff_stage_features = pd.DataFrame(diff_stage_features)
    diff_stage_features = diff_stage_features.sort_values(by='coef', key=abs, ascending=False)
    
    # 3. 获取物种名列名称
    species_col = species_abund.columns[0]

    #去除tax_前缀
    diff_stage_features['feature'] = diff_stage_features['feature'].str.replace('tax_', '')
    
    # 4. 筛选物种丰度数据
    # 使用部分匹配（因为物种名可能不完全一致）
    mask = species_abund[species_col].apply(
        lambda x: any(feature in str(x) for feature in diff_stage_features['feature'])
    )
    
    filtered_data = species_abund[mask]
    
    # 5. 保存结果
    filtered_data.to_csv(f"dataset/{output_prefix}_abundance.txt", sep='\t', index=False)
    diff_stage_features.to_csv(f"dataset/{output_prefix}_species_list.csv", sep='\t', index=False)
    
    print(f"找到 {len(diff_stage_features)} 个diff_stage相关菌种")
    print(f"筛选出 {len(filtered_data)} 行丰度数据")
    print(f"结果已保存到 {output_prefix}_abundance.txt 和 {output_prefix}_species_list.csv")
    
    return filtered_data

# 使用示例
if __name__ == "__main__":
    # 直接指定文件路径
    all_results_path = "results/maaslin2/all_results.tsv"
    species_abund_path = "dataset/taxonomy_Species_abund.txt"
    
    if os.path.exists(all_results_path) and os.path.exists(species_abund_path):
        simple_filter_diff_stage_species(all_results_path, species_abund_path)
    else:
        print("请指定正确的文件路径")