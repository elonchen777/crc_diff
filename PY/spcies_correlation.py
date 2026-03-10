#!/usr/bin/env python3
"""
分组并行计算宏基因组相关性脚本
分别计算CRC早期/晚期吸烟/非吸烟四组的相关性
使用dataset.py导入数据
"""

import pandas as pd
import numpy as np
from scipy.stats import spearmanr, pearsonr
from typing import List, Tuple, Optional, Dict
import multiprocessing as mp
from concurrent.futures import ProcessPoolExecutor, as_completed
import warnings
import time
from pathlib import Path
import os

# 导入项目中的模块
from dataset import BioSmokeDataset
from split_group import prepare_early_late_groups

warnings.filterwarnings('ignore')

class GroupParallelCorrelation:
    """
    分组并行计算宏基因组相关性的类
    """
    
    def __init__(self, 
                 n_jobs: int = -1,
                 correlation_method: str = 'spearman',
                 output_dir: str = 'results/group_correlations'):
        """
        初始化分组并行相关性计算器
        
        Args:
            n_jobs: 并行进程数，-1表示使用所有可用CPU
            correlation_method: 相关性计算方法 ('spearman' 或 'pearson')
            output_dir: 结果输出目录
        """
        self.n_jobs = mp.cpu_count() if n_jobs == -1 else min(n_jobs, mp.cpu_count())
        self.correlation_method = correlation_method
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        print(f"初始化分组并行计算器: {self.n_jobs} 个进程, 方法: {correlation_method}")
    
    def load_group_data(self) -> Dict[str, pd.DataFrame]:
        """
        使用dataset.py加载数据并按组分割
        
        Returns:
            Dict: 包含四个分组数据的字典
        """
        print("加载宏基因组数据...")
        start_time = time.time()
        
        # 使用BioSmokeDataset加载数据
        ds = BioSmokeDataset()
        df = ds.merge_to_dataframe()
        
        # 准备分组
        grouped_data = prepare_early_late_groups(df)
        
        # 按组分割数据
        groups = ['CRC_early_nonsmoking', 'CRC_early_smoking', 
                 'CRC_late_nonsmoking', 'CRC_late_smoking']
        
        group_data = {}
        for group in groups:
            group_samples = grouped_data[grouped_data['group'] == group]
            if len(group_samples) > 0:
                # 只保留物种数据列（以'tax_'开头的列）
                species_cols = [col for col in group_samples.columns if col.startswith('tax_')]
                group_data[group] = group_samples[species_cols]
                print(f"{group}: {len(group_samples)} 个样本, {len(species_cols)} 个物种")
            else:
                print(f"警告: {group} 组没有样本数据")
                group_data[group] = pd.DataFrame()
        
        elapsed_time = time.time() - start_time
        print(f"数据加载和分组完成，耗时: {elapsed_time:.2f} 秒")
        
        return group_data
    
    def _calculate_correlation_pair(self, args: Tuple[str, str, pd.Series, pd.Series]) -> Tuple[str, str, float, float]:
        """
        计算两个特征之间的相关性（工作进程函数）
        
        Args:
            args: 包含特征1名称、特征2名称、特征1数据、特征2数据的元组
            
        Returns:
            Tuple: (特征1, 特征2, 相关性系数, p值)
        """
        feature1, feature2, data1, data2 = args
        
        # 移除NaN值
        valid_mask = ~(np.isnan(data1) | np.isnan(data2))
        data1_valid = data1[valid_mask]
        data2_valid = data2[valid_mask]
        
        if len(data1_valid) < 3:  # 至少需要3个有效数据点
            return feature1, feature2, np.nan, np.nan
        
        if self.correlation_method == 'spearman':
            corr, p_value = spearmanr(data1_valid, data2_valid)
        else:
            corr, p_value = pearsonr(data1_valid, data2_valid)
        
        return feature1, feature2, corr, p_value
    
    def calculate_group_correlations(self, 
                                   group_data: Dict[str, pd.DataFrame],
                                   top_k: int = 200,
                                   correlation_threshold: float = 0.3) -> Dict[str, pd.DataFrame]:
        """
        为每个组计算相关性网络
        
        Args:
            group_data: 包含各组数据的字典
            top_k: 每个组计算前k个最重要的物种
            correlation_threshold: 相关性阈值，只保留大于此值的相关性
            
        Returns:
            Dict: 包含各组相关性网络的字典
        """
        print(f"\n=== 开始计算各组相关性网络 ===")
        print(f"每个组计算前 {top_k} 个重要物种")
        print(f"相关性阈值: {correlation_threshold}")
        
        group_correlations = {}
        
        for group_name, data in group_data.items():
            if data.empty:
                print(f"跳过 {group_name} (无数据)")
                group_correlations[group_name] = pd.DataFrame()
                continue
            
            print(f"\n--- 计算 {group_name} 组相关性 ---")
            print(f"样本数: {data.shape[0]}, 物种数: {data.shape[1]}")
            
            # 选择最重要的物种（基于方差）
            if data.shape[1] > top_k:
                variances = data.var()
                important_species = variances.nlargest(top_k).index.tolist()
                data = data[important_species]
                print(f"选择前 {top_k} 个高方差物种")
            
            species_list = data.columns.tolist()
            n_species = len(species_list)
            
            # 准备所有需要计算的相关性对
            correlation_pairs = []
            for i in range(n_species):
                for j in range(i + 1, n_species):  # 只计算上三角，避免重复
                    correlation_pairs.append((
                        species_list[i], 
                        species_list[j], 
                        data[species_list[i]], 
                        data[species_list[j]]
                    ))
            
            total_pairs = len(correlation_pairs)
            print(f"需要计算的相关性对数量: {total_pairs}")
            
            # 并行计算相关性
            results = []
            with ProcessPoolExecutor(max_workers=self.n_jobs) as executor:
                # 分批提交任务以避免内存问题
                batch_size = min(5000, total_pairs // self.n_jobs + 1)
                
                for batch_start in range(0, total_pairs, batch_size):
                    batch_end = min(batch_start + batch_size, total_pairs)
                    batch_pairs = correlation_pairs[batch_start:batch_end]
                    
                    futures = [
                        executor.submit(self._calculate_correlation_pair, pair) 
                        for pair in batch_pairs
                    ]
                    
                    for future in as_completed(futures):
                        feature1, feature2, corr, p_value = future.result()
                        # 只采信p值<0.05的相关性结果
                        if not np.isnan(corr) and abs(corr) >= correlation_threshold and p_value < 0.05:
                            results.append({
                                'Species1': feature1,
                                'Species2': feature2,
                                'Correlation': corr,
                                'PValue': p_value,
                                'AbsCorrelation': abs(corr),
                                'Significant': True  # 所有结果都是显著的
                            })
                    
                    print(f"{group_name}: 完成批次 {batch_end}/{total_pairs}")
            
            # 转换为DataFrame
            if results:
                correlation_network = pd.DataFrame(results)
                correlation_network = correlation_network.sort_values('AbsCorrelation', ascending=False)
                group_correlations[group_name] = correlation_network
                print(f"{group_name}: 发现 {len(correlation_network)} 个强相关性对")
            else:
                group_correlations[group_name] = pd.DataFrame()
                print(f"{group_name}: 未发现强相关性对")
        
        return group_correlations
    
    def analyze_group_correlations(self, group_correlations: Dict[str, pd.DataFrame]):
        """
        分析各组的相关性结果
        
        Args:
            group_correlations: 包含各组相关性网络的字典
        """
        print(f"\n=== 各组相关性分析结果 ===")
        
        for group_name, correlations in group_correlations.items():
            if correlations.empty:
                print(f"\n{group_name}: 无相关性数据")
                continue
            
            print(f"\n{group_name}:")
            print(f"  显著相关性边数 (p < 0.05): {len(correlations)}")
            print(f"  平均相关性: {correlations['Correlation'].mean():.4f}")
            print(f"  最强正相关性: {correlations['Correlation'].max():.4f}")
            print(f"  最强负相关性: {correlations['Correlation'].min():.4f}")
            
            # 正负相关性统计
            positive = correlations[correlations['Correlation'] > 0]
            negative = correlations[correlations['Correlation'] < 0]
            print(f"  正相关性比例: {len(positive)/len(correlations):.2%}")
            
            # p值统计（所有结果都是显著的）
            if 'PValue' in correlations.columns:
                highly_significant = correlations[correlations['PValue'] < 0.01]
                very_significant = correlations[correlations['PValue'] < 0.001]
                
                print(f"  高度显著性 (p < 0.01): {len(highly_significant)} ({len(highly_significant)/len(correlations):.2%})")
                print(f"  极高度显著性 (p < 0.001): {len(very_significant)} ({len(very_significant)/len(correlations):.2%})")
                print(f"  所有相关性均显著 (p < 0.05)")
                
                # 显示前5个最强相关性（所有都是显著的）
                top_5 = correlations.head(5)
                print(f"  前5个最强相关性 (p < 0.05):")
                for idx, row in top_5.iterrows():
                    significance = "***" if row['PValue'] < 0.001 else "**" if row['PValue'] < 0.01 else "*"
                    print(f"    {row['Species1']} <-> {row['Species2']}: {row['Correlation']:.4f} (p={row['PValue']:.4f}) {significance}")
            else:
                print(f"  警告: 未找到PValue列")
                
                # 显示前5个最强相关性
                top_5 = correlations.head(5)
                print(f"  前5个最强相关性:")
                for idx, row in top_5.iterrows():
                    print(f"    {row['Species1']} <-> {row['Species2']}: {row['Correlation']:.4f}")
    
    def save_group_results(self, group_correlations: Dict[str, pd.DataFrame]):
        """
        保存各组的相关性结果
        
        Args:
            group_correlations: 包含各组相关性网络的字典
        """
        print(f"\n=== 保存结果 ===")
        
        for group_name, correlations in group_correlations.items():
            if correlations.empty:
                print(f"跳过 {group_name} (无数据)")
                continue
            
            # 保存完整的相关性网络
            filename = f"{group_name}_correlations.csv"
            filepath = self.output_dir / filename
            correlations.to_csv(filepath, index=False)
            print(f"保存: {filepath}")
            
            # 保存前50个最强相关性
            top_50 = correlations.head(50)
            top_filename = f"{group_name}_top50_correlations.csv"
            top_filepath = self.output_dir / top_filename
            top_50.to_csv(top_filepath, index=False)
            print(f"保存: {top_filepath}")
    
    def compare_groups(self, group_correlations: Dict[str, pd.DataFrame], top_n: int = 50):
        """
        比较不同组的相关性模式，找出差异最大的相关性对
        
        Args:
            group_correlations: 包含各组相关性网络的字典
            top_n: 返回差异最大的前n组
        """
        print(f"\n=== 组间相关性差异分析 ===")
        
        # 收集各组的基本统计信息
        group_stats = []
        
        for group_name, correlations in group_correlations.items():
            if correlations.empty:
                continue
            
            stats = {
                'Group': group_name,
                'TotalEdges': len(correlations),
                'MeanCorrelation': correlations['Correlation'].mean(),
                'MaxCorrelation': correlations['Correlation'].max(),
                'MinCorrelation': correlations['Correlation'].min(),
                'PositiveRatio': len(correlations[correlations['Correlation'] > 0]) / len(correlations)
            }
            group_stats.append(stats)
        
        if group_stats:
            stats_df = pd.DataFrame(group_stats)
            print("各组统计信息:")
            print(stats_df.to_string(index=False))
            
            # 保存比较结果
            compare_file = self.output_dir / "group_comparison.csv"
            stats_df.to_csv(compare_file, index=False)
            print(f"组间比较结果保存至: {compare_file}")
        
        # 分析相关性差异
        self._analyze_correlation_differences(group_correlations, top_n)
        
        # 专门分析吸烟/非吸烟差异
        # self._analyze_smoking_vs_nonsmoking_differences(group_correlations, top_n)
        
        # 分析四组都高度显著的相关性 (p < 0.0005)
        self._analyze_highly_significant_correlations(group_correlations, p_threshold=0.0005)
    
    def _analyze_correlation_differences(self, group_correlations: Dict[str, pd.DataFrame], top_n: int = 50):
        """
        分析不同组之间相关性差异
        
        Args:
            group_correlations: 包含各组相关性网络的字典
            top_n: 返回差异最大的前n组
        """
        print(f"\n--- 分析相关性差异 (前{top_n}个最大差异) ---")
        
        # 获取所有组的名称
        group_names = [name for name, data in group_correlations.items() if not data.empty]
        
        if len(group_names) < 2:
            print("需要至少2个有效组才能进行差异分析")
            return
        
        # 创建所有组的合并相关性表（只使用p值<0.05的数据）
        all_correlations = {}
        all_pvalues = {}
        
        for group_name in group_names:
            correlations = group_correlations[group_name]
            # 为每个相关性对创建唯一标识符
            for _, row in correlations.iterrows():
                # 只使用p值<0.05的数据
                if row['PValue'] < 0.05:
                    species_pair = tuple(sorted([row['Species1'], row['Species2']]))
                    if species_pair not in all_correlations:
                        all_correlations[species_pair] = {}
                        all_pvalues[species_pair] = {}
                    all_correlations[species_pair][group_name] = row['Correlation']
                    all_pvalues[species_pair][group_name] = row['PValue']
        
        # 计算差异
        difference_results = []
        
        for species_pair, group_corrs in all_correlations.items():
            # 确保在所有组中都有数据，并且所有组的p值都<0.05
            if len(group_corrs) == len(group_names):
                # 检查所有组的p值是否都<0.05
                all_significant = all(p < 0.05 for p in all_pvalues[species_pair].values())
                
                if all_significant:
                    # 计算最大差异
                    correlations = list(group_corrs.values())
                    max_diff = max(correlations) - min(correlations)
                    
                    # 计算标准差
                    std_dev = np.std(correlations)
                    
                    # 计算组间差异（所有可能的两两组合差异）
                    pairwise_diffs = []
                    for i in range(len(group_names)):
                        for j in range(i+1, len(group_names)):
                            diff = abs(group_corrs[group_names[i]] - group_corrs[group_names[j]])
                            pairwise_diffs.append(diff)
                    
                    avg_pairwise_diff = np.mean(pairwise_diffs)
                    
                    # 收集结果（包括p值信息）
                    result = {
                        'Species1': species_pair[0],
                        'Species2': species_pair[1],
                        'MaxDifference': max_diff,
                        'StdDev': std_dev,
                        'AvgPairwiseDifference': avg_pairwise_diff,
                        'AllGroupsSignificant': True  # 标记所有组都显著
                    }
                    
                    # 添加每个组的相关性值和p值
                    for group_name in group_names:
                        result[f'{group_name}_Correlation'] = group_corrs.get(group_name, np.nan)
                        result[f'{group_name}_PValue'] = all_pvalues[species_pair].get(group_name, np.nan)
                    
                    # 计算平均p值和显著性比例
                    pvalues = [all_pvalues[species_pair].get(group_name, np.nan) for group_name in group_names]
                    valid_pvalues = [p for p in pvalues if not np.isnan(p)]
                    if valid_pvalues:
                        result['AvgPValue'] = np.mean(valid_pvalues)
                        result['MinPValue'] = np.min(valid_pvalues)
                        result['MaxPValue'] = np.max(valid_pvalues)
                        result['SignificantRatio'] = 1.0  # 所有组都显著
                    else:
                        result['AvgPValue'] = np.nan
                        result['MinPValue'] = np.nan
                        result['MaxPValue'] = np.nan
                        result['SignificantRatio'] = np.nan
                    
                    difference_results.append(result)
        
        if not difference_results:
            print("没有找到在所有组中都存在的相关性对")
            return
        
        # 转换为DataFrame并排序
        diff_df = pd.DataFrame(difference_results)
        
        # 按最大差异排序
        top_diff_max = diff_df.nlargest(top_n, 'MaxDifference')
        
        # 按平均两两差异排序
        top_diff_avg = diff_df.nlargest(top_n, 'AvgPairwiseDifference')
        
        # 按标准差排序（变异性）
        top_diff_std = diff_df.nlargest(top_n, 'StdDev')
        
        print(f"\n1. 按最大差异排序 (前{top_n}个):")
        print("物种对 | 最大差异 | 各组相关性")
        for idx, row in top_diff_max.head(10).iterrows():
            group_corrs_str = ", ".join([f"{group}:{row[f'{group}_Correlation']:.3f}" 
                                       for group in group_names])
            print(f"{row['Species1']} <-> {row['Species2']} | {row['MaxDifference']:.3f} | {group_corrs_str}")
        
        print(f"\n2. 按平均两两差异排序 (前{top_n}个):")
        print("物种对 | 平均差异 | 各组相关性")
        for idx, row in top_diff_avg.head(10).iterrows():
            group_corrs_str = ", ".join([f"{group}:{row[f'{group}_Correlation']:.3f}" 
                                       for group in group_names])
            print(f"{row['Species1']} <-> {row['Species2']} | {row['AvgPairwiseDifference']:.3f} | {group_corrs_str}")
        
        # 保存差异分析结果
        self._save_difference_analysis(diff_df, top_diff_max, top_diff_avg, top_diff_std)
    
    def _analyze_smoking_vs_nonsmoking_differences(self, group_correlations: Dict[str, pd.DataFrame], top_n: int = 50):
        """
        专门分析吸烟与非吸烟组之间的相关性模式差异
        
        Args:
            group_correlations: 包含各组相关性网络的字典
            top_n: 返回差异最大的前n组
        """
        print(f"\n--- 吸烟 vs 非吸烟相关性差异分析 (前{top_n}个最大差异) ---")
        
        # 定义吸烟组和非吸烟组
        smoking_groups = ['CRC_early_smoking', 'CRC_late_smoking']
        nonsmoking_groups = ['CRC_early_nonsmoking', 'CRC_late_nonsmoking']
        
        # 检查各组数据是否存在
        available_smoking = [g for g in smoking_groups if g in group_correlations and not group_correlations[g].empty]
        available_nonsmoking = [g for g in nonsmoking_groups if g in group_correlations and not group_correlations[g].empty]
        
        if not available_smoking or not available_nonsmoking:
            print("吸烟组或非吸烟组数据不足，无法进行比较")
            return
        
        print(f"吸烟组: {available_smoking}")
        print(f"非吸烟组: {available_nonsmoking}")
        
        # 创建吸烟组和非吸烟组的合并相关性表（包含p值）
        smoking_correlations = {}
        nonsmoking_correlations = {}
        smoking_pvalues = {}
        nonsmoking_pvalues = {}
        
        # 合并吸烟组数据（只使用p值<0.05的数据）
        for group in available_smoking:
            correlations = group_correlations[group]
            for _, row in correlations.iterrows():
                # 只使用p值<0.05的数据
                if row['PValue'] < 0.05:
                    species_pair = tuple(sorted([row['Species1'], row['Species2']]))
                    if species_pair not in smoking_correlations:
                        smoking_correlations[species_pair] = []
                        smoking_pvalues[species_pair] = []
                    smoking_correlations[species_pair].append(row['Correlation'])
                    smoking_pvalues[species_pair].append(row['PValue'])
        
        # 合并非吸烟组数据（只使用p值<0.05的数据）
        for group in available_nonsmoking:
            correlations = group_correlations[group]
            for _, row in correlations.iterrows():
                # 只使用p值<0.05的数据
                if row['PValue'] < 0.05:
                    species_pair = tuple(sorted([row['Species1'], row['Species2']]))
                    if species_pair not in nonsmoking_correlations:
                        nonsmoking_correlations[species_pair] = []
                        nonsmoking_pvalues[species_pair] = []
                    nonsmoking_correlations[species_pair].append(row['Correlation'])
                    nonsmoking_pvalues[species_pair].append(row['PValue'])
        
        # 计算吸烟/非吸烟差异
        difference_results = []
        
        # 找出在两个集合中都存在的物种对
        common_pairs = set(smoking_correlations.keys()) & set(nonsmoking_correlations.keys())
        
        print(f"在吸烟组和非吸烟组中都存在的物种对数量: {len(common_pairs)}")
        
        for species_pair in common_pairs:
            # 检查吸烟组和非吸烟组的所有p值是否都<0.05
            smoking_all_significant = all(p < 0.05 for p in smoking_pvalues.get(species_pair, []))
            nonsmoking_all_significant = all(p < 0.05 for p in nonsmoking_pvalues.get(species_pair, []))
            
            # 只有当吸烟组和非吸烟组的所有p值都<0.05时才计入比较
            if smoking_all_significant and nonsmoking_all_significant:
                # 计算吸烟组的平均相关性
                smoking_avg = np.mean(smoking_correlations[species_pair])
                
                # 计算非吸烟组的平均相关性
                nonsmoking_avg = np.mean(nonsmoking_correlations[species_pair])
                
                # 计算差异
                difference = smoking_avg - nonsmoking_avg
                abs_difference = abs(difference)
                
                # 计算效应大小（标准化差异）
                smoking_std = np.std(smoking_correlations[species_pair])
                nonsmoking_std = np.std(nonsmoking_correlations[species_pair])
                pooled_std = np.sqrt((smoking_std**2 + nonsmoking_std**2) / 2)
                effect_size = difference / pooled_std if pooled_std > 0 else 0
                
                # 计算p值相关的统计信息
                smoking_pvalue_avg = np.mean(smoking_pvalues.get(species_pair, [np.nan]))
                nonsmoking_pvalue_avg = np.mean(nonsmoking_pvalues.get(species_pair, [np.nan]))
                
                # 计算显著性比例（应该都是1.0，因为所有p值都<0.05）
                smoking_sig_ratio = 1.0 if smoking_all_significant else 0.0
                nonsmoking_sig_ratio = 1.0 if nonsmoking_all_significant else 0.0
                
                # 使用t检验比较相关性差异的显著性
                try:
                    if len(smoking_correlations[species_pair]) >= 2 and len(nonsmoking_correlations[species_pair]) >= 2:
                        from scipy.stats import ttest_ind
                        t_stat, p_value_diff = ttest_ind(smoking_correlations[species_pair], 
                                                        nonsmoking_correlations[species_pair], 
                                                        equal_var=False)
                    else:
                        t_stat, p_value_diff = np.nan, np.nan
                except:
                    t_stat, p_value_diff = np.nan, np.nan
                
                # 收集结果
                result = {
                    'Species1': species_pair[0],
                    'Species2': species_pair[1],
                    'Smoking_AvgCorrelation': smoking_avg,
                    'Nonsmoking_AvgCorrelation': nonsmoking_avg,
                    'Difference': difference,
                    'AbsDifference': abs_difference,
                    'EffectSize': effect_size,
                    'T_Statistic': t_stat,
                    'PValue_Difference': p_value_diff,
                    'Smoking_PValue_Avg': smoking_pvalue_avg,
                    'Nonsmoking_PValue_Avg': nonsmoking_pvalue_avg,
                    'Smoking_Significant_Ratio': smoking_sig_ratio,
                    'Nonsmoking_Significant_Ratio': nonsmoking_sig_ratio,
                    'Smoking_Groups': len(smoking_correlations[species_pair]),
                    'Nonsmoking_Groups': len(nonsmoking_correlations[species_pair]),
                    'Significant_Difference': p_value_diff < 0.05 if not np.isnan(p_value_diff) else False,
                    'AllGroupsSignificant': True  # 标记所有组都显著
                }
                
                difference_results.append(result)
        
        if not difference_results:
            print("没有找到在吸烟组和非吸烟组中都存在的相关性对")
            return
        
        # 转换为DataFrame并排序
        diff_df = pd.DataFrame(difference_results)
        
        # 按绝对差异排序（最大差异）
        top_diff_abs = diff_df.nlargest(top_n, 'AbsDifference')
        
        # 按效应大小排序（标准化差异）
        top_diff_effect = diff_df.nlargest(top_n, 'EffectSize')
        
        # 按差异显著性排序（p值最小）
        top_diff_significant = diff_df[~diff_df['PValue_Difference'].isna()].nsmallest(top_n, 'PValue_Difference')
        
        # 按吸烟组相关性增强排序（吸烟组相关性更高的）
        smoking_enhanced = diff_df[diff_df['Difference'] > 0].nlargest(top_n, 'Difference')
        
        # 按非吸烟组相关性增强排序（非吸烟组相关性更高的）
        nonsmoking_enhanced = diff_df[diff_df['Difference'] < 0].nsmallest(top_n, 'Difference')
        
        # 统计显著性差异
        significant_differences = diff_df[diff_df['Significant_Difference'] == True]
        print(f"显著性差异的物种对数量 (p < 0.05): {len(significant_differences)}")
        
        print(f"\n1. 吸烟 vs 非吸烟最大差异 (前{top_n}个):")
        print("物种对 | 吸烟组相关性 | 非吸烟组相关性 | 差异 | 效应大小 | p值")
        for idx, row in top_diff_abs.head(10).iterrows():
            direction = "吸烟增强" if row['Difference'] > 0 else "非吸烟增强"
            p_value_info = f"{row['PValue_Difference']:.4f}" if not np.isnan(row['PValue_Difference']) else "N/A"
            significance = "*" if row['Significant_Difference'] else ""
            print(f"{row['Species1']} <-> {row['Species2']} | {row['Smoking_AvgCorrelation']:.3f} | {row['Nonsmoking_AvgCorrelation']:.3f} | {row['Difference']:.3f} ({direction}) | {row['EffectSize']:.3f} | {p_value_info}{significance}")
        
        print(f"\n2. 吸烟 vs 非吸烟效应大小最大 (前{top_n}个):")
        print("物种对 | 吸烟组相关性 | 非吸烟组相关性 | 差异 | 效应大小 | p值")
        for idx, row in top_diff_effect.head(10).iterrows():
            direction = "吸烟增强" if row['Difference'] > 0 else "非吸烟增强"
            p_value_info = f"{row['PValue_Difference']:.4f}" if not np.isnan(row['PValue_Difference']) else "N/A"
            significance = "*" if row['Significant_Difference'] else ""
            print(f"{row['Species1']} <-> {row['Species2']} | {row['Smoking_AvgCorrelation']:.3f} | {row['Nonsmoking_AvgCorrelation']:.3f} | {row['Difference']:.3f} ({direction}) | {row['EffectSize']:.3f} | {p_value_info}{significance}")
        
        print(f"\n3. 吸烟 vs 非吸烟显著性差异最大 (前{top_n}个):")
        print("物种对 | 吸烟组相关性 | 非吸烟组相关性 | 差异 | p值 | 显著性")
        for idx, row in top_diff_significant.head(10).iterrows():
            direction = "吸烟增强" if row['Difference'] > 0 else "非吸烟增强"
            significance = "***" if row['PValue_Difference'] < 0.001 else "**" if row['PValue_Difference'] < 0.01 else "*"
            print(f"{row['Species1']} <-> {row['Species2']} | {row['Smoking_AvgCorrelation']:.3f} | {row['Nonsmoking_AvgCorrelation']:.3f} | {row['Difference']:.3f} ({direction}) | {row['PValue_Difference']:.4f} | {significance}")
        
        # 保存吸烟/非吸烟差异分析结果
        self._save_smoking_difference_analysis(diff_df, top_diff_abs, top_diff_effect, 
                                              smoking_enhanced, nonsmoking_enhanced)
    
    
    def _analyze_highly_significant_correlations(self, group_correlations: Dict[str, pd.DataFrame], 
                                                p_threshold: float = 0.0005):
        """
        分析四组都高度显著的相关性 (p < threshold)
        
        Args:
            group_correlations: 包含各组相关性网络的字典
            p_threshold: p值阈值
        """
        print(f"\n--- 分析四组都高度显著的相关性 (p < {p_threshold}) ---")
        
        # 定义四个分组
        groups = ['CRC_early_smoking', 'CRC_early_nonsmoking', 
                 'CRC_late_smoking', 'CRC_late_nonsmoking']
        
        # 检查所有组是否都有数据
        available_groups = [g for g in groups if g in group_correlations and not group_correlations[g].empty]
        
        if len(available_groups) < 4:
            print(f"警告: 只有 {len(available_groups)} 个组有数据，需要4个组才能进行分析")
            return
        
        print(f"分析 {len(available_groups)} 个组: {available_groups}")
        
        # 创建所有组的合并相关性表（只使用p值 < threshold的数据）
        all_correlations = {}
        all_pvalues = {}
        
        for group_name in available_groups:
            correlations = group_correlations[group_name]
            for _, row in correlations.iterrows():
                # 只使用p值 < threshold的数据
                if row['PValue'] < p_threshold:
                    species_pair = tuple(sorted([row['Species1'], row['Species2']]))
                    if species_pair not in all_correlations:
                        all_correlations[species_pair] = {}
                        all_pvalues[species_pair] = {}
                    all_correlations[species_pair][group_name] = row['Correlation']
                    all_pvalues[species_pair][group_name] = row['PValue']
        
        # 找出在所有四个组中都存在的物种对
        highly_significant_pairs = []
        
        for species_pair, group_corrs in all_correlations.items():
            # 确保在所有四个组中都有数据
            if len(group_corrs) == 4:
                # 检查所有组的p值都 < threshold
                all_significant = all(p < p_threshold for p in all_pvalues[species_pair].values())
                
                if all_significant:
                    # 收集结果
                    result = {
                        'Species1': species_pair[0],
                        'Species2': species_pair[1]
                    }
                    
                    # 添加每个组的相关性值和p值
                    for group_name in available_groups:
                        result[f'{group_name}_Correlation'] = group_corrs.get(group_name, np.nan)
                        result[f'{group_name}_PValue'] = all_pvalues[species_pair].get(group_name, np.nan)
                    
                    # 计算统计信息
                    correlations = list(group_corrs.values())
                    pvalues = list(all_pvalues[species_pair].values())
                    
                    result['AvgCorrelation'] = np.mean(correlations)
                    result['MaxCorrelation'] = np.max(correlations)
                    result['MinCorrelation'] = np.min(correlations)
                    result['CorrelationRange'] = np.max(correlations) - np.min(correlations)
                    result['AvgPValue'] = np.mean(pvalues)
                    result['MinPValue'] = np.min(pvalues)
                    result['MaxPValue'] = np.max(pvalues)
                    
                    highly_significant_pairs.append(result)
        
        if not highly_significant_pairs:
            print(f"没有找到在所有四个组中都高度显著 (p < {p_threshold}) 的相关性对")
            return
        
        # 转换为DataFrame
        highly_significant_df = pd.DataFrame(highly_significant_pairs)
        
        # 按平均相关性排序
        highly_significant_df = highly_significant_df.sort_values('AvgCorrelation', ascending=False)
        
        print(f"\n找到 {len(highly_significant_df)} 个在所有四个组中都高度显著 (p < {p_threshold}) 的相关性对")
        
        # 显示前20个结果
        print(f"\n前20个高度显著的相关性对 (按平均相关性排序):")
        print("=" * 120)
        print("物种对 | 平均相关性 | 相关性范围 | 平均p值 | 最小p值 | 四组相关性值")
        print("-" * 120)
        
        for idx, row in highly_significant_df.head(20).iterrows():
            # 简化物种名称
            species1_short = row['Species1'].replace('tax_s__', '')[:20]
            species2_short = row['Species2'].replace('tax_s__', '')[:20]
            
            # 收集四组相关性值
            group_correlations_str = []
            for group in available_groups:
                corr = row[f'{group}_Correlation']
                group_correlations_str.append(f"{group[-2:]}:{corr:.3f}")
            
            print(f"{species1_short} <-> {species2_short} | {row['AvgCorrelation']:.3f} | {row['CorrelationRange']:.3f} | "
                  f"{row['AvgPValue']:.6f} | {row['MinPValue']:.6f} | {' '.join(group_correlations_str)}")
        
        print("-" * 120)
        
        # 保存结果
        if self.output_dir:
            # 保存完整表格
            filename = f"highly_significant_correlations_p{p_threshold}.csv"
            filepath = self.output_dir / filename
            highly_significant_df.to_csv(filepath, index=False)
            print(f"\n高度显著相关性表格已保存至: {filepath}")
            
            # 保存简化表格（只包含关键信息）
            simplified_cols = ['Species1', 'Species2', 'AvgCorrelation', 'CorrelationRange', 
                             'AvgPValue', 'MinPValue'] + [f'{g}_Correlation' for g in available_groups]
            simplified_df = highly_significant_df[simplified_cols]
            
            simplified_filename = f"highly_significant_correlations_simplified_p{p_threshold}.csv"
            simplified_filepath = self.output_dir / simplified_filename
            simplified_df.to_csv(simplified_filepath, index=False)
            print(f"简化表格已保存至: {simplified_filepath}")
        
        # 生成统计报告
        self._generate_highly_significant_report(highly_significant_df, p_threshold, available_groups)
    
    def _generate_highly_significant_report(self, df: pd.DataFrame, p_threshold: float, groups: List[str]):
        """
        生成高度显著相关性的统计报告
        """
        if self.output_dir:
            report_file = self.output_dir / f"highly_significant_report_p{p_threshold}.txt"
            
            with open(report_file, 'w', encoding='utf-8') as f:
                f.write(f"高度显著相关性分析报告 (p < {p_threshold})\n")
                f.write("=" * 80 + "\n\n")
                
                f.write(f"总相关性对数量: {len(df)}\n")
                f.write(f"分析组数: {len(groups)}\n")
                f.write(f"分析组: {', '.join(groups)}\n\n")
                
                # 基本统计
                f.write("基本统计信息:\n")
                f.write(f"  平均相关性: {df['AvgCorrelation'].mean():.4f} ± {df['AvgCorrelation'].std():.4f}\n")
                f.write(f"  相关性范围: {df['CorrelationRange'].mean():.4f} ± {df['CorrelationRange'].std():.4f}\n")
                f.write(f"  平均p值: {df['AvgPValue'].mean():.6f} ± {df['AvgPValue'].std():.6f}\n")
                f.write(f"  最小p值: {df['MinPValue'].min():.6f}\n\n")
                
                # 相关性强度分布
                f.write("相关性强度分布:\n")
                strong_corr = df[df['AvgCorrelation'] >= 0.7]
                moderate_corr = df[(df['AvgCorrelation'] >= 0.5) & (df['AvgCorrelation'] < 0.7)]
                weak_corr = df[df['AvgCorrelation'] < 0.5]
                
                f.write(f"  强相关性 (≥0.7): {len(strong_corr)} ({len(strong_corr)/len(df):.1%})\n")
                f.write(f"  中等相关性 (0.5-0.7): {len(moderate_corr)} ({len(moderate_corr)/len(df):.1%})\n")
                f.write(f"  弱相关性 (<0.5): {len(weak_corr)} ({len(weak_corr)/len(df):.1%})\n\n")
                
                # 前10个最强相关性
                f.write("前10个最强相关性对:\n")
                top_10 = df.nlargest(10, 'AvgCorrelation')
                for idx, row in top_10.iterrows():
                    species1_short = row['Species1'].replace('tax_s__', '')[:25]
                    species2_short = row['Species2'].replace('tax_s__', '')[:25]
                    f.write(f"  {species1_short} <-> {species2_short}: {row['AvgCorrelation']:.3f} (p={row['AvgPValue']:.6f})\n")
                
                # 组间一致性分析
                f.write("\n组间一致性分析:\n")
                consistent_pairs = df[df['CorrelationRange'] <= 0.2]
                variable_pairs = df[df['CorrelationRange'] > 0.2]
                
                f.write(f"  一致性高 (范围≤0.2): {len(consistent_pairs)} ({len(consistent_pairs)/len(df):.1%})\n")
                f.write(f"  变异性大 (范围>0.2): {len(variable_pairs)} ({len(variable_pairs)/len(df):.1%})\n")
            
            print(f"统计报告已保存至: {report_file}")
    
    def _save_smoking_difference_analysis(self, diff_df: pd.DataFrame, 
                                        top_diff_abs: pd.DataFrame, 
                                        top_diff_effect: pd.DataFrame,
                                        smoking_enhanced: pd.DataFrame,
                                        nonsmoking_enhanced: pd.DataFrame):
        """
        保存吸烟/非吸烟差异分析结果
        
        Args:
            diff_df: 完整的差异数据框
            top_diff_abs: 按绝对差异排序的前n个
            top_diff_effect: 按效应大小排序的前n个
            smoking_enhanced: 吸烟组相关性增强的前n个
            nonsmoking_enhanced: 非吸烟组相关性增强的前n个
        """
        # 保存完整的吸烟/非吸烟差异表
        smoking_diff_file = self.output_dir / "smoking_vs_nonsmoking_differences_full.csv"
        diff_df.to_csv(smoking_diff_file, index=False)
        print(f"\n吸烟/非吸烟完整差异表保存至: {smoking_diff_file}")
        
        # 保存按绝对差异排序的结果
        abs_diff_file = self.output_dir / "top_smoking_differences_by_abs.csv"
        top_diff_abs.to_csv(abs_diff_file, index=False)
        print(f"按绝对差异排序结果保存至: {abs_diff_file}")
        
        # 保存按效应大小排序的结果
        effect_diff_file = self.output_dir / "top_smoking_differences_by_effect.csv"
        top_diff_effect.to_csv(effect_diff_file, index=False)
        print(f"按效应大小排序结果保存至: {effect_diff_file}")
        
        # 保存吸烟组相关性增强的结果
        smoking_enhanced_file = self.output_dir / "smoking_enhanced_correlations.csv"
        smoking_enhanced.to_csv(smoking_enhanced_file, index=False)
        print(f"吸烟组相关性增强结果保存至: {smoking_enhanced_file}")
        
        # 保存非吸烟组相关性增强的结果
        nonsmoking_enhanced_file = self.output_dir / "nonsmoking_enhanced_correlations.csv"
        nonsmoking_enhanced.to_csv(nonsmoking_enhanced_file, index=False)
        print(f"非吸烟组相关性增强结果保存至: {nonsmoking_enhanced_file}")
        
        # 创建吸烟/非吸烟差异汇总表
        summary_data = []
        
        # 合并前50个绝对差异的结果
        for idx, row in top_diff_abs.iterrows():
            summary_data.append({
                'Species1': row['Species1'],
                'Species2': row['Species2'],
                'DifferenceType': 'AbsoluteDifference',
                'DifferenceValue': row['AbsDifference'],
                'EffectSize': row['EffectSize'],
                'Direction': 'SmokingEnhanced' if row['Difference'] > 0 else 'NonsmokingEnhanced',
                'Rank': idx + 1
            })
        
        # 合并前50个效应大小的结果
        for idx, row in top_diff_effect.iterrows():
            summary_data.append({
                'Species1': row['Species1'],
                'Species2': row['Species2'],
                'DifferenceType': 'EffectSize',
                'DifferenceValue': row['EffectSize'],
                'EffectSize': row['EffectSize'],
                'Direction': 'SmokingEnhanced' if row['Difference'] > 0 else 'NonsmokingEnhanced',
                'Rank': idx + 1
            })
        
        summary_df = pd.DataFrame(summary_data)
        summary_file = self.output_dir / "smoking_difference_summary.csv"
        summary_df.to_csv(summary_file, index=False)
        print(f"吸烟/非吸烟差异汇总表保存至: {summary_file}")
        
        # 生成吸烟/非吸烟差异统计报告
        self._generate_smoking_difference_report(diff_df, smoking_enhanced, nonsmoking_enhanced)
    
    def _generate_smoking_difference_report(self, diff_df: pd.DataFrame,
                                           smoking_enhanced: pd.DataFrame,
                                           nonsmoking_enhanced: pd.DataFrame):
        """
        生成吸烟/非吸烟差异统计报告
        
        Args:
            diff_df: 完整的差异数据框
            smoking_enhanced: 吸烟组相关性增强的数据
            nonsmoking_enhanced: 非吸烟组相关性增强的数据
        """
        report_file = self.output_dir / "smoking_difference_report.txt"
        
        with open(report_file, 'w', encoding='utf-8') as f:
            f.write("吸烟 vs 非吸烟相关性差异分析报告\n")
            f.write("=" * 50 + "\n\n")
            
            f.write(f"总分析物种对数量: {len(diff_df)}\n")
            f.write(f"吸烟组相关性增强的物种对数量: {len(smoking_enhanced)}\n")
            f.write(f"非吸烟组相关性增强的物种对数量: {len(nonsmoking_enhanced)}\n")
            f.write(f"吸烟组增强比例: {len(smoking_enhanced)/len(diff_df):.2%}\n")
            f.write(f"非吸烟组增强比例: {len(nonsmoking_enhanced)/len(diff_df):.2%}\n\n")
            
            # 基本统计
            f.write("差异统计:\n")
            f.write(f"  平均差异: {diff_df['Difference'].mean():.4f}\n")
            f.write(f"  平均绝对差异: {diff_df['AbsDifference'].mean():.4f}\n")
            f.write(f"  平均效应大小: {diff_df['EffectSize'].mean():.4f}\n")
            f.write(f"  最大正差异: {diff_df['Difference'].max():.4f}\n")
            f.write(f"  最大负差异: {diff_df['Difference'].min():.4f}\n\n")
            
            # 吸烟组增强的top物种对
            f.write("吸烟组相关性增强的Top 10物种对:\n")
            for idx, row in smoking_enhanced.head(10).iterrows():
                f.write(f"  {row['Species1']} <-> {row['Species2']}: ")
                f.write(f"吸烟组={row['Smoking_AvgCorrelation']:.3f}, ")
                f.write(f"非吸烟组={row['Nonsmoking_AvgCorrelation']:.3f}, ")
                f.write(f"差异={row['Difference']:.3f}, ")
                f.write(f"效应大小={row['EffectSize']:.3f}\n")
            
            f.write("\n非吸烟组相关性增强的Top 10物种对:\n")
            for idx, row in nonsmoking_enhanced.head(10).iterrows():
                f.write(f"  {row['Species1']} <-> {row['Species2']}: ")
                f.write(f"吸烟组={row['Smoking_AvgCorrelation']:.3f}, ")
                f.write(f"非吸烟组={row['Nonsmoking_AvgCorrelation']:.3f}, ")
                f.write(f"差异={row['Difference']:.3f}, ")
                f.write(f"效应大小={row['EffectSize']:.3f}\n")
        
        print(f"吸烟/非吸烟差异分析报告保存至: {report_file}")
    
    def _save_difference_analysis(self, diff_df: pd.DataFrame, 
                                 top_diff_max: pd.DataFrame, 
                                 top_diff_avg: pd.DataFrame, 
                                 top_diff_std: pd.DataFrame):
        """
        保存差异分析结果
        
        Args:
            diff_df: 完整的差异数据框
            top_diff_max: 按最大差异排序的前n个
            top_diff_avg: 按平均差异排序的前n个
            top_diff_std: 按标准差排序的前n个
        """
        # 保存完整的差异表
        diff_file = self.output_dir / "correlation_differences_full.csv"
        diff_df.to_csv(diff_file, index=False)
        print(f"\n完整差异表保存至: {diff_file}")
        
        # 保存按最大差异排序的结果
        max_diff_file = self.output_dir / "top_differences_by_max.csv"
        top_diff_max.to_csv(max_diff_file, index=False)
        print(f"按最大差异排序结果保存至: {max_diff_file}")
        
        # 保存按平均差异排序的结果
        avg_diff_file = self.output_dir / "top_differences_by_avg.csv"
        top_diff_avg.to_csv(avg_diff_file, index=False)
        print(f"按平均差异排序结果保存至: {avg_diff_file}")
        
        # 保存按标准差排序的结果
        std_diff_file = self.output_dir / "top_differences_by_std.csv"
        top_diff_std.to_csv(std_diff_file, index=False)
        print(f"按变异性排序结果保存至: {std_diff_file}")
        
        # 创建汇总表
        summary_data = []
        
        # 合并前50个最大差异的结果
        for idx, row in top_diff_max.iterrows():
            summary_data.append({
                'Species1': row['Species1'],
                'Species2': row['Species2'],
                'DifferenceType': 'MaxDifference',
                'DifferenceValue': row['MaxDifference'],
                'Rank': idx + 1
            })
        
        # 合并前50个平均差异的结果
        for idx, row in top_diff_avg.iterrows():
            summary_data.append({
                'Species1': row['Species1'],
                'Species2': row['Species2'],
                'DifferenceType': 'AvgPairwiseDifference',
                'DifferenceValue': row['AvgPairwiseDifference'],
                'Rank': idx + 1
            })
        
        summary_df = pd.DataFrame(summary_data)
        summary_file = self.output_dir / "difference_summary.csv"
        summary_df.to_csv(summary_file, index=False)
        print(f"差异汇总表保存至: {summary_file}")


def main():
    """
    主函数：分组并行相关性计算
    """
    # 初始化分组并行计算器
    calculator = GroupParallelCorrelation(
        n_jobs=-1,  # 使用所有CPU
        correlation_method='spearman',
        output_dir='results/group_correlations'
    )
    
    # 加载分组数据
    print("=" * 60)
    group_data = calculator.load_group_data()
    
    # 检查各组样本数量
    valid_groups = {name: data for name, data in group_data.items() if not data.empty}
    if not valid_groups:
        print("错误: 没有找到有效的分组数据")
        return
    
    print(f"\n有效分组: {list(valid_groups.keys())}")
    
    # 计算各组相关性
    print("=" * 60)
    group_correlations = calculator.calculate_group_correlations(
        valid_groups,
        top_k=200,  # 每个组计算前150个重要物种
        correlation_threshold=0.4  # 只保留强相关性
    )
    
    # 分析结果
    print("=" * 60)
    calculator.analyze_group_correlations(group_correlations)
    
    # 比较各组并分析差异
    calculator.compare_groups(group_correlations, top_n=50)
        
    # 保存结果
    print("=" * 60)
    calculator.save_group_results(group_correlations)
    
    print("\n" + "=" * 60)
    print("分组相关性分析完成!")
    print(f"结果保存在: {calculator.output_dir}")


if __name__ == "__main__":
    main()