import os
import pandas as pd
import numpy as np
from dataset import BioSmokeDataset
from sklearn.preprocessing import StandardScaler
from sklearn.cross_decomposition import PLSRegression
from sklearn.model_selection import cross_val_predict, LeaveOneOut
from sklearn.metrics import accuracy_score, confusion_matrix, classification_report
from scipy import stats
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import Ellipse
import warnings
warnings.filterwarnings('ignore')
from split_group import prepare_crc_diff_groups


def oplsda_analysis(df: pd.DataFrame, prefix: str, group1: str, group2: str, 
                   output_dir: str = 'oplsda_figures', n_components: int = 2):
    """
    执行OPLS-DA分析
    
    参数:
        df: 包含特征和分组的数据框
        prefix: 特征前缀 ('tax_' for species, 'met_' for metabolites)
        group1: 第一组名称
        group2: 第二组名称
        output_dir: 输出目录
        n_components: PLS成分数
    
    返回:
        OPLS-DA分析结果字典
    """
    # 筛选两组数据
    df_pair = df[df['group'].isin([group1, group2])].copy()
    
    if df_pair.empty:
        print(f"警告: 没有找到 {group1} 或 {group2} 的数据")
        return None
    
    print(f"\n分析 {group1} vs {group2}: {len(df_pair)} 个样本")
    
    # 提取特征
    # 使用列索引来选择特征，避免列名重复导致的问题
    prefix_column_indices = [i for i, c in enumerate(df_pair.columns) if c.startswith(prefix)]
    # 确保索引唯一
    unique_indices = list(dict.fromkeys(prefix_column_indices))
    
    if len(unique_indices) == 0:
        print(f'没有找到以 "{prefix}" 开头的特征列')
        return None
    
    # 选择特征列
    selected_df = df_pair.iloc[:, unique_indices]
    # 更新features列表为实际选择的列名
    features = list(selected_df.columns)
    
    print(f"  特征数: {len(features)}")
    
    # 准备数据
    X = selected_df.values.astype(float)
    y = df_pair['group'].values
    
    # 将组标签转换为数值（0和1）
    y_numeric = np.where(y == group1, 0, 1)
    
    # 标准化数据
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    
    # 执行PLS-DA（OPLS-DA的简化版本）
    print("  执行PLS-DA分析...")
    plsda = PLSRegression(n_components=n_components)
    plsda.fit(X_scaled, y_numeric)
    
    # 获取得分
    X_scores = plsda.x_scores_
    X_loadings = plsda.x_loadings_
    
    # 交叉验证
    print("  执行留一法交叉验证...")
    loo = LeaveOneOut()
    y_pred_cv = cross_val_predict(plsda, X_scaled, y_numeric, cv=loo)
    
    # 计算预测准确率
    accuracy = accuracy_score(y_numeric, np.round(y_pred_cv))
    
    # 计算VIP值（变量重要性投影）
    print("  计算VIP值...")
    vip_scores = calculate_vip(plsda, X_scaled, y_numeric)
    
    # 执行置换检验
    print("  执行置换检验（1000次置换）...")
    permutation_results = permutation_test(plsda, X_scaled, y_numeric, n_permutations=1000)
    
    # 创建结果字典
    results = {
        'groups': (group1, group2),
        'n_samples': len(df_pair),
        'n_features': len(features),
        'accuracy': accuracy,
        'X_scores': X_scores,
        'X_loadings': X_loadings,
        'vip_scores': vip_scores,
        'feature_names': features,
        'y_true': y_numeric,
        'y_pred': y_pred_cv,
        'permutation_results': permutation_results
    }
    
    # 绘制OPLS-DA得分图
    plot_oplsda_scores(results, prefix, group1, group2, output_dir)
    
    # 绘制VIP图
    plot_vip_scores(results, prefix, group1, group2, output_dir)
    
    # 保存重要特征
    save_important_features(results, prefix, group1, group2, output_dir)
    
    return results


def calculate_vip(plsda, X, y):
    """
    计算VIP（变量重要性投影）值
    
    参数:
        plsda: 训练好的PLS模型
        X: 特征矩阵
        y: 目标变量
    
    返回:
        VIP值数组
    """
    t = plsda.x_scores_  # 得分矩阵
    w = plsda.x_weights_  # 权重矩阵
    q = plsda.y_loadings_  # Y载荷
    
    # 计算每个成分的VIP
    vip = np.zeros((X.shape[1],))
    
    # 计算总方差解释
    s = np.diag(t.T @ t @ q.T @ q).reshape(-1, 1)
    
    # 计算VIP值
    for i in range(X.shape[1]):
        weight = np.array([(w[i, j] / np.linalg.norm(w[:, j])) for j in range(w.shape[1])])
        vip[i] = np.sqrt(X.shape[1] * np.sum(s * weight**2) / np.sum(s))
    
    return vip

def permutation_test(plsda, X_scaled, y_numeric, n_permutations: int = 1000, random_state: int = 42):
    """
    执行置换检验来评估OPLS-DA模型的显著性
    
    参数:
        plsda: 训练好的PLS模型
        X_scaled: 标准化后的特征矩阵
        y_numeric: 数值型目标变量
        n_permutations: 置换次数（默认1000）
        random_state: 随机种子
    
    返回:
        包含置换检验结果的字典
    """
    np.random.seed(random_state)
    
    original_score = plsda.score(X_scaled, y_numeric)
    
    permuted_scores = []
    
    for i in range(n_permutations):
        y_permuted = np.random.permutation(y_numeric)
        
        plsda_perm = PLSRegression(n_components=plsda.n_components)
        plsda_perm.fit(X_scaled, y_permuted)
        
        permuted_score = plsda_perm.score(X_scaled, y_permuted)
        permuted_scores.append(permuted_score)
        
        if (i + 1) % 100 == 0:
            print(f"    置换进度: {i + 1}/{n_permutations}")
    
    permuted_scores = np.array(permuted_scores)
    
    p_value = (np.sum(permuted_scores >= original_score) + 1) / (n_permutations + 1)
    
    return {
        'original_score': original_score,
        'permuted_scores': permuted_scores,
        'p_value': p_value,
        'n_permutations': n_permutations
    }

def plot_oplsda_scores(results, prefix: str, group1: str, group2: str, output_dir: str):
    """
    绘制OPLS-DA得分图
    
    参数:
        results: OPLS-DA结果字典
        prefix: 特征前缀
        group1: 第一组名称
        group2: 第二组名称
        output_dir: 输出目录
    """
    X_scores = results['X_scores']
    y_true = results['y_true']
    
    # 创建图形
    plt.figure(figsize=(12, 10))
    
    # 定义配色方案（与pcoa_plot.py相同）
    colors = {
            'CRC_poordiff': '#FF6B6B',
            'CRC_welldiff': '#FFD166',
            'CTRL': '#06D6A0'
        }
    
    # 绘制得分散点图
    for i, group_idx in enumerate([0, 1]):
        group_mask = y_true == group_idx
        group_label = group1 if group_idx == 0 else group2
        # 根据组名选择颜色
        color = colors.get(group_label, '#999999')  # 默认灰色
        plt.scatter(X_scores[group_mask, 0], X_scores[group_mask, 1],
                   color=color, marker='o', s=100,
                   label=group_label, alpha=0.8, linewidth=0)
    
    # 添加置信椭圆（95%）
    for i, group_idx in enumerate([0, 1]):
        group_mask = y_true == group_idx
        if np.sum(group_mask) >= 2:  # 至少需要2个点才能计算椭圆
            group_scores = X_scores[group_mask, :2]
            mean = group_scores.mean(axis=0)
            cov = np.cov(group_scores.T)
            
            # 计算椭圆参数
            try:
                vals, vecs = np.linalg.eigh(cov)
                order = vals.argsort()[::-1]
                vals = vals[order]
                vecs = vecs[:, order]
                
                # 95%置信椭圆（卡方分布，2个自由度）
                chi2_val = 5.991  # chi2.ppf(0.95, 2)
                angle = np.degrees(np.arctan2(vecs[1, 0], vecs[0, 0]))
                width, height = 2 * np.sqrt(vals * chi2_val)
                
                # 根据组名选择颜色
                group_label = group1 if group_idx == 0 else group2
                color = colors.get(group_label, '#999999')  # 默认灰色
                
                ellipse = Ellipse(xy=mean, width=width, height=height, angle=angle,
                                  edgecolor='none', facecolor=color,
                                  lw=0, alpha=0.3)
                plt.gca().add_patch(ellipse)
            except:
                pass
    
    # 设置图表属性
    plt.xlabel(f'OPLS-DA Component 1', fontsize=14)
    plt.ylabel(f'OPLS-DA Component 2', fontsize=14)
    
    title = f'OPLS-DA Scores Plot: {group1} vs {group2}'
    if prefix == 'tax_':
        title += ' (Species)'
    else:
        title += ' (Metabolites)'
    
    plt.title(title, fontsize=16, fontweight='bold')
    
    # 添加准确率信息
    # accuracy = results['accuracy']
    # plt.text(0.02, 0.98, f'Accuracy: {accuracy:.3f}', transform=plt.gca().transAxes,
    #          fontsize=12, fontweight='bold', va='top', ha='left',
    #          bbox=dict(facecolor='white', alpha=0.8, edgecolor='black'))
    
    # 添加置换检验p值信息
    if 'permutation_results' in results:
        p_value = results['permutation_results']['p_value']
        p_text = f'Permutation p-value: {p_value:.4f}' if p_value >= 0.0001 else 'Permutation p-value: < 0.0001'
        plt.text(0.02, 0.98, p_text, transform=plt.gca().transAxes,
                 fontsize=12, fontweight='bold', va='top', ha='left',
                 bbox=dict(facecolor='white', alpha=0.8, edgecolor='black'))
    
    # 添加图例
    plt.legend(loc='upper right', fontsize=12, framealpha=0.9)
    
    # 添加网格
    # plt.grid(True, alpha=0.3, linestyle='--')
    
    # 添加零线
    # plt.axhline(y=0, color='k', linestyle='-', alpha=0.3)
    # plt.axvline(x=0, color='k', linestyle='-', alpha=0.3)
    
    plt.tight_layout()
    
    # 保存图片
    safe_group1 = group1.replace('CRC_', '').replace('_', '-')
    safe_group2 = group2.replace('CRC_', '').replace('_', '-')
    filename_suffix = 'species' if prefix == 'tax_' else 'metabolites'
    out_png = os.path.join(output_dir, f'oplsda_{filename_suffix}_{safe_group1}_vs_{safe_group2}_scores.png')
    
    plt.savefig(out_png, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"  已保存得分图: {out_png}")


def plot_vip_scores(results, prefix: str, group1: str, group2: str, output_dir: str, top_n: int = 20):
    """
    绘制VIP（变量重要性投影）图
    
    参数:
        results: OPLS-DA结果字典
        prefix: 特征前缀
        group1: 第一组名称
        group2: 第二组名称
        output_dir: 输出目录
        top_n: 显示前N个重要特征
    """
    vip_scores = results['vip_scores']
    feature_names = results['feature_names']
    
    # 确保feature_names和vip_scores长度一致
    min_len = min(len(feature_names), len(vip_scores))
    if min_len < len(feature_names) or min_len < len(vip_scores):
        print(f"  警告: feature_names和vip_scores长度不一致，仅使用前{min_len}个元素")
        feature_names = feature_names[:min_len]
        vip_scores = vip_scores[:min_len]
    
    # 创建VIP值数据框
    vip_df = pd.DataFrame({
        'Feature': feature_names,
        'VIP': vip_scores
    })
    
    # 按VIP值降序排序
    vip_df = vip_df.sort_values('VIP', ascending=False)
    
    # 选择前N个特征
    top_features = vip_df.head(top_n)
    
    # 创建图形
    plt.figure(figsize=(14, 10))
    
    # 绘制水平条形图
    colors = plt.cm.viridis(np.linspace(0.3, 0.9, len(top_features)))
    bars = plt.barh(range(len(top_features)), top_features['VIP'].values, color=colors)
    
    # 设置y轴标签
    plt.yticks(range(len(top_features)), top_features['Feature'].values, fontsize=10)
    plt.xlabel('VIP Score', fontsize=14)
    plt.ylabel('Features', fontsize=14)
    
    # 添加标题
    title = f'Top {top_n} VIP Scores: {group1} vs {group2}'
    if prefix == 'tax_':
        title += ' (Species)'
    else:
        title += ' (Metabolites)'
    
    plt.title(title, fontsize=16, fontweight='bold')
    
    # 添加VIP=1的参考线（通常认为VIP>1的特征是重要的）
    plt.axvline(x=1.0, color='r', linestyle='--', alpha=0.7, linewidth=1.5)
    plt.text(1.02, len(top_features) * 0.95, 'VIP = 1.0', color='r', 
             fontsize=11, fontweight='bold', va='top')
    
    # 添加数值标签
    for i, (vip, bar) in enumerate(zip(top_features['VIP'].values, bars)):
        plt.text(vip + 0.02, bar.get_y() + bar.get_height()/2, 
                f'{vip:.3f}', va='center', fontsize=9)
    
    plt.gca().invert_yaxis()  # 反转y轴使最高VIP值在顶部
    plt.tight_layout()
    
    # 保存图片
    safe_group1 = group1.replace('CRC_', '').replace('_', '-')
    safe_group2 = group2.replace('CRC_', '').replace('_', '-')
    filename_suffix = 'species' if prefix == 'tax_' else 'metabolites'
    out_png = os.path.join(output_dir, f'oplsda_{filename_suffix}_{safe_group1}_vs_{safe_group2}_vip.png')
    
    plt.savefig(out_png, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"  已保存VIP图: {out_png}")


def save_important_features(results, prefix: str, group1: str, group2: str, output_dir: str, vip_threshold: float = 1.0):
    """
    保存重要特征信息
    
    参数:
        results: OPLS-DA结果字典
        prefix: 特征前缀
        group1: 第一组名称
        group2: 第二组名称
        output_dir: 输出目录
        vip_threshold: VIP阈值（通常VIP>1被认为是重要的）
    """
    vip_scores = results['vip_scores']
    feature_names = results['feature_names']
    
    # 确保feature_names和vip_scores长度一致
    min_len = min(len(feature_names), len(vip_scores))
    if min_len < len(feature_names) or min_len < len(vip_scores):
        print(f"  警告: feature_names和vip_scores长度不一致，仅使用前{min_len}个元素")
        feature_names = feature_names[:min_len]
        vip_scores = vip_scores[:min_len]
    
    # 创建特征重要性数据框
    feature_df = pd.DataFrame({
        'Feature': feature_names,
        'VIP_Score': vip_scores,
        'Important': vip_scores > vip_threshold
    })
    
    # 按VIP值降序排序
    feature_df = feature_df.sort_values('VIP_Score', ascending=False)
    
    # 保存到CSV文件
    safe_group1 = group1.replace('CRC_', '').replace('_', '-')
    safe_group2 = group2.replace('CRC_', '').replace('_', '-')
    filename_suffix = 'species' if prefix == 'tax_' else 'metabolites'
    out_csv = os.path.join(output_dir, f'oplsda_{filename_suffix}_{safe_group1}_vs_{safe_group2}_features.csv')
    
    feature_df.to_csv(out_csv, index=False)
    
    # 统计重要特征数量
    n_important = feature_df['Important'].sum()
    print(f"  重要特征数 (VIP > {vip_threshold}): {n_important}/{len(feature_df)}")
    print(f"  特征数据已保存: {out_csv}")
    
    return feature_df

def calculate_volcano_data(df: pd.DataFrame, prefix: str, group1: str, group2: str):
    """
    计算火山图所需数据（log2 fold change和调整p值）
    
    参数:
        df: 包含特征和分组的数据框
        prefix: 特征前缀 ('met_' for metabolites)
        group1: 第一组名称
        group2: 第二组名称
    
    返回:
        包含log2FC和调整p值的数据框
    """
    df_pair = df[df['group'].isin([group1, group2])].copy()
    
    if df_pair.empty:
        return None
    
    prefix_column_indices = [i for i, c in enumerate(df_pair.columns) if c.startswith(prefix)]
    unique_indices = list(dict.fromkeys(prefix_column_indices))
    
    if len(unique_indices) == 0:
        return None
    
    selected_df = df_pair.iloc[:, unique_indices]
    features = list(selected_df.columns)
    
    X = selected_df.values.astype(float)
    y = df_pair['group'].values
    
    group1_mask = y == group1
    group2_mask = y == group2
    
    volcano_data = []
    
    for i, feature in enumerate(features):
        group1_values = X[group1_mask, i]
        group2_values = X[group2_mask, i]
        
        mean1 = np.mean(group1_values)
        mean2 = np.mean(group2_values)
        
        log2fc = np.log2(mean2 / mean1 + 1e-10)
        
        t_stat, p_value = stats.ttest_ind(group1_values, group2_values)
        
        volcano_data.append({
            'Feature': feature,
            'log2FC': log2fc,
            'p_value': p_value
        })
    
    volcano_df = pd.DataFrame(volcano_data)
    
    volcano_df['adj_p_value'] = stats.false_discovery_control(volcano_df['p_value'], method='bh')
    
    return volcano_df

def plot_single_volcano(volcano_df: pd.DataFrame, vip_scores: np.ndarray, feature_names: list, 
                       group1: str, group2: str, ax=None, top_n: int = 10, 
                       show_yaxis: bool = True, show_spines: bool = True):
    """
    绘制单个火山图，标注VIP前10代谢物
    
    参数:
        volcano_df: 包含火山图数据的数据框
        vip_scores: VIP分数数组
        feature_names: 特征名称列表
        group1: 第一组名称
        group2: 第二组名称
        ax: matplotlib轴对象
        top_n: 标注前N个VIP特征
        show_yaxis: 是否显示Y轴
        show_spines: 是否显示坐标轴方框
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 10))
    
    volcano_df['neg_log10_adj_p'] = -np.log10(volcano_df['adj_p_value'] + 1e-10)
    
    volcano_df['color'] = '#999999'
    volcano_df.loc[(volcano_df['log2FC'] > 0) & (volcano_df['adj_p_value'] >= 0.05), 'color'] = '#FFB6B6'
    volcano_df.loc[(volcano_df['log2FC'] > 1) & (volcano_df['adj_p_value'] < 0.05), 'color'] = '#DC143C'
    volcano_df.loc[(volcano_df['log2FC'] < 0) & (volcano_df['adj_p_value'] >= 0.05), 'color'] = '#ADD8E6'
    volcano_df.loc[(volcano_df['log2FC'] < -1) & (volcano_df['adj_p_value'] < 0.05), 'color'] = '#00008B'
    
    ax.scatter(volcano_df['neg_log10_adj_p'], volcano_df['log2FC'], 
               c=volcano_df['color'], alpha=0.6, s=50, edgecolors='none')
    
    # ax.axvline(x=-np.log10(0.05), color='gray', linestyle='--', alpha=0.5)
    ax.axhline(y=1, color='gray', linestyle='--', alpha=0.5)
    ax.axhline(y=-1, color='gray', linestyle='--', alpha=0.5)
    
    feature_vip_dict = dict(zip(feature_names, vip_scores))
    volcano_df['VIP'] = volcano_df['Feature'].map(feature_vip_dict)
    top_vip_features = volcano_df.nlargest(top_n, 'VIP')
    
    for _, row in top_vip_features.iterrows():
        feature_name = row['Feature'].replace('met_', '')
        ax.annotate(feature_name, 
                   (row['neg_log10_adj_p'], row['log2FC']),
                   xytext=(5, 5), textcoords='offset points',
                   fontsize=7, alpha=0.8,
                   arrowprops=dict(arrowstyle='->', color='black', alpha=0.3))
    
    # ax.set_xlabel('-log10(Adjusted P-value)', fontsize=12)
    # ax.set_ylabel('log2 Fold Change', fontsize=12)
    ax.set_ylim(-2, 2)
    ax.set_xlim(0,0.75)
    
    title = f'{group1.replace("CRC_", "")} vs {group2.replace("CRC_", "")}'
    ax.set_title(title, fontsize=14, fontweight='bold')
    
    if not show_spines:
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.spines['left'].set_visible(False)
    
    if not show_yaxis:
        ax.set_ylabel('')
        ax.set_yticklabels([])
        ax.tick_params(left=False)
    
    return ax

def main():
    # 加载数据
    ds = BioSmokeDataset()
    ds.preprocess_metabolomics_data(pqn_normalization=True, transform=True, transform_method='log')
    merged = ds.merge_to_dataframe()
    
    print(f"总样本数: {merged.shape[0]}")
    print(f"总特征数: {merged.shape[1]}")
    
    # 准备分组
    merged = prepare_crc_diff_groups(merged)
    
    # 打印各组样本数
    print("\n各组样本数:")
    groups = ['CTRL', 'CRC_welldiff', 'CRC_poordiff']
    group_counts = merged['group'].value_counts()
    for group in groups:
        if group in group_counts:
            print(f"  {group}: {group_counts[group]} 样本")
        else:
            print(f"  {group}: 0 样本")
    
    # 确保输出目录存在
    output_dir = './results/oplsda_plots'
    os.makedirs(output_dir, exist_ok=True)
    
    # 定义比较组
    comparisons = [
        ('CRC_poordiff', 'CRC_welldiff'),
        ('CRC_poordiff', 'CTRL'),
        ('CRC_welldiff', 'CTRL'),
    ]
    
    results_dict = {}
    
    # 对每个比较组执行分析
    for group1, group2 in comparisons:
        print(f"\n{'='*60}")
        print(f"分析比较: {group1} vs {group2}")
        print(f"{'='*60}")
        
        # 执行代谢物分析
        print("\n代谢物分析:")
        met_results = oplsda_analysis(
            merged, 
            prefix='met_', 
            group1=group1, 
            group2=group2, 
            output_dir=output_dir
        )
        
        if met_results is not None:
            key = f'{group1}_vs_{group2}'
            results_dict[key] = met_results
    
    print(f"\n{'='*60}")
    print("分析完成！")
    print(f"{'='*60}")

if __name__ == "__main__":
    main()
