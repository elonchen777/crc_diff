import pandas as pd
import numpy as np
import re
from typing import Dict, List, Optional, Tuple, Set
from scipy import stats
from sklearn.preprocessing import StandardScaler, MinMaxScaler


class BioSmokeDataset:
    def __init__(self, 
                 sample_file: str = 'dataset/id_sample.xlsx',
                 taxonomy_file: str = 'dataset/taxonomy_Species_abund.txt',
                 metabolomics_file: str = 'dataset/metabolome_data.csv',
                 kegg_file: str = 'dataset/KEGG/4_KOEntry/KEGG_KOEntry_abund.txt',
                 load_kegg: bool = False):   
        """
        初始化数据集
        
        Args:
            sample_file: 包含样本信息的Excel文件路径
            taxonomy_file: 宏基因组数据文件路径
            metabolomics_file: 代谢组数据文件路径
            kegg_file: KEGG KO丰度数据文件路径
            load_kegg: 是否读取KEGG数据（默认False）
        """
        self.sample_file = sample_file
        self.taxonomy_file = taxonomy_file
        self.metabolomics_file = metabolomics_file
        self.kegg_file = kegg_file
        self.load_kegg = load_kegg
        
        # 存储数据
        self.sample_ids: List[str] = []
        self.taxonomy_data: Optional[pd.DataFrame] = None
        self.metabolomics_data: Optional[pd.DataFrame] = None
        self.kegg_data: Optional[pd.DataFrame] = None
        self.smoking_labels: Dict[str, int] = {}
        self.gender_labels: Dict[str, int] = {}
        self.tnm_labels: Dict[str, int] = {}
        self.age_labels: Dict[str, float] = {}
        self.differentiation_labels: Dict[str, int] = {}
        
        # 加载数据
        self._load_data()
    
    def _load_data(self) -> None:
        """加载所有数据"""
        self._load_sample_ids()
        self._load_taxonomy_data()
        self._load_metabolomics_data()
        if self.load_kegg:
            self._load_kegg_data()
        self._extract_smoking_labels()
        self._extract_gender_labels()
        self._extract_tnm_labels()
        self._extract_age_labels()
        self._extract_differentiation_labels()
    
    def preprocess_taxonomy_data(self, remove_low_expression: bool = True,
                                  min_abund_threshold: float = 0.01,
                                  min_prevalence_threshold: float = 0.1,
                                  remove_outliers: bool = False,
                                  outlier_method: str = 'iqr',
                                  outlier_threshold: float = 3.0,
                                  relative_abund: bool = True,
                                  transform: bool = False,
                                  transform_method: str = 'log') -> None:
        """预处理宏基因组数据"""
        if self.taxonomy_data is None:
            return
            
        print("\n预处理宏基因组数据...")
        original_shape = self.taxonomy_data.shape
        
        # 去除重复的物种
        species_col = self.taxonomy_data.columns[0]
        self.taxonomy_data = self.taxonomy_data.drop_duplicates(subset=[species_col], keep='first')
        sample_cols = self.taxonomy_data.columns[1:]
        
        # 去除低表达物种（在大多数样本中丰度都很低的物种）
        if remove_low_expression:
            # 获取样本列（除第一列外的所有列）
            sample_cols = self.taxonomy_data.columns[1:]
            self.taxonomy_data[sample_cols] = self.taxonomy_data[sample_cols].apply(pd.to_numeric, errors='coerce')
            
            # 计算每个物种在多少样本中表达
            expression_counts = (self.taxonomy_data[sample_cols] > min_abund_threshold).sum(axis=1)
            
            # 保留在至少min_prevalence_threshold个样本中表达的物种
            self.taxonomy_data = self.taxonomy_data[expression_counts >= min_prevalence_threshold * len(sample_cols)]
            
            print(f"去除低表达物种: 保留 {self.taxonomy_data.shape[0]} 个物种")
        
        # 去除异常值
        if remove_outliers:
            sample_cols = self.taxonomy_data.columns[1:]
            self.taxonomy_data[sample_cols] = self.taxonomy_data[sample_cols].apply(pd.to_numeric, errors='coerce')
            outlier_count = 0
            
            if outlier_method == 'iqr':
                for col in sample_cols:
                    Q1 = self.taxonomy_data[col].quantile(0.25)
                    Q3 = self.taxonomy_data[col].quantile(0.75)
                    IQR = Q3 - Q1
                    lower_bound = Q1 - outlier_threshold * IQR
                    upper_bound = Q3 + outlier_threshold * IQR
                    
                    outliers = (self.taxonomy_data[col] < lower_bound) | (self.taxonomy_data[col] > upper_bound)
                    self.taxonomy_data.loc[outliers, col] = self.taxonomy_data.loc[outliers, col].clip(lower_bound, upper_bound)
                    outlier_count += outliers.sum()
                    
            elif outlier_method == 'zscore':
                for col in sample_cols:
                    z_scores = np.abs(stats.zscore(self.taxonomy_data[col], nan_policy='omit'))
                    outliers = z_scores > outlier_threshold
                    
                    median_val = self.taxonomy_data[col].median()
                    self.taxonomy_data.loc[outliers, col] = median_val
                    outlier_count += outliers.sum()
            
            print(f"去除异常值 ({outlier_method}方法): 处理了 {outlier_count} 个异常值")
        
        if relative_abund:
            # 相对丰度转换（按列/样本归一化）
            self.taxonomy_data[sample_cols] = self.taxonomy_data[sample_cols].apply(pd.to_numeric, errors='coerce')
            col_sums = self.taxonomy_data[sample_cols].sum(axis=0)
            col_sums = col_sums.replace(0, 1)
            self.taxonomy_data[sample_cols] = self.taxonomy_data[sample_cols].div(col_sums, axis=1)*100
            print(f"宏基因组数据归一化 (相对丰度转换)")
        
        if transform:
            sample_cols = self.taxonomy_data.columns[1:]
            self.taxonomy_data[sample_cols] = self.taxonomy_data[sample_cols].apply(pd.to_numeric, errors='coerce')
            
            if transform_method == 'clr':
                # 真正的CLR变换：按样本列计算 log(x + pseudocount) - mean(log(x + pseudocount))
                pseudo_count = 1e-9
                clipped = self.taxonomy_data[sample_cols].clip(lower=0)
                logged = np.log(clipped + pseudo_count)
                self.taxonomy_data[sample_cols] = logged.sub(logged.mean(axis=0), axis=1)
                print(f"宏基因组数据变换 (CLR)")
                
            elif transform_method == 'log':
                # 对数变换
                self.taxonomy_data[sample_cols] = np.log1p(self.taxonomy_data[sample_cols].clip(lower=0))
        
        print(f"宏基因组数据预处理完成: {original_shape} -> {self.taxonomy_data.shape}\n")
    
    def preprocess_metabolomics_data(self, remove_high_missing: bool = True,
                                     missing_threshold: float = 0.5,
                                     min_samples: int = None,
                                     remove_outliers: bool = False,
                                     outlier_method: str = 'iqr',
                                     outlier_threshold: float = 3.0,
                                     relative_abund: bool = False,
                                     pqn_normalization: bool = True,
                                     transform: bool = False,
                                     transform_method: str = 'log',
                                     scale: bool = False) -> None:
        """预处理代谢组数据"""
        if self.metabolomics_data is None:
            return
            
        print("\n预处理代谢组数据...")
        
        meta_data = self.metabolomics_data.copy()
        if meta_data.shape[0] >= 1 and meta_data.shape[1] >= 1:
            first_cell = str(meta_data.iloc[0, 0]).strip()
            if first_cell.lower() in ('lable', 'label'):
                meta_data = meta_data.iloc[1:].copy()
        
        original_shape = meta_data.shape
        
        metab_col = meta_data.columns[0]
        meta_data = meta_data.drop_duplicates(subset=[metab_col], keep='first')
        sample_cols = meta_data.columns[1:]
        
        if remove_high_missing:
            meta_data[sample_cols] = meta_data[sample_cols].apply(pd.to_numeric, errors='coerce')
            n_samples = len(sample_cols)
            
            missing_counts = meta_data[sample_cols].isna().sum(axis=1)
            missing_ratio = missing_counts / n_samples
            
            if min_samples is not None:
                kept_rows = (n_samples - missing_counts) >= min_samples
            else:
                kept_rows = missing_ratio < missing_threshold
            
            meta_data = meta_data[kept_rows]
            
            print(f"去除高缺失代谢物: 保留 {meta_data.shape[0]} 个代谢物 (缺失率阈值: {missing_threshold*100:.1f}%)")
        
        if remove_outliers:
            sample_cols = meta_data.columns[1:]
            meta_data[sample_cols] = meta_data[sample_cols].apply(pd.to_numeric, errors='coerce')
            outlier_count = 0
            
            if outlier_method == 'iqr':
                for col in sample_cols:
                    Q1 = meta_data[col].quantile(0.25)
                    Q3 = meta_data[col].quantile(0.75)
                    IQR = Q3 - Q1
                    lower_bound = Q1 - outlier_threshold * IQR
                    upper_bound = Q3 + outlier_threshold * IQR
                    
                    outliers = (meta_data[col] < lower_bound) | (meta_data[col] > upper_bound)
                    meta_data.loc[outliers, col] = meta_data.loc[outliers, col].clip(lower_bound, upper_bound)
                    outlier_count += outliers.sum()
                    
            elif outlier_method == 'zscore':
                for col in sample_cols:
                    z_scores = np.abs(stats.zscore(meta_data[col], nan_policy='omit'))
                    outliers = z_scores > outlier_threshold
                    
                    median_val = meta_data[col].median()
                    meta_data.loc[outliers, col] = median_val
                    outlier_count += outliers.sum()
            
            print(f"去除异常值 ({outlier_method}方法): 处理了 {outlier_count} 个异常值")
        
        if relative_abund:
            # 相对丰度转换（按列/样本归一化）
            meta_data[sample_cols] = meta_data[sample_cols].apply(pd.to_numeric, errors='coerce')
            col_sums = meta_data[sample_cols].sum(axis=0)
            col_sums = col_sums.replace(0, 1)
            meta_data[sample_cols] = meta_data[sample_cols].div(col_sums, axis=1)*100
            print(f"代谢组数据归一化 (相对丰度转换)")

        if pqn_normalization:
            # 1. 确保数据为数值型
            meta_data[sample_cols] = meta_data[sample_cols].apply(pd.to_numeric, errors='coerce')
            # 2. 计算参考样本（所有样本的中位数）
            reference_sample = meta_data[sample_cols].median(axis=1)
            # 避免除以 0
            min_pos = reference_sample[reference_sample > 0].min() if (reference_sample > 0).any() else 1e-6
            reference_sample = reference_sample.replace(0, min_pos)
            
            # 3. 计算每个样本相对于参考样本的比率
            quotients = meta_data[sample_cols].div(reference_sample, axis=0)
            
            # 4. 计算每个样本比率的中位数（稀释因子）
            dilution_factors = quotients.median(axis=0)
            dilution_factors = dilution_factors.replace(0, 1) # 防止除以 0
            
            # 5. 用稀释因子对原始数据进行归一化
            meta_data[sample_cols] = meta_data[sample_cols].div(dilution_factors, axis=1)
            print(f"代谢组数据归一化 (PQN)")
        
        if transform:
            sample_cols = meta_data.columns[1:]
            meta_data[sample_cols] = meta_data[sample_cols].apply(pd.to_numeric, errors='coerce')
            
            meta_data[sample_cols] = meta_data[sample_cols].fillna(meta_data[sample_cols].median())
            
            if transform_method == 'log':
                meta_data[sample_cols] = np.log1p(meta_data[sample_cols].clip(lower=0))
                print(f"代谢组数据变换 (log1p)")
                
            elif transform_method == 'log2':
                positive_vals = meta_data[sample_cols].to_numpy(dtype=float)
                positive_vals = positive_vals[positive_vals > 0]
                if positive_vals.size > 0:
                    eps = min(1e-6, float(np.nanmin(positive_vals)) / 2.0)
                    eps = max(eps, 1e-12)
                else:
                    eps = 1e-6
                meta_data[sample_cols] = np.log2(meta_data[sample_cols].clip(lower=0) + eps)
                print(f"代谢组数据变换 (log2, eps={eps:.2e})")
                
            elif transform_method == 'sqrt':
                meta_data[sample_cols] = np.sqrt(meta_data[sample_cols].clip(lower=0))
                print(f"代谢组数据变换 (sqrt)")
        
        if scale:
            sample_cols = meta_data.columns[1:]
            
            scaler = StandardScaler()
            meta_data[sample_cols] = scaler.fit_transform(meta_data[sample_cols])
            print(f"代谢组数据标准化 (StandardScaler)")
        
        if self.metabolomics_data.shape[0] >= 1 and self.metabolomics_data.shape[1] >= 1:
            first_cell = str(self.metabolomics_data.iloc[0, 0]).strip()
            if first_cell.lower() in ('lable', 'label'):
                self.metabolomics_data = pd.concat([self.metabolomics_data.iloc[[0]], meta_data], ignore_index=True)
            else:
                self.metabolomics_data = meta_data
        else:
            self.metabolomics_data = meta_data
        
        print(f"代谢组数据预处理完成: {original_shape} -> {self.metabolomics_data.shape}\n")

    def preprocess_kegg_data(self, remove_low_expression: bool = True,
                             min_abund_threshold: float = 0.01,
                             min_prevalence_threshold: float = 0.1,
                             remove_outliers: bool = False,
                             outlier_method: str = 'iqr',
                             outlier_threshold: float = 3.0,
                             relative_abund: bool = True,
                             transform: bool = False,
                             transform_method: str = 'log') -> None:
        """预处理KEGG KO丰度数据"""
        if self.kegg_data is None:
            return

        print("\n预处理KEGG KO丰度数据...")
        original_shape = self.kegg_data.shape

        ko_col = self.kegg_data.columns[0]
        self.kegg_data = self.kegg_data.drop_duplicates(subset=[ko_col], keep='first')
        sample_cols = self.kegg_data.columns[1:]

        if remove_low_expression:
            self.kegg_data[sample_cols] = self.kegg_data[sample_cols].apply(pd.to_numeric, errors='coerce')

            expression_counts = (self.kegg_data[sample_cols] > min_abund_threshold).sum(axis=1)
            self.kegg_data = self.kegg_data[expression_counts >= min_prevalence_threshold * len(sample_cols)]

            print(f"去除低表达KO: 保留 {self.kegg_data.shape[0]} 个KO")

        if remove_outliers:
            sample_cols = self.kegg_data.columns[1:]
            self.kegg_data[sample_cols] = self.kegg_data[sample_cols].apply(pd.to_numeric, errors='coerce')
            outlier_count = 0

            if outlier_method == 'iqr':
                for col in sample_cols:
                    Q1 = self.kegg_data[col].quantile(0.25)
                    Q3 = self.kegg_data[col].quantile(0.75)
                    IQR = Q3 - Q1
                    lower_bound = Q1 - outlier_threshold * IQR
                    upper_bound = Q3 + outlier_threshold * IQR

                    outliers = (self.kegg_data[col] < lower_bound) | (self.kegg_data[col] > upper_bound)
                    self.kegg_data.loc[outliers, col] = self.kegg_data.loc[outliers, col].clip(lower_bound, upper_bound)
                    outlier_count += outliers.sum()

            elif outlier_method == 'zscore':
                for col in sample_cols:
                    z_scores = np.abs(stats.zscore(self.kegg_data[col], nan_policy='omit'))
                    outliers = z_scores > outlier_threshold

                    median_val = self.kegg_data[col].median()
                    self.kegg_data.loc[outliers, col] = median_val
                    outlier_count += outliers.sum()

            print(f"去除异常值 ({outlier_method}方法): 处理了 {outlier_count} 个异常值")

        if relative_abund:
            self.kegg_data[sample_cols] = self.kegg_data[sample_cols].apply(pd.to_numeric, errors='coerce')
            col_sums = self.kegg_data[sample_cols].sum(axis=0)
            col_sums = col_sums.replace(0, 1)
            self.kegg_data[sample_cols] = self.kegg_data[sample_cols].div(col_sums, axis=1) * 100
            print("KEGG数据归一化 (相对丰度转换)")

        if transform:
            sample_cols = self.kegg_data.columns[1:]
            self.kegg_data[sample_cols] = self.kegg_data[sample_cols].apply(pd.to_numeric, errors='coerce')

            if transform_method == 'clr':
                pseudo_count = 1e-9
                clipped = self.kegg_data[sample_cols].clip(lower=0)
                logged = np.log(clipped + pseudo_count)
                self.kegg_data[sample_cols] = logged.sub(logged.mean(axis=0), axis=1)
                print("KEGG数据变换 (CLR)")

            elif transform_method == 'log':
                self.kegg_data[sample_cols] = np.log1p(self.kegg_data[sample_cols].clip(lower=0))
                print("KEGG数据变换 (log1p)")

        print(f"KEGG数据预处理完成: {original_shape} -> {self.kegg_data.shape}\n")
        
    def _load_sample_ids(self) -> None:
        """从Excel文件加载样本ID"""
        try:
            # 读取两个sheet
            df_patients = pd.read_excel(self.sample_file, sheet_name=0)  # 患者
            df_controls = pd.read_excel(self.sample_file, sheet_name=1)  # 健康对照组
            
            # 提取SAMPLE_ID列，并去除空格
            patient_ids = df_patients['SAMPLE_ID'].dropna().astype(str).str.strip().tolist()
            control_ids = df_controls['SAMPLE_ID'].dropna().astype(str).str.strip().tolist()
            
            # 合并所有样本ID
            self.sample_ids = patient_ids + control_ids
            
            # 存储DataFrame供后续使用
            self.df_patients = df_patients
            self.df_controls = df_controls
            
            print(f"加载了 {len(patient_ids)} 个患者样本和 {len(control_ids)} 个健康对照样本")
            
        except Exception as e:
            raise ValueError(f"加载样本ID失败: {e}")
    
    def _load_taxonomy_data(self) -> None:
        """加载宏基因组数据"""
        try:
            # 读取宏基因组数据
            self.taxonomy_data = pd.read_csv(self.taxonomy_file, sep='\t')
            
            # 去除列名中的空格
            self.taxonomy_data.columns = self.taxonomy_data.columns.str.strip()

            # 去除unclassified
            self.taxonomy_data = self.taxonomy_data[~self.taxonomy_data.iloc[:,0].str.contains('unclassified', na=False)]

            # self.taxonomy_data = self.taxonomy_data[~self.taxonomy_data.iloc[:,0].str.contains('CAG', na=False)]

            self.taxonomy_data = self.taxonomy_data[~self.taxonomy_data.iloc[:,0].str.contains('_sp', na=False)]

            self.taxonomy_data = self.taxonomy_data[~self.taxonomy_data.iloc[:,0].str.contains('\\[', na=False)]

            # 第一列是Species，其他列是样本
            print(f"宏基因组数据形状: {self.taxonomy_data.shape}")
            # print(f"宏基因组样本列数: {len(self.taxonomy_data.columns) - 1}")
            
        except Exception as e:
            raise ValueError(f"加载宏基因组数据失败: {e}")
    
    def _load_metabolomics_data(self) -> None:
        """加载代谢组数据"""
        try:
            # 读取代谢组数据
            self.metabolomics_data = pd.read_csv(self.metabolomics_file, sep=',', low_memory=False)
            
            # 去除列名中的空格
            self.metabolomics_data.columns = self.metabolomics_data.columns.str.strip()
            
            # 第一行是标签行（分组信息），第二行开始是实际数据
            # 第一列是代谢物标识
            print(f"代谢组数据形状: {self.metabolomics_data.shape}")
            
            # 检查第一行是否为标签行
            if self.metabolomics_data.iloc[0, 0] == 'Lable':
                print("检测到标签行，将在处理时跳过")
            
            # 尝试将数值列转换为浮点数
            # 跳过第一列（代谢物标识）
            # if len(self.metabolomics_data.columns) > 1:
            #     # 获取数值列（从第二列开始）
            #     numeric_cols = self.metabolomics_data.columns[1:]
            #     for col in numeric_cols:
            #         try:
            #             # 尝试转换为数值，错误时强制转换为NaN
            #             self.metabolomics_data[col] = pd.to_numeric(self.metabolomics_data[col], errors='coerce')
            #         except Exception as e:
            #             print(f"警告: 无法将列 {col} 转换为数值: {e}")
            
        except Exception as e:
            raise ValueError(f"加载代谢组数据失败: {e}")

    def _load_kegg_data(self) -> None:
        """加载KEGG KO丰度数据"""
        try:
            self.kegg_data = pd.read_csv(self.kegg_file, sep='\t')

            self.kegg_data.columns = self.kegg_data.columns.str.strip()

            print(f"KEGG数据形状: {self.kegg_data.shape}")

        except Exception as e:
            raise ValueError(f"加载KEGG数据失败: {e}")
    
    def _extract_smoking_labels(self) -> None:
        """提取吸烟标签"""
        try:
            # 对于患者样本，从sheet1获取吸烟信息
            if hasattr(self, 'df_patients'):
                # 查找包含'吸烟'的列
                smoking_cols = [col for col in self.df_patients.columns if '吸烟' in str(col)]
                
                if smoking_cols:
                    # 使用第一个包含'吸烟'的列
                    smoking_col = smoking_cols[0]
                    print(f"使用列 '{smoking_col}' 获取吸烟信息")
                    
                    # 提取患者样本的吸烟标签
                    for _, row in self.df_patients.iterrows():
                        sample_id = str(row['SAMPLE_ID']).strip()  # 去除空格
                        smoking_value = row[smoking_col]
                        
                        # 处理缺失值
                        if pd.isna(smoking_value):
                            self.smoking_labels[sample_id] = 0  # 默认不吸烟
                        else:
                            # 转换为整数
                            try:
                                self.smoking_labels[sample_id] = int(float(smoking_value))
                            except:
                                self.smoking_labels[sample_id] = 0
                else:
                    print("警告: 未找到吸烟信息列，所有患者样本标记为0")
                    for sample_id in self.sample_ids:
                        if sample_id in [str(id).strip() for id in self.df_patients['SAMPLE_ID'].dropna().tolist()]:
                            self.smoking_labels[sample_id] = 0
            
            # 对于健康对照组，默认不吸烟（0）
            if hasattr(self, 'df_controls'):
                for sample_id in self.df_controls['SAMPLE_ID'].dropna().astype(str).str.strip():
                    self.smoking_labels[sample_id] = 0
                    
            print(f"提取了 {len(self.smoking_labels)} 个样本的吸烟标签")
            
        except Exception as e:
            print(f"提取吸烟标签时出错: {e}")
            # 如果出错，将所有样本标记为0
            for sample_id in self.sample_ids:
                self.smoking_labels[sample_id] = 0
    
    def _extract_gender_labels(self) -> None:
        """提取性别标签"""
        try:
            # 对于患者样本，从sheet1获取性别信息
            if hasattr(self, 'df_patients'):
                # 查找包含'性别'的列
                gender_cols = [col for col in self.df_patients.columns if '性别' in str(col)]
                
                if gender_cols:
                    # 使用第一个包含'性别'的列
                    gender_col = gender_cols[0]
                    print(f"使用列 '{gender_col}' 获取性别信息")
                    
                    # 提取患者样本的性别标签
                    for _, row in self.df_patients.iterrows():
                        sample_id = str(row['SAMPLE_ID']).strip()  # 去除空格
                        gender_value = row[gender_col]
                        
                        # 处理缺失值
                        if pd.isna(gender_value):
                            self.gender_labels[sample_id] = 0  # 默认女性
                        else:
                            # 转换为整数
                            try:
                                self.gender_labels[sample_id] = int(float(gender_value))
                            except:
                                self.gender_labels[sample_id] = 0
                else:
                    print("警告: 未找到性别信息列，所有患者样本标记为0")
                    for sample_id in self.sample_ids:
                        if sample_id in [str(id).strip() for id in self.df_patients['SAMPLE_ID'].dropna().tolist()]:
                            self.gender_labels[sample_id] = 0
            
            # 对于健康对照组，从sheet2获取性别信息
            if hasattr(self, 'df_controls'):
                # 查找包含'性别'的列
                gender_cols = [col for col in self.df_controls.columns if '性别' in str(col)]
                
                if gender_cols:
                    # 使用第一个包含'性别'的列
                    gender_col = gender_cols[0]
                    print(f"使用列 '{gender_col}' 获取健康对照组性别信息")
                    
                    # 提取健康对照样本的性别标签
                    for _, row in self.df_controls.iterrows():
                        sample_id = str(row['SAMPLE_ID']).strip()  # 去除空格
                        gender_value = row[gender_col]
                        
                        # 处理缺失值
                        if pd.isna(gender_value):
                            self.gender_labels[sample_id] = 0  # 默认女性
                        else:
                            # 转换为整数：0=female, 1=male
                            try:
                                # 如果性别是字符串，转换为数字
                                if isinstance(gender_value, str):
                                    if '男' in gender_value or 'male' in gender_value.lower() or gender_value == '1':
                                        self.gender_labels[sample_id] = 1
                                    else:
                                        self.gender_labels[sample_id] = 0
                                else:
                                    self.gender_labels[sample_id] = int(float(gender_value))
                            except:
                                self.gender_labels[sample_id] = 0
                else:
                    print("警告: 未找到健康对照组性别信息列，所有健康对照样本标记为0")
                    for sample_id in self.df_controls['SAMPLE_ID'].dropna().astype(str).str.strip():
                        self.gender_labels[sample_id] = 0
                    
            print(f"提取了 {len(self.gender_labels)} 个样本的性别标签")
            
        except Exception as e:
            print(f"提取性别标签时出错: {e}")
            # 如果出错，将所有样本标记为0
            for sample_id in self.sample_ids:
                self.gender_labels[sample_id] = 0
    
    def _extract_tnm_labels(self) -> None:
        """提取TNM分期标签"""
        try:
            # 对于患者样本，从sheet1获取TNM分期信息
            if hasattr(self, 'df_patients'):
                # 查找包含'TNM'或'分期'的列
                tnm_cols = [col for col in self.df_patients.columns if 'TNM' in str(col) or '分期' in str(col)]
                
                if tnm_cols:
                    # 使用第一个包含'TNM'或'分期'的列
                    tnm_col = tnm_cols[0]
                    print(f"使用列 '{tnm_col}' 获取TNM分期信息")
                    
                    # 提取患者样本的TNM分期标签
                    for _, row in self.df_patients.iterrows():
                        sample_id = str(row['SAMPLE_ID']).strip()  # 去除空格
                        tnm_value = row[tnm_col]
                        
                        # 处理缺失值
                        if pd.isna(tnm_value):
                            self.tnm_labels[sample_id] = 0
                        else:
                            try:
                                self.tnm_labels[sample_id] = int(tnm_value) 
                            except:
                                self.tnm_labels[sample_id] = 0
                            
                else:
                    print("警告: 未找到TNM分期信息列，所有患者样本标记为Unknown")
                    for sample_id in self.sample_ids:
                        if sample_id in [str(id).strip() for id in self.df_patients['SAMPLE_ID'].dropna().tolist()]:
                            self.tnm_labels[sample_id] = 0
            
            # 对于健康对照组，标记为'Healthy'
            if hasattr(self, 'df_controls'):
                for sample_id in self.df_controls['SAMPLE_ID'].dropna().astype(str).str.strip():
                    self.tnm_labels[sample_id] = 0
                    
            print(f"提取了 {len(self.tnm_labels)} 个样本的TNM分期标签")
            
        except Exception as e:
            print(f"提取TNM分期标签时出错: {e}")
            # 如果出错，将所有样本标记为Unknown
            for sample_id in self.sample_ids:
                self.tnm_labels[sample_id] = 'Unknown'

    def _extract_age_labels(self) -> None:
        """提取年龄标签"""
        try:
            # 患者样本年龄列
            if hasattr(self, 'df_patients'):
                age_cols = [col for col in self.df_patients.columns if '龄' in str(col) or 'age' in str(col).lower()]
                if age_cols:
                    age_col = age_cols[0]
                    print(f"使用列 '{age_col}' 获取年龄信息")
                    for _, row in self.df_patients.iterrows():
                        sample_id = str(row['SAMPLE_ID']).strip()
                        val = row[age_col]
                        if pd.isna(val):
                            self.age_labels[sample_id] = 0
                        else:
                            try:
                                self.age_labels[sample_id] = int(val)
                            except Exception:
                                m = re.search(r"\d+", str(val))
                                if m:
                                    self.age_labels[sample_id] = int(m.group())
                                else:
                                    self.age_labels[sample_id] = 0
                else:
                    print("警告: 未找到年龄列，所有患者样本标记为0")
                    for sample_id in self.df_patients['SAMPLE_ID'].dropna().astype(str).str.strip():
                        self.age_labels[sample_id] = 0

            # 健康对照组年龄
            if hasattr(self, 'df_controls'):
                age_cols = [col for col in self.df_controls.columns if '龄' in str(col) or 'age' in str(col).lower()]
                if age_cols:
                    age_col = age_cols[0]
                    print(f"使用列 '{age_col}' 获取健康对照组年龄信息")
                    for _, row in self.df_controls.iterrows():
                        sample_id = str(row['SAMPLE_ID']).strip()
                        val = row[age_col]
                        if pd.isna(val):
                            self.age_labels[sample_id] = 0
                        else:
                            try:
                                self.age_labels[sample_id] = int(val)
                            except Exception:
                                m = re.search(r"\d+", str(val))
                                if m:
                                    self.age_labels[sample_id] = int(m.group())
                                else:
                                    self.age_labels[sample_id] = 0
                else:
                    print("警告: 未找到健康对照组年龄列，所有健康对照样本标记为0")
                    for sample_id in self.df_controls['SAMPLE_ID'].dropna().astype(str).str.strip():
                        self.age_labels[sample_id] = 0

            print(f"提取了 {len(self.age_labels)} 个样本的年龄标签")
        except Exception as e:
            print(f"提取年龄标签时出错: {e}")
            for sample_id in self.sample_ids:
                self.age_labels[sample_id] = 0

    def _extract_differentiation_labels(self) -> None:
        """提取分化标签: 0=中-高分化, 1=低分化"""
        try:
            if hasattr(self, 'df_patients'):
                diff_cols = [col for col in self.df_patients.columns if '分化' in str(col)]
                
                if diff_cols:
                    diff_col = diff_cols[0]
                    print(f"使用列 '{diff_col}' 获取分化信息")
                    
                    for _, row in self.df_patients.iterrows():
                        sample_id = str(row['SAMPLE_ID']).strip()
                        val = row[diff_col]
                        
                        if pd.isna(val):
                            self.differentiation_labels[sample_id] = 0
                        else:
                            val_str = str(val).strip()
                            if val_str in ['1', '低分化', '低', 'poor', 'Poor']:
                                self.differentiation_labels[sample_id] = 1
                            else:
                                self.differentiation_labels[sample_id] = 0
                else:
                    print("警告: 未找到分化信息列，所有患者样本标记为0")
                    for sample_id in self.df_patients['SAMPLE_ID'].dropna().astype(str).str.strip():
                        self.differentiation_labels[sample_id] = 0
            
            if hasattr(self, 'df_controls'):
                for sample_id in self.df_controls['SAMPLE_ID'].dropna().astype(str).str.strip():
                    self.differentiation_labels[sample_id] = 0
                    
            print(f"提取了 {len(self.differentiation_labels)} 个样本的分化标签")
            
        except Exception as e:
            print(f"提取分化标签时出错: {e}")
            for sample_id in self.sample_ids:
                self.differentiation_labels[sample_id] = 0

    def merge_to_dataframe(self) -> pd.DataFrame:
        """将 taxonomy、metabolomics 和 smoking_labels 合并为以样本为行的 DataFrame

        返回:
            pd.DataFrame: 行是样本ID，列为 `tax_<Species>`、`met_<Metabolite>`，以及 `smoking_label`
        """
        # 检查数据加载
        if self.taxonomy_data is None or self.metabolomics_data is None:
            raise ValueError("taxonomy_data 或 metabolomics_data 未加载")

        # 处理宏基因组数据：Species列为索引，转置为样本为行
        tax_df = self.taxonomy_data.copy()
        species_col = tax_df.columns[0]
        tax_df[species_col] = tax_df[species_col].astype(str).str.strip()
        tax_t = tax_df.set_index(species_col).T
        tax_t.index = tax_t.index.astype(str).str.strip()
        tax_t.index.name = 'SAMPLE_ID'
        tax_t.columns = [f"tax_{c}" for c in tax_t.columns]

        # 处理代谢组数据：如果第一行为标签则跳过，再将代谢物列设为索引并转置
        meta = self.metabolomics_data.copy()
        meta_data = meta
        if meta.shape[0] >= 1 and meta.shape[1] >= 1:
            first_cell = str(meta.iloc[0, 0]).strip()
            if first_cell.lower() in ('lable', 'label'):
                meta_data = meta.iloc[1:].copy()

        metab_col = meta_data.columns[0]
        meta_data[metab_col] = meta_data[metab_col].astype(str).str.strip()
        meta_t = meta_data.set_index(metab_col).T
        meta_t.index = meta_t.index.astype(str).str.strip()
        meta_t.index.name = 'SAMPLE_ID'
        meta_t.columns = [f"met_{c}" for c in meta_t.columns]

        # 按样本ID合并（外连接），保持标签尽可能完整
        merge_parts = [tax_t, meta_t]
        if self.kegg_data is not None:
            # 处理KEGG数据：KO列为索引并转置
            kegg_df = self.kegg_data.copy()
            ko_col = kegg_df.columns[0]
            kegg_df[ko_col] = kegg_df[ko_col].astype(str).str.strip()
            kegg_t = kegg_df.set_index(ko_col).T
            kegg_t.index = kegg_t.index.astype(str).str.strip()
            kegg_t.index.name = 'SAMPLE_ID'
            kegg_t.columns = [f"kegg_{c}" for c in kegg_t.columns]
            merge_parts.append(kegg_t)

        merged = pd.concat(merge_parts, axis=1, join='outer')

        # 如果 sample_ids 存在且在 merged 中，按 sample_ids 顺序排列
        ordered_idx = [str(s).strip() for s in self.sample_ids if str(s).strip() in merged.index]
        if ordered_idx:
            merged = merged.reindex(ordered_idx)

        # 添加 CRC/CTRL 标签：1 表示 CRC（来自患者表），0 表示 CTRL（来自对照表）
        patient_set: Set[str] = set()
        control_set: Set[str] = set()
        if hasattr(self, 'df_patients'):
            patient_set = set(self.df_patients['SAMPLE_ID'].dropna().astype(str).str.strip().tolist())
        if hasattr(self, 'df_controls'):
            control_set = set(self.df_controls['SAMPLE_ID'].dropna().astype(str).str.strip().tolist())

        def _crc_label_for(sid: str) -> int:
            sid = str(sid).strip()
            if sid in patient_set:
                return 1
            if sid in control_set:
                return 0
            return 0

        crc_labels = [_crc_label_for(sid) for sid in merged.index]
        merged['crc_label'] = pd.Series(crc_labels, index=merged.index).astype(int)

        # 添加吸烟标签列，缺失的标记为0
        smoking_series = pd.Series({str(k).strip(): int(v) for k, v in self.smoking_labels.items()})
        smoking_series = smoking_series.reindex(merged.index).fillna(0).astype(int)
        merged['smoking_label'] = smoking_series.values

        # 添加性别标签列，缺失的标记为0
        gender_series = pd.Series({str(k).strip(): int(v) for k, v in self.gender_labels.items()})
        gender_series = gender_series.reindex(merged.index).fillna(0).astype(int)
        merged['gender_label'] = gender_series.values

        # 添加TNM分期标签列，缺失的标记为'0'
        tnm_series = pd.Series({str(k).strip(): v for k, v in self.tnm_labels.items()})
        tnm_series = tnm_series.reindex(merged.index).fillna(0).astype(int)
        merged['tnm_stage'] = tnm_series.values

        # 添加年龄列，缺失时填充0
        age_series = pd.Series({str(k).strip(): float(v) for k, v in self.age_labels.items()})
        age_series = age_series.reindex(merged.index).fillna(0.0).astype(float)
        merged['age'] = age_series.values

        # 添加分化标签列：0=中-高分化，1=低分化
        diff_series = pd.Series({str(k).strip(): int(v) for k, v in self.differentiation_labels.items()})
        diff_series = diff_series.reindex(merged.index).fillna(0).astype(int)
        merged['differentiation'] = diff_series.values

        # 将缺失值替换为0（数值数据）
        merged = merged.fillna(0)

        print(f'合并后数据形状: {merged.shape}')

        return merged
    
    def intersection_to_dataframe(self) -> pd.DataFrame:
        # 检查数据加载
        if self.taxonomy_data is None or self.metabolomics_data is None:
            raise ValueError("taxonomy_data 或 metabolomics_data 未加载")

        # 处理宏基因组数据：Species列为索引，转置为样本为行
        tax_df = self.taxonomy_data.copy()
        species_col = tax_df.columns[0]
        tax_df[species_col] = tax_df[species_col].astype(str).str.strip()
        tax_t = tax_df.set_index(species_col).T
        tax_t.index = tax_t.index.astype(str).str.strip()
        tax_t.index.name = 'SAMPLE_ID'
        tax_t.columns = [f"tax_{c}" for c in tax_t.columns]

        # 处理代谢组数据：如果第一行为标签则跳过，再将代谢物列设为索引并转置
        meta = self.metabolomics_data.copy()
        meta_data = meta
        if meta.shape[0] >= 1 and meta.shape[1] >= 1:
            first_cell = str(meta.iloc[0, 0]).strip()
            if first_cell.lower() in ('lable', 'label'):
                meta_data = meta.iloc[1:].copy()

        metab_col = meta_data.columns[0]
        meta_data[metab_col] = meta_data[metab_col].astype(str).str.strip()
        meta_t = meta_data.set_index(metab_col).T
        meta_t.index = meta_t.index.astype(str).str.strip()
        meta_t.index.name = 'SAMPLE_ID'
        meta_t.columns = [f"met_{c}" for c in meta_t.columns]

        # 按样本ID合并（内连接），只保留两者都有的样本
        merged = pd.merge(tax_t, meta_t, left_index=True, right_index=True, how='inner')
        if self.kegg_data is not None:
            # 处理KEGG数据：KO列为索引并转置
            kegg_df = self.kegg_data.copy()
            ko_col = kegg_df.columns[0]
            kegg_df[ko_col] = kegg_df[ko_col].astype(str).str.strip()
            kegg_t = kegg_df.set_index(ko_col).T
            kegg_t.index = kegg_t.index.astype(str).str.strip()
            kegg_t.index.name = 'SAMPLE_ID'
            kegg_t.columns = [f"kegg_{c}" for c in kegg_t.columns]
            merged = pd.merge(merged, kegg_t, left_index=True, right_index=True, how='inner')

        # 如果 sample_ids 存在且在 merged 中，按 sample_ids 顺序排列
        ordered_idx = [str(s).strip() for s in self.sample_ids if str(s).strip() in merged.index]
        if ordered_idx:
            merged = merged.reindex(ordered_idx)

        # 添加 CRC/CTRL 标签：1 表示 CRC（来自患者表），0 表示 CTRL（来自对照表）
        patient_set: Set[str] = set()
        control_set: Set[str] = set()
        if hasattr(self, 'df_patients'):
            patient_set = set(self.df_patients['SAMPLE_ID'].dropna().astype(str).str.strip().tolist())
        if hasattr(self, 'df_controls'):
            control_set = set(self.df_controls['SAMPLE_ID'].dropna().astype(str).str.strip().tolist())

        def _crc_label_for(sid: str) -> int:
            sid = str(sid).strip()
            if sid in patient_set:
                return 1
            if sid in control_set:
                return 0
            return 0

        crc_labels = [_crc_label_for(sid) for sid in merged.index]
        merged['crc_label'] = pd.Series(crc_labels, index=merged.index).astype(int)

        # 添加吸烟标签列，缺失的标记为0
        smoking_series = pd.Series({str(k).strip(): int(v) for k, v in self.smoking_labels.items()})
        smoking_series = smoking_series.reindex(merged.index).fillna(0).astype(int)
        merged['smoking_label'] = smoking_series.values

        # 添加性别标签列，缺失的标记为0
        gender_series = pd.Series({str(k).strip(): int(v) for k, v in self.gender_labels.items()})
        gender_series = gender_series.reindex(merged.index).fillna(0).astype(int)
        merged['gender_label'] = gender_series.values

        # 添加TNM分期标签列，缺失的标记为'0'
        tnm_series = pd.Series({str(k).strip(): v for k, v in self.tnm_labels.items()})
        tnm_series = tnm_series.reindex(merged.index).fillna(0).astype(int)
        merged['tnm_stage'] = tnm_series.values

        # 添加年龄列，缺失时填充0
        age_series = pd.Series({str(k).strip(): float(v) for k, v in self.age_labels.items()})
        age_series = age_series.reindex(merged.index).fillna(0).astype(int)
        merged['age'] = age_series.values

        # 添加分化标签列：0=中-高分化，1=低分化
        diff_series = pd.Series({str(k).strip(): int(v) for k, v in self.differentiation_labels.items()})
        diff_series = diff_series.reindex(merged.index).fillna(0).astype(int)
        merged['differentiation'] = diff_series.values

        # 将缺失值替换为0（数值数据）
        merged = merged.fillna(0)

        print(f"合并后的数据形状: {merged.shape}")

        return merged

if __name__ == '__main__':
    # 加载数据
    ds = BioSmokeDataset(load_kegg=True)
    ds.preprocess_taxonomy_data()
    ds.preprocess_metabolomics_data()
    if ds.kegg_data is not None:
        ds.preprocess_kegg_data()
    merged = ds.merge_to_dataframe()
    merged = ds.intersection_to_dataframe()