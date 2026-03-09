import pandas as pd
import numpy as np

# 读取TSV文件（文件扩展名是xls但实际是TSV格式）
df = pd.read_csv('dataset/rawdata/combine.intensity.xls', sep='\t')

# 筛选出有注释的代谢物（IsAnnotated=1）
annotated_df = df[df['IsAnnotated'] == 1].copy()

# 提取注释信息，优先使用MS2Metabolite，其次是MS1hmdbName，最后是MS1keggName
def get_annotation(row):
    if pd.notna(row['MS2Metabolite']) and row['MS2Metabolite'] != '':
        return row['MS2Metabolite']
    elif pd.notna(row['MS1hmdbName']) and row['MS1hmdbName'] != '':
        return row['MS1hmdbName']
    elif pd.notna(row['MS1keggName']) and row['MS1keggName'] != '':
        return row['MS1keggName']
    else:
        return None

annotated_df['Annotation'] = annotated_df.apply(get_annotation, axis=1)

# 移除没有注释的代谢物
annotated_df = annotated_df[annotated_df['Annotation'].notna()].copy()

# 提取样本列（所有以oCRC_、oCTRL_、yCRC_、yCTRL_开头的列）
sample_columns = [col for col in df.columns if any(prefix in col for prefix in ['oCRC_', 'oCTRL_', 'yCRC_', 'yCTRL_'])]

# 构建结果数据框
result_df = annotated_df[['Annotation'] + sample_columns].copy()

# 去重，保留每个注释的第一行
result_df = result_df.drop_duplicates(subset=['Annotation'], keep='first')

# 重置索引
result_df = result_df.reset_index(drop=True)

# 保存为CSV文件
result_df.to_csv('metabolome_annotated_samples.csv', index=False, encoding='utf-8-sig')

print(f"处理完成！")
print(f"原始数据行数: {len(df)}")
print(f"有注释的代谢物行数: {len(annotated_df)}")
print(f"去重后唯一注释数: {len(result_df)}")
print(f"保存文件: metabolome_annotated_samples.csv")
