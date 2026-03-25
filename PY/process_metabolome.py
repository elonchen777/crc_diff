import pandas as pd
import numpy as np

# 读取TSV文件（文件扩展名是xls但实际是TSV格式）
df = pd.read_csv('dataset/rawdata/combine.intensity.xls', sep='\t')

# 直接处理全部代谢物，未注释条目用 ID 兜底命名
processed_df = df.copy()

# 提取注释信息，优先使用MS2Metabolite，其次是MS1hmdbName，最后是MS1keggName
def _clean_text(value):
    if pd.isna(value):
        return None
    text = str(value).strip()
    if text == '' or text == '-':
        return None
    return text


def _clean_ms1_text(value):
    text = _clean_text(value)
    if text is None:
        return None
    return f'ms1-{text}'


def get_annotation(row):
    value = _clean_text(row['MS2Metabolite'])
    if value is not None:
        return value

    for col in ['MS1hmdbName', 'MS1keggName']:
        value = _clean_ms1_text(row[col])
        if value is not None:
            return value

    # 保留未注释代谢物，用 ID 作为后缀避免被去重合并
    raw_id = _clean_text(row.get('ID'))
    if raw_id is None:
        raw_id = str(row.name)
    return f'unkown-{raw_id}'

processed_df['Annotation'] = processed_df.apply(get_annotation, axis=1)

# 提取样本列（所有以oCRC_、oCTRL_、yCRC_、yCTRL_开头的列）
sample_columns = [
    col for col in df.columns
    if col.startswith(('oCRC_', 'oCTRL_', 'yCRC_', 'yCTRL_'))
]

# 构建结果数据框
result_df = processed_df[['Annotation'] + sample_columns].copy()

# 重复项按注释聚合丰度：可选 'max' 或 'mean'
aggregate_method = 'max'
if aggregate_method == 'mean':
    result_df = result_df.groupby('Annotation', as_index=False)[sample_columns].mean()
else:
    result_df = result_df.groupby('Annotation', as_index=False)[sample_columns].max()

# 重置索引
result_df = result_df.reset_index(drop=True)

# 只保留 MS2 标注的代谢物作为可选项
result_annotated_df = result_df[
    (~result_df['Annotation'].str.startswith('ms1-')) &
    (~result_df['Annotation'].str.startswith('unkown-'))
].copy()
result_df = result_df.reset_index(drop=True)

# 统计标注类型数量
ms2_count = len(result_annotated_df)
ms1_count = len(result_df[result_df['Annotation'].str.startswith('ms1-')])
unannotated_count = len(result_df[result_df['Annotation'].str.startswith('unkown-')])

# 保存为CSV文件
result_df.to_csv('dataset/metabolome_all_data.csv', index=False, encoding='utf-8-sig')
result_annotated_df.to_csv('dataset/metabolome_annotated_data.csv', index=False, encoding='utf-8-sig')

print(f"处理完成！")
print(f"原始数据行数: {len(df)}")
print(f"处理后代谢物行数: {len(processed_df)}")
print(f"去重后唯一注释数: {len(result_df)}")
print(f"保留的MS2标注数量: {ms2_count}")
print(f"MS1标注数量: {ms1_count}")
print(f"未标注数量: {unannotated_count}")
print(f"重复项聚合方式: {aggregate_method}")
print(f"保存文件: dataset/metabolome_all_data.csv")
print(f"保存文件: dataset/metabolome_annotated_data.csv")
