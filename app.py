import streamlit as st
import pandas as pd
import numpy as np
import joblib

# 加载模型
model = joblib.load("voting_clf.pkl")

# 加载特征表格
exon_df = pd.read_excel("特征表格.xlsx", sheet_name="第一张表", engine="openpyxl")
domain_df = pd.read_excel("特征表格.xlsx", sheet_name="第二张表", engine="openpyxl")

# 定义获取 exon 的函数
def get_exon(mutation_position_start):
    row = exon_df[(exon_df['start'] <= mutation_position_start) & (exon_df['end'] >= mutation_position_start)]
    if not row.empty:
        return row['exon'].values[0]
    return None

# 定义获取 functional_area 的函数
def get_functional_area(exon):
    row = exon_df[exon_df['exon'] == exon]
    if not row.empty:
        return row['exon'].values[0]  # 这里改为适当的列名，如果需要
    return None

# 定义获取 domain_order 的函数
def get_domain_order(mutation_position_start):
    row = domain_df[(domain_df['Start_Position'] <= mutation_position_start) & (domain_df['End_Position'] >= mutation_position_start)]
    if not row.empty:
        return row['Domain_order'].values[0]
    return None

# Streamlit 页面标题
st.title("DMD Mutation Prediction App")

# 用户输入
mutation_position_start = st.number_input("Enter mutation position start", min_value=0)
mutation_position_stop = st.number_input("Enter mutation position stop", min_value=0)
mutation_type = st.selectbox("Mutation Type", options=[1, 2, 3, 4, 5], format_func=lambda x: ["Nonsense", "Frameshift", "Missense", "Synonymous", "Non-frameshift small indel"][x-1])

# 自动计算特征
exon = get_exon(mutation_position_start)
functional_area = get_functional_area(exon) if exon is not None else -999
domain_order = get_domain_order(mutation_position_start) if mutation_position_start is not None else -999

# 组装输入数据
input_data = {
    'Mutation_position_start': mutation_position_start,
    'Mutation_position_stop': mutation_position_stop,
    'Mutation_types': mutation_type,
    'Functional_area': functional_area,
    'Domain_order': domain_order,
    'exon': exon if exon is not None else -999,
    'frame_of_exons': 1 if exon in [1, 2, 6, 7, 8, 11, 12, 17, 18, 19, 20, 21, 22, 43, 44, 45, 46, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 61, 62, 63, 65, 66, 67, 68, 69, 70, 75, 76, 78, 79] else 0,
    'skipping_of_in_frame_exons': 1 if exon in [9, 25, 27, 29, 31, 37, 38, 39, 41, 72, 74] else 0,
    'frame': 1 if mutation_type in [1, 2] else 0,
    'Amino_acid_properties_changed': -999,
    'SpliceAI_pred_DS_DL': -999,
    'CADD_PHRED': -999,
    'CADD_RAW': -999,
    'GERP++_NR': -999,
    'GERP++_RS': -999,
    'GERP++_RS_rankscore': -999,
    'BayesDel_addAF_score': -999,
    'BayesDel_noAF_rankscore': -999,
    'BayesDel_noAF_score': -999,
    'DANN_rankscore': -999,
    'DANN_score': -999,
    'PrimateAI_score': -999,
    'MetaLR_rankscore': -999
}

# 转换为 DataFrame
input_data_df = pd.DataFrame([input_data])

# 显示自动推算的特征值
st.write("Exon:", exon)
st.write("Functional Area:", functional_area)
st.write("Domain Order:", domain_order)

# 预测
if st.button("Predict"):
    try:
        prediction = model.predict(input_data_df)
        st.write("Prediction (1=DMD, 0=BMD):", prediction[0])
    except Exception as e:
        st.error(f"Error in prediction: {e}")
