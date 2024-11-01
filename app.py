import streamlit as st
import pandas as pd
import numpy as np
import joblib
import requests

# 加载模型
model = joblib.load("voting_clf.pkl")

# 加载特征表格
exon_df = pd.read_excel("小程序自行计算的特征.xlsx", sheet_name="第一张表", engine="openpyxl")
domain_df = pd.read_excel("小程序自行计算的特征.xlsx", sheet_name="第二张表", engine="openpyxl")

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
        return row['Functional_area'].values[0]  # 假设列名为 Functional_area
    return -999

# 定义获取 domain_order 的函数
def get_domain_order(mutation_position_start):
    row = domain_df[(domain_df['Start_Position'] <= mutation_position_start) & (domain_df['End_Position'] >= mutation_position_start)]
    if not row.empty:
        return int(row['Domain_order'].values[0])  # 确保返回整数
    return -999  # 缺失值时用 -999 表示，且是整数

# 检查氨基酸是否在同一组
def check_amino_acid_group(amino_acid_before, amino_acid_after):
    hydrophobic = {'A', 'V', 'I', 'L', 'M', 'F', 'Y', 'W'}
    polar = {'S', 'T', 'N', 'Q'}
    positive = {'K', 'R', 'H'}
    negative = {'D', 'E'}
    groups = [hydrophobic, polar, positive, negative]

    for group in groups:
        if amino_acid_before in group and amino_acid_after in group:
            return 1  # 同组
    return 0  # 不同组

# 获取 Ensembl API 中的变异位置（仅获取指定转录本 ENST00000357033 的位置）
def fetch_ensembl_data(variant_suffix):
    url = f"https://rest.ensembl.org/vep/human/hgvs/NM_004006.3:{variant_suffix}"
    headers = {"Content-Type": "application/json"}
    response = requests.get(url, headers=headers, proxies={"http": None, "https": None})
    
    if response.status_code == 200:
        data = response.json()
        for transcript in data[0].get('transcript_consequences', []):
            if transcript.get('transcript_id') == 'ENST00000357033':
                start = transcript.get('cdna_start', -999)
                end = transcript.get('cdna_end', -999)
                amino_acids = transcript.get('amino_acids', '/').split('/')
                amino_acid_before, amino_acid_after = amino_acids if len(amino_acids) == 2 else (None, None)
                return start, end, amino_acid_before, amino_acid_after
    return -999, -999, None, None

# Streamlit 页面标题
st.title("DMD Mutation Prediction App")

# 用户输入
variant_suffix = st.text_input("Enter HGVS suffix (e.g., c.1399del)")
mutation_type = st.selectbox("Mutation Type", options=[1, 2, 3, 4], format_func=lambda x: ["Nonsense", "Frameshift", "Missense", "Synonymous"][x-1])

# 从 Ensembl API 获取变异位置
if variant_suffix:
    mutation_position_start, mutation_position_stop, amino_acid_before, amino_acid_after = fetch_ensembl_data(variant_suffix)
    
    # 自动设置 Amino_acid_properties_changed
    amino_acid_properties_changed = -999
    if mutation_type == 4:  # Synonymous
        amino_acid_properties_changed = 1
    elif mutation_type == 1 or mutation_type == 2:  # Nonsense or Frameshift
        amino_acid_properties_changed = -999
    elif mutation_type == 3:  # Missense
        if amino_acid_before and amino_acid_after:
            amino_acid_properties_changed = check_amino_acid_group(amino_acid_before, amino_acid_after)
    
    # 获取其他特征
    exon = get_exon(mutation_position_start)
    functional_area = get_functional_area(exon) if exon is not None else -999
    domain_order = get_domain_order(mutation_position_start)

    # 显示获取的自动填充数据
    st.write("Mutation Position Start:", mutation_position_start)
    st.write("Mutation Position Stop:", mutation_position_stop)
    st.write("Amino Acid Before:", amino_acid_before)
    st.write("Amino Acid After:", amino_acid_after)
    st.write("Amino Acid Properties Changed:", amino_acid_properties_changed)
    st.write("Mutation Type:", mutation_type)

    # 组装输入数据，确保特征名和顺序与训练时一致
    input_data = {
        'Functional_area': functional_area,
        'frame_of_exons': 1 if exon in [1, 2, 6, 7, 8, 11, 12, 17, 18, 19, 20, 21, 22, 43, 44, 45, 46, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 61, 62, 63, 65, 66, 67, 68, 69, 70, 75, 76, 78, 79] else 0,
        'Amino_acid_properties_changed': amino_acid_properties_changed,
        'exon': exon if exon is not None else -999,
        'Mutation_position_start': mutation_position_start,
        'Mutation_position_stop': mutation_position_stop,
        'frame': 1 if mutation_type in [1, 2] else 0,
        'Mutation_types': mutation_type,
        'Domain_order': domain_order,
        'skipping_of_in_frame_exons': 1 if exon in [9, 25, 27, 29, 31, 37, 38, 39, 41, 72, 74] else 0,
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

    # 转换为 DataFrame，并确保列顺序一致
    input_data_df = pd.DataFrame([input_data], columns=[
        'Functional_area', 'frame_of_exons', 'Amino_acid_properties_changed', 'exon', 
        'Mutation_position_start', 'Mutation_position_stop', 'frame', 'Mutation_types', 
        'Domain_order', 'skipping_of_in_frame_exons', 'SpliceAI_pred_DS_DL', 
        'CADD_PHRED', 'CADD_RAW', 'GERP++_NR', 'GERP++_RS', 'GERP++_RS_rankscore', 
        'BayesDel_addAF_score', 'BayesDel_noAF_rankscore', 'BayesDel_noAF_score', 
        'DANN_rankscore', 'DANN_score', 'PrimateAI_score', 'MetaLR_rankscore'
    ])

    # 预测
    if st.button("Predict"):
        try:
            # 获取预测类别
            prediction = model.predict(input_data_df)
            
            # 获取预测概率
            prediction_proba = model.predict_proba(input_data_df)
            
            # 显示预测结果
            result = "DMD" if prediction[0] == 1 else "BMD"
            probability = prediction_proba[0][1] if prediction[0] == 1 else prediction_proba[0][0]
            
            st.write(f"Prediction: {result}")
            st.write(f"Prediction Probability: {probability:.2f}")
            
        except Exception as e:
            st.error(f"Error in prediction: {e}")
