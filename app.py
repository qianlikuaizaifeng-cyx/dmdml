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

# 固定的转录本 ID
TRANSCRIPT_ID = "NM_004006.3"

# 定义获取 exon 的函数
def get_exon(mutation_position_start):
    row = exon_df[(exon_df['start'] <= mutation_position_start) & (exon_df['end'] >= mutation_position_start)]
    if not row.empty:
        return row['exon'].values[0]
    return None

# 定义获取 domain_order 的函数
def get_domain_order(mutation_position_start):
    row = domain_df[(domain_df['Start_Position'] <= mutation_position_start) & (domain_df['End_Position'] >= mutation_position_start)]
    if not row.empty:
        return int(row['Domain_order'].values[0])
    return -999

# 获取 Ensembl API 中的氨基酸变化和位置信息
def fetch_variant_info(hgvs_notation):
    url = f"https://rest.ensembl.org/vep/human/hgvs/{TRANSCRIPT_ID}:{hgvs_notation}"
    headers = {"Content-Type": "application/json"}
    response = requests.get(url, headers=headers, proxies={"http": None, "https": None})
    
    if response.status_code == 200:
        data = response.json()
        variant_info = data[0]
        
        start = variant_info.get('start')
        end = variant_info.get('end')
        consequence_terms = variant_info.get('most_severe_consequence')
        
        # 提取氨基酸变化信息
        transcript = variant_info.get('transcript_consequences', [{}])[0]
        amino_acids = transcript.get('amino_acids')
        amino_acid_before, amino_acid_after = (amino_acids.split("/") if amino_acids else (None, None))
        
        return start, end, consequence_terms, amino_acid_before, amino_acid_after
    return None, None, None, None, None

# 检查氨基酸是否在同一组
def check_amino_acid_group(amino_acid_before, amino_acid_after):
    hydrophobic = {'A', 'V', 'I', 'L', 'M', 'F', 'Y', 'W'}
    polar = {'S', 'T', 'N', 'Q'}
    positive = {'K', 'R', 'H'}
    negative = {'D', 'E'}
    groups = [hydrophobic, polar, positive, negative]

    for group in groups:
        if amino_acid_before in group and amino_acid_after in group:
            return 1  # Same group
    return 0  # Different groups

# Streamlit 页面标题
st.title("DMD Mutation Prediction App")

# 用户输入 HGVS 表达式的变异部分
hgvs_notation = st.text_input("Enter HGVS notation (e.g., c.1399A>T)")

# 自动获取和计算特征
if hgvs_notation:
    start, end, consequence_terms, amino_acid_before, amino_acid_after = fetch_variant_info(hgvs_notation)
    
    if start and end:
        st.write("Mutation Position Start:", start)
        st.write("Mutation Position Stop:", end)
        st.write("Consequence Terms:", consequence_terms)
        st.write("Amino Acid Before:", amino_acid_before)
        st.write("Amino Acid After:", amino_acid_after)
        
        # 自动填充特征
        exon = get_exon(start)
        functional_area = exon if exon is not None else -999
        domain_order = get_domain_order(start)
        
        # 根据变异类型设置 Amino_acid_properties_changed
        amino_acid_properties_changed = -999  # 默认值
        if consequence_terms == "synonymous_variant":
            amino_acid_properties_changed = 1
        elif consequence_terms in ["nonsense", "frameshift_variant"]:
            amino_acid_properties_changed = -999
        elif consequence_terms == "missense_variant":
            if amino_acid_before and amino_acid_after:
                amino_acid_properties_changed = check_amino_acid_group(amino_acid_before, amino_acid_after)

        # 组装输入数据
        input_data = {
            'Functional_area': functional_area,
            'frame_of_exons': 1 if exon in [1, 2, 6, 7, 8, 11, 12, 17, 18, 19, 20, 21, 22, 43, 44, 45, 46, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 61, 62, 63, 65, 66, 67, 68, 69, 70, 75, 76, 78, 79] else 0,
            'Amino_acid_properties_changed': amino_acid_properties_changed,
            'exon': exon if exon is not None else -999,
            'Mutation_position_start': start,
            'Mutation_position_stop': end,
            'frame': 1 if consequence_terms in ["nonsense", "frameshift_variant"] else 0,
            'Mutation_types': 3 if consequence_terms == "missense_variant" else (1 if consequence_terms == "nonsense" else (2 if consequence_terms == "frameshift_variant" else 4)),
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

        # 转换为 DataFrame 并确保列顺序一致
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
                prediction = model.predict(input_data_df)
                prediction_proba = model.predict_proba(input_data_df)
                
                result = "DMD" if prediction[0] == 1 else "BMD"
                probability = prediction_proba[0][1] if prediction[0] == 1 else prediction_proba[0][0]
                
                st.write(f"Prediction: {result}")
                st.write(f"Prediction Probability: {probability:.2f}")
                
            except Exception as e:
                st.error(f"Error in prediction: {e}")
