import streamlit as st
import pandas as pd
import joblib

# 加载模型
model_path = "voting_clf.pkl"
model = joblib.load(model_path)

# 加载 exon 和 Domain_order 数据表
# 请确保 '小程序自行计算的特征.xlsx' 文件存在，并且包含以下工作表
exon_df = pd.read_excel('小程序自行计算的特征.xlsx', sheet_name='第一张表')  # 第一张表格应包含 exon 的位置信息
domain_order_df = pd.read_excel('小程序自行计算的特征.xlsx', sheet_name='第二张表')  # 第二张表格应包含 Domain_order 的位置信息

# 定义获取 exon 的函数
def get_exon(mutation_position_start):
    # 查询起始位置是否在 exon 区间内
    row = exon_df[(exon_df['start'] <= mutation_position_start) & (exon_df['end'] >= mutation_position_start)]
    if not row.empty:
        return row['exon'].values[0]
    return None

# 根据 exon 获取 Functional_area
def get_functional_area(exon):
    if 1 <= exon <= 8:
        return 1
    elif 9 <= exon <= 64:
        return 2
    elif 65 <= exon <= 70:
        return 3
    elif exon >= 71:
        return 4
    return None

# 根据 mutation_position_start 获取 Domain_order
def get_domain_order(mutation_position_start):
    row = domain_order_df[
        (domain_order_df['Start_Position'] <= mutation_position_start) &
        (domain_order_df['End_Position'] >= mutation_position_start)
    ]
    if not row.empty:
        return row['Domain_order'].values[0]
    return None

# Streamlit 用户输入
st.title("DMD/BMD Prediction Model")

mutation_position_start = st.number_input("Mutation Position Start", min_value=1, max_value=12000, step=1)
mutation_position_stop = st.number_input("Mutation Position Stop", min_value=1, max_value=12000, step=1)
mutation_type = st.selectbox("Mutation Type", options=[1, 2, 3, 4, 5], format_func=lambda x: ["Nonsense", "Frameshift", "Missense", "Synonymous", "Non-frameshift"][x-1])

# 计算特征
exon = get_exon(mutation_position_start)
functional_area = get_functional_area(exon) if exon is not None else None
domain_order = get_domain_order(mutation_position_start)

# 显示自动生成的特征
st.write(f"Exon: {exon}")
st.write(f"Functional Area: {functional_area}")
st.write(f"Domain Order: {domain_order}")

# 将生成的特征组合到输入数据中并预测
if exon is not None and functional_area is not None and domain_order is not None:
    input_data = {
        'mutation_position_start': mutation_position_start,
        'mutation_position_stop': mutation_position_stop,
        'mutation_type': mutation_type,
        'exon': exon,
        'Functional_area': functional_area,
        'Domain_order': domain_order,
        # 默认值或其他自动生成的特征
        'frame_of_exons': 1 if exon in [1, 2, 6, 7, 8, 11, 12, 17, 18, 19, 20, 21, 22, 43, 44, 45, 46, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 61, 62, 63, 65, 66, 67, 68, 69, 70, 75, 76, 78, 79] else 0,
        'skipping_of_in_frame_exons': 1 if exon in [9, 25, 27, 29, 31, 37, 38, 39, 41, 72, 74] else 0,
        'frame': 1 if mutation_type in [1, 2] else 0,
        # 添加计算机预测特征的默认值 -999
        'CADD_PHRED': -999, 'CADD_RAW': -999, 'GERP++_NR': -999, 'GERP++_RS': -999,
        'GERP++_RS_rankscore': -999, 'BayesDel_addAF_score': -999, 'BayesDel_noAF_rankscore': -999,
        'BayesDel_noAF_score': -999, 'DANN_rankscore': -999, 'DANN_score': -999,
        'PrimateAI_score': -999, 'MetaLR_rankscore': -999
    }
    input_df = pd.DataFrame([input_data])

    # 预测结果
    if st.button("Predict"):
        prediction = model.predict(input_df)
        st.write(f"The prediction result is: {'DMD' if prediction[0] == 1 else 'BMD'}")
else:
    st.write("Please ensure all features are calculated correctly.")
