import streamlit as st
import pandas as pd
import joblib

# 加载 exon 和 domain_order 数据表
exon_df = pd.read_excel('小程序自行计算的特征.xlsx', sheet_name='第一张')
domain_order_df = pd.read_excel('小程序自行计算的特征.xlsx', sheet_name='第二张')

# 定义获取 exon 的函数
def get_exon(mutation_position_start):
    # 查找 mutation_position_start 是否在 exon 区间内
    row = exon_df[(exon_df['Mutation_position_start'] <= mutation_position_start) &
                  (exon_df['Mutation_position_stop'] >= mutation_position_start)]
    if not row.empty:
        return row['exon'].values[0]
    return None

# 定义获取 functional_area 的函数
def get_functional_area(exon):
    # 根据 exon 获取功能区域信息
    row = exon_df[exon_df['exon'] == exon]
    if not row.empty:
        return row['functional_area'].values[0]
    return None

# 定义获取 domain_order 的函数
def get_domain_order(mutation_position_start):
    # 查找 mutation_position_start 是否在 domain 区间内
    row = domain_order_df[(domain_order_df['Domain_start'] <= mutation_position_start) &
                          (domain_order_df['Domain_stop'] >= mutation_position_start)]
    if not row.empty:
        return row['Domain_order'].values[0]
    return None

# 加载模型
model = joblib.load('model.pkl')

# Streamlit 应用程序
st.title("DMD Mutation Prediction App")

# 输入
mutation_position_start = st.number_input("Enter mutation position start", min_value=0, step=1)
mutation_position_stop = st.number_input("Enter mutation position stop", min_value=0, step=1)
mutation_type = st.selectbox("Mutation Type", options=[1, 2, 3, 4, 5])

# 计算特征
if st.button("Calculate Features"):
    exon = get_exon(mutation_position_start)
    functional_area = get_functional_area(exon) if exon is not None else None
    domain_order = get_domain_order(mutation_position_start)

    st.write(f"Exon: {exon}")
    st.write(f"Functional Area: {functional_area}")
    st.write(f"Domain Order: {domain_order}")

    # 创建输入数据框
    input_data = {
        'Mutation_position_start': mutation_position_start,
        'Mutation_position_stop': mutation_position_stop,
        'Mutation_types': mutation_type,
        'exon': exon,
        'functional_area': functional_area,
        'Domain_order': domain_order
    }
    input_df = pd.DataFrame([input_data])

    # 预测
    if st.button("Predict"):
        prediction = model.predict(input_df)
        st.write(f"The prediction result is: {'DMD' if prediction[0] == 1 else 'Non-DMD'}")
else:
    st.write("Please calculate features first before predicting.")
