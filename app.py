import streamlit as st
import pandas as pd
import joblib

# 加载 exon 和 domain_order 数据表
exon_df = pd.read_excel('小程序自行计算的特征.xlsx', sheet_name='第一张表')
domain_order_df = pd.read_excel('小程序自行计算的特征.xlsx', sheet_name='第二张表')

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
        return row['exon'].values[0]  # 修改为合适的列名以匹配functional_area
    return None

# 定义获取 domain_order 的函数
def get_domain_order(mutation_position_start):
    row = domain_order_df[(domain_order_df['Start_Position'] <= mutation_position_start) & (domain_order_df['End_Position'] >= mutation_position_start)]
    if not row.empty:
        return row['Domain_order'].values[0]
    return None

# 加载模型
model = joblib.load('voting_clf.pkl')

# Streamlit 应用程序
st.title("DMD Mutation Prediction App")

# 输入字段
mutation_position_start = st.number_input("Enter mutation position start", min_value=0, step=1)
mutation_position_stop = st.number_input("Enter mutation position stop", min_value=0, step=1)
mutation_type = st.selectbox("Mutation Type", options=[1, 2, 3, 4, 5])

# 计算特征
if st.button("Calculate Features"):
    exon = get_exon(mutation_position_start)
    functional_area = get_functional_area(exon) if exon is not None else None
    domain_order = get_domain_order(mutation_position_start)
    
    # 显示计算结果
    st.write("Exon:", exon)
    st.write("Functional Area:", functional_area)
    st.write("Domain Order:", domain_order)
    
    # 确保所有特征已计算
    if exon is not None and functional_area is not None and domain_order is not None:
        # 准备模型输入
        input_data = [[mutation_position_start, mutation_position_stop, mutation_type]]
        prediction = model.predict(input_data)
        
        # 显示预测结果
        st.write("Prediction:", prediction[0])
    else:
        st.write("Error: One or more features could not be calculated. Please check the input values.")
