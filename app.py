# 计算特征并更新特征名称以匹配模型的期望
exon = get_exon(mutation_position_start)
functional_area = get_functional_area(exon) if exon is not None else None
domain_order = get_domain_order(mutation_position_start)

# 确保输入数据包含模型期望的特征名称
if exon is not None and functional_area is not None and domain_order is not None:
    input_data = {
        'Mutation_position_start': mutation_position_start,  # 修改名称
        'Mutation_position_stop': mutation_position_stop,    # 修改名称
        'Mutation_types': mutation_type,                     # 修改名称
        'exon': exon,
        'Functional_area': functional_area,
        'Domain_order': domain_order,
        'frame_of_exons': 1 if exon in [1, 2, 6, 7, 8, 11, 12, 17, 18, 19, 20, 21, 22, 43, 44, 45, 46, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 61, 62, 63, 65, 66, 67, 68, 69, 70, 75, 76, 78, 79] else 0,
        'skipping_of_in_frame_exons': 1 if exon in [9, 25, 27, 29, 31, 37, 38, 39, 41, 72, 74] else 0,
        'frame': 1 if mutation_type in [1, 2] else 0,
        # 模型期望的其他特征名称和默认值
        'Amino_acid_properties_changed': -999,   # 添加缺失的特征
        'SpliceAI_pred_DS_DL': -999,             # 添加缺失的特征
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
