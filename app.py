import streamlit as st
import joblib
import pandas as pd
import numpy as np

# Path to your model file
model_path = "voting_clf.pkl"  # 使用相对路径
model = joblib.load(model_path)  # Load the model


# Set Streamlit page title
st.title("DMD/BMD Prediction Model")

# Input features form
st.write("Please enter the following feature values for prediction:")
features = {
    'Functional_area': st.number_input("Functional Area", value=0.0),
    'frame_of_exons': st.number_input("Frame of Exons", value=0.0),
    'Amino_acid_properties_changed': st.number_input("Amino Acid Properties Changed", value=0.0),
    'exon': st.number_input("Exon", value=0),
    'Mutation_position_start': st.number_input("Mutation Position Start", value=0),
    'Mutation_position_stop': st.number_input("Mutation Position Stop", value=0),
    'frame': st.number_input("Frame", value=0.0),
    'Mutation_types': st.number_input("Mutation Types", value=0.0),
    'Domain_order': st.number_input("Domain Order", value=0.0),
    'skipping_of_in_frame_exons': st.number_input("Skipping of In-frame Exons", value=0.0),
    'SpliceAI_pred_DS_DL': st.number_input("SpliceAI Prediction DS_DL", value=0.0),
    'CADD_PHRED': st.number_input("CADD PHRED", value=0.0),
    'CADD_RAW': st.number_input("CADD RAW", value=0.0),
    'GERP++_NR': st.number_input("GERP++ NR", value=0.0),
    'GERP++_RS': st.number_input("GERP++ RS", value=0.0),
    'GERP++_RS_rankscore': st.number_input("GERP++ RS Rankscore", value=0.0),
    'BayesDel_addAF_score': st.number_input("BayesDel addAF Score", value=0.0),
    'BayesDel_noAF_rankscore': st.number_input("BayesDel noAF Rankscore", value=0.0),
    'BayesDel_noAF_score': st.number_input("BayesDel noAF Score", value=0.0),
    'DANN_rankscore': st.number_input("DANN Rankscore", value=0.0),
    'DANN_score': st.number_input("DANN Score", value=0.0),
    'PrimateAI_score': st.number_input("PrimateAI Score", value=0.0),
    'MetaLR_rankscore': st.number_input("MetaLR Rankscore", value=0.0),
}

# Predict button
if st.button("Predict"):
    # Convert input features to DataFrame
    input_data = pd.DataFrame([features])
    
    # Predict using the model
    prediction_proba = model.predict_proba(input_data)[:, 1][0]  # Get probability value
    prediction = model.predict(input_data)[0]  # Get classification result
    
    # Display prediction result
    st.write(f"Predicted Phenotype: {'DMD' if prediction == 1 else 'BMD'}")
    st.write(f"Prediction Probability (DMD): {prediction_proba:.2f}")
