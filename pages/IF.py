import streamlit as st
import pandas as pd
import joblib
import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors

# Load the trained model, scaler, and dataset
model = joblib.load('if_predictor_model.pkl')
scaler = joblib.load('scaler.pkl')
dataset = pd.read_csv('your_dataset.csv')  # Replace with your dataset file path

# Clean and inspect dataset columns
dataset.columns = dataset.columns.str.strip()  # Removes leading/trailing spaces

# Verify updated column names
expected_columns = [
    "Compound name", "M.wt", "tempConc(mmol)", "FUNCTIONAL MONOMER",
    "M.wt", "MonomerConc(mmol)", "crosslinker", "M.wt", "CROSSLINKER CONC(mmol)",
    "Solubility", "Xlogp", "H bond acceptor count", "H bond donor count",
    "IF", "MIP", "template smile", "crosslinker smile", "momomer smile"
]

# Inject custom CSS for styling
st.markdown("""
    <style>
        body {
            font-family: 'Arial', sans-serif;
            background-color: #f4f4f9;
            color: #333;
        }
        .header {
            text-align: center;
            color: #4CAF50;
            font-size: 3rem;
            font-weight: bold;
            margin-bottom: 0rem;
        }
        .sub-header {
            text-align: center;
            color: #777;
            font-size: 1.2rem;
            margin-top: 0rem;
            margin-bottom: 0rem;
        }
        .input-section {
            padding: 2rem;
            background-color: #ffffff;
            border-radius: 10px;
            box-shadow: 0px 4px 6px rgba(0, 0, 0, 0.1);
            margin-bottom: 2rem;
        }
        .button-section {
            display: flex;
            justify-content: center;
            margin-top: 1rem;
        }
        .stButton button {
            background-color: #FFC107;
            color: white;
            border-radius: 5px;
            font-size: 18px;
            width: 200px;
        }
        .stButton button:hover {
            background-color: #FF9800;
        }
        .stTable {
            margin-top: 2rem;
        }
        .footer {
            text-align: center;
            font-size: 1rem;
            color: #aaa;
            margin-top: 2rem;
        }
        .if-value {
            font-size: 2rem;
            font-weight: bold;
            color: #FF6347;  /* Tomato color for IF value */
        }
    </style>
""", unsafe_allow_html=True)

# App title and description
st.markdown('<div class="header">Imprinting Factor (IF) Prediction</div>', unsafe_allow_html=True)
st.markdown('<div class="sub-header">🔬 Provide a compound name to fetch details or predict the Imprinting Factor (IF).</div>', unsafe_allow_html=True)

# Input section for compound name
compound_name = st.text_input("🔎 Compound Name", placeholder="Enter compound name to search")
submit_button = st.button("🔍 Submit")  # Added Submit Button

# Remove the unwanted box if no content is displayed
if submit_button and compound_name:
    compound_data = dataset[dataset['Compound name'].str.lower() == compound_name.lower()]
    if not compound_data.empty:
        if len(compound_data) > 3:
            st.write(f"⚠️ Multiple entries ({len(compound_data)}) found for '{compound_name}'. Please select one:")

            # Allow user to select specific row
            selected_index = st.selectbox(
                "Select an entry:",
                compound_data.index.tolist(),
                format_func=lambda idx: f"Entry {idx + 1} - MIP: {compound_data.loc[idx, 'MIP']}"
            )

            selected_data = compound_data.loc[selected_index]

            st.write("### Selected Compound Details")
            st.table(selected_data)

            # Allow additional inputs for concentrations
            st.write("🔧 Adjust additional details:")
            template_conc = st.number_input("💧 Template Concentration (tempConc in mmol)", value=selected_data["tempConc(mmol)"])
            monomer_conc = st.number_input("💧 Monomer Concentration (MonomerConc in mmol)", value=selected_data["MonomerConc(mmol)"])
            crosslinker_conc = st.number_input("💧 Crosslinker Concentration (CROSSLINKER CONC in mmol)", value=selected_data["CROSSLINKER CONC(mmol)"])

            # Submit button to trigger IF prediction
            if st.button("🔮 submit"):
                # Create the input features for the model (including the SMILES and concentrations)
                features = [
                    selected_data['M.wt'], template_conc, selected_data['FUNCTIONAL MONOMER'], 
                    selected_data['MonomerConc(mmol)'], crosslinker_conc, selected_data['Solubility'], 
                    selected_data['Xlogp'], selected_data['H bond acceptor count'], 
                    selected_data['H bond donor count'], selected_data['MIP']
                ]
                
                # Normalize the features using the loaded scaler
                scaled_features = scaler.transform([features])

                # Make prediction using the model
                predicted_if = model.predict(scaled_features)[0]

                # Show the prediction result in a separate section
                st.markdown("### Predicted Imprinting Factor (IF)")
                st.markdown(f'<div class="if-value">IF Value: {predicted_if:.2f}</div>', unsafe_allow_html=True)
        else:
            # Display up to 3 matching rows
            st.write("### Compound Details")
            st.table(compound_data)
    else:
        st.write("⚠️ Compound not found in the dataset  CHOOSE  [ ML IF ] to predict using ML")
        

