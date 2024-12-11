import streamlit as st
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
import py3Dmol

# Load dataset
dataset = pd.read_csv('your_dataset.csv')  # Replace with your dataset file path
dataset.columns = dataset.columns.str.strip()  # Clean column names

# Custom CSS for styling
st.markdown("""
    <style>
        body {
            background-color: #f9f9fb;
            color: #333;
            font-family: 'Arial', sans-serif;
        }
        .header {
            font-size: 2.5rem;
            font-weight: bold;
            text-align: center;
            color: #4CAF50;
            margin-bottom: 1rem;
        }
        .sub-header {
            text-align: center;
            color: #555;
            font-size: 1.2rem;
            margin-bottom: 2rem;
        }
        .if-value {
            font-size: 1.5rem;
            font-weight: bold;
            color: #FF6347;  
            background-color: #FFE4E1; 
            padding: 7px;
            border-radius: 4px;
            text-align: center;
            margin-bottom: 0.5rem;
        }
        .details-title {
            font-size: 1.5rem;
            color: #333;
            margin-bottom: 10px;
        }
        .stTable {
            margin-top: 1rem;
        }
    </style>
""", unsafe_allow_html=True)

# App title and description
st.markdown('<div class="header">🔬 Imprinting Factor (IF) Prediction</div>', unsafe_allow_html=True)
st.markdown('<div class="sub-header">🔍 Search by Template (Compound Name) or SMILES to view IF and details.</div>', unsafe_allow_html=True)

# Input section
input_query = st.text_input(
    "Enter Template (Compound Name) or SMILES",
    placeholder="E.g., compound name or CC(=O)OC1=CC=CC=C1C(=O)O"
)

def makeblock(smi):
    """Convert SMILES to a 3D MolBlock."""
    mol = Chem.MolFromSmiles(smi)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol)
    mblock = Chem.MolToMolBlock(mol)
    return mblock

def render_mol(xyz):
    """Render a molecule in 3D using py3Dmol."""
    viewer = py3Dmol.view(width=800, height=400)
    viewer.addModel(xyz, "mol")
    viewer.setStyle({'stick': {}})
    viewer.setBackgroundColor("white")
    viewer.zoomTo()
    return viewer._make_html()

if input_query:
    # Filter dataset for matches
    matches = dataset[
        (dataset['Compound name'].str.lower() == input_query.lower()) |
        (dataset['template smile'].str.lower() == input_query.lower()) |
        (dataset['momomer smile'].str.lower() == input_query.lower()) |
        (dataset['crosslinker smile'].str.lower() == input_query.lower())
    ]

    if matches.empty:
        st.warning("⚠️ No matching Template or SMILES found.")
    else:
        # Create a dropdown to select the IF value from the matches
        if_values = matches[['Compound name', 'IF']].drop_duplicates()
        selected_if = st.selectbox("Select an IF value", if_values['IF'].unique())

        # Filter the dataset for the selected IF value
        selected_data = matches[matches['IF'] == selected_if].iloc[0]

        # Display IF Value
        st.markdown("###  Retrieved Imprinting Factor (IF) Value")
        st.markdown(f'<div class="if-value">IF Value: {selected_data["IF"]:.2f}</div>', unsafe_allow_html=True)

        # Display Template (Compound Name), Crosslinker, and Functional Monomer
        st.markdown('<div class="details-section">', unsafe_allow_html=True)
        st.markdown('<div class="details-title" style="color: #FFD700;">📋 Compound name and Related Details</div>', unsafe_allow_html=True)
        st.write(f"**🧪 Template:** {selected_data['Compound name']}")
        st.write(f"**🧬 Functional Monomer:** {selected_data['FUNCTIONAL MONOMER']}")
        st.write(f"**🔗 Crosslinker:** {selected_data['crosslinker']}")
        st.write(f"**⚗️ Template Concentration:** {selected_data['tempConc(mmol)']} mmol")
        st.write(f"**⚗️ Monomer Concentration:** {selected_data['MonomerConc(mmol)']} mmol")
        st.write(f"**⚗️ Crosslinker Concentration:** {selected_data['CROSSLINKER CONC(mmol)']} mmol")
        st.markdown('</div>', unsafe_allow_html=True)

        # Display SMILES and generate 3D models
        smiles_data = {
            "SMILES Type": ["Template SMILES", "Monomer SMILES", "Crosslinker SMILES"],
            "SMILES": [
                selected_data["template smile"],
                selected_data["momomer smile"],
                selected_data["crosslinker smile"]
            ]
        }
        smiles_df = pd.DataFrame(smiles_data)
        st.table(smiles_df)

        # Render 3D models for each SMILES
        st.markdown("### 🧪 3D Molecular Structures")
        for i, row in smiles_df.iterrows():
            smi_type = row["SMILES Type"]
            smi = row["SMILES"]
            st.subheader(f"{smi_type}")
            try:
                mol_block = makeblock(smi)
                mol_html = render_mol(mol_block)
                st.components.v1.html(mol_html, height=450, scrolling=False)
            except Exception as e:
                st.warning(f"Could not generate 3D model for {smi_type}: {e}")
