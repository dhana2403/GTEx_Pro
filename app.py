import streamlit as st
import pandas as pd
import plotly.express as px

# ---------------------------
# Load Preprocessed GTEx Data
# ---------------------------
# For demo purposes, we'll simulate with a toy dataframe.
# Replace this with loading your preprocessed GTEx-Pro output (e.g., Parquet/HDF5).

data = {
    "gene": ["TP53", "TP53", "TP53", "BRCA1", "BRCA1", "BRCA1"],
    "tissue": ["Brain", "Liver", "Heart", "Brain", "Liver", "Heart"],
    "expression": [1.2, 0.8, 0.5, 0.6, 1.1, 0.9]
}

df = pd.DataFrame(data)

# ---------------------------
# Streamlit App Layout
# ---------------------------

st.set_page_config(page_title="GTEx Explorer", layout="centered")
st.title("GTEx Explorer")

# Gene input box
gene_id = st.text_input("Enter Gene ID", "TP53")

# Plot type options
plot_type = st.radio("Select plot type", ["Boxplot", "Violin"], horizontal=True)

if gene_id:
    subset = df[df["gene"] == gene_id]

    if not subset.empty:
        if plot_type == "Boxplot":
            fig = px.box(subset, x="tissue", y="expression", points="all",
                         title=f"Expression of {gene_id} across tissues",
                         labels={"expression": "Expression level", "tissue": "Tissue"})
        else:
            fig = px.violin(subset, x="tissue", y="expression", box=True, points="all",
                            title=f"Expression of {gene_id} across tissues",
                            labels={"expression": "Expression level", "tissue": "Tissue"})

        st.plotly_chart(fig, use_container_width=True)
    else:
        st.warning(f"No data available for {gene_id}.")

st.caption("Data source: Preprocessed GTEx")
