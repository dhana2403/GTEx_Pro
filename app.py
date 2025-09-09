import streamlit as st
import pandas as pd
import requests
import json

st.title("GTEx Expression Viewer")

# Load JSON mapping of tissue -> CSV link
json_url = "https://www.dropbox.com/scl/fi/mjk0737mb30c2tkw7y7y3/files.json?rlkey=e6jy92mfmy61asx51t9cax8f3&st=6ra8e95i&dl=1"
response = requests.get(json_url)
tissue_map = response.json()

# Dropdown for tissues
tissue = st.selectbox("Select a tissue:", list(tissue_map.keys()))

# Load expression data
df = pd.read_csv(tissue_map[tissue], index_col=0)

# Select gene(s)
genes = st.multiselect("Select gene(s):", df.index)

if genes:
    sub_df = df.loc[genes].T  # samples = rows
    st.write(sub_df.head())   # preview

    # Example visualization
    import plotly.express as px
    for gene in genes:
        fig = px.box(sub_df, y=gene, title=f"{gene} expression in {tissue}")
        st.plotly_chart(fig)

