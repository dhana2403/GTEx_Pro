import streamlit as st
import pandas as pd
import requests
import json
import plotly.express as px
import plotly.figure_factory as ff

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

# Checkbox: show correlation matrix
show_corr = st.checkbox("Show gene-gene correlation heatmap")

if genes:
    # Filter and transpose: samples as rows, genes as columns
    sub_df = df.loc[genes].T
    st.write("Expression Table (samples x genes):")
    st.write(sub_df.head())

    # Reshape to long format for violin plot
    melted = sub_df.reset_index().melt(id_vars="index", 
                                      var_name="Gene", 
                                      value_name="Expression")
    melted = melted.rename(columns={"index": "Sample"})

    # Violin plot
    fig_violin = px.violin(
        melted,
        x="Gene",
        y="Expression",
        box=True,
        points="all",
        title=f"Expression of Selected Genes in {tissue}"
    )
    for trace in fig_violin.data:
        if trace.type == 'violin':
            trace.pointpos = 0
            trace.jitter = 0.05

    st.plotly_chart(fig_violin, use_container_width=True)

    # Correlation heatmap
    if show_corr and len(genes) > 1:
        corr_matrix = sub_df.corr(method='spearman')  # or 'spearman'
        st.write("Gene-Gene Correlation (Spearman):")
        fig_corr = ff.create_annotated_heatmap(
            z=corr_matrix.values,
            x=corr_matrix.columns.tolist(),
            y=corr_matrix.index.tolist(),
            colorscale='Viridis',
            showscale=True,
            reversescale=False,
            annotation_text=corr_matrix.round(2).values
        )
        fig_corr.update_layout(title_text=f"Gene-Gene Correlation in {tissue}")
        st.plotly_chart(fig_corr, use_container_width=True)


