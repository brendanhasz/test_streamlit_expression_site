from collections import namedtuple
import altair as alt
import math
import pandas as pd
import plotly.graph_objects as go
import streamlit as st


"""
# Page title

Introduction to project
"""


@st.experimental_singleton
def get_data_male():
    """Load the male data, and cache it."""
    return pd.read_csv(
        "https://raw.githubusercontent.com/brendanhasz/test_streamlit_expression_site/main/data/deseq2_male.csv"
    )


@st.experimental_singleton
def get_data_female():
    """Load the female data, and cache it."""
    return pd.read_csv(
        "https://raw.githubusercontent.com/brendanhasz/test_streamlit_expression_site/main/data/deseq2_female.csv"
    )


df_male = get_data_male()
df_female = get_data_female()

comparisons = [
    "WT Fracture vs WT Control",
    "DTR Fracture vs WT Fracture",
    "DTR Fracture vs DTR Control",
    "DTR Control vs WT Control",
]

# Input boxes
with st.container():
    col1, col2 = st.columns(2)
    with col1:
        gene_select = st.selectbox(label="Select a gene:", options=df_male["Gene Symbol"])
    with col2:
        comparison_select = st.selectbox(label="Comparison:", options=comparisons)

# Data for selected gene
gene_data_male = df_male[df_male["Gene Symbol"] == gene_select].iloc[0, :]
gene_data_female = df_female[df_female["Gene Symbol"] == gene_select].iloc[0, :]

# Get columns to plot
col1, col2 = comparison_select.split(" vs ")

# Bar chart
fig = go.Figure()
fig.add_trace(
    go.Bar(
        name=col1,
        x=["Male", "Female"],
        y=[gene_data_male[col1], gene_data_female[col1]], 
        error_y=dict(type='data', array=[gene_data_male["SEM "+col1], gene_data_female["SEM "+col1]])
    )
)
fig.add_trace(
    go.Bar(
        name=col2,
        x=["Male", "Female"],
        y=[gene_data_male[col2], gene_data_female[col2]], 
        error_y=dict(type='data', array=[gene_data_male["SEM "+col2], gene_data_female["SEM "+col2]])
    )
)
fig.update_layout(barmode='group')
st.plotly_chart(fig, use_container_width=True)

# Log fold change
fold_change_map = {
    "WT Fracture vs WT Control": "LogFC WT Fx vs WT Ctrl",
    "DTR Fracture vs WT Fracture": "LogFC DTR Fx vs WT FX",
    "DTR Fracture vs DTR Control": "LogFC DTR Fx vs DTR Ctrl",
    "DTR Control vs WT Control": "LogFC DTR Ctrl vs WT Ctrl",
}
st.write("Log2 of Transcript Count Fold Change:")
st.write(f"Male: *{gene_data_male[fold_change_map[comparison_select]]}*")
st.write(f"Female: *{gene_data_female[fold_change_map[comparison_select]]}*")

# P-values
pvalue_map = {
    "WT Fracture vs WT Control": "Adj-p WT Fx vs WT Ctrl",
    "DTR Fracture vs WT Fracture": "Adj-p DTR Fx vs WT FX",
    "DTR Fracture vs DTR Control": "Adj-p DTR Fx vs DTR Ctrl",
    "DTR Control vs WT Control": "Adj-p DTR Ctrl vs WT Ctrl",
}
st.write("P values:")
st.write(f"Male: p=*{gene_data_male[pvalue_map[comparison_select]]}*")
st.write(f"Female: p=*{gene_data_female[pvalue_map[comparison_select]]}*")
