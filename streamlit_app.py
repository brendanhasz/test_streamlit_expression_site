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

data_url = "https://raw.githubusercontent.com/brendanhasz/test_expression_website/main/data/data.csv"
df = pd.read_csv(data_url)

comparisons = [
    "WT Fracture vs WT Control",
    "DTR Fracture vs WT Fracture",
    "DTR Fracture vs DTR Control",
    "DTR Control vs WT Control",
]

# Input boxes
with st.container():
    col1, col2, col3 = st.columns(3)
    with col1:
        gene_select = st.selectbox(label="Select a gene:", options=df["Gene Symbol"])
    with col2:
        comparison_select = st.selectbox(label="Comparison:", options=comparisons)
    with col3:
        gender_select = st.selectbox(label="Gender:", options=["Male", "Female"])

# Data for selected gene
gene_data = df[df["Gene Symbol"] == gene_select].iloc[0, :]

# Get columns to plot
col1, col2 = comparison_select.split(" vs ")

# Bar chart
fig = go.Figure()
fig.add_trace(go.Bar(
    name="Male",
    x=[col1, col2],
    y=[gene_data[col1], gene_data[col2]], 
    error_y=dict(type='data', array=[gene_data["SEM "+col1], gene_data["SEM "+col2]])
))
fig.update_layout(barmode='group')
st.plotly_chart(fig, use_container_width=True)
