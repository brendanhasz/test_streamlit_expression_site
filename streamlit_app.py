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
    "WT Fracture vs WT Control",
    "WT Fracture vs WT Control",
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

# Bar chart
fig = go.Figure()
fig.add_trace(go.Bar(
    name='Control',
    x=['Trial 1', 'Trial 2', 'Trial 3'], y=[3, 6, 4],
    error_y=dict(type='data', array=[1, 0.5, 1.5])
))
fig.add_trace(go.Bar(
    name='Experimental',
    x=['Trial 1', 'Trial 2', 'Trial 3'], y=[4, 7, 3],
    error_y=dict(type='data', array=[0.5, 1, 2])
))
fig.update_layout(barmode='group')
st.plotly_chart(fig, use_container_width=True)
