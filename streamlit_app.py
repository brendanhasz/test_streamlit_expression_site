from collections import namedtuple
import altair as alt
import math
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import streamlit as st


"""
# Repopulated microglia exhibit a unique transcriptome and contribute to sex-independent pain resolution

This data originates from a study in a mouse model of Complex Regional Pain Syndrome that showed microglial
depletion after injury attenuates pain and improves peripheral symtoms. Lumbar spinal cord microglia were
FACS-Isolated and RNA-Sequenced to study the transcriptome of the repopulating microglia. Sequencing data was
analyzed using [Deseq2](https://doi.org/doi:10.18129/B9.bioc.DESeq2).
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
    "Injured-Resident vs Uninjured-Resident",
    "Injured-Repopulated vs Injured-Resident",
    "Injured-Repopulated vs Uninjured-Repopulated",
    "Uninjured-Repopulated vs Uninjured-Resident",
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

# P-values
p_value_threshold = 0.05
st.write(comparison_select)
st.write(gene_data_male.columns)
male_p_value = gene_data_male["Adj-p " + comparison_select]
female_p_value = gene_data_female["Adj-p " + comparison_select]

# Bar chart
fig = go.Figure()
fig.add_trace(
    go.Bar(
        name=col2,
        x=["Male", "Female"],
        y=[gene_data_male[col2], gene_data_female[col2]], 
        error_y=dict(type='data', array=[gene_data_male["SEM "+col2], gene_data_female["SEM "+col2]])
    )
)
fig.add_trace(
    go.Bar(
        name=col1,
        x=["Male", "Female"],
        y=[gene_data_male[col1], gene_data_female[col1]], 
        error_y=dict(type='data', array=[gene_data_male["SEM "+col1], gene_data_female["SEM "+col1]])
    )
)
star_y = 1.05 * (
    np.max([gene_data_male[col1], gene_data_female[col1], gene_data_male[col2], gene_data_female[col2]])
    + np.max([gene_data_male["SEM "+col1], gene_data_female["SEM "+col1], gene_data_male["SEM "+col2], gene_data_female["SEM "+col2]])
)
if male_p_value < p_value_threshold:
    fig.add_trace(
        go.Scatter(
            name=f"M p={male_p_value}",
            x=["Male"],
            y=[star_y],
            mode="markers",  # or "lines"
            marker_symbol="asterisk",  # or "star"
            marker_line_color="black",
            marker_color="black",
            marker_line_width=1,
            marker_size=10,
        )
    )
else:
    fig.add_annotation(
        text="n.s.",
        x="Male",
        y=star_y,
        showarrow=False,
        align="center",
        xanchor="center",
        font=dict(size=11, color="black"),
    )
if female_p_value < p_value_threshold:
    fig.add_trace(
        go.Scatter(
            name=f"F p={male_p_value}",
            x=["Female"],
            y=[star_y],
            mode="markers",  # or "lines"
            marker_symbol="asterisk",  # or "star"
            marker_line_color="black",
            marker_color="black",
            marker_line_width=1,
            marker_size=10,
        )
    )
else:
    fig.add_annotation(
        text="n.s.",
        x="Female",
        y=star_y,
        showarrow=False,
        align="center",
        xanchor="center",
        font=dict(size=11, color="black"),
    )
fig.update_layout(
    barmode='group',
    title_text=gene_select,
    yaxis_title="Normalized Transcript Counts",
    font=dict(size=12),
)
for trace in fig['data']: 
    if(trace['name'] not in [col1, col2]):
        trace['showlegend'] = False
st.plotly_chart(fig, use_container_width=True)

col1, col2 = st.columns(2)

# Log fold change
with col1:
    st.write("Log2 of Transcript Count Fold Change:")
    st.write(f"Male: *{gene_data_male['LogFC ' + comparison_select]}*")
    st.write(f"Female: *{gene_data_female['LogFC ' + comparison_select]}*")

# P-values
with col2:
    st.write("P values:")
    st.write(f"Male: p=*{male_p_value}*")
    st.write(f"Female: p=*{female_p_value}*")
