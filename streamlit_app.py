from collections import namedtuple
import altair as alt
import math
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import streamlit as st


"""
# Newly repopulated spinal cord microglia exhibit a unique transcriptome and coincide with sex-independent pain resolution

These data originate from a study in the tibial fracture/casting mouse model of Complex Regional Pain Syndrome (CRPS) 
which showed microglial depletion (using the Cx3CR1-creER;iDTR-LSL transgenic mouse) after injury attenuates pain and
improves peripheral symptoms (hind paw edema and warmth). Lumbar spinal cord microglia were FACS-Isolated
(CD45mid CD11+ CX3CR1+) and RNA-Sequenced to study the transcriptome of resident and repopulating microglia.
Sequencing data were analyzed using Deseq2.

The following groups were analyzed and compared:

* Uninjured-Resident = microglia isolated from uninjured mice with resident microglia
* Injured-Resident = microglia isolated from mice 5 weeks after tibial fracture/casting with resident microglia
* Uninjured-Repopulated = microglia isolated from uninjured mice with repopulated microglia (2 weeks after depletion/repopulation)
* Injured-Repopulated = microglia isolated from mice 5 weeks after tibial fracture/casting with repopulated microglia (2 weeks after depletion/repopulation)

[https://www.biorxiv.org/content/10.1101/2022.12.20.521295v1](https://www.biorxiv.org/content/10.1101/2022.12.20.521295v1)
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
        gene_select = st.selectbox(label="Select a gene:", options=df_male["Gene Symbol"].unique())
    with col2:
        comparison_select = st.selectbox(label="Comparison:", options=comparisons)

# Data for selected gene
all_gene_data_male = df_male[df_male["Gene Symbol"] == gene_select]
all_gene_data_female = df_female[df_female["Gene Symbol"] == gene_select]

# Get avg p value by transcript_id
overlapping_transcript_ids = set(all_gene_data_male["Transcript ID"]).union(set(all_gene_data_female["Transcript ID"]))
pvalue_col = "Adj-p " + comparison_select
pvalues = pd.DataFrame()
pvalues["Transcript ID"] = list(overlapping_transcript_ids)
pvalues = pvalues.merge(all_gene_data_male[["Transcript ID", pvalue_col]], on="Transcript ID")
pvalues = pvalues.merge(all_gene_data_female[["Transcript ID", pvalue_col]], on="Transcript ID")
pvalues["avg pvalue"] = (pvalues[pvalue_col+"_x"] + pvalues[pvalue_col+"_y"]) / 2
pvalues.sort_values("avg pvalue", inplace=True)
sorted_transcript_ids = pvalues["Transcript ID"].tolist()


def plot_single_transcript(transcript_id, gene_data_male, gene_data_female):
    """Plot gene data for a single transcript."""
    
    # Get columns to plot
    col1, col2 = comparison_select.split(" vs ")
    col1_full = f"FC {comparison_select}: {col1}" 
    col2_full = f"FC {comparison_select}: {col2}" 

    # P-values
    p_value_threshold = 0.05
    male_p_value = gene_data_male["Adj-p " + comparison_select]
    female_p_value = gene_data_female["Adj-p " + comparison_select]

    # Bar chart
    fig = go.Figure()
    fig.add_trace(
        go.Bar(
            name=col2,
            x=["Male", "Female"],
            y=[gene_data_male[col2_full], gene_data_female[col2_full]], 
        )
    )
    fig.add_trace(
        go.Bar(
            name=col1,
            x=["Male", "Female"],
            y=[gene_data_male[col1_full], gene_data_female[col1_full]], 
        )
    )
    star_y = 1.05 * np.max([gene_data_male[col1_full], gene_data_female[col1_full], gene_data_male[col2_full], gene_data_female[col2_full]])
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
        title_text=f"{gene_select} ({transcript_id})",
        yaxis_title="Transcript Count Fold Change",
        font=dict(size=12),
    )
    for trace in fig['data']: 
        if(trace['name'] not in [col1, col2]):
            trace['showlegend'] = False
    st.plotly_chart(fig, use_container_width=True)

    # Show exact p-values
    st.write("P values:")
    st.write(f"Male: p=*{male_p_value}*")
    st.write(f"Female: p=*{female_p_value}*")


# Show a plot for each transcript with this gene symbol
for transcript_id in sorted_transcript_ids:
    gene_data_male = all_gene_data_male.loc[all_gene_data_male["Transcript ID"] == transcript_id, :].iloc[0, :]
    gene_data_female = all_gene_data_female.loc[all_gene_data_female["Transcript ID"] == transcript_id, :].iloc[0, :]
    plot_single_transcript(
        transcript_id,
        gene_data_male,
        gene_data_female,
    )
    if len(overlapping_transcript_ids) > 1:
        st.markdown("""---""")
