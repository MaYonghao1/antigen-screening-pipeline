import streamlit as st
import os
import zipfile
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
from antigen_evaluator import run_cascading_pipeline

st.set_page_config(page_title="Cascading Antigen Screening Pipeline", layout="wide")
st.title("🧬 High-Throughput Cascading Antigen Screening System")

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
mapped_csv = os.path.join(BASE_DIR, "local_antigen_db", "parent_mapped_antigens.csv")
mapped_zip = os.path.join(BASE_DIR, "local_antigen_db", "parent_mapped_antigens.zip")

@st.cache_data(ttl=3600)
def load_mapped_db(csv_path, zip_path):
    # 优先直接读取 CSV；若不存在则自动解压 ZIP
    if os.path.exists(csv_path):
        return pd.read_csv(csv_path)
    elif os.path.exists(zip_path):
        with zipfile.ZipFile(zip_path, 'r') as z:
            z.extractall(os.path.dirname(zip_path))
        return pd.read_csv(csv_path)
    return pd.DataFrame()

df_raw_full = load_mapped_db(mapped_csv, mapped_zip)

if len(df_raw_full) == 0:
    st.error("❌ Warning: Database not found! Please ensure `parent_mapped_antigens.zip` exists.")
else:
    st.toast(f"✅ Successfully loaded ALL {len(df_raw_full)} parent antigen entries.", icon="📦")

st.sidebar.header("🎛️ Cascading Filter Thresholds")

tax_option = st.sidebar.selectbox(
    "1. Taxonomy Filtering Strategy",
    ["Xeno-Antigens (Homo/Mus Excluded)", "All Species", "Non-Mouse (Mus Excluded)", "Non-Human (Homo Excluded)", "Human Only (Homo sapiens)", "Mouse Only (Mus musculus)", "Bacterial Only", "Viral Only", "Others"]
)

min_length = st.sidebar.slider("2. Min Sequence Length (AA)", min_value=50, max_value=1000, value=100)
min_mw, max_mw = st.sidebar.slider("3. Molecular Weight Range (kDa)", min_value=10, max_value=200, value=(20, 100))
max_instability = st.sidebar.slider("4. Max Instability Index", min_value=10, max_value=100, value=45)
min_density = st.sidebar.slider("5. Min Epitope Density (%)", min_value=0.0, max_value=10.0, value=0.1, step=0.1)

filters = {
    "tax_option": tax_option,
    "min_length": min_length,
    "min_mw": min_mw, "max_mw": max_mw,
    "max_instability": max_instability,
    "min_density": min_density
}

res = run_cascading_pipeline(df_raw_full, filters)
counts = res['counts']
df_final = res['df5']

plotly_config = {
    'toImageButtonOptions': {'format': 'svg', 'filename': 'Antigen_Screening_Plot', 'height': 600, 'width': 1000, 'scale': 2},
    'displayModeBar': True
}

st.subheader("🔻 Cascading Filtration Funnel Chart")
fig_funnel = go.Figure(go.Funnel(
    y = ['Taxonomy Filter', 'Length Cutoff', 'MW Range Cutoff', 'Stability Cutoff', 'Final Candidates'],
    x = [counts['Taxonomy'], counts['Length'], counts['MW'], counts['Stability'], counts['Final']],
    textinfo = "value+percent initial",
    marker = {"color": ["#2b5c8f", "#4682b4", "#e67e22", "#d35400", "#27ae60"]}
))
fig_funnel.update_layout(margin=dict(l=20, r=20, t=20, b=20))
st.plotly_chart(fig_funnel, use_container_width=True, config=plotly_config)

st.subheader("📊 Candidate Antigen Distribution Scatter Plot")
total_final = len(df_final)

if total_final == 0:
    st.warning("⚠️ Thresholds are too strict! No candidate antigens matched the criteria.")
else:
    top_display_limit = min(50, total_final) if total_final > 1 else total_final
    df_scatter_data = df_final.head(top_display_limit)

    fig_scatter = px.scatter(
        df_scatter_data, 
        x="MW_kDa", 
        y="Instability_Index", 
        size="Epitope_Density_%" if df_scatter_data["Epitope_Density_%"].sum() > 0 else None,
        color="Composite_Score",
        text="Antigen_Name",
        hover_data=["Organism", "Seq_Length", "Antigen_ID"],
        title=f"Top {len(df_scatter_data)} Selected Antigens (Click camera button to export Vector SVG/PDF)"
    )
    fig_scatter.update_traces(textposition='top center', marker=dict(line=dict(width=1, color='DarkSlateGrey')))
    fig_scatter.update_layout(height=700, margin=dict(l=20, r=20, t=50, b=20))
    st.plotly_chart(fig_scatter, use_container_width=True, config=plotly_config)

st.subheader("🏆 Ranked Candidate Antigens Leaderboard")
display_cols = ["Antigen_Name", "Organism", "Seq_Length", "MW_kDa", "Instability_Index", "Epitope_Density_%", "Composite_Score", "Antigen_ID"]
st.dataframe(df_final[[c for c in display_cols if c in df_final.columns]], use_container_width=True)

with st.expander("🔍 Detailed Data Inspection at Each Stage"):
    t1, t2, t3, t4 = st.tabs(["Post Taxonomy", "Post Length", "Post MW", "Post Stability"])
    safe_cols = [c for c in display_cols if c in res['df1'].columns]
    t1.dataframe(res['df1'][safe_cols].head(100), use_container_width=True)
    t2.dataframe(res['df2'][safe_cols].head(100), use_container_width=True)
    t3.dataframe(res['df3'][safe_cols].head(100), use_container_width=True)
    t4.dataframe(res['df4'][safe_cols].head(100), use_container_width=True)
