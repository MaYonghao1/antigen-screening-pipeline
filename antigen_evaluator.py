import os
import re
import pandas as pd
import numpy as np

def run_cascading_pipeline(df_raw, filters, save_dir="./output_results"):
    df0 = df_raw.copy()

    # Ensure required columns exist
    for col in ['Organism', 'Epitope_Density_%', 'Instability_Index', 'Seq_Length', 'MW_kDa']:
        if col not in df0.columns:
            if col == 'Organism': df0[col] = "Unknown"
            elif col == 'MW_kDa': df0[col] = np.round(df0['Seq_Length'] * 0.11, 2)
            else: df0[col] = 0.0

    # Calculate Composite Score
    df0['Composite_Score'] = np.round(
        df0['Epitope_Density_%'] * 40.0 + 
        (100.0 - df0['Instability_Index']) * 0.3 + 
        (df0['Seq_Length'] / 10.0), 2
    )

    # Strict Taxonomy Filtering
    tax_opt = filters['tax_option']
    org_clean = df0['Organism'].astype(str).str.lower().str.replace('_', ' ')

    is_mouse = org_clean.str.contains(r'mus musculus|mouse|10090', regex=True, na=False)
    is_human = org_clean.str.contains(r'homo sapiens|human|9606', regex=True, na=False)
    is_bact = org_clean.str.contains(r'bacteri|coli|staph|strep|salmonella|mycobacteri|bacillus|clostridium', regex=True, na=False)
    is_viral = org_clean.str.contains(r'virus|viral|phage|sars|cov|hiv|influenza', regex=True, na=False)

    if tax_opt == "Non-Mouse (Mus Excluded)": mask1 = ~is_mouse
    elif tax_opt == "Non-Human (Homo Excluded)": mask1 = ~is_human
    elif tax_opt == "Xeno-Antigens (Homo/Mus Excluded)": mask1 = (~is_human) & (~is_mouse)
    elif tax_opt == "Human Only (Homo sapiens)": mask1 = is_human
    elif tax_opt == "Mouse Only (Mus musculus)": mask1 = is_mouse
    elif tax_opt == "Bacterial Only": mask1 = is_bact
    elif tax_opt == "Viral Only": mask1 = is_viral
    elif tax_opt == "Others": mask1 = (~is_human) & (~is_mouse) & (~is_bact) & (~is_viral)
    else: mask1 = np.ones(len(df0), dtype=bool)

    df1 = df0[mask1].copy()

    # Step-by-step Cascading Filters
    df2 = df1[df1['Seq_Length'] >= filters['min_length']].copy()
    df3 = df2[(df2['MW_kDa'] >= filters['min_mw']) & (df2['MW_kDa'] <= filters['max_mw'])].copy()
    df4 = df3[df3['Instability_Index'] <= filters['max_instability']].copy()
    df5 = df4[df4['Epitope_Density_%'] >= filters['min_density']].copy()

    df_final = df5.sort_values(by="Composite_Score", ascending=False).reset_index(drop=True)

    stage_counts = {
        "Taxonomy": len(df1),
        "Length": len(df2),
        "MW": len(df3),
        "Stability": len(df4),
        "Final": len(df5)
    }

    if not os.path.exists(save_dir):
        try: os.makedirs(save_dir, exist_ok=True)
        except Exception: save_dir = "."
    out_path = os.path.join(save_dir, "Cascading_Filtered_Antigens.csv")
    try:
        df_final.to_csv(out_path, index=False, encoding="utf-8-sig")
    except Exception:
        pass

    return {
        "df1": df1, "df2": df2, "df3": df3, "df4": df4, "df5": df_final,
        "counts": stage_counts, "out_path": out_path
    }