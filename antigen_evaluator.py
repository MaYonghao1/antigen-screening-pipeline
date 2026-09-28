import os
import pandas as pd
import numpy as np
import gc 
def run_cascading_pipeline(df_raw, filters, save_dir="./output_results"):
    df0 = df_raw 
    
    if 'Composite_Score' not in df0.columns:
        df0['Composite_Score'] = np.round(
            df0['Epitope_Density_%'] * 40.0 + 
            (100.0 - df0['Instability_Index']) * 0.3 + 
            (df0['Seq_Length'] / 10.0), 2
        )
    
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
    else: mask1 = pd.Series(True, index=df0.index)
        
    mask2 = mask1 & (df0['Seq_Length'] >= filters['min_length'])
    mask3 = mask2 & (df0['MW_kDa'] >= filters['min_mw']) & (df0['MW_kDa'] <= filters['max_mw'])
    mask4 = mask3 & (df0['Instability_Index'] <= filters['max_instability'])
    mask5 = mask4 & (df0['Epitope_Density_%'] >= filters['min_density'])
    
    df_final = df0[mask5].sort_values(by="Composite_Score", ascending=False).reset_index(drop=True)
    
    stage_counts = {
        "Taxonomy": int(mask1.sum()),
        "Length": int(mask2.sum()),
        "MW": int(mask3.sum()),
        "Stability": int(mask4.sum()),
        "Final": int(mask5.sum())
    }
    
    df1_head = df0[mask1].head(100)
    df2_head = df0[mask2].head(100)
    df3_head = df0[mask3].head(100)
    df4_head = df0[mask4].head(100)
    
    del org_clean, mask1, mask2, mask3, mask4, mask5
    gc.collect()
    
    return {
        "df1": df1_head, "df2": df2_head, "df3": df3_head, "df4": df4_head, "df5": df_final,
        "counts": stage_counts
    }
