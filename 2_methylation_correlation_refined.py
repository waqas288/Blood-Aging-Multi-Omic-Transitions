import pandas as pd
import numpy as np
import os, gc
from scipy.stats import t
from statsmodels.stats.multitest import multipletests

def run_refined_methylation_discovery():
    print("--- AGING PIPELINE: STEP 2b (REFINED METHYLATION DISCOVERY) ---")
    
    geo_age_file = "C:/Users/WAQAS/Desktop/quick_paper/datasets/aging/GEO_GSE87571_Ages.csv"
    meth_folder = "C:/Users/WAQAS/Desktop/quick_paper/datasets/aging/GSE87571_matrix1of2 (1).txt/GSE87571_matrix1of2.txt"
    meth_file = os.path.join(meth_folder, "GSE87571_matrix1of2.txt")
    OUTPUT_FILE = "Significant_Aging_CpGs_Refined.csv"
    
    if not (os.path.exists(geo_age_file) and os.path.exists(meth_file)):
        print("Missing GEO files. Check GSE87571 data.")
        return

    geo_ages = pd.read_csv(geo_age_file)
    geo_ages['SampleID'] = geo_ages['SampleID'].astype(str)
    geo_ages = geo_ages.set_index('SampleID')
    geo_ages['Age'] = pd.to_numeric(geo_ages['Age'], errors='coerce')
    geo_ages = geo_ages.dropna()
    
    temp_df = pd.read_csv(meth_file, sep='\t', nrows=0)
    clean_cols = [str(x).replace('"', '').strip() for x in temp_df.columns]
    meth_samples = clean_cols[1:]
    valid_geo = [s for s in meth_samples if s in geo_ages.index]
    
    col_mapping = {clean: orig for clean, orig in zip(clean_cols, temp_df.columns)}
    pandas_valid_cols = [col_mapping[s] for s in valid_geo]
    
    age_vector = geo_ages.loc[valid_geo, 'Age'].values
    n = len(age_vector)
    age_centered = age_vector - np.mean(age_vector)
    age_std = np.sqrt(np.sum(age_centered**2) + 1e-8)
    
    print(f"Calculating Correlation and P-values for {n} patients...")
    
    chunksize = 20000
    meth_iter = pd.read_csv(meth_file, sep='\t', index_col=0, chunksize=chunksize, low_memory=False)
    
    all_results = []
    for i, chunk in enumerate(meth_iter):
        chunk_sub = chunk[pandas_valid_cols].apply(pd.to_numeric, errors='coerce').astype(np.float32)
        
        # Handle NaNs by filling with row mean (standard imputation for correlation)
        # This prevents NaN poisoning of the matrix multiplication
        chunk_sub = chunk_sub.T.fillna(chunk_sub.mean(axis=1)).T
        
        cpg_mean = chunk_sub.mean(axis=1).values[:, np.newaxis]
        cpg_centered = chunk_sub.values - cpg_mean
        
        cov = np.dot(cpg_centered, age_centered)
        cpg_std = np.sqrt(np.sum(cpg_centered**2, axis=1) + 1e-8)
        
        r_vals = cov / (cpg_std * age_std + 1e-8)
        
        if i == 0:
            print(f"    DEBUG: n={n}, age_std={age_std:.2f}, max_r={np.nanmax(np.abs(r_vals)):.4f}")
        
        # Calculate p-values using the t-distribution
        # t = r * sqrt((n-2)/(1-r^2))
        t_vals = r_vals * np.sqrt((n - 2) / (1 - r_vals**2 + 1e-8))
        p_vals = 2 * t.sf(np.abs(t_vals), n - 2)
        
        for idx in range(len(r_vals)):
            all_results.append({
                'CpG': chunk_sub.index[idx],
                'Pearson_R': r_vals[idx],
                'p_value': p_vals[idx]
            })
        
        if (i+1) % 5 == 0:
            print(f"  Processed {(i+1) * chunksize} CpGs...")
        gc.collect()

    print("\nApplying FDR Multiple Hypothesis Correction (Benjamini-Hochberg)...")
    res_df = pd.DataFrame(all_results).dropna(subset=['p_value'])
    
    if len(res_df) > 0:
        reject, p_adj, _, _ = multipletests(res_df['p_value'], method='fdr_bh')
        res_df['FDR_adj_p_value'] = p_adj
        
        # Final Filter: FDR < 0.05 AND |R| > 0.3
        sig_df = res_df[(res_df['FDR_adj_p_value'] < 0.05) & (np.abs(res_df['Pearson_R']) > 0.3)]
        
        sig_df.to_csv(OUTPUT_FILE, index=False)
        print(f"SUCCESS! Found {len(sig_df)} high-confidence Age-Associated CpGs.")
        print(f"Saved results to '{OUTPUT_FILE}'")
    else:
        print("Error: No valid p-values calculated.")

if __name__ == "__main__":
    run_refined_methylation_discovery()
