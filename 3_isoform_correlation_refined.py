import pandas as pd
import numpy as np
import os, gc
from scipy.stats import t
from statsmodels.stats.multitest import multipletests

def run_refined_isoform_discovery():
    print("--- AGING PIPELINE: STEP 6b (REFINED ISOFORM DISCOVERY) ---")
    
    iso_file = "filtered_blood_isoforms_tpm.csv"
    age_file = "GTEx_Blood_Ages.csv"
    OUTPUT_FILE = "Significant_Aging_Isoforms_Refined.csv"
    
    if not (os.path.exists(iso_file) and os.path.exists(age_file)):
        print("Missing GTEx files. Check preprocessing outputs.")
        return

    # Load and parse ages
    gtex_ages = pd.read_csv(age_file)
    def get_midpoint(age_str):
        try:
            parts = str(age_str).split('-')
            return (int(parts[0]) + int(parts[1])) / 2.0
        except: return np.nan
    
    gtex_ages['Age_Midpoint'] = gtex_ages['AgeBracket'].apply(get_midpoint)
    gtex_ages = gtex_ages.dropna().set_index('SampleID')
    
    print("  Loading Isoform matrix...")
    iso_df = pd.read_csv(iso_file, index_col=0)
    iso_samples = iso_df.columns[1:]
    
    valid_samples = gtex_ages.index.intersection(iso_samples)
    age_vector = gtex_ages.loc[valid_samples, 'Age_Midpoint'].values
    n = len(age_vector)
    iso_data = iso_df[valid_samples].apply(pd.to_numeric, errors='coerce').astype(np.float32)
    
    print(f"Calculating Correlation and P-values for {len(iso_data)} isoforms across {n} patients...")
    
    # Vectorized calculation
    age_centered = age_vector - np.mean(age_vector)
    iso_mean = iso_data.mean(axis=1).values[:, np.newaxis]
    iso_centered = iso_data.values - iso_mean
    
    cov = np.dot(iso_centered, age_centered)
    age_std = np.sqrt(np.sum(age_centered**2) + 1e-8)
    iso_std = np.sqrt(np.sum(iso_centered**2, axis=1) + 1e-8)
    
    r_vals = cov / (iso_std * age_std + 1e-8)
    
    # Calculate p-values
    t_vals = r_vals * np.sqrt((n - 2) / (1 - r_vals**2 + 1e-8))
    p_vals = 2 * t.sf(np.abs(t_vals), n - 2)
    
    all_results = pd.DataFrame({
        'Isoform': iso_data.index,
        'Gene': iso_df.iloc[:, 0].values,
        'Pearson_R': r_vals,
        'p_value': p_vals
    })
    
    print("\nApplying FDR Multiple Hypothesis Correction...")
    reject, p_adj, _, _ = multipletests(all_results['p_value'], method='fdr_bh')
    all_results['FDR_adj_p_value'] = p_adj
    
    # Final Filter: FDR < 0.05 AND |R| > 0.3
    sig_df = all_results[(all_results['FDR_adj_p_value'] < 0.05) & (np.abs(all_results['Pearson_R']) > 0.3)]
    
    sig_df.to_csv(OUTPUT_FILE, index=False)
    print(f"SUCCESS! Found {len(sig_df)} high-confidence Age-Associated Isoforms.")
    print(f"Saved results to '{OUTPUT_FILE}'")

if __name__ == "__main__":
    run_refined_isoform_discovery()
