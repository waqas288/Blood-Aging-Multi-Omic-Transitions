import pandas as pd
import numpy as np
import os
from statsmodels.duration.hazard_regression import PHReg
from statsmodels.stats.multitest import multipletests

# Set working directory
os.chdir(r'c:\Users\WAQAS\Desktop\quick_paper\datasets\aging')

def run_cell_adjusted_survival():
    print("--- AGING PIPELINE: STEP 9e (CELL-TYPE ADJUSTED SURVIVAL) ---")
    
    ROBUST_FILE = "Final_Robust_Elite_Survival.csv"
    ISO_FILE = "filtered_blood_isoforms_tpm.csv"
    GENE_FILE = "filtered_blood_genes_tpm.csv" # For marker genes
    AGE_FILE = "GTEx_Blood_Ages.csv"
    PHENO_FILE = "GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt"
    ATTR_FILE = "GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt"
    OUTPUT_FILE = "Final_Bulletproof_Survival_Results.csv"
    
    if not os.path.exists(ROBUST_FILE):
        print("Error: Previous survival results missing.")
        return

    print("[1/4] Loading Data and Defining Cell-Type Proxies...")
    elite_df = pd.read_csv(ROBUST_FILE)
    target_isoforms = elite_df['Isoform'].unique()
    
    # Manual Ensembl IDs for Blood Cell Type Markers (Gencode v26/v38 compatible)
    # Neutrophils: ELANE, MPO
    # T-cells: CD3D, CD8A
    marker_map = {
        'ELANE': 'ENSG00000197561',
        'MPO': 'ENSG00000105711',
        'CD3D': 'ENSG00000167286',
        'CD8A': 'ENSG00000153563'
    }
    
    # Load gene matrix to get markers
    gene_df = pd.read_csv(GENE_FILE, index_col=0)
    
    # Strip versions from gene_df index if present (e.g., ENSG000.1 -> ENSG000)
    gene_df.index = [str(i).split('.')[0] for i in gene_df.index]
    
    # Extract marker expression
    valid_ensgs = [e for e in marker_map.values() if e in gene_df.index]
    marker_data = gene_df.loc[valid_ensgs].T
    # Rename columns back to symbols for clarity
    inv_map = {v: k for k, v in marker_map.items()}
    marker_data = marker_data.rename(columns=inv_map)
    
    print("[2/4] Merging Technical and Biological Covariates...")
    # Load phenotypes and technical attrs (from Step 9d logic)
    ages_df = pd.read_csv(AGE_FILE)
    pheno_df = pd.read_csv(PHENO_FILE, sep='\t')[['SUBJID', 'SEX']]
    ages_df['SUBJID'] = ages_df['SampleID'].apply(lambda x: "-".join(x.split("-")[:2]))
    ages_df = ages_df.merge(pheno_df, on='SUBJID', how='left')
    
    attr_df = pd.read_csv(ATTR_FILE, sep='\t', usecols=['SAMPID', 'SMRIN', 'SMTSISCH'])
    ages_df = ages_df.merge(attr_df, left_on='SampleID', right_on='SAMPID', how='left')
    
    # Add Cell Type markers
    ages_df = ages_df.merge(marker_data, left_on='SampleID', right_index=True, how='left')
    
    def get_midpoint(a):
        try:
            p = str(a).split('-')
            return (int(p[0]) + int(p[1])) / 2.0
        except: return np.nan
        
    ages_df['time'] = ages_df['AgeBracket'].apply(get_midpoint)
    ages_df['event'] = 1
    # Final data cleaning
    covariates = ['SEX', 'SMRIN', 'SMTSISCH', 'ELANE', 'CD3D'] # Using top markers to avoid overfitting
    ages_df = ages_df.dropna(subset=['time'] + covariates).set_index('SampleID')
    
    print(f"  Cohort Size for Bulletproof Model: {len(ages_df)} patients.")

    print("[3/4] Loading Target Isoforms...")
    iso_df = pd.read_csv(ISO_FILE, index_col=0)
    iso_df = iso_df[iso_df.index.isin(target_isoforms)]
    valid_samples = ages_df.index.intersection(iso_df.columns[1:])
    iso_data = iso_df[valid_samples].astype(np.float32)
    ages_sub = ages_df.loc[valid_samples].copy()
    
    # Normalize
    np.random.seed(42)
    ages_sub['time_jittered'] = ages_sub['time'] + np.random.uniform(-4.5, 4.5, size=len(ages_sub))
    
    iso_mean = iso_data.mean(axis=1).values[:, np.newaxis]
    iso_std = iso_data.std(axis=1).values[:, np.newaxis]
    iso_norm = (iso_data.values - iso_mean) / (iso_std + 1e-8)
    iso_data = pd.DataFrame(iso_norm, index=iso_data.index, columns=iso_data.columns)
    
    print("[4/4] Running Bulletproof Cox Models (Adjusting for Cell Types)...")
    results = []
    
    for iso_id, row in iso_data.iterrows():
        model_data = pd.DataFrame({
            'time': ages_sub['time_jittered'].values,
            'event': ages_sub['event'].values,
            'expr': row.values,
            'sex': ages_sub['SEX'].values,
            'rin': ages_sub['SMRIN'].values,
            'isch': ages_sub['SMTSISCH'].values,
            'neu': ages_sub['ELANE'].values,
            'tcell': ages_sub['CD3D'].values
        }).apply(pd.to_numeric, errors='coerce').dropna()
        
        try:
            # Model: Survival ~ Expression + Sex + RIN + Isch + Neutrophils + T-cells
            mod = PHReg(model_data['time'], model_data[['expr', 'sex', 'rin', 'isch', 'neu', 'tcell']], status=model_data['event'])
            res = mod.fit()
            
            results.append({
                'GeneSymbol': elite_df[elite_df['Isoform'] == iso_id].iloc[0]['GeneSymbol'],
                'Isoform': iso_id,
                'Hazard_Ratio_Iso': np.exp(res.params[0]),
                'p_value_Iso': res.pvalues[0],
                'p_value_Sex': res.pvalues[1],
                'p_value_RIN': res.pvalues[2],
                'p_value_Isch': res.pvalues[3],
                'p_value_Neutrophils': res.pvalues[4],
                'p_value_Tcells': res.pvalues[5]
            })
        except Exception as e: 
            print(f"Model failed for {iso_id}: {e}")
            continue
            
    final_df = pd.DataFrame(results)
    if len(final_df) > 0:
        reject, p_adj, _, _ = multipletests(final_df['p_value_Iso'], method='fdr_bh')
        final_df['FDR_adj_p_value'] = p_adj
        final_df = final_df.sort_values('p_value_Iso')
        final_df.to_csv(OUTPUT_FILE, index=False)
        
        sig = final_df[final_df['FDR_adj_p_value'] < 0.05]
        print(f"SUCCESS! {len(sig)} Elite genes are BULLETPROOF (Significant after Cell-Type + Technical adjustment).")
        print(f"Results saved to '{OUTPUT_FILE}'")
        print("\nTop Bulletproof Hits:")
        print(final_df.head(10).to_string(index=False))
    else:
        print("Models failed.")

if __name__ == "__main__":
    run_cell_adjusted_survival()
