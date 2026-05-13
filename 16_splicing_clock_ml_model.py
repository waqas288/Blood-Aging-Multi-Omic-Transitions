import pandas as pd
import numpy as np
import os
from sklearn.linear_model import ElasticNetCV
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split
from sklearn.metrics import r2_score, mean_absolute_error
import matplotlib.pyplot as plt
import seaborn as sns

# Set working directory
os.chdir(r'c:\Users\WAQAS\Desktop\quick_paper\datasets\aging')

def run_splicing_clock():
    print("--- AGING PIPELINE: STEP 27 (THE SPLICING CLOCK ML MODEL) ---")
    
    # 1. LOAD DATA
    ELITE_FILE = "Final_Bulletproof_Survival_Results.csv"
    ISO_FILE = "filtered_blood_isoforms_tpm.csv"
    AGE_FILE = "GTEx_Blood_Ages.csv"
    
    if not os.path.exists(ELITE_FILE):
        print("Error: Results missing.")
        return

    # Load Elite targets
    df_elite = pd.read_csv(ELITE_FILE)
    target_isoforms = df_elite['Isoform'].unique()
    
    # Load Expression data
    iso_df = pd.read_csv(ISO_FILE, index_col=0)
    X = iso_df.loc[target_isoforms].T
    
    # Load Metadata (Sex, RIN, Cells)
    GENE_FILE = "filtered_blood_genes_tpm.csv"
    gene_df = pd.read_csv(GENE_FILE, index_col=0) # Use the gene file for markers
    gene_df.index = [str(i).split('.')[0] for i in gene_df.index]
    marker_map = {'ELANE': 'ENSG00000197561', 'CD3D': 'ENSG00000167286', 'CD14': 'ENSG00000170458'}
    marker_data = gene_df.loc[[e for e in marker_map.values() if e in gene_df.index]].T
    
    # Load Age and Phenotype data (Sex)
    PHENO_FILE = "GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt"
    ATTR_FILE = "GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt"
    ages_df = pd.read_csv(AGE_FILE)
    pheno_df = pd.read_csv(PHENO_FILE, sep='\t')[['SUBJID', 'SEX']]
    attr_df = pd.read_csv(ATTR_FILE, sep='\t', usecols=['SAMPID', 'SMRIN', 'SMTSISCH'])
    
    ages_df['SUBJID'] = ages_df['SampleID'].apply(lambda x: "-".join(x.split("-")[:2]))
    ages_df = ages_df.merge(pheno_df, on='SUBJID', how='left')
    ages_df = ages_df.merge(attr_df, left_on='SampleID', right_on='SAMPID', how='left')
    
    def get_midpoint(a):
        try:
            p = str(a).split('-')
            return (int(p[0]) + int(p[1])) / 2.0
        except: return np.nan
    ages_df['Age'] = ages_df['AgeBracket'].apply(get_midpoint)
    ages_df = ages_df.set_index('SampleID')
    
    # Align samples
    common_samples = X.index.intersection(ages_df.index).intersection(marker_data.index)
    X = X.loc[common_samples].astype(float)
    y = ages_df.loc[common_samples, 'Age'].astype(float)
    
    # 1. Log transform expression (Standard for RNA-seq)
    X = np.log2(X + 1)
    
    # 2. ADD COVARIATES (After Log Transform)
    X['Sex'] = ages_df.loc[common_samples, 'SEX']
    X['RIN'] = ages_df.loc[common_samples, 'SMRIN']
    for m_name, m_ensg in marker_map.items():
        if m_ensg in marker_data.columns:
            X[m_name] = marker_data.loc[common_samples, m_ensg]
    
    X = X.dropna()
    y = y.loc[X.index]
    
    print(f"Dataset Size: {len(X)} samples, {len(X.columns)} features.")

    # 2. TRAIN/TEST SPLIT
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=42)
    
    # 3. SCALE FEATURES
    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train)
    X_test_scaled = scaler.transform(X_test)
    
    # 4. TRAIN ELASTIC NET CV
    print("Training Elastic Net Splicing Clock (with 5-fold Cross-Validation)...")
    model = ElasticNetCV(l1_ratio=[.1, .5, .7, .9, .95, .99, 1], cv=5, random_state=42, max_iter=10000)
    model.fit(X_train_scaled, y_train)
    
    # 5. EVALUATE
    y_pred = model.predict(X_test_scaled)
    r2 = r2_score(y_test, y_pred)
    mae = mean_absolute_error(y_test, y_pred)
    
    print("\n--- MODEL PERFORMANCE ---")
    print(f"R-squared: {r2:.3f}")
    print(f"Median Absolute Error: {mae:.2f} years")

    # 6. PLOT RESULTS
    plt.figure(figsize=(10, 8))
    sns.regplot(x=y_test, y=y_pred, scatter_kws={'alpha':0.5}, line_kws={'color':'red'})
    plt.title(f'Splicing Clock: Chronological vs. Predicted Age\n(R2 = {r2:.3f}, MAE = {mae:.2f} yrs)', fontsize=14)
    plt.xlabel('Chronological Age (Years)', fontsize=12)
    plt.ylabel('Predicted "Splicing" Age (Years)', fontsize=12)
    plt.savefig('figures/Figure5_Splicing_Clock_Accuracy.png', dpi=300)
    plt.close()
    
    # 7. FEATURE IMPORTANCE (WEIGHTS)
    coefs = pd.Series(model.coef_, index=X.columns).sort_values(key=abs, ascending=False).head(15)
    # Map back to Gene Symbols
    mapping = df_elite.set_index('Isoform')['GeneSymbol'].to_dict()
    coef_df = pd.DataFrame({'Weight': coefs.values, 'Feature': [f"{mapping.get(idx, idx)}" for idx in coefs.index]})
    
    plt.figure(figsize=(10, 8))
    sns.barplot(data=coef_df, x='Weight', y='Feature', palette='vlag')
    plt.axvline(0, color='black', linewidth=1)
    plt.title('Top 15 Predictors (Splicing Clock Weights)', fontsize=14)
    plt.xlabel('Elastic Net Coefficient (Magnitude & Direction)', fontsize=12)
    plt.savefig('figures/Figure6_Clock_Weights.png', dpi=300)
    plt.close()

    print("\nSUCCESS: ML Clock figures saved to 'figures/' folder.")

if __name__ == "__main__":
    run_splicing_clock()
