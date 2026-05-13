import pandas as pd
import numpy as np
import os, gc

def calculate_stepminer(y_sorted):
    """
    Computes the optimal StepMiner binary split for 1D array.
    Returns: (best_index, min_sse, step_magnitude, direction)
    """
    n = len(y_sorted)
    if n < 4:
        return 0, np.inf, 0, "None"
        
    # We must skip the very extreme edges so we don't pick single-outlier boundaries
    min_size = int(n * 0.1) # 10% minimum group size
    
    # Calculate cumulative sums for speed (O(N) instead of O(N^2))
    cumsum_y = np.cumsum(y_sorted)
    cumsum_y2 = np.cumsum(y_sorted**2)
    
    best_sse = np.inf
    best_idx = -1
    best_mag = 0
    direction = "None"
    
    total_y = cumsum_y[-1]
    total_y2 = cumsum_y2[-1]
    
    for i in range(min_size, n - min_size):
        # Left group (0 to i-1)
        mean_L = cumsum_y[i-1] / i
        sse_L = cumsum_y2[i-1] - (cumsum_y[i-1]**2) / i
        
        # Right group (i to n-1)
        n_R = n - i
        sum_R = total_y - cumsum_y[i-1]
        sum_R2 = total_y2 - cumsum_y2[i-1]
        mean_R = sum_R / n_R
        sse_R = sum_R2 - (sum_R**2) / n_R
        
        total_sse = sse_L + sse_R
        
        if total_sse < best_sse:
            best_sse = total_sse
            best_idx = i
            best_mag = abs(mean_R - mean_L)
            direction = "Up" if mean_R > mean_L else "Down"
            
    return best_idx, best_sse, best_mag, direction

def run_age_stepminer():
    print("--- AGING PIPELINE: STEP 4 (STEPMINER AGE-SWITCHES) ---")
    
    PAIRS_FILE = "Verified_Age_Switch_Pairs.csv"
    ISO_FILE = "filtered_blood_isoforms_tpm.csv"
    AGE_FILE = "GTEx_Blood_Ages.csv"
    OUTPUT_FILE = "Age_Switch_Thresholds.csv"
    
    if not os.path.exists(PAIRS_FILE):
        print(f"Error: {PAIRS_FILE} not found. Please run Step 3 first.")
        return
        
    # 1. Load Golden Pairs
    print("\n[1/3] Loading 969 Verified Epigenetic-Splicing Pairs...")
    pairs = pd.read_csv(PAIRS_FILE)
    target_isoforms = pairs['Isoform'].unique()
    
    # 2. Load Ages and Calculate Midpoints
    print("\n[2/3] Loading and sorting Patient Ages...")
    ages_df = pd.read_csv(AGE_FILE)
    def get_midpoint(a):
        try:
            p = str(a).split('-')
            return (int(p[0]) + int(p[1])) / 2.0
        except: return np.nan
        
    ages_df['Age'] = ages_df['AgeBracket'].apply(get_midpoint)
    ages_df = ages_df.dropna().set_index('SampleID')
    
    # 3. Load Isoform TPM (Only for target isoforms)
    print(f"\n[3/3] Loading expression data for {len(target_isoforms)} selected isoforms...")
    iso_df = pd.read_csv(ISO_FILE, index_col=0)
    # Filter for only our gold list
    iso_df = iso_df[iso_df.index.isin(target_isoforms)]
    
    valid_samples = ages_df.index.intersection(iso_df.columns[1:])
    print(f"  Validating {len(valid_samples)} clinical patients.")
    
    # Sort samples strictly by Chronological Age
    sorted_ages = ages_df.loc[valid_samples].sort_values(by='Age')
    sorted_samples = sorted_ages.index.tolist()
    sorted_age_vector = sorted_ages['Age'].values
    
    # Align matrix to sorted patient ages
    iso_matrix = iso_df[sorted_samples].astype(np.float32)
    
    print("\n--- RUNNING STEPMINER ON SORTED AGES ---")
    results = []
    
    total = len(iso_matrix)
    for idx, (iso_id, row) in enumerate(iso_matrix.iterrows()):
        expr_vals = row.values
        
        # StepMiner strictly searches for the mathematical "cliff" where expression permanent shifts
        best_split_idx, min_sse, magnitude, direction = calculate_stepminer(expr_vals)
        
        switch_age = sorted_age_vector[best_split_idx]
        
        # Match back to the Gene Name from pairs
        gene_name = pairs[pairs['Isoform'] == iso_id]['GeneSymbol'].values[0]
        
        results.append({
            'GeneSymbol': gene_name,
            'Isoform': iso_id,
            'Switch_Age': switch_age,
            'Switch_Direction': direction,
            'Switch_Magnitude': magnitude
        })
        
        if (idx+1) % 100 == 0:
            print(f"  Processed {idx+1}/{total} isoforms...")
            
    final_df = pd.DataFrame(results)
    # Sort by magnitude of the biological shift
    final_df = final_df.sort_values(by='Switch_Magnitude', ascending=False)
    
    final_df.to_csv(OUTPUT_FILE, index=False)
    print(f"\nSUCCESS! Found precise universal Age-Switches for {len(final_df)} validated mechanisms.")
    print(f"Saved to '{OUTPUT_FILE}'\n")
    
    print("Top 5 Largest Age Switches:")
    print(final_df.head(5).to_string(index=False))

if __name__ == "__main__":
    run_age_stepminer()
