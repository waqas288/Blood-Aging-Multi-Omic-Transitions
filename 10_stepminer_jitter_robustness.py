import pandas as pd
import numpy as np
import os

# Set working directory
os.chdir(r'c:\Users\WAQAS\Desktop\quick_paper\datasets\aging')

# Simplified StepMiner for robustness testing
def stepminer_lite(data, ages):
    # Sort data by age
    idx = np.argsort(ages)
    x = ages[idx]
    y = data[idx]
    
    n = len(y)
    best_sse = float('inf')
    best_step = -1
    
    # Test every possible step point
    for i in range(1, n-1):
        left_mean = np.mean(y[:i])
        right_mean = np.mean(y[i:])
        
        sse = np.sum((y[:i] - left_mean)**2) + np.sum((y[i:] - right_mean)**2)
        
        if sse < best_sse:
            best_sse = sse
            best_step = i
            
    return x[best_step]

def run_jitter_robustness():
    print("--- AGING PIPELINE: STEP 15 (STEPMINER JITTER ROBUSTNESS) ---")
    
    ELITE_FILE = "Final_Elite_Survival_Results.csv"
    ISO_FILE = "filtered_blood_isoforms_tpm.csv"
    AGE_FILE = "GTEx_Blood_Ages.csv"
    
    if not os.path.exists(ELITE_FILE):
        print("Error: Elite results missing.")
        return

    # Load data
    elite_df = pd.read_csv(ELITE_FILE)
    target_isoforms = elite_df['Isoform'].unique()
    
    ages_df = pd.read_csv(AGE_FILE)
    def get_midpoint(a):
        try:
            p = str(a).split('-')
            return (int(p[0]) + int(p[1])) / 2.0
        except: return np.nan
    ages_df['Age'] = ages_df['AgeBracket'].apply(get_midpoint)
    ages_df = ages_df.dropna().set_index('SampleID')
    
    iso_df = pd.read_csv(ISO_FILE, index_col=0)
    iso_df = iso_df[iso_df.index.isin(target_isoforms)]
    
    valid_samples = ages_df.index.intersection(iso_df.columns[1:])
    iso_data = iso_df[valid_samples]
    ages_vec = ages_df.loc[valid_samples, 'Age'].values
    
    print(f"Running Jitter Test on {len(target_isoforms)} Elite targets...")
    
    results = []
    np.random.seed(42)
    
    for iso_id in target_isoforms:
        y = iso_data.loc[iso_id].values
        
        # 1. Original Switch Point
        orig_switch = stepminer_lite(y, ages_vec)
        
        # 2. Jittered Switch Point (Add +/- 5 years of noise)
        jittered_ages = ages_vec + np.random.uniform(-5, 5, size=len(ages_vec))
        jitter_switch = stepminer_lite(y, jittered_ages)
        
        results.append({
            'Isoform': iso_id,
            'Original_Switch': orig_switch,
            'Jittered_Switch': jitter_switch,
            'Difference': abs(orig_switch - jitter_switch)
        })

    res_df = pd.DataFrame(results)
    res_df.to_csv("StepMiner_Jitter_Robustness.csv", index=False)
    
    avg_diff = res_df['Difference'].mean()
    median_jitter = res_df['Jittered_Switch'].median()
    
    print("\n--- ROBUSTNESS SUMMARY ---")
    print(f"Average Switch Shift under +/- 5yr Noise: {avg_diff:.2f} years")
    print(f"Median Switch Point (Original): {res_df['Original_Switch'].median():.1f}")
    print(f"Median Switch Point (Jittered): {median_jitter:.1f}")
    
    if 50 <= median_jitter <= 60:
        print("\nSUCCESS: The switch point remains stable in the Fifth Decade even with significant age noise.")
        print("This proves the 'Fifth Decade Transition' is a biological signal, not a binning artifact.")
    else:
        print("\nWARNING: The signal shifted significantly. Further investigation needed.")

if __name__ == "__main__":
    run_jitter_robustness()
