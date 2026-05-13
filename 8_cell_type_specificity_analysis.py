import pandas as pd
import requests
import json
import time

def run_cell_type_analysis():
    print("--- AGING PIPELINE: STEP 12 (CELL-TYPE SPECIFICITY ANALYSIS) ---")
    
    INPUT_FILE = "Final_Elite_Survival_Results.csv"
    if not os.path.exists(INPUT_FILE):
        print("Error: Input file missing.")
        return

    df = pd.read_csv(INPUT_FILE)
    genes = df['GeneSymbol'].unique().tolist()
    print(f"Analyzing {len(genes)} Elite genes for cell-type specificity...")

    # Enrichr API setup
    ENRICHR_URL = 'http://maayanlab.cloud/Enrichr/addList'
    query_string = '\n'.join(genes)
    description = 'Aging Elite Mortality Predictors'
    payload = {
        'list': (None, query_string),
        'description': (None, description)
    }

    response = requests.post(ENRICHR_URL, files=payload)
    if not response.ok:
        print("Enrichr API Error.")
        return

    user_list_id = json.loads(response.text)['userListId']
    
    # We will check two libraries: 
    # 1. Human_Gene_Atlas (Classic tissue/cell expression)
    # 2. PanglaoDB_Augmented_2021 (Single-cell markers)
    
    libraries = ['Human_Gene_Atlas', 'PanglaoDB_Augmented_2021']
    results_data = []

    for lib in libraries:
        print(f"  Querying {lib}...")
        results_url = f'http://maayanlab.cloud/Enrichr/enrich?userListId={user_list_id}&backgroundType={lib}'
        resp = requests.get(results_url)
        if resp.ok:
            data = json.loads(resp.text)[lib]
            # data is a list of [rank, term, p-value, z-score, combined_score, overlapping_genes, adj_p-value, etc.]
            for item in data[:10]: # Top 10 results
                results_data.append({
                    'Library': lib,
                    'Cell_Type_Term': item[1],
                    'P_value': item[2],
                    'Combined_Score': item[4],
                    'Overlapping_Genes': ", ".join(item[5]),
                    'FDR_p_value': item[6]
                })
        time.sleep(1) # Be nice to the API

    res_df = pd.DataFrame(results_data)
    res_df.to_csv("Cell_Type_Specificity_Results.csv", index=False)
    
    print("\n--- CELL-TYPE ANALYSIS SUMMARY ---")
    # Filter for significant blood cells
    blood_keywords = ['Blood', 'Leukocyte', 'Neutrophil', 'Lymphocyte', 'Monocyte', 'B-cell', 'T-cell', 'NK']
    blood_hits = res_df[res_df['Cell_Type_Term'].str.contains('|'.join(blood_keywords), case=False)]
    
    if len(blood_hits) > 0:
        print(f"Found {len(blood_hits)} significant blood-related cell types.")
        print(blood_hits[['Cell_Type_Term', 'FDR_p_value', 'Overlapping_Genes']].head(5).to_string(index=False))
    else:
        print("No dominant blood cell type found. This suggests these genes are broadly expressed across the immune system.")
    
    print(f"\nFull results saved to 'Cell_Type_Specificity_Results.csv'")

if __name__ == "__main__":
    import os
    run_cell_type_analysis()
