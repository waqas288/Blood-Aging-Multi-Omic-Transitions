import pandas as pd
import requests
import json
import os
import time

def run_gtex_validation():
    print("--- AGING PIPELINE: STEP 14 (GTEx PORTAL sQTL VALIDATION) ---")
    
    INPUT_FILE = "Final_Elite_Survival_Results.csv"
    HIGH_CONF_FILE = "High_Confidence_Proximity_Pairs_Refined.csv"
    
    if not (os.path.exists(INPUT_FILE) and os.path.exists(HIGH_CONF_FILE)):
        print("Error: Input files missing.")
        return

    # Load genes and their Ensembl IDs
    df_elite = pd.read_csv(INPUT_FILE)
    df_proximity = pd.read_csv(HIGH_CONF_FILE)
    
    # Merge to get ENSG IDs for the elite symbols
    elite_targets = df_elite.merge(df_proximity[['GeneSymbol', 'Gene']].drop_duplicates(), on='GeneSymbol', how='left')
    
    genes_to_test = elite_targets[['GeneSymbol', 'Gene']].drop_duplicates().dropna()
    print(f"Querying GTEx Portal for {len(genes_to_test)} Elite Genes...")

    # GTEx API v2 Endpoint for sQTLs
    # https://gtexportal.org/api/v2/association/singleTissueSqtl?gencodeId=ENSG00000103194.15&tissueSiteDetailId=Whole_Blood
    
    GTEX_API_URL = "https://gtexportal.org/api/v2/association/singleTissueSqtl"
    TISSUE = "Whole_Blood"
    
    validation_results = []

    for idx, row in genes_to_test.iterrows():
        symbol = row['GeneSymbol']
        ensg = row['Gene']
        
        print(f"  Checking {symbol} ({ensg})...")
        
        params = {
            'gencodeId': ensg,
            'tissueSiteDetailId': TISSUE,
            'datasetId': 'gtex_v8'
        }
        
        try:
            response = requests.get(GTEX_API_URL, params=params)
            if response.ok:
                data = response.json()
                # data is usually a list of sQTL associations
                if 'data' in data and len(data['data']) > 0:
                    top_sqtl = data['data'][0] # Get the strongest association
                    validation_results.append({
                        'GeneSymbol': symbol,
                        'GTEx_Validated': 'YES',
                        'Top_sQTL_Variant': top_sqtl.get('variantId'),
                        'GTEx_p_value': top_sqtl.get('pValue'),
                        'GTEx_NES': top_sqtl.get('nes') # Normalized Effect Size
                    })
                else:
                    validation_results.append({
                        'GeneSymbol': symbol,
                        'GTEx_Validated': 'NO',
                        'Top_sQTL_Variant': None,
                        'GTEx_p_value': None,
                        'GTEx_NES': None
                    })
            else:
                print(f"    API Error for {symbol}")
        except Exception as e:
            print(f"    Connection error for {symbol}: {str(e)}")
        
        time.sleep(0.5) # Be gentle with the GTEx API

    val_df = pd.DataFrame(validation_results)
    val_df.to_csv("GTEx_sQTL_Validation_Results.csv", index=False)
    
    # Summary
    validated_count = len(val_df[val_df['GTEx_Validated'] == 'YES'])
    print(f"\n--- VALIDATION SUMMARY ---")
    print(f"Total Genes Tested: {len(val_df)}")
    print(f"GTEx-Validated sQTLs: {validated_count} ({int(validated_count/len(val_df)*100)}%)")
    
    if validated_count > 0:
        print("\nTop Validated Genes:")
        print(val_df[val_df['GTEx_Validated'] == 'YES'].head(10).to_string(index=False))
    
    print(f"\nFull results saved to 'GTEx_sQTL_Validation_Results.csv'")

if __name__ == "__main__":
    run_gtex_validation()
