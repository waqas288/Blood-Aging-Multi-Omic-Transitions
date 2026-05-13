import pandas as pd
import mygene
import os

# Set working directory
os.chdir(r'c:\Users\WAQAS\Desktop\quick_paper\datasets\aging')

def run_refined_proximity():
    print("--- AGING PIPELINE: STEP 7d (REFINED PROXIMITY CHECK) ---")
    
    INPUT_FILE = "Verified_Age_Switch_Pairs_Refined.csv"
    ANNO_FILE = "illumina_450k_annotation.csv.gz"
    
    if not os.path.exists(INPUT_FILE):
        print("Error: Input missing. Run 7c first.")
        return

    df_pairs = pd.read_csv(INPUT_FILE)
    unique_isoforms = df_pairs['Isoform'].unique().tolist()

    print(f"Fetching coordinates for {len(unique_isoforms)} unique isoforms...")
    mg = mygene.MyGeneInfo()

    def get_transcript_coords(isoform_list):
        clean_ids = [i.split('.')[0] for i in isoform_list]
        results = mg.querymany(clean_ids, scopes='ensembl.transcript', fields='genomic_pos', species='human')
        coord_map = {}
        for res in results:
            query_id = res['query']
            if 'genomic_pos' in res:
                pos = res['genomic_pos']
                if isinstance(pos, list): pos = pos[0]
                coord_map[query_id] = {
                    'chr': str(pos.get('chr')),
                    'start': int(pos.get('start')),
                    'end': int(pos.get('end'))
                }
        return coord_map

    isoform_coords = get_transcript_coords(unique_isoforms)

    print("Loading CpG annotations...")
    cpg_anno = pd.read_csv(ANNO_FILE, compression='gzip', skiprows=7, usecols=['IlmnID', 'CHR', 'MAPINFO'])
    cpg_anno = cpg_anno[cpg_anno['IlmnID'].isin(df_pairs['CpG'].unique())]

    print("Calculating proximity...")
    results = []
    for idx, row in df_pairs.iterrows():
        iso_id = row['Isoform'].split('.')[0]
        cpg_id = row['CpG']
        
        if iso_id in isoform_coords:
            iso_pos = isoform_coords[iso_id]
            cpg_info = cpg_anno[cpg_anno['IlmnID'] == cpg_id]
            
            if not cpg_info.empty:
                c_chr = str(cpg_info.iloc[0]['CHR'])
                c_pos = int(cpg_info.iloc[0]['MAPINFO'])
                
                if c_chr != iso_pos['chr']:
                    dist = float('inf')
                else:
                    dist = min(abs(c_pos - iso_pos['start']), abs(c_pos - iso_pos['end']))
                    if iso_pos['start'] <= c_pos <= iso_pos['end']: dist = 0
                
                row_dict = row.to_dict()
                row_dict['Distance_bp'] = dist
                results.append(row_dict)

    df_results = pd.DataFrame(results)
    proximity_threshold = 5000

    high_conf = df_results[df_results['Distance_bp'] <= proximity_threshold]
    distal = df_results[df_results['Distance_bp'] > proximity_threshold]

    print(f"High-Confidence Pairs (<=5kb): {len(high_conf)}")
    print(f"Distal Regulatory Pairs (>5kb): {len(distal)}")

    high_conf.to_csv('High_Confidence_Proximity_Pairs_Refined.csv', index=False)
    distal.to_csv('Distal_Regulatory_Pairs_Refined.csv', index=False)
    print("Files saved successfully.")

if __name__ == "__main__":
    run_refined_proximity()
