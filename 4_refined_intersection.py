import pandas as pd
import os

def run_refined_intersection():
    print("--- AGING PIPELINE: STEP 7c (REFINED INTERSECTION) ---")
    
    CPG_FILE = "Significant_Aging_CpGs_Refined.csv"
    ISO_FILE = "Significant_Aging_Isoforms_Refined.csv"
    ANNO_FILE = "illumina_450k_annotation.csv.gz"
    OUTPUT_FILE = "Verified_Age_Switch_Pairs_Refined.csv"
    
    if not (os.path.exists(CPG_FILE) and os.path.exists(ISO_FILE)):
        print("Missing refined discovery files. Run 2b and 6b first.")
        return

    print("[1/3] Loading Refined Datasets...")
    cpg_df = pd.read_csv(CPG_FILE)
    iso_df = pd.read_csv(ISO_FILE)
    
    # Map CpGs to Genes
    print("[2/3] Mapping CpGs to Genes...")
    anno = pd.read_csv(ANNO_FILE, compression='gzip', skiprows=7, usecols=['IlmnID', 'UCSC_RefGene_Name'])
    anno = anno[anno['IlmnID'].isin(cpg_df['CpG'])]
    
    # Split Gene Names (some probes map to multiple genes like 'GeneA;GeneB')
    anno['GeneSymbol'] = anno['UCSC_RefGene_Name'].str.split(';')
    anno = anno.explode('GeneSymbol').dropna(subset=['GeneSymbol'])
    
    cpg_mapped = cpg_df.merge(anno, left_on='CpG', right_on='IlmnID').drop(columns=['IlmnID', 'UCSC_RefGene_Name'])
    
    print("[3/3] Performing Multi-Omics Intersection...")
    # Join on GeneSymbol
    # Note: iso_df has 'Isoform', 'Gene' (Ensembl ID), 'Pearson_R'
    # We need to map Ensembl IDs to symbols or use the 'GeneSymbol' if it exists.
    # Looking at Step 6, 'Gene' is actually the Ensembl ID.
    
    # Wait, the original Step 7 used the MyGene mapping.
    # Let's ensure we have Gene Symbols in the Isoform DF.
    # If not, we fetch them.
    
    import mygene
    mg = mygene.MyGeneInfo()
    ensembl_ids = iso_df['Gene'].unique().tolist()
    # Strip versions
    clean_ids = [str(x).split('.')[0] for x in ensembl_ids]
    results = mg.querymany(clean_ids, scopes='ensembl.gene', fields='symbol', species='human')
    
    symbol_map = {res['query']: res.get('symbol') for res in results if 'symbol' in res}
    iso_df['GeneSymbol'] = iso_df['Gene'].apply(lambda x: symbol_map.get(str(x).split('.')[0]))
    iso_df = iso_df.dropna(subset=['GeneSymbol'])
    
    # Final Join
    merged = iso_df.merge(cpg_mapped, on='GeneSymbol', suffixes=('_Isoform', '_CpG'))
    
    merged.to_csv(OUTPUT_FILE, index=False)
    print(f"SUCCESS! Identified {len(merged)} Refined Regulatory Pairs across {merged['GeneSymbol'].nunique()} genes.")
    print(f"Results saved to '{OUTPUT_FILE}'")

if __name__ == "__main__":
    run_refined_intersection()
