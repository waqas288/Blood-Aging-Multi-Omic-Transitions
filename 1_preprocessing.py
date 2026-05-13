import pandas as pd
import gc
import os

print("--- AGING MULTI-OMICS PIPELINE: STEP 1 (PREPROCESSING) ---")

# -------------------------------------------------------------
# 1. GTEx PREPROCESSING (Isoform level, Memory Efficient)
# -------------------------------------------------------------
print("\n[GTEx] Loading Clinical Metadata...")

# A. Load Annotations
sample_df = pd.read_csv("C:/Users/WAQAS/Desktop/quick_paper/datasets/aging/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt", sep='\t', low_memory=False)
blood_samples = sample_df[sample_df['SMTSD'] == 'Whole Blood']['SAMPID'].tolist()
print(f"[GTEx] Found {len(blood_samples)} 'Whole Blood' samples.")

pheno_df = pd.read_csv("C:/Users/WAQAS/Desktop/quick_paper/datasets/aging/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt", sep='\t')
# Map GTEX-1117F -> "60-69"
age_map = dict(zip(pheno_df['SUBJID'], pheno_df['AGE']))

blood_sample_ages = {}
for samp_id in blood_samples:
    subj_id = "-".join(samp_id.split("-")[:2])
    if subj_id in age_map:
        blood_sample_ages[samp_id] = age_map[subj_id]

valid_blood_samples = set(blood_sample_ages.keys())
print(f"[GTEx] Mapped Age brackets for {len(valid_blood_samples)} GTEx Blood samples.")

# Save mapping to file
pd.DataFrame(list(blood_sample_ages.items()), columns=['SampleID', 'AgeBracket']).to_csv("GTEx_Blood_Ages.csv", index=False)

print("\n[GTEx] Processing 11GB GTEx Isoform Matrix (Line-by-line)...")
gct_file = "C:/Users/WAQAS/Desktop/quick_paper/datasets/aging/GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_transcript_tpm.gct"
out_file = "filtered_blood_isoforms_tpm.csv"

if not os.path.exists(out_file):
    with open(gct_file, 'r') as f_in, open(out_file, 'w') as f_out:
        f_in.readline(); f_in.readline() # Skip headers
        header = f_in.readline().strip('\n').split('\t')
        
        col_indices = [0, 1] # transcript_id, gene_id
        col_names = ["transcript_id", "gene_id"]
        
        for i in range(2, len(header)):
            if header[i] in valid_blood_samples:
                col_indices.append(i)
                col_names.append(header[i])
        
        f_out.write(",".join(col_names) + "\n")
        print(f"[GTEx Isoforms] Extracting {len(col_indices) - 2} matched columns...")
        
        line_count = 0
        for line in f_in:
            parts = line.strip('\n').split('\t')
            filtered_parts = [parts[i] for i in col_indices]
            f_out.write(",".join(filtered_parts) + "\n")
            line_count += 1
            if line_count % 50000 == 0:
                print(f"  ...processed {line_count} transcripts")
else:
    print(f"[GTEx] {out_file} already exists. Skipping extraction.")


print("\n[GTEx] Processing GTEx Whole-Gene Matrix (Line-by-line)...")
gene_gct_file = "mGTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct/mGTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct"
# wait, the directory is mGTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct but the file is usually GTEx_Analysis...
# To be safe, let's just let the user know we'll generate the scripts first.
# Wait, I can just use os.listdir to find the exact file inside.
import glob
gene_gct_files = glob.glob("mGTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct/*.gct")
if gene_gct_files:
    gene_gct_file = gene_gct_files[0]
    gene_out_file = "filtered_blood_genes_tpm.csv"
    
    if not os.path.exists(gene_out_file):
        with open(gene_gct_file, 'r') as f_in, open(gene_out_file, 'w') as f_out:
            f_in.readline(); f_in.readline()
            header = f_in.readline().strip('\n').split('\t')
            
            col_indices = [0, 1] # Name, Description
            col_names = ["gene_id", "gene_name"]
            
            for i in range(2, len(header)):
                if header[i] in valid_blood_samples:
                    col_indices.append(i)
                    col_names.append(header[i])
            
            f_out.write(",".join(col_names) + "\n")
            print(f"[GTEx Genes] Extracting {len(col_indices) - 2} matched columns. Writing to {gene_out_file}...")
            
            line_count = 0
            for line in f_in:
                parts = line.strip('\n').split('\t')
                filtered_parts = [parts[i] for i in col_indices]
                f_out.write(",".join(filtered_parts) + "\n")
                line_count += 1
                if line_count % 10000 == 0:
                    print(f"  ...processed {line_count} genes")
    else:
        print(f"[GTEx] {gene_out_file} already exists. Skipping extraction.")
else:
    print("[GTEx] Failed to locate the gene.gct file inside the mGTEx folder.")
# -------------------------------------------------------------
# 2. GEO PREPROCESSING (Methylation Metadata)
# -------------------------------------------------------------
print("\n[GEO] Checking GSE87571 Methylation Metadata...")

geo_metadata = "GSE87571_series_matrix.txt"
geo_sample_ids = []
geo_ages = []

if os.path.exists(geo_metadata):
    with open(geo_metadata, 'r', encoding='utf-8', errors='ignore') as f:
        for line in f:
            if line.startswith('!Sample_title'):
                # Extract clean sample ID: e.g. "X1 genomic DNA from whole blood" -> "X1"
                titles = line.strip().split('\t')[1:]
                geo_sample_ids = [t.replace('"', '').split(' ')[0] for t in titles]
            
            elif line.startswith('!Sample_characteristics_ch1') and 'age' in line.lower():
                # Extract age value: e.g. "age: 44" or "Age (yrs): 44"
                chars = line.strip().split('\t')[1:]
                geo_ages = [c.split(':')[1].strip().replace('"', '') for c in chars if ':' in c]
    
    if len(geo_sample_ids) > 0 and len(geo_ages) > 0:
        # Save mapping to file
        df_geo = pd.DataFrame({'SampleID': geo_sample_ids, 'Age': geo_ages})
        df_geo.to_csv("GEO_GSE87571_Ages.csv", index=False)
        print(f"[GEO] Saved Ages for {len(df_geo)} GSE87571 Methylation samples.")
    else:
        print("[GEO] Warning: Could not automatically locate 'age:' row in the series matrix.")
else:
    print(f"[GEO] File {geo_metadata} not found.")

print("\n--- PREPROCESSING SCRIPT FINISHED ---")