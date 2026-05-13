# Multi-Omic Aging & Splicing Clock Pipeline

This repository contains the source code for the "Multi-Omic Aging & Splicing Clock" project. It includes a comprehensive pipeline for analyzing isoform-level transcriptomic data, identifying aging-associated signals, and building predictive machine learning models (Splicing Clocks).

## Repository Overview

The core of the analysis is broken down into modular scripts (1-16) and a unified Master Pipeline for production runs.

### Core Scripts
1.  **Preprocessing**: Data cleaning and normalization.
2.  **Correlation Analysis**: Identifying methylation and isoform correlations with age.
3.  **Survival Modeling**: Running Cox Proportional Hazards models adjusted for cell types, RIN, and sex.
4.  **StepMiner Analysis**: Detecting "switches" in isoform expression across age.
5.  **Validation**: Independent validation using GTEx and other datasets.
6.  **Machine Learning**: Building and evaluating the Splicing Clock using Elastic Net.

### Master Pipeline
The `MASTER_AGING_PIPELINE.py` script aggregates the discovery, survival, and ML stages into a single execution flow for reproducibility.

## How to Run

### Prerequisites
- Python 3.8+
- Required Packages: `pandas`, `numpy`, `scikit-learn`, `statsmodels`, `scipy`, `matplotlib`, `seaborn`

### Execution
To run the full production pipeline:
```bash
python MASTER_AGING_PIPELINE.py
```
Outputs (results and figures) will be saved in the `production_results/` folder.

## Data Availability
This pipeline is designed to work with GTEx and multi-tissue transcriptomic data. Metadata requirements include:
- `GTEx_Blood_Ages.csv`
- `filtered_blood_isoforms_tpm.csv`
- Subject phenotypes and sample attributes.

## License
This project is licensed under the **MIT License** - see the [LICENSE](LICENSE) file for details.

## Citation
If you use this code in your research, please cite:
*(Insert Paper Citation Here)*
