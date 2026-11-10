# BLOOMPARK
Official code repository for the analysis associated with the paper titled “Integrative Blood Transcriptomics Identifies U1-snRNP Gene Repression and miRNA–mRNA Regulatory Hubs in Pre-clinical Parkinson’s Disease.”

Paper aviaible at https://doi.org/10.1038/s41531-026-01526-1

## Repository contents

| File | Description |
|---|---|
| `README.md` | Project description and link to the associated article (doi: 10.1038/s41531-026-01526-1) |
| `BLOOMPARK_Pipeline.R` | Full R pipeline, from raw count assembly to miRNA–mRNA correlation analysis |
| `Input_Data_Directory_Structure.txt` | Expected layout of the input data directory |

The input data (PPMI featureCounts and miRNA count files, the curated PPMI clinical spreadsheet, and the blood chemistry/hematology table) are not included in this repository. The script uses absolute paths and must be adapted to the local environment before running.

## Analysis summary

The pipeline analyses peripheral blood transcriptomes from the Parkinson's Disease Progression Markers Initiative (PPMI) cohort to find genes and miRNAs associated with disease progression, and to test whether differentially expressed miRNAs correlate with their predicted target genes.

**1. Count matrix assembly.** RNA-seq counts (featureCounts) and miRNA-seq counts are read per sample and indexed by `PATNO_EVENT_ID`.

**2. Metadata and blood chemistry.** Samples are assigned to four severity stages based on the NSD stage: Pre-clinical (NSD 1a, 1b, 2a, and the prodromal cohort), Mild (2b), Moderate (3), and Severe (4–6). Healthy controls form the reference group. White blood cell count, platelets, neutrophils, lymphocytes, and BMI are imputed hierarchically (within-patient mean, then within-stage mean, then global mean). The neutrophil-to-lymphocyte ratio (NLR) is derived from the neutrophil and lymphocyte counts. For healthy controls, missing levodopa equivalent daily dose (LEDD) is set to 0.

**3. Sample harmonization.** Only samples present in the RNA-seq data, the miRNA data, and the metadata are retained.

**4. Gene annotation.** Gene identifiers are mapped to symbols and annotated.

**5. Differential expression.** Each stage is compared with Healthy controls using edgeR. The design model is `~ disease + age + SEX + BMI + WBC + Platelets + NLR + LEDD`, which adjusts for age, sex, BMI, blood counts, NLR, and medication. Genes and miRNAs are considered significant at adjusted p < 0.05 and |log2FC| ≥ 0.58 (about 1.5-fold). The outputs include DEG and DEmiRNA tables per stage, volcano plots, MDS plots, and an UpSet plot showing the overlap of DEGs across stages.

**6. miRNA targets and integration.** Predicted and validated targets of the DEmiRNAs are retrieved with MIRit. For each stage, the miRNA and gene data are batch-corrected (patient and sex as batch factors, with covariates and surrogate variables) and each DEmiRNA is tested for Spearman correlation with its target genes. The Moderate stage is run serially to limit memory use; the others run in parallel on 16 cores.

Summary counts (sample numbers, number of DEGs and DEmiRNAs per stage, number of target pairs, and correlation pairs) are written to `General_statistics.txt` during the run.
