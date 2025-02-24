# VEP Annotation Analysis Scripts Version Control
We created a pipeline, for analysing VEP annotated VCF files till three results txt files, for uploading to cloud and generate PDF/APP contents automatically.

# Major functions
### cell 1 
- **Path setting**: config of storing all pathes of input/outputs.
- **Input files/databases Handling**:  
  - VCF/GZ files: `clinvar_20240611.vcf.gz`
  - Tabular databases: `GeneDB_GenCC.txt`, `diseaseDB_1115_3.txt`,`genedb.ontology.all0307.csv`,`Merged_GWAS_vcf_2024.txt`,`Merged_Pharma_vcf_2024.txt`,`Reports_genome_databases_traits_merged_2.txt`,`Clingen-variants-2024-12-09.txt`,`pheno_OMIM_all.txt`

### cell 2
- **Whole genome variants annotation and Categorization**:
VCF Reading → ClinVar Updates if annotated by a earlier version → Multi-criteria Analysis → Variants filtration,Categorization → Initial Report Generation
**Five-Dimensional Classification**

| Flag | Criteria | Report |
|------|----------|--------|
| `saveFlag1` | Pathogenic/Likely Pathogenic, Conflicting Pathogenicity with P/LP & VUS submissions | A |
| `saveFlag2` | Low Frequency (MAF<5%) + Any in -sillico Predictive Models over the recommended threshold + Non-known benign from ClinVar | B | 
| `saveFlag3` | GWAS-associated (p <= 5e-8 and iGenotype == risk_genotype ) | C |
| `saveFlag4` | Pharmacogenomics associated variants | D |
| `saveFlag5` | Wellness Trait-related | E |

#### 2. **Runtime Parameter Management**
**Current Production Version**: `v4 (2025-02-20)`

## Version History

| Version | Date       | Key Changes                              |
|---------|------------|------------------------------------------|
| v3.1      | 2025-02-20 | Including v3 features, while adding APOE e2,e3,e4 special handling,and clinvar updates from py script and keep P/LP/CPLP based on all same version |
| v3.0      | 2025-01-10 | Fast version, running tim around 90 sec.nodup4 file refinement, known benign variants removed, adding new db integration (PanelApp, GenCC, ClinGen, geneBe), removing unneeded middle file generation process|
| v2.1      | 2024-01-30 | Fixing bugs, such as PGx variants in parallel with risk variants (rs6025, which should be  a famous PGx genes for F5 gene, and tratis section keep homozygous reference genotype for patients ; vep.gz file used ; Hemizygous on X chromosome for male/female distinguishment         |
| v2.0      | 2023-09-10 | Second established version of a script for analyzing vep annotated file, giving scores, detection of both known phenotypic-associated variants and predicted variants; multiple gene panels, as well as completed GWAS,PGx,nodup4 file and traits outputs.         |
| v1.0      | 2022-09-10 | Initial established version of a script for detecting ClinVar Pathogenic/likely pathogenic variants only from a non-annotated VCF single sample file. The variants-genes-disease association is used as input database.    |

### Key Improvements on v3.1
- **Muchfaster processing** compared to v1/v2: **15 min** to **1.5 min**.
- **Multiple new databases intergrated** Added gene curation db: PanelApp, and GenCC; Added variant curation db: ClinGen, geneBe (waiting for server permission fixing), disease categories.
- **Review Star**: over or equal to 2 should be considered as solid known findings.
- **pipeline outputs**:
  ```bash
  output/
  ├── _sp_Inheritance_4_nodup.txt  # Known pathogenic and predicted risk variants
  ├── _GWAS.txt      # GWAS-specific findings (single SNP matching)
  └── _pharmaco.txt      # PGx specific findings (single SNP matching, haplotype waiting for improvements)
  ```

## Environment Requirements


## Usage Example
```bash
# Run analysis pipeline
python pythonpipeline.py {input_vep_annotated_file} {output_file_prefix}
```
