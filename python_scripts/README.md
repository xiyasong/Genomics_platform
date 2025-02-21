# VEP Annotation Analysis Scripts Version Control

**Current Production Version**: `v4 (2025-02-20)`

## Version History

| Version | Date       | Key Changes                              |
|---------|------------|------------------------------------------|
| v4      | 2025-02-20 | Including v3 features, while adding APOE e2,e3,e4 special handling,and clinvar updates from py script and keep P/LP/CPLP based on all same version |
| v3      | 2025-01-10 | Fast version, running tim around 90 sec.nodup4 file refinement, known benign variants removed, adding new db integration (PanelApp, GenCC, ClinGen, geneBe), removing unneeded middle file generation process|
| v2.1      | 2024-01-30 | Fixing bugs, such as PGx variants in parallel with risk variants (rs6025, which should be  a famous PGx genes for F5 gene, and tratis section keep homozygous reference genotype for patients ; vep.gz file used ; Hemizygous on X chromosome for male/female distinguishment         |
| v2      | 2023-09-10 | Second established version of a script for analyzing vep annotated file, giving scores, detection of both known phenotypic-associated variants and predicted variants; multiple gene panels, as well as completed GWAS,PGx,nodup4 file and traits outputs.         |
| v1      | 2022-09-10 | Initial established version of a script for detecting ClinVar Pathogenic/likely pathogenic variants only from a non-annotated VCF single sample file. The variants-genes-disease association is used as input database.    |

### Key Improvements on v4
- 🚀 **Muchfaster processing** compared to v1/v2: **15 min** to **1.5 min**.
- 🛡️ **Multiple new databases intergrated** Added gene curation db: PanelApp, and GenCC; Added variant curation db: ClinGen, geneBe (waiting for server permission fixing), disease categories.
- ✨ **Review Star**: over or equal to 2 should be considered as solid known findings.
- 📊 **pipeline outputs**:
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
