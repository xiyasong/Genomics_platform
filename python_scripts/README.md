# VEP Annotation Analysis Scripts Version Control

**Current Production Version**: `v4 (2025-02-20)`

## Version History

| Version | Date       | Key Changes                              |
|---------|------------|------------------------------------------|
| v4      | 2024-02-20 | Including v3 features, while adding APOE e2,e3,e4 special handling,and clinvar updates from py script and keep P/LP/CPLP based on all same version |
| v3      | 2024-01-10 | nodup4 file refinement, known benign variants removed, adding new db integration (PanelApp, GenCC, ClinGen, geneBe), removing unneeded middle file generation process|
| v2      | 2023-11-30 | Fixing bugs, such as PGx variants in parallel with risk variants, and tratis section keep homozygous reference genotype for patients ; vep.gz file used ; Hemizygous on X chromosome for male/female distinguishment         |
| v1      | 2023-09-10 | Initial fixed version of a script for analyzing vep annotated file, giving scores, detection of both known phenotypic-associated variants and predicted variants; multiple gene panels, as well as completed GWAS,PGx,nodup4 file and traits outputs.         |

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
