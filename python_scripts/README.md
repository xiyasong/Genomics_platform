# VEP Annotation Analysis Scripts Version Control
We created a pipeline, for analysing VEP annotated VCF files till three results txt files, for uploading to cloud and generate PDF/APP contents automatically.

**Current Production Version**: `v3.1 (2025-02-20)`

## Environment Requirements
pandas>=1.3.4 
genebe 
paramiko>=2.11.0 
scp>=0.14.4 
numpy>=1.21.2 
python-dateutil>=2.8 
csv>=1.0 
re>=2.2.1 
gzip 
datetime
## Usage Example
```bash
# Run analysis pipeline
python pythonpipeline.py {input_vep_annotated_file} {output_file_prefix}
```
### Key Improvements on v 3.1
- **Much faster processing** compared to v1/v2: **15 min** to **1.5 min**.
- **Multiple new databases intergrated** Added gene curation db: PanelApp, and GenCC; Added variant curation db: ClinGen, geneBe (waiting for server permission fixing), disease categories.
- **Review Star**: over or equal to 2 should be considered as solid known findings.
- **pipeline outputs**:
  ```bash
  output/
  ├── _sp_Inheritance_4_nodup.txt  # Known pathogenic and predicted risk variants
  ├── _GWAS.txt      # GWAS-specific findings (single SNP matching)
  └── _pharmaco.txt      # PGx specific findings (single SNP matching, haplotype waiting for improvements)
  ```
---
## Version History

| Version | Date       | Key Changes                              |
|---------|------------|------------------------------------------|
| v3.1      | 2025-02-24 | Including v3 features, while 1) adding APOE e2,e3,e4 haplotype special handling for AD; 2) ClinVar updates from py script and keep P/LP/CPLP based on all same/newest version |
| v3.0      | 2025-01-10 | Fast version, running tim around 90 sec. The nodup4 file refinement, known benign variants removed, adding new db integration (PanelApp, GenCC, ClinGen, geneBe); removing unneeded middle file generation |
| v2.1      | 2024-01-30 | Fixing bugs, such as PGx variants in parallel with risk variants (rs6025, which should be  a famous PGx genes for F5 gene, and tratis section keep homozygous reference genotype for patients ; vep.gz file used ; Hemizygous on X chromosome for male/female distinguishment         |
| v2.0      | 2023-09-10 | Second established version of a script for analyzing vep annotated file, giving scores, detection of both known phenotypic-associated variants and predicted variants; multiple gene panels, as well as completed GWAS,PGx,nodup4 file and traits outputs.         |
| v1.0      | 2022-09-10 | Initial established version of a script for detecting ClinVar Pathogenic/likely pathogenic variants only from a non-annotated VCF single sample file. The variants-genes-disease association is used as input database.    |

---
## Major functions

### Cell 1 
- **Path setting**: config of storing all pathes of input/outputs.
- **Input files/databases Handling**:  
  - VCF/GZ files: `clinvar_20240611.vcf.gz`
  - Tabular databases: `GeneDB_GenCC.txt`, `diseaseDB_1115_3.txt`,`genedb.ontology.all0307.csv`,`Merged_GWAS_vcf_2024.txt`,`Merged_Pharma_vcf_2024.txt`,`Reports_genome_databases_traits_merged_2.txt`,`Clingen-variants-2024-12-09.txt`,`pheno_OMIM_all.txt`

---
### Cell 2
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

---
### Cell 3
To prepare of mapping variants from ClinVar and predicted parts to genes, then mapping genes to diseases.

- **Gene Database Initialization**
- Parses TSV from `geneBaseFile`
- Populates core lists: 
  - `ConfidenceLevel`, `TargetGroup`, `geneBasedRef`, `Disease`, `SZAdiseaseID_GDB`

- **Disease Database Loading**
- Dual storage strategy:
  - Raw lists: `DiseaseID_DSDB`, `DiseaseName_DSDB`
  - Lookup dictionary(for speed up): `disease_lookup[standardized_name] = (id, name, index)` 
- Implements case-insensitive search prep via disease names standardization

---
### Cell 4
Core Functionality: ClinVar Pathogenic Variant Reporting

#### 1. **Data Preparation**
- Extracts ClinVar annotations from VCF INFO/CSQ field
- Handles genotype parsing with sex-aware logic:
  ```python
  gpf.parse_genotype(..., check_male=check_male_flag)
  ```
#### 2. **Disease Matching**
- Standardizes disease names using `disease_lookup`
- Implements two-level matching:
  1. Direct ID match: `SZAdiseaseID_GDB`
  2. Semantic match: Word-set comparison between ClinVar and GeneDB terms

#### 3. **Integer fixed Scoring**
- `scoreFlag` values: variant-centric
  - 5: Pathogenic
  - 4: Likely Pathogenic
  - 2: Conflicting interpretations with P & VUS
  - 1: Conflicting interpretations with LP & VUS

- `scoreFlag` values: gene-centric
 - ConfidenceLevel == 'High': 10 (Genes including in major screening projects)
 - ConfidenceLevel == 'Moderate': 5 (Genes defined by ClinVar that is casual genes for certain diseases)
 - ConfidenceLevel == 'Low': 0 (Genes defined by ClinVar that is just one of associated genes for certain diseases)

#### 4. **Three-Tier Case Handling**
| Case | Gene/Disease Matching Condition | What's being kept for records |
|------|-----------|-----------|
| 1 | Gene in GeneDB + Disease match | Full scoring + GeneDB metadata |
| 2 | Novel gene + Disease match | Low-confidence scoring |
| 3 | No gene/disease match | Basic flagging |

---
### Cell 8 
Core Functionality: Variants that are from samples, without known pathogenicity neither known benign (VUS) Risk Assessment

#### 1. **Basic Variant Filter**
- MAF < 5% (`MAX_AF` check)
- Excludes ClinVar benign/likely_benign variants

#### 2. **Multi-in-sillico scores filter**
- High Impact (VEP)
- LoF (Loftee) == 'HC'
- ada_score>0.6
- rf_score >0.6
- REVEL>0.75
- Any 4 of the Splice_DS_score >0.5
- BayesDel_addAF_score >0.0692655
- BayesDel_noAF_score > -0.0570105
- Alphamissense am_classs == likely_pathogenic or am_pathogenicity>0.564

#### 3. **Gene-Disease Mapping**
- Matches genes to GeneDB entries with that gene's associated diseases
- Calculates composite score:
  ```
  SZAscore = Base(2) + Impact(1 if HIGH) + ConfidenceLevel(0/5/10)
  ```
- `jSZADiseaseID`: Mapped disease identifier
- `predicted_impact`: Boolean aggregation of prediction models
-  Output Generation

#### Score Flag Values
| Confidence Level | Base Score | Final Range |
|------------------|------------|-------------|
| High             | +10        | 12-13       |
| Moderate         | +5         | 7-8         | 
| Low              | +0         | 2-3         |

---
### Cell 9 

- **CSQ Field Expansion**: Splits VEP annotations into discrete columns
- **Inheritance Mapping**: Links variants to OMIM inheritance patterns
- **Disease Context Enrichment**: Combines disease descriptions from multiple sources,adding diseases ontology information
- **Deduplication**: `drop_duplicates()` on 5 key columns ["Target.group", "Disease", "Genes", "SZAID", "SZAreportCategory"]
- **Other Clinical Evidence databases Integration**
`gnb.annotate() # ACMG classification from geneBe
pd.merge(clingen) # ClinGen curated variants evidence
pd.merge(genecc) # GenCC gene-diseases confidence level`
- **Output Optimization**
 - Generates HGVS standardized names ['HGVS_Naming']
 - Adds ClinVar review star ratings
 - Simplifies gene lists (first gene only)

---
### Cell 5 
**GWAS Filter & Reporting Pipeline**
VCF Input → Genotype Parsing → GWAS Matching → Filtering by genotype match → TSV Export → Cloud Sync
- **Mandatory Criteria**:
  - p-value ≤ 5e-8 
  - Risk genotype match
- **Output Fields**:
  - Odds Ratio (OR)
  - PubMed ID
  - Mapped Gene
  - Variant Context
 
---
### Cell 6
**PGx Filter & Reporting Pipeline**
VCF Input → Genotype Parsing → DrugDB Matching → Filtering by genotype match → TSV Export → Cloud Sync
- Requires exact genotype match in PharmGKB
- Captures evidence metrics from PharmGKB:
  - Level of Evidence
  - Affected Drugs
  - Phenotype Category
  - Clinical Annotation Score
  - PharmGKB URL
  - Genotype-based annotation.texts

---


