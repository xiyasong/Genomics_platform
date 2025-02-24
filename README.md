# Genomics_platform

   
## Overview

A repository for the preparing manuscript "A platform/webserver for genomic data analysis, reporting, and interpretation"

------------------------------------------------------------------------


why PGS scores are hard to generate based on singel-sample vcf of WGS files: see discussion https://github.com/PGScatalog/pgsc_calc/discussions/98

```bash
my_pipeline/
├── Snakefile             # Snakemake 
├── config.yaml           # config
├── envs/
│   ├── vep.yaml          # vep conda env
│   ├── python.yaml       # python conda env
│   └── ...
├── scripts/
│   └── pythonpipelinefix_AWSconnection.py  # vep annotated analysis py script
├── vcf/                  # original vcf
│   └── sample1.vcf
├── vcf_pass/             # vcf file after filtering "PASS"
├── vep/                  # vep.annotated
└── python_results/              # results files from py script
