# Genomics_platform

   
## Overview

A repository for the preparing manuscript "A platform/webserver for genomic data analysis, reporting, and interpretation"

------------------------------------------------------------------------


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
