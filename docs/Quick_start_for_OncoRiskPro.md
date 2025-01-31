
# **Quick start**
**OncoRiskPro: A Web Server for Oncology Therapy Reports from Tumor Sequencing Data**  


## Introduction

OncoRiskPro is a specialised web server that streamlines the analysis and interpretation of tumour sequencing data and generates clinical reports from Illumina's established TruSight Oncology 500 (TSO500) assay (1) and other validated oncology sequencing workflows. The platform processes standard TSO500 outputs, including annotated variants in JSON format, raw sequencing data in VCF format, and DNA/RNA-combined variant outputs in TSV format, implementing clinically validated filtering strategies based on ClinVar (2) classifications and standard variant annotations. The server employs a systematic variant prioritization approach incorporating well-established clinical resources and databases. The analysis pipeline integrates 54 FDA-approved biomarkers from FoundationOne CDx (Foundation Medicine, Inc.)(3), providing a validated framework for therapeutic associations. OncoRiskPro evaluates both DNA alterations (SNVs/Indels and structural variants) and RNA changes (gene fusions and splice variants), along with standardized assessments of tumour mutation burden (TMB ≥10 mutations/megabase) and microsatellite instability (MSI ≥30%). The platform generates structured clinical reports prioritising FDA-approved therapy-associated biomarkers, followed by variants in biomarker-associated genes, and provides comprehensive variant data in an appendix format. Each positive biomarker result links to validated clinical resources, including FDA-approved drugs and drugs tested in clinical trials from clinicaltrials.gov(4). This integrated platform will be particularly valuable for clinical oncologists, pathologists, and researchers needing standardized analysis of tumour sequencing data using established biomarkers and reporting frameworks. The web server is freely accessible, with login requirements only for sensitive data analysis, and will be maintained for at least ten years.


## **Submit jobs**

Users are required to sign up with a personal email address to start analysis for sensitive genome data. BabyGenePro currently only accepts files in the vcf format. Four types of data are available to be analysed in BabyGenePro: 1) single nucleotide polymorphism (SNPs), 2) insertions and deletions (INDELs); 3) copy number variants (CNVs); 4) other structural variants (SVs).

Users can also play with the example file without login. 

![image](https://github.com/user-attachments/assets/8c758e1f-cb59-4c77-b3d5-b202b3abd10c)


**References:**


1.	Zhao C, Jiang T, Ju JH, Zhang S, Tao J, Fu Y, et al. TruSight Oncology 500: Enabling Comprehensive Genomic Profiling and Biomarker Reporting with Targeted Sequencing [Internet]. bioRxiv; 2020 [cited 2024 Dec 20]. p. 2020.10.21.349100. Available from: https://www.biorxiv.org/content/10.1101/2020.10.21.349100v1
2.	Landrum MJ, Lee JM, Riley GR, Jang W, Rubinstein WS, Church DM, et al. ClinVar: public archive of relationships among sequence variation and human phenotype. Nucleic Acids Res. 2014 Jan;42(Database issue):D980-985. 
3.	Health C for D and R. List of Cleared or Approved Companion Diagnostic Devices (In Vitro and Imaging Tools). FDA [Internet]. 2024 Nov 15 [cited 2024 Dec 20]; Available from: https://www.fda.gov/medical-devices/in-vitro-diagnostics/list-cleared-or-approved-companion-diagnostic-devices-in-vitro-and-imaging-tools
4.	Home | ClinicalTrials.gov [Internet]. [cited 2024 Dec 20]. Available from: https://clinicaltrials.gov/



        
