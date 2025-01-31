# **OncoRiskPro: A Web Server for Oncology Therapy Reports from Tumor Sequencing Data** 

## **Introduction**

OncoRiskPro is a specialised web server that streamlines the analysis and interpretation of tumour sequencing data and generates clinical reports from Illumina's established TruSight Oncology 500 (TSO500) assay (1) and other validated oncology sequencing workflows. The platform processes standard TSO500 outputs, including annotated variants in JSON format, raw sequencing data in VCF format, and DNA/RNA-combined variant outputs in TSV format, implementing clinically validated filtering strategies based on ClinVar (2) classifications and standard variant annotations. The server employs a systematic variant prioritization approach incorporating well-established clinical resources and databases. The analysis pipeline integrates 54 FDA-approved biomarkers from FoundationOne CDx (Foundation Medicine, Inc.)(3), providing a validated framework for therapeutic associations. OncoRiskPro evaluates both DNA alterations (SNVs/Indels and structural variants) and RNA changes (gene fusions and splice variants), along with standardized assessments of tumour mutation burden (TMB ≥10 mutations/megabase) and microsatellite instability (MSI ≥30%). The platform generates structured clinical reports prioritising FDA-approved therapy-associated biomarkers, followed by variants in biomarker-associated genes, and provides comprehensive variant data in an appendix format. Each positive biomarker result links to validated clinical resources, including FDA-approved drugs and drugs tested in clinical trials from clinicaltrials.gov(4). This integrated platform will be particularly valuable for clinical oncologists, pathologists, and researchers needing standardized analysis of tumour sequencing data using established biomarkers and reporting frameworks. The web server is freely accessible, with login requirements only for sensitive data analysis, and will be maintained for at least ten years.


## **Submit jobs**

Users are required to sign up with a personal email address to start analysis for sensitive genome data.
![image](https://github.com/user-attachments/assets/ba621a51-84a8-439b-93f5-cbe3ef99a819)

---
## Background
<img width="1210" alt="image" src="https://github.com/user-attachments/assets/92ef7884-3cae-4a8a-a116-35e9af6b0e56" />

<img width="780" alt="image" src="https://github.com/user-attachments/assets/da6a4c50-4804-4879-be8c-e9e337903a84" />

## **Methodology Explainations**

   1. **Data Upload**  
   2. **Automated Filtering & Annotation or adjustable filtering**  
   3. **Result Prioritization**  
   4. **Reporting & Clinical/Therapy/medicine Association and generate PDF format**  

---

## **Quick start**

### 1. Data Upload

Options of accepted input files:
- TSO500 JSON (annotated variants file) 
- Combined DNA TSV (Outputs from Illumina TSO500)
- Combined DNA TSV (Outputs from Illumina TSO500)
- Raw VCF files for SNV/INdels from other sequencing facilities.

### 2. Automated Filtering & Annotation

- Remove ClinVar-known benign variants.
- Discard synonymous, UTR, and intronic variants lacking any pathogenic evidence.
- Adjustable thresholds for VAF, global MAF, and variant consequences.


### 3. Result Prioritization

**1)Therapy-Associated Biomarkers**    
Matches current FDA-approved companion diagnostics list (54 biomarkers, including SNVs/Indels, CNVs, fusions, splice variants), TMB ≥10 mut/Mb, and MSI ≥30%.

**2) Other Prioritized Variants**    
Genes related to known biomarkers but not directly FDA-approved.

**3) Appendix**    
All remaining variants not covered in the above categories.

- **Tech Note:**  
  > “OncoRisk maps detected variants to the latest FDA-approved biomarker list, highlights other relevant oncogenes, and retains an appendix of all filtered variants.”


### 4. Reporting & Clinical Association

- **Tech Note:**  
  > “A comprehensive, structured report is generated, highlighting actionable biomarkers, corresponding therapies, and clinical trial information. Users can export the final report as a PDF and perform manual review or modifications as needed.”

---

Users can also play with the example file without login. 

![image](https://github.com/user-attachments/assets/8c758e1f-cb59-4c77-b3d5-b202b3abd10c)


**References:**


1.	Zhao C, Jiang T, Ju JH, Zhang S, Tao J, Fu Y, et al. TruSight Oncology 500: Enabling Comprehensive Genomic Profiling and Biomarker Reporting with Targeted Sequencing [Internet]. bioRxiv; 2020 [cited 2024 Dec 20]. p. 2020.10.21.349100. Available from: https://www.biorxiv.org/content/10.1101/2020.10.21.349100v1
2.	Landrum MJ, Lee JM, Riley GR, Jang W, Rubinstein WS, Church DM, et al. ClinVar: public archive of relationships among sequence variation and human phenotype. Nucleic Acids Res. 2014 Jan;42(Database issue):D980-985. 
3.	Health C for D and R. List of Cleared or Approved Companion Diagnostic Devices (In Vitro and Imaging Tools). FDA [Internet]. 2024 Nov 15 [cited 2024 Dec 20]; Available from: https://www.fda.gov/medical-devices/in-vitro-diagnostics/list-cleared-or-approved-companion-diagnostic-devices-in-vitro-and-imaging-tools
4.	Home | ClinicalTrials.gov [Internet]. [cited 2024 Dec 20]. Available from: https://clinicaltrials.gov/



        
