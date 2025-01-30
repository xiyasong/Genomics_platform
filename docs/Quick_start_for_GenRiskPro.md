# **Quick start**
**GenRiskPro: A Webserver for Generating Clinical and Wellness Reports from NGS Data**  


## Introduction

GenRiskPro is a comprehensive web server for genomic variant annotation, interpretation, and reporting that implements established clinical genomics pipelines and guidelines. The platform systematically applies the American College of Medical Genetics and Genomics (ACMG) variant classification criteria by integrating widely validated genomic resources, including ClinVar(1), OMIM(2), PharmGKB(3), and the GWAS Catalog(4). The analysis begins with standard variant call format (VCF) files and employs curated gene panels from established NGS screening projects, with flexibility for user-defined panels.
The server generates detailed variant annotations by leveraging extensively validated databases and computational tools. Outputs include prioritized clinically relevant variants with comprehensive annotations of pathogenicity, population frequencies, inheritance patterns, and phenotypic associations. The system incorporates established resources to provide broader insights, including GWAS-derived phenotypic associations, evidence-based wellness and trait information across 15 complex trait categories, and pharmacogenetic annotations from PharmGKB with validated genotype-phenotype correlations. GenRiskPro implements standard statistical methods for population-level analysis to calculate allele frequencies, as demonstrated through analysis of whole genome sequencing data from Turkish (n=275) and Swedish (n=101) cohorts. This integrated platform will be particularly valuable for clinical laboratories, genetic counsellors, and researchers requiring standardized variant interpretation and reporting workflows based on established genomic resources. The web server is freely accessible with login requirements only for sensitive data analysis and will be maintained for at least ten years.
<img width="1018" alt="image" src="https://github.com/user-attachments/assets/01ad80da-ce90-4ad3-abdf-7be746aaed07" />


## **Submit jobs**

Users are required to sign up with a personal email address to start analysis for sensitive genome data. BabyGenePro currently only accepts files in the vcf format. Four types of data are available to be analysed in BabyGenePro: 1) single nucleotide polymorphism (SNPs), 2) insertions and deletions (INDELs); 3) copy number variants (CNVs); 4) other structural variants (SVs).

Users can also play with the example file without login. 

![image](https://github.com/user-attachments/assets/8c758e1f-cb59-4c77-b3d5-b202b3abd10c)


**References:**

1.	Landrum MJ, Lee JM, Riley GR, Jang W, Rubinstein WS, Church DM, et al. ClinVar: public archive of relationships among sequence variation and human phenotype. Nucleic Acids Res. 2014 Jan;42(Database issue):D980-985. 
2.	Amberger JS, Bocchini CA, Schiettecatte F, Scott AF, Hamosh A. OMIM.org: Online Mendelian Inheritance in Man (OMIM®), an online catalog of human genes and genetic disorders. Nucleic Acids Res. 2015 Jan;43(Database issue):D789-798. 
3.	Whirl-Carrillo M, Huddart R, Gong L, Sangkuhl K, Thorn CF, Whaley R, et al. An Evidence-Based Framework for Evaluating Pharmacogenomics Knowledge for Personalized Medicine. Clin Pharmacol Ther. 2021 Sep;110(3):563–72. 
4.	Sollis E, Mosaku A, Abid A, Buniello A, Cerezo M, Gil L, et al. The NHGRI-EBI GWAS Catalog: knowledgebase and deposition resource. Nucleic Acids Res. 2023 Jan 6;51(D1):D977–85. 

        
