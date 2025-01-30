![image](https://github.com/user-attachments/assets/d0b12757-23b0-4249-8bf4-a0f2e26f60f0)# Allele Frequency Tools  

## Introduction of the tool  

This tool is used as a "blueprint" of all variants that found in the in-house uploaded sequencing samples. The webpage mainly including two sections:
- The **SNPs** in-house allele frequencies that dynamically calculated from all the uploaded samples.
- The **CNV/SVs** searching functions, that allow users to find a specific genomic regions, checking the CNVs/SVs records that existed in the uploaded samples.

## Quick start

### SNP allele frequencies
<img width="1051" alt="image" src="https://github.com/user-attachments/assets/7f217088-84e6-4584-9ded-e3184c170105" />

1. **View files**  
By clicking this button, a list of all uploaded files (by sample name) are shown.  
  
2. **Export all calculation results of in-house allele frequencies**  
The calculation is based on the formula:
Allele Frequency = Total allele count of the SNP / 2* (Total Number of Individuals) whereas Total allele count is calculated referring the heterozygous/homozygous status of the variant.

4. **Search / Scrolling tab**
Here, by searching (input the chromosome location of the variants that interested) and scrolling the user is able to view all variants detected.

5.**Export variant-specific results** 
This button allows downloading of a CSV file including all the sample IDs that carrying the specific variant being searched. 


<img width="1051" alt="image" src="https://github.com/user-attachments/assets/ba737f00-7deb-44e1-8065-e7d706c4881b" />

<img width="1055" alt="image" src="https://github.com/user-attachments/assets/60511624-2f52-4d1d-80b8-fea41673a0f4" />


All of SNPs are named as the format chr:pos:Ref/Alt, the percentage indicated the newest calculated allele frequency. 

---
### CNVs/SVs investigation based on DECIPHER database
<img width="1197" alt="image" src="https://github.com/user-attachments/assets/4458a96a-8810-4613-9ce5-a4effa740842" />

<img width="1186" alt="image" src="https://github.com/user-attachments/assets/de2768dd-0fcf-4742-9f87-5320efedadea" />
