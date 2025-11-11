# Analysis of CNV and SV files for Patient P00422

The translation from Russian report that saying the patient's diseases:
- 1. Chromosome 9p24.3–p23 Microdeletion
	•	Genomic location (GRCh38): chr9: 0–13,587,692
	•	Type: Heterozygous microdeletion
	•	Size: ~13.6 Mb
	•	Copy number: 1 (loss)
	•	Affected genes include: DOCK8, KANK1, SMARCA2, VLDLR, KCNV2, GLIS3, SLC1A1, JAK2, CD274, RIC1, GLDC
	•	Associated syndrome: Chromosome 9p deletion syndrome (OMIM 158170)
	•	Clinical description: This syndrome is typically associated with global developmental delay, intellectual disability, craniofacial dysmorphism, and sometimes structural brain or heart abnormalities.
	•	Interpretation: Classified as pathogenic and considered the most likely cause of the patient’s clinical presentation.


- 2. Chromosome 10q25.1–q26.3 Microduplication
	•	Genomic location (GRCh38): chr10: 107,558,240–133,797,422
	•	Type: Heterozygous microduplication
	•	Size: ~26.2 Mb
	•	Copy number: 3 (gain)
	•	Affected genes include: ACADSB, FGFR2, HTRA1, SHOC2, TCF7L2, WDR11, BAG3, EMX2, GFRA1, MXI1, VAX1, and others
	•	Associated syndrome: Distal duplication 10q syndrome
	•	Clinical description: Usually linked to developmental delay, growth retardation, abnormal muscle tone, facial and organ malformations, and neurodevelopmental abnormalities.
	•	Interpretation: Classified as pathogenic and also considered a likely cause of the patient’s phenotype.
"


Here, I summarizes the query of the `.cnv.vcf` and `.sv.vcf` files to verify the presence of two reported structural variants.

## 1. Chromosome 9p24.3–p23 Microdeletion

### Grep Command

To find deletions on chromosome 9, the following `grep` command was used on both the `cnv` and `sv` files:

```bash
grep 'chr9' /Users/xiyas/Downloads/szaomics_p00422_24_6496_2025-10-10_1106/P00422_24_6496.cnv.vcf | grep 'DEL'
grep 'chr9' /Users/xiyas/Downloads/szaomics_p00422_24_6496_2025-10-10_1106/P00422_24_6496.sv.vcf | grep 'DEL'
```

### Findings

A large, heterozygous deletion was identified on chromosome 9 in both files, consistent with the clinical report.

-   **In `P00422_24_6496.cnv.vcf`:**
    -   A deletion spanning from position **114,000 to 11,495,204** (approximately 11.4 Mb) was found.
    -   The genotype (`GT`) is `0/1` (heterozygous) and the copy number (`CN`) is `1`.
    -   **VCF Line:**
      ```
        chr9	114000	DRAGEN:LOSS:chr9:114001-11495204	N	<DEL>	104	PASS	SVLEN=-11381204;SVTYPE=CNV;END=11495204;REFLEN=11381204	GT:SM:CN:BC:GC:CT:AC:PE	0/1:0.507353:1:10057:0.386923:0.499821:0.500456:0,18
      ```
    -   **Annotation with AnnotSV**
    	-   ACMG prediction: pathogenic (score = 5)
     	-  	Overlapped with genes:
      		- 		DMRT2, DMRT3, DOCK8, DOCK8-AS1, DOCK8-AS2, ERMP1, FOXD4, GLDC, GLIS3, GLIS3-AS1, GLIS3-AS2, IL33, INSL4, INSL6, JAK2, KANK1, KCNV2, KDM4C, MIR101-2, MIR4665, MLANA, PDCD1LG2, PLGRKT, PLPP6, PTPRD, PTPRD-AS1, PTPRD-DT, PUM3, RANBP6, RCL1, RFX3, RFX3-DT, RIC1, RLN1, RLN2, SLC1A1, SMARCA2, SPATA6L, TPD52L3, UHRF2, VLDLR, VLDLR-AS1, ZNG1A
    	- Overlapped pathogenic clinical evidence: (See attachments)
       	- Overlapped pathogenic phenotypes:
        	- 	Blepharophimosis-impaired intellectual development syndrome, 619293 (3) AD;Nicolaides-Baraitser syndrome, 601358 (3) AD;CATIFA syndrome, 618761 (3) AR;Cerebellar hypoplasia, impaired intellectual development, and dysequilibrium syndrome 1, 224050 (3) AR;Cerebral palsy, spastic quadriplegic, 2, 612900 (3);Diabetes mellitus, neonatal, with congenital hypothyroidism, 610199 (3) AR;Dicarboxylic aminoaciduria, 222730 (3) AR;?Schizophrenia susceptibility 18, 615232 (3);Erythrocytosis, somatic, 133100 (3);Leukemia, acute myeloid, somatic, 601626 (3);Myelofibrosis, somatic, 254450 (3);Polycythemia vera, somatic, 263300 (3);Thrombocythemia 3, 614521 (3) Somatic mutation,AD;Budd-Chiari syndrome, somatic, 600880 (3);Glycine encephalopathy1, 605899 (3) AR;Hyper-IgE syndrome 2, AR, with recurrent infections, 243700 (3) AR;Non-ketotic_hyperglycinemia;Retinal cone dystrophy 3B, 610356 (3) AR

-   **In `P00422_24_6496.sv.vcf`:**
    -   A large deletion from position **6,995,965 to 21,038,487** (approximately 14 Mb) was identified.
    -   The genotype (`GT`) is `0/1` (heterozygous).
    -   **VCF Line:**
        ```
        chr9	6995965	DRAGEN:DEL:135161:0:1:0:0:0	G	<DEL>	659	NoPairSupport	END=21038487;SVTYPE=DEL;SVLEN=-14042522;CIPOS=0,7;CIEND=0,7;HOMLEN=7;HOMSEQ=CCCTGCA	GT:FT:GQ:PL:PR:SR	0/1:PASS:342:709,0,339:15,0:16,19
        ```
    -   **Annotation with AnnotSV**
    	-   ACMG prediction: pathogenic (score = 5)
     	-  	Overlapped with genes within 0–13,587,692:
      		- 		KDM4C, DMAC1, PTPRD, PTPRD-AS1, PTPRD-DT, TYRP1, LURAP1L-AS1, LURAP1L, SNORD137, MPDZ
    	- Overlapped pathogenic clinical evidence: dbVar:nssv18786580;morbid:TYRP1; dbVar:nssv17171485;dbVar:nssv18788314;morbid:MPDZ
       	- Overlapped pathogenic phenotypes:
        	- 	Albinism, oculocutaneous, type III, 203290 (3) AR;Skin;hair;eye pigmentation, variation in, 11 (Melanesian blond hair), 612271 (3); Hydrocephalus, congenital, 2, with or without brain or eye anomalies, 615219 (3) AR


## 2. Chromosome 10q25.1–q26.3 Microduplication

### Grep Command

To find duplications on chromosome 10, the following `grep` command was used:

```bash
grep 'chr10' /Users/xiyas/Downloads/szaomics_p00422_24_6496_2025-10-10_1106/P00422_24_6496.cnv.vcf | grep 'DUP'
grep 'chr10' /Users/xiyas/Downloads/szaomics_p00422_24_6496_2025-10-10_1106/P00422_24_6496.sv.vcf | grep 'DUP'
```

### Findings

A large, heterozygous duplication was identified on chromosome 10 in the `cnv` file. It was represented as five contiguous events.

-   **In `P00422_24_6496.cnv.vcf`:**
    -   The combined duplication spans from position **107,558,632 to 133,621,535** (approximately 26 Mb).
    -   The copy number (`CN`) for this region is `3`.
    -   **VCF Lines:**
        ```
        chr10	107558632	DRAGEN:GAIN:chr10:107558633-112352736...
        chr10	112357126	DRAGEN:GAIN:chr10:112357127-120467592...
        chr10	120469915	DRAGEN:GAIN:chr10:120469916-131622969...
        chr10	131624948	DRAGEN:GAIN:chr10:131624949-133333049...
        chr10	133339658	DRAGEN:GAIN:chr10:133339659-133621535...
        ```
    -   **Annotation with AnnotSV**
    	-	Five duplications were identified:
			- 10_107558632_112352736: ACMG Uncertain significance (score = 3)
			- 10_112357126_120467592: ACMG Pathogenic (score = 5)
			- 10_120469915_131622969: ACMG Pathogenic (score = 5)
 			- 10_131624948_133333049: ACMG Uncertain significance (score = 3)
 			- 10_133339658_133621535: ACMG Uncertain significance (score = 3)
     	- 10_112357126_120467592:
  			- Overlapped genes: ABLIM1, ACSL5, ADRB1, AFAP1L2, ATRNL1, BAG3, CACUL1, CASC2, CASP7, CCDC172, CCDC186, DCLRE1A, DENND10, EIF3A, EMX2, EMX2OS, ENO4, FAM204A, FHIP2A, GFRA1, GRK5, GRK5-IT1, HABP2, HSPA12A, HSPA12A-AS1, INPP5F, KCNK18, MCMBP, MIR2110, MIR3663, MIR3663HG, MIR4295, MIR4483, MIR4681, MIR4682, MIR9851, NANOS1, NHLRC2, NRAP, PDZD8, PLEKHS1, PLPP4, PNLIP, PNLIPRP1, PNLIPRP2, PNLIPRP3, PRDX3, PRLHR, RAB11FIP2, RGS10, SEC23IP, SFXN4, SHTN1, SLC18A2, SLC18A2-AS1, SNORA19, SNORA87, SNORD158, SPMIP5, TCF7L2, TDRD1, TIAL1, TRUB1, VAX1, VTI1A, VWA2, ZDHHC6
      		- Overlapped clinical evidence: dbVar:nssv15135427;dbVar:nssv18790444
        	- Overlapped phenotypes: Cardiomyopathy, dilated, 1HH, 613881 (3) AD;Myopathy, myofibrillar, 6, 612954 (3) AD;Combined oxidative phosphorylation deficiency 18, 615578 (3) AR;Corneal dystrophy, punctiform and polychromatic pre-Descemet, 619871 (3) AD;Spinocerebellar ataxia, AR 32, 619862 (3) AR;Dilated_cardiomyopathy_1HH;FINCA syndrome, 618278 (3) AR;Intellectual developmental disorder with autism and dysmorphic facies, 620021 (3) AR;Parkinsonism-dystonia, infantile, 2, 618049 (3) AR;Renal hypodysplasia;aplasia 4, 619887 (3) AR;Schizencephaly, 269160 (3);Spermatogenic failure 12, 615413 (3) AD 
      	 - 10_120469915_131622969:
			- Overlapped genes: ABRAXAS2, ACADSB, ADAM12, ARMS2, AS-PTPRE, ATE1, ATE1-AS1, BCCIP, BTBD16, BUB3, C10orf120, C10orf143, C10orf88, C10orf88B, C10orf90, CHST15, CLRN3, CPXM2, CTAGE7P, CTBP2, CUZD1, DHX32, DMBT1, DMBT1L1, DOCK1, EBF3, EBF3-AS1, EDRF1, EDRF1-AS1, EDRF1-DT, EEF1AKMT2, FAM24A, FAM24B, FAM24B-CUZD1, FAM53B, FAM53B-AS1, FANK1, FANK1-AS1, FGFR2, FOXI2, GLRX3, GPR26, HMX2, HMX3, HTRA1, IKZF5, INSYN2A, LHPP, MGMT, MIR378C, MIR3941, MIR4296, MIR4297, MIR4484, MKI67, MMP21, NKX1-2, NPS, NSMCE4A, OAT, PLEKHA1, PLPP4, PSTK, PTPRE, SPADH, TACC2, TCERG1L, TCERG1L-AS1, TEX36, TEX36-AS1, UROS, WDR11, WDR11-DT, ZRANB1
        	- Overlapped clinical evidence: dbVar:nssv15121959;dbVar:nssv15148104;dbVar:nssv18792228
         	- Overalpped pheontypes: (See attachment)


-   **In `P00422_24_6496.sv.vcf`:**
    -   No large duplications corresponding to the reported event were found.

	 -   **Annotation with AnnotSV**
    	-	Two breakend duplications were found within this region:
    		- 1) 10_125502228_125512653, Breakend: chr10:125512646_C
            - 2) 10_125501872_125509034, Breakend: chr10:125508664_T	
      	-   ACMG prediction: Uncertain significance (score = 3)
     	-  	These two mutations are not overlapped with any gene. 

## Conclusion

Both the pathogenic microdeletion on chromosome 9 and the microduplication on chromosome 10 were confirmed to be present in the provided sequencing files.
