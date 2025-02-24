# PGS_TEDS

```mermaid
flowchart TD
  A[Prepare Environment] -->|"Open Interactive Session<br> Add R Module<br> Install R Packages<br> Set up logging file<br> Define input/output paths<br>"| B[Start Script]

  B -->|"Check SumStats format<br> Rename SumStats file <br> Load SumStats"| C[GWAS Summary Statistics]
  C -->|Save/Replace/Update| C2[scores_to_create.csv]

  B -->|HapMap3+ reference map| E[LDpred2 Reference Data]
  
  C -->|"Read Summary Statistics"| F[QC on GWAS Data]
  F -->|"-Check Column names <br> -Uppercase alleles <br> -Effective Sample Size<br> -check/remove MAFs<br> -SNPs with INFO <0.6 <br> -Remove duplicate SNPs <br>"| G[Match SNPs with LD Reference]
  E --> G 

  G --> |"-Track flipped SNPs<br> -Track reversed SNPs<br> -Remove duplicates<br>"| H[Compute GWAS Statistics]
  H --> |"-Calculate Chi-square<br> -Calculate SD<br>"| H1[Compute SNP Standard Deviation]
  
  H1 --> |"-Calculate SD using UKBB<br>-Check for Over/Under SD<br>-Flag Bad SNPs"| H2[Filter Bad SNPs]
  
  H2 --> |"-Remove SNPs with Unreliable SD"| I[Match SNPs with Test Data]
  H2 -->|"If >50% fail QC then:<br>Generate SD Plot<br>"| V2[Save SD Plot - sd.medianN.png]

  I -->|"Load Test Data (Genotypes + Metadata)"| T1[Load Test Genotype Data]
  T1 -->|"Extract SNP Metadata<br>(Chromosome, Position, Alleles)"| T2[Extract Test SNP Map]
  T2 -->|"Check + Remove Duplicates"| J[Prepare Test Data for Analysis]

  J -->|"Match GWAS SNPs with Test Data"| T3[Match SNPs in Test Data]
  T3 -->|"Align Test SNPs with GWAS Data<br>(Track Flipped/Reversed SNPs)"| T4[Matched SNPs Ready for PRS]

  T4 -->|"Estimate Heritability (LDSC)"| K[Estimate Heritability LDSC]
  G --> K

  K -->|Convert LD Matrix| L[Sparse Matrix Format]
  E --> L
  L -->|Run Model| M[LDpred2-inf]
  L -->|Run Model| N[LDpred2-auto]
  
  M --> O[Compute PRS LDpred2-inf]
  N --> P[Compute PRS LDpred2-auto]
   
  O -->|"PRS Computed Using Genotypes + Effect Sizes (big_prodVec)"| T5[Final PRS for Test Data]
  P --> T5
  T5 --> Q[Save PRS Results]

  Q --> R[End Script]

  %% Completely Separate Highlighted Boxes
  V1(["📈 Save sd.png"]):::highlight
  V3(["📈 Save auto_chains.png"]):::highlight

  H1 -->|SD Comparison Plot| V1
  N -->|LDpred2-auto Plot| V3

  %% Styling for Highlighted Boxes
  classDef highlight fill:#FFD700,stroke:#000,stroke-width:2px;
