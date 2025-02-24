# PGS_TEDS

```mermaid
flowchart TD
  %% Prepare Environment
  A[Prepare Environment] -->|"Open Interactive Session<br> Add R Module<br> Install R Packages<br> Set up logging file<br> Define input/output paths"| B[Start Script]

  %% Move "Load Test Genotype Data" to the Far Left
  TL1["Load Test Genotype Data"]
  B -->|"Load Test Data (Genotypes + Metadata)"| TL1
  TL1 -->|"Extract SNP Metadata<br>(Chromosome, Position, Alleles)"| T2[Extract Test SNP Map]
  T2 -->|"Check + Remove Duplicates"| T3[Prepare Test Data for Analysis]

  %% GWAS Data Processing
  B -->|"Check SumStats format<br> Rename SumStats file <br> Load SumStats"| C[GWAS Summary Statistics]
  C -->|Save/Replace/Update| C2[scores_to_create.csv]

  B -->|"HapMap3+ reference map"| E[LDpred2 Reference Data]

  C -->|"Read Summary Statistics"| F[QC on GWAS Data]
  F -->|"-Check Column names <br> -Uppercase alleles <br> -Effective Sample Size<br> -Check/remove MAFs<br> -SNPs with INFO <0.6 <br> -Remove duplicate SNPs"| G[Match SNPs with LD Reference]
  E --> G

  G --> |"-Track flipped SNPs<br> -Track reversed SNPs<br> -Remove duplicates"| H[Compute GWAS Statistics]
  H --> |"-Calculate Chi-square<br> -Calculate SD"| H1[Compute SNP Standard Deviation]

  H1 --> |"-Calculate SD using UKBB<br>-Check for Over/Under SD<br>-Flag Bad SNPs"| H2[Filter Bad SNPs]
  
  H2 --> |"-Remove SNPs with Unreliable SD"| I[Prepare Data for PRS]
  H2 -->|"If >50% fail QC:<br>Generate SD Plot (sd.medianN.png)"| V2[Save Updated SD Plot]

  %% Connect Test Data to Matching Step
  I -->|"Match GWAS SNPs with Test Data"| T4[Match SNPs in Test Data]
  T3 --> T4
  T4 -->|"Align Test SNPs with GWAS Data<br>(Track Flipped/Reversed SNPs)"| T5[Matched SNPs Ready for PRS]

  T5 -->|"Estimate Heritability (LDSC)"| K[Estimate Heritability LDSC]
  G --> K

  K -->|Convert LD Matrix| L[Sparse Matrix Format]
  E --> L
  L -->|Run Model| M[LDpred2-inf]
  L -->|Run Model| N[LDpred2-auto]

  M -->|"Compute PRS (LDpred2-inf)"| O[Compute PGS LDpred2-inf]
  N -->|"Compute PRS (LDpred2-auto)"| P[Compute PGS LDpred2-auto]

  O -->|"PGS Computed Using Genotypes + Effect Sizes (big_prodVec)"| T6[Final PGS for Test Data]
  P --> T6
  T6 --> Q[Save PGS Results in log]:::highlight

  Q --> R[End Script]

  %% Separate Highlighted Boxes for Outputs
  V1(["📈 Save SD Comparison Plot (sd.png)"]):::highlight
  V3(["📈 Save LDpred2-auto Plot (auto_chains.png)"]):::highlight

  H1 -->|Generates| V1
  N -->|Generates| V3

  %% Styling for Highlighted Boxes
  classDef highlight fill:#FFD700,stroke:#000,stroke-width:2px;
