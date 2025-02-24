# PGS_TEDS

For case/control traits: use cases and controls
For continuous traits: use provided N
Calculate chi-square statistics


Match SNPs to Reference Panel

Read HapMap3+ reference map
Match SNPs to reference
Identify flipped and reversed SNPs
Remove duplicates


Standard Deviation Processing

Calculate reference and summary statistics SDs
Flag bad SNPs based on SD criteria
If >50% SNPs are bad, impute effective N
Generate SD plots


Process Test Data

Load test genotype data
Match with variants in test data
Prepare for matrix multiplication


Run LDSC (Linkage Disequilibrium Score Regression)

Calculate heritability estimates
Process LD scores


Create Sparse Matrix

Process chromosomes 1-22
Create correlation matrices
Save as sparse format


LDpred2 Analysis
a. LDpred2-inf:

Compute infinitesimal model
Generate weights
Calculate scores

b. LDpred2-auto:

Run auto model with multiple chains
Check convergence
Filter outlier betas
Calculate final scores


Results Processing

Save weights and scores
Generate plots
Clean up temporary files
Complete logging

```mermaid
flowchart TD
  A[Prepare Environment] -->|"Open Interactive Session<br> Add R Module<br> Install R Packages<br> Install R Packages<br> Set up logging file<br> Define input/output paths<br>"| B[Start Script]

  B -->|"Check SumStats format<br> Rename SumStats file <br> Load SumStats"| C[GWAS Summary Statistics]
  C -->|Save/Replace/Update| C2[scores_to_create.csv]

  B -->|Load HapMap3+ .rds file| E[ldpred2 Reference Data]
  
  C --> |"Read Summary Statistics"| F[QC on GWAS Data]
  F -->|"Check Column names <br> Uppercase alleles <br> Effective Sample Size<br> check/remove MAFs<br> SNPs with INFO <0.6 <br> Remove duplicate SNPs <br>"| G[Match SNPs with LD Reference]
  E --> G

  G --> H[Compute GWAS Statistics]
  H --> I[Match SNPs with Test Data]

  I -->|Check & Remove Duplicates| J[Prepare Data for Analysis]
  J --> K[Estimate Heritability LDSC]
  G --> K

  K -->|Convert LD Matrix| L[Sparse Matrix Format]
  E --> L
  L -->|Run Model| M[LDpred2-inf]
  L -->|Run Model| N[LDpred2-auto]
  
  M --> O[Compute PRS LDpred2-inf]
  N --> P[Compute PRS LDpred2-auto]
   
  O --> Q[Save PRS Results]
  P --> Q

  Q --> R[End Script]
