# PGS_TEDS
Creating PolyGenic Scores

Input Processing

Load required libraries
Parse command-line arguments
Set up logging
Define input/output paths


Read Summary Statistics

Load GWAS summary statistics
Check required columns (CHR, BP, A1, A2, BETA, SE)
Log dimensions and headers


Quality Control

Check and filter Minor Allele Frequency (MAF)
Filter on INFO scores if present
Convert alleles to uppercase
Remove duplicate SNPs


Calculate Effective Sample Size

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
    A[1 - Start and Setup] --> B[2 - Input Processing]
    B --> C[3 - Read Summary Statistics]
    
    C --> D{4 - Quality Control}
    D --> D1[Check MAF]
    D --> D2[Check INFO scores]
    D --> D3[Convert alleles to uppercase]
    
    D1 & D2 & D3 --> E[5 - Calculate Effective Sample Size]
    
    E --> F[6 - Match SNPs to Reference Panel]
    F --> G[7 - Calculate Standard Deviations]
    
    G --> H{8 - SD Check}
    H -->|Good| I[9 - Process Test Data]
    H -->|Bad > 50%| J[10 - Impute Effective N]
    J --> I
    
    I --> K[11 - Run LDSC]
    K --> L[12 - Create Sparse Matrix]
    
    L --> M1[13 - Run LDpred2-inf]
    L --> M2[14 - Run LDpred2-auto]
    
    M1 --> N1[15 - Generate infinitesimal weights]
    M2 --> N2[15 - Generate auto weights]
    
    N1 --> O1[16 - Calculate infinitesimal scores]
    N2 --> O2[16 - Calculate auto scores]
    
    O1 & O2 --> P[17 - Save Results]
    P --> Q[End]

