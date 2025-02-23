# PGS_TEDS
Creating PolyGenic Scores

```mermaid
flowchart TD
    A[1. Start] --> B[2. Input Processing]
    B --> C[3. Read Summary Statistics]
    
    C --> D{4. Quality Control}
    D --> D1[Check MAF]
    D --> D2[Check INFO scores]
    D --> D3[Convert alleles to uppercase]
    
    D1 & D2 & D3 --> E[5. Calculate Effective Sample Size]
    
    E --> F[6. Match SNPs to Reference Panel]
    F --> G[7. Calculate Standard Deviations]
    
    G --> H{8. SD Check}
    H -->|Good| I[9. Process Test Data]
    H -->|Bad > 50%| J[10. Impute Effective N]
    J --> I
    
    I --> K[11. Run LDSC]
    K --> L[12. Create Sparse Matrix]
    
    L --> M1[13. Run LDpred2-inf]
    L --> M2[14. Run LDpred2-auto]
    
    M1 --> N1[15. Generate infinitesimal weights]
    M2 --> N2[15. Generate auto weights]
    
    N1 --> O1[16. Calculate infinitesimal scores]
    N2 --> O2[16. Calculate auto scores]
    
    O1 & O2 --> P[17. Save Results]
    P --> Q[End]

