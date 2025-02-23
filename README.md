# PGS_TEDS
Creating PolyGenic Scores

```mermaid
flowchart TD
    A[Start] --> B[Input Processing]
    B --> C[Read Summary Statistics]
    
    C --> D{Quality Control}
    D --> D1[Check MAF]
    D --> D2[Check INFO scores]
    D --> D3[Convert alleles to uppercase]
    
    D1 & D2 & D3 --> E[Calculate Effective Sample Size]
    
    E --> F[Match SNPs to Reference Panel]
    F --> G[Calculate Standard Deviations]
    
    G --> H{SD Check}
    H -->|Good| I[Process Test Data]
    H -->|Bad > 50%| J[Impute Effective N]
    J --> I
    
    I --> K[Run LDSC]
    K --> L[Create Sparse Matrix]
    
    L --> M1[Run LDpred2-inf]
    L --> M2[Run LDpred2-auto]
    
    M1 --> N1[Generate infinitesimal weights]
    M2 --> N2[Generate auto weights]
    
    N1 --> O1[Calculate infinitesimal scores]
    N2 --> O2[Calculate auto scores]
    
    O1 & O2 --> P[Save Results]
    P --> Q[End]


%% Annotations
    classDef note fill:#f9f,stroke:#333,stroke-width:2px;
    class A note;

    %% Adding text annotations
    A:::note -->|Note| A [This is the start point]
