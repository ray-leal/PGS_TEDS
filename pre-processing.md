```mermaid
flowchart TD

A[Start] --> B(Load Libraries and Arguments);
B --> C{Check for CHR, ORIGBP or RSID};
C -- CHR & ORIGBP --> D[Determine Genome Build];
C -- RSID --> E[Harmonize with Reference by RSID];
D --> F[Harmonize with Reference by CHR and ORIGBP];
F --> G;
E --> G;
G --> H{Check INFO};
H -- INFO present --> I[Filter by INFO];
H -- INFO absent --> J[Log: INFO absent];
I --> K{Check FREQ};
J --> K{Check FREQ};
K -- FREQ present --> L[Filter by Reported MAF];
L --> M[Filter by Reference MAF];
M --> N[Filter by MAF Discordance];
K -- FREQ absent --> O[Log: FREQ absent];
N --> P[Filter by P-value];
O --> P[Filter by P-value];
P --> Q[Filter Duplicates];
Q --> R[Filter by N];
R --> S[Check Genomic Control];
S -- Genomic Control Detected --> T[Recalculate P];
S -- Genomic Control Not Detected --> U[Log: No Genomic Control];
T --> V{Check for SE};
U --> V{Check for SE};
V -- SE present --> W[Filter SE == 0];
V -- SE absent --> X[Calculate SE];
X --> W[Filter SE == 0];
W --> Y[Uppercase A1/A2];
Y --> Z{Gzip Output?};
Z -- Yes --> A1[Write Gzipped Output];
Z -- No --> A2[Write Output];
A1 --> A3[End];
A2 --> A3[End];
