# TEDS - Creating new polygenics scores
This document describes how you should creat a new polygenic score in the TEDS twin dataset, using our standard LDpred2 pipeline. This pipeline came into effect in November 2022; the old LDpred1 pipeline is no longer functioning.

The links here document the TEDS genotypic dataset, our QC procedures, and how the TEDS dataset has been used to prepare our polygenic scores. The TEDS LDpred2 pipeline will generate a file containing a single polygenic score variable. 
This is equivalent to the threshold-1 score generated in the earlier TEDS LDpred1 pipeline. This new score will be made available for more widespread sharing, by inclusion in the list in the data dictionary:
https://www.teds.ac.uk/datadictionary/studies/measures/polygenic_scores.htm

Every new score will immediately be made available within the TEDS team for analysis. Polygenic scores are included in shared TEDS datasets following our usual data request procedure (see http://www.teds.ac.uk/researchers/teds-data-access-policy).

# Steps
[1] Basic information
- Published paper details 
- Readme file
- GWAS data

[2] Data check 
- Var name checker/rename
- Check columns vs required
- Create file for MetaData

[3] Pre-clean (SumState cleaner)
