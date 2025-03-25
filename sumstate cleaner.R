```R
#!/usr/bin/Rscript
# This script was written by Oliver Pain whilst at King's College London University. 
# Now it is maintaining by Wangjingyi Liao at Queen Mary University of London.

# START CLEANING
# ---------------------------------------------------------------------------------
# record start date and time for logging the script's execution duration.
start.time <- Sys.time()
# Loads the `optparse` library to handle command-line arguments without displaying loading messages.
suppressMessages(library("optparse"))

### Command-Line Options Handling ###
# Defines the list of command-line options the script accepts.
option_list = list(
  # Option for the path to the input summary statistics file. Mandatory.
  make_option("--sumstats", action="store", default=NA, type='character',
              help="Path to summary statistics file [required]"),
  # Option for the path to the reference panel files (split by chromosome, in .rds format). Mandatory.
  make_option("--ref_chr", action="store", default=NA, type='character',
              help="Path to per chromosome reference .rds files [required]"),
  # Option for the minimum imputation INFO score threshold. Optional, default is 0.6.
  make_option("--info", action="store", default=0.6, type='numeric',
              help="INFO threshold [optional]"),
  # Option for the minimum Minor Allele Frequency (MAF) threshold. Optional, default is 0.01.
  make_option("--maf", action="store", default=0.01, type='numeric',
              help="MAF threshold [optional]"),
  # Option for the maximum allowed absolute difference between reference MAF and reported MAF. Optional, default is 0.2.
  make_option("--maf_diff", action="store", default=0.2, type='numeric',
              help="Difference between reference and reported MAF threshold [optional]"),
  # Option to specify whether the output summary statistics should be gzipped. Optional, default is TRUE.
  make_option("--gz", action="store", default=T, type='logical',
              help="Set to T to gzip summary statistics [optional]"),
  # Option for the base path/filename for output files. Optional, default is './Output'.
  make_option("--output", action="store", default='./Output', type='character',
              help="Path for output files [optional]")
)
# Parses the provided command-line arguments based on the defined option_list.
opt = parse_args(OptionParser(option_list=option_list))

# Prints the input filename being processed to the console.
cat("Cleaning file:", opt$sumstats, "\n")

# Loads the `data.table` library for efficient data manipulation, especially for large files.
library(data.table)
# Creates the full output directory path based on the --output option.
opt$output_dir<-paste0(dirname(opt$output),'/')
# Creates the output directory if it doesn't already exist, using a system command. '-p' prevents errors if it exists.
system(paste0('mkdir -p ',opt$output_dir))
# Redirects all subsequent console output (stdout and stderr) to a log file named based on the --output option. Append=F overwrites existing log.
sink(file = paste(opt$output,'.log',sep=''), append = F)

# Writes a header and metadata (script name, contact info) to the log file.
cat(
  '#################################################################
# sumstat_cleaner.R
# For questions contact Oliver Pain (oliver.pain@kcl.ac.uk)
#################################################################
Analysis started at',as.character(start.time),'
Options are:\n')

# Prints the parsed command-line options to the console (which is currently redirected to the log file).
cat('Options are:\n')
# Prints the R object containing the parsed options to the log file.
print(opt)
# Logs the analysis start time again.
cat('Analysis started at',as.character(start.time),'\n')
# Stops redirecting output to the log file, restoring output to the console for now.
sink()

#####
# Reads the GWAS summary statistics file and ensures the `P` column is numeric.
#####

# Redirects console output back to the log file, appending to the existing content.
sink(file = paste(opt$output,'.log',sep=''), append = T)
# Logs a message indicating that the GWAS summary statistics file is being read.
cat('Reading in GWAS sumstats.\n')
# Stops redirecting output to the log file.
sink()

# Reads the summary statistics file specified by --sumstats using data.table's fread for speed.
GWAS<-fread(opt$sumstats)

# Ensures the 'P' (p-value) column is treated as numeric, handling potential read-in issues.
GWAS$P <- as.numeric(GWAS$P)

# Redirects console output to the log file (appending).
sink(file = paste(opt$output,'.log',sep=''), append = T)
# Logs the total number of variants (rows) initially read from the GWAS file.
cat('GWAS contains',dim(GWAS)[1],'variants.\n')
# Stops redirecting output to the log file.
sink()


# IUPAC - Defines a function to generate IUPAC nucleotide ambiguity codes based on two alleles.
# Link provided for IUPAC codes reference.
#------------------------------------------------------------------------------

# Defines a function `snp_iupac` that takes two vectors of alleles (x, y).
snp_iupac<-function(x=NA, y=NA){
  # Checks if the input allele vectors have the same length.
  if(length(x) != length(y)){
    # Prints an error message if lengths differ.
    print('x and y are different lengths')
  } else {
    # Initializes an empty vector 'iupac' with the same length as the input alleles.
    iupac<-rep(NA, length(x))
    # Assigns 'W' for A/T polymorphisms.
    iupac[x == 'A' & y =='T' | x == 'T' & y =='A' ]<-'W'
    # Assigns 'S' for C/G polymorphisms.
    iupac[x == 'C' & y =='G' | x == 'G' & y =='C' ]<-'S'
    # Assigns 'R' for A/G polymorphisms (purines).
    iupac[x == 'A' & y =='G' | x == 'G' & y =='A' ]<-'R'
    # Assigns 'Y' for C/T polymorphisms (pyrimidines).
    iupac[x == 'C' & y =='T' | x == 'T' & y =='C' ]<-'Y'
    # Assigns 'K' for G/T polymorphisms.
    iupac[x == 'G' & y =='T' | x == 'T' & y =='G' ]<-'K'
    # Assigns 'M' for A/C polymorphisms.
    iupac[x == 'A' & y =='C' | x == 'C' & y =='A' ]<-'M'
    # Returns the vector of IUPAC codes.
    return(iupac)
  }
}
# Convert alleles to upper case
#------------------------------------------------------------------------------
# Ensures the A1 allele column contains only uppercase letters.
GWAS$A1<-toupper(GWAS$A1)
# Ensures the A2 allele column contains only uppercase letters.
GWAS$A2<-toupper(GWAS$A2)

# Insert IUPAC codes into target
#------------------------------------------------------------------------------
# Calls the `snp_iupac` function to create a new column 'IUPAC' based on the A1 and A2 alleles.
GWAS$IUPAC<-snp_iupac(GWAS$A1, GWAS$A2)

# Retain only non-ambiguous SNPs (non AT/GC, i.e., where strand cannot be uniquely determined just from alleles)
# Note: The comment says "non-ambiguous", but the codes kept (R, Y, K, M) represent transversions/transitions,
# while W and S (AT/GC pairs) are the ambiguous ones often removed. This code KEEPS R, Y, K, M.
# Let's assume the goal is to keep SNPs where strand flips can be detected/corrected reliably.
#------------------------------------------------------------------------------
# Filters the GWAS data table, keeping only rows where the IUPAC code is R, Y, K, or M.
GWAS<-GWAS[(GWAS$IUPAC %in% c('R', 'Y', 'K', 'M')),]
# Redirects output to the log file (appending).
sink(file = paste(opt$output,'.log',sep=''), append = T)
# Logs the number of variants remaining after removing ambiguous (or potentially non-SNP) variants.
cat('After removal of variants that are not SNPs or are ambiguous,',dim(GWAS)[1],'variants remain.\n')
# Stops redirecting output.
sink()


# Harmonise per chromosome with reference panel
#------------------------------------------------------------------------------
# Checks if both 'CHR' (chromosome) and 'ORIGBP' (original base pair position) columns exist in the GWAS data. Result is TRUE or FALSE.
chr_bp_avail<-sum(c('CHR','ORIGBP') %in% names(GWAS)) == 2
# Checks if the 'SNP' column exists and if more than 90% of the entries look like RSIDs (start with 'rs'). Result is TRUE or FALSE.
rsid_avail<-(sum(grepl('rs', GWAS$SNP)) > 0.9*length(GWAS$SNP))

# Defines a function to get the complement of a DNA base.
snp_allele_comp<-function(x=NA){
  # Creates a copy of the input vector.
  x_new<-x
  # Replaces 'A' with 'T'.
  x_new[x == 'A']<-'T'
  # Replaces 'T' with 'A'.
  x_new[x == 'T']<-'A'
  # Replaces 'G' with 'C'.
  x_new[x == 'G']<-'C'
  # Replaces 'C' with 'G'.
  x_new[x == 'C']<-'G'
  # Sets any non-standard base characters to NA.
  x_new[!(x %in% c('A','T','G','C'))]<-NA
  # Returns the vector of complementary bases.
  return(x_new)
}

# Initializes a variable to store the detected genome build (e.g., GRCh37). Will be NA initially.
# Comments providing rough details about different genome builds.
# ???? GRCh36 (March 3 2006)  ? genes (~? million SNPs) ? billion base pairs
# GRCh37 (February 2009) 20,000 genes (~28 million SNPs) 3.1 billion base pairs
# GRCh38 (December 2013) 19,000 genes (~80 million SNPs) 3.2 billion base pairs
target_build<-NA

# Determines the genome build and matches the reference data per chromosome.
# This block executes only if CHR and ORIGBP columns are available in the GWAS data.
if(chr_bp_avail){
  # Redirect output to log file (appending).
  sink(file = paste(opt$output,'.log',sep=''), append = T)
  # Log message indicating that CHR/BP will be used for merging.
  cat('Using CHR, BP, A1 and A2 to merge with the reference.\n')
  # Stop redirecting output.
  sink()

  ### Section to Determine the Genome Build ###
  # Read in reference data specifically for chromosome 22 to test build concordance.
  i<-22 # Sets chromosome number to 22.

  # Reads the reference panel RDS file for the selected chromosome (i=22).
  tmp<-readRDS(file = paste0(opt$ref_chr,i,'.rds'))

  # Prepare data for merging: ensure CHR columns are character type for safe merging.
  tmp$CHR<-as.character(tmp$CHR)
  GWAS$CHR<-as.character(GWAS$CHR)
  # Create an empty list to store merge results for each build.
  matched<-list()
  # Merge GWAS (chr 22 subset implicitly) with reference (chr 22) using CHR, ORIGBP, and IUPAC, trying GRCh36 BP column.
  matched[['GRCh36']]<-merge(GWAS, tmp, by.x=c('CHR','ORIGBP','IUPAC'), by.y=c('CHR','BP_GRCh36','IUPAC'))
  # Merge GWAS with reference using CHR, ORIGBP, and IUPAC, trying GRCh37 BP column.
  matched[['GRCh37']]<-merge(GWAS, tmp, by.x=c('CHR','ORIGBP','IUPAC'), by.y=c('CHR','BP_GRCh37','IUPAC'))
  # Merge GWAS with reference using CHR, ORIGBP, and IUPAC, trying GRCh38 BP column.
  matched[['GRCh38']]<-merge(GWAS, tmp, by.x=c('CHR','ORIGBP','IUPAC'), by.y=c('CHR','BP_GRCh38','IUPAC'))

  # Redirect output to log file (appending).
  sink(file = paste(opt$output,'.log',sep=''), append = T)
  # Calculate and log the percentage match for GRCh36 based on chromosome 22 SNPs.
  cat('GRCh36 match: ',round(nrow(matched[['GRCh36']])/sum(GWAS$CHR == i)*100, 2),'%\n',sep='')
  # Calculate and log the percentage match for GRCh37 based on chromosome 22 SNPs.
  cat('GRCh37 match: ',round(nrow(matched[['GRCh37']])/sum(GWAS$CHR == i)*100, 2),'%\n',sep='')
  # Calculate and log the percentage match for GRCh38 based on chromosome 22 SNPs.
  cat('GRCh38 match: ',round(nrow(matched[['GRCh38']])/sum(GWAS$CHR == i)*100, 2),'%\n',sep='')
  # Stop redirecting output.
  sink()

  # If the match percentage for GRCh36 on chr 22 is > 70%, assume the GWAS data is on GRCh36.
  if((nrow(matched[['GRCh36']])/sum(GWAS$CHR == i)) > 0.7){
    target_build<-'GRCh36'
  }
  # If the match percentage for GRCh37 on chr 22 is > 70%, assume the GWAS data is on GRCh37. (Overwrites previous if necessary)
  if((nrow(matched[['GRCh37']])/sum(GWAS$CHR == i)) > 0.7){
    target_build<-'GRCh37'
  }
  # If the match percentage for GRCh38 on chr 22 is > 70%, assume the GWAS data is on GRCh38. (Overwrites previous if necessary)
  if((nrow(matched[['GRCh38']])/sum(GWAS$CHR == i)) > 0.7){
    target_build<-'GRCh38'
  }

  # Remove the temporary matched list and reference data to free memory.
  rm(matched,tmp)

  # Proceed with harmonization only if a likely target build was identified.
  if(!is.na(target_build)){

    # Log that the detected build will be used.
    sink(file = paste(opt$output,'.log',sep=''), append = T)
    cat('Target build detected as:', target_build, '\nProceeding with reference harmonisation.\n')
    sink()

    # Define the autosomes to iterate over.
    chrs<-1:22
    # Initialize an empty variable to store the harmonized GWAS results across all chromosomes.
    GWAS_matched<-NULL
    # Loop through each chromosome (1 to 22).
    for(i in chrs){
      # Print the current chromosome number being processed to the console.
      print(i)

      # Read the reference panel RDS file for the current chromosome `i`.
      tmp<-readRDS(file = paste0(opt$ref_chr,i,'.rds'))

      # Subset reference data: keep only SNPs (allele lengths = 1).
      tmp<-tmp[nchar(tmp$A1) == 1 & nchar(tmp$A2) == 1,]
      # Rename reference columns to avoid conflicts after merging, prefixing with 'REF.'.
      names(tmp)[names(tmp) == 'CHR']<-'REF.CHR'
      names(tmp)[names(tmp) == 'SNP']<-'REF.SNP'
      names(tmp)[names(tmp) == 'BP_GRCh36']<-'REF.BP_GRCh36'
      names(tmp)[names(tmp) == 'BP_GRCh37']<-'REF.BP_GRCh37'
      names(tmp)[names(tmp) == 'BP_GRCh38']<-'REF.BP_GRCh38'
      # Rename reference frequency column. Assuming input might be REF.FRQ
      names(tmp)[names(tmp) == 'REF.FRQ']<-'REF.FREQ' # Correcting potential typo in original script comment? Assuming it should be FREQ or FRQ consistently.

      # Subset the GWAS data for the current chromosome `i`.
      GWAS_chr<-GWAS[GWAS$CHR == i,]
      # Merge the GWAS subset with the reference panel based on the base pair position of the detected `target_build`.
      # Note: Merge happens on 'ORIGBP' from GWAS and 'REF.BP_build' from reference.
      ref_target<-merge(GWAS_chr, tmp, by.x='ORIGBP', by.y=paste0('REF.BP_',target_build))

      # Identify SNPs that are likely on opposite strands based on IUPAC codes
      # (e.g., R (A/G) in GWAS matches Y (C/T) in reference, suggesting a strand flip).
      flipped<-ref_target[(ref_target$IUPAC.x == 'R' & ref_target$IUPAC.y == 'Y') |
                            (ref_target$IUPAC.x == 'Y' & ref_target$IUPAC.y == 'R') |
                            (ref_target$IUPAC.x == 'K' & ref_target$IUPAC.y == 'M') |
                            (ref_target$IUPAC.x == 'M' & ref_target$IUPAC.y == 'K'),]

      # If any flipped SNPs were identified:
      if(nrow(flipped) > 0) {
        # Change the alleles (A1.x, A2.x from GWAS) in the 'flipped' subset to their complements.
        flipped$A1.x<-snp_allele_comp(flipped$A1.x)
        flipped$A2.x<-snp_allele_comp(flipped$A2.x)

        # Update the IUPAC codes (IUPAC.x from GWAS) for the flipped SNPs based on the complemented alleles.
        flipped$IUPAC.x<-snp_iupac(flipped$A1.x, flipped$A2.x)
      }

      # Identify SNPs where the (potentially flipped) GWAS IUPAC code now matches the reference IUPAC code.
      matched<-ref_target[ref_target$IUPAC.x == ref_target$IUPAC.y,]
      # Combine the directly matched SNPs with the strand-corrected (flipped) SNPs.
      matched<-rbind(matched, flipped)

      # If the GWAS A1 allele (A1.x) does not match the reference A1 allele (A1.y),
      # it means the alleles are swapped relative to the reference frequency.
      # In this case, flip the reference allele frequency (REF.FREQ) to represent the frequency of A1.x.
      # Note: This assumes REF.FREQ corresponds to A1.y in the reference data.
      matched$REF.FREQ[matched$A1.x != matched$A1.y]<-1-matched$REF.FREQ[matched$A1.x != matched$A1.y]

      # Clean up columns after harmonization:
      # Keep the (potentially flipped) GWAS alleles as the final A1 and A2.
      matched$A1<-matched$A1.x
      matched$A1.y<-NULL # Remove reference A1
      matched$A1.x<-NULL # Remove original GWAS A1
      matched$A2<-matched$A2.x
      matched$A2.y<-NULL # Remove reference A2
      matched$A2.x<-NULL # Remove original GWAS A2
      # Remove IUPAC columns used for matching.
      matched$IUPAC.y<-NULL
      matched$IUPAC.x<-NULL
      # Use the reference SNP identifier as the primary SNP ID.
      matched$SNP<-matched$REF.SNP
      matched$REF.SNP<-NULL # Remove the prefixed reference SNP column.
      # Use the reference chromosome identifier.
      matched$CHR<-matched$REF.CHR
      matched$REF.CHR<-NULL # Remove the prefixed reference CHR column.
      # Keep the original GWAS base pair position as the final BP.
      matched$BP<-matched$ORIGBP
      matched$ORIGBP<-NULL # Remove the original GWAS BP column name.
      # Remove all reference BP columns using pattern matching (data.table syntax with=F).
      matched<-matched[,!grepl('REF.BP_GRCh',names(matched)), with=F]

      # Assign the processed data for this chromosome (though variable name is confusing).
      GWAS_matched_chr<-rbind(matched) # rbind seems redundant here if 'matched' already contains all needed rows.
      # Append the harmonized data for the current chromosome to the overall results.
      GWAS_matched<-rbind(GWAS_matched, GWAS_matched_chr)
    } # End of chromosome loop
   } else { # If !is.na(target_build) is FALSE
      # Log a warning if the genome build could not be determined using CHR/BP.
      sink(file = paste(opt$output,'.log',sep=''), append = T)
      cat('WARNING: Could not determine genome build based on CHR/BP matching. Attempting RSID matching if available.\n')
      sink()
   }
} # End of if(chr_bp_avail) block

# This block executes if CHR/BP were not available OR if build detection failed using CHR/BP,
# BUT if RSIDs are available.
if(is.na(target_build) & rsid_avail){
  # Redirect output to log file (appending).
  sink(file = paste(opt$output,'.log',sep=''), append = T)
  # Log message indicating that SNP IDs (RSIDs) will be used for merging.
  cat('Using SNP, A1 and A2 to merge with the reference.\n')
  # Stop redirecting output.
  sink()

  # Define the autosomes to iterate over.
  chrs<-c(1:22)

  # Initialize an empty variable to store the harmonized GWAS results.
  GWAS_matched<-NULL
  # Loop through each chromosome (1 to 22).
  for(i in chrs){
    # Print the current chromosome number being processed.
    print(i)

    # Read the reference panel RDS file for the current chromosome `i`.
    tmp<-readRDS(file = paste0(opt$ref_chr,i,'.rds'))
    # Alternative path commented out: #tmp<-readRDS(file = paste0('reference_panels/forTEDS.chr',i,'.rds'))

    # Subset reference data: keep only SNPs (allele lengths = 1).
    tmp<-tmp[nchar(tmp$A1) == 1 & nchar(tmp$A2) == 1,]
    # Rename reference columns to avoid conflicts, prefixing with 'REF.'.
    names(tmp)[names(tmp) == 'CHR']<-'REF.CHR'
    names(tmp)[names(tmp) == 'BP_GRCh36']<-'REF.BP_GRCh36'
    names(tmp)[names(tmp) == 'BP_GRCh37']<-'REF.BP_GRCh37'
    names(tmp)[names(tmp) == 'BP_GRCh38']<-'REF.BP_GRCh38'
    # Rename reference frequency column.
    names(tmp)[names(tmp) == 'REF.FRQ']<-'REF.FREQ' # Assuming consistency needed.

    # Merge the GWAS data (entire table, not subsetted by chr here) with the reference panel based on the 'SNP' (RSID) column.
    # Note: This merges the *entire* GWAS table with *each* chromosome's reference. This might be inefficient
    # if the GWAS table is large and not pre-filtered by chromosome. Assumes RSIDs are unique across chromosomes.
    ref_target<-merge(GWAS, tmp, by='SNP')

    # Identify SNPs that are likely on opposite strands based on IUPAC codes (same logic as before).
    flipped<-ref_target[(ref_target$IUPAC.x == 'R' & ref_target$IUPAC.y == 'Y') |
                          (ref_target$IUPAC.x == 'Y' & ref_target$IUPAC.y == 'R') |
                          (ref_target$IUPAC.x == 'K' & ref_target$IUPAC.y == 'M') |
                          (ref_target$IUPAC.x == 'M' & ref_target$IUPAC.y == 'K'),]

    # If any flipped SNPs were identified:
    if(nrow(flipped) > 0) {
       # Change the GWAS alleles (A1.x, A2.x) in the 'flipped' subset to their complements.
       flipped$A1.x<-snp_allele_comp(flipped$A1.x)
       flipped$A2.x<-snp_allele_comp(flipped$A2.x)

       # Update the GWAS IUPAC codes (IUPAC.x) for the flipped SNPs.
       flipped$IUPAC.x<-snp_iupac(flipped$A1.x, flipped$A2.x)
    }

    # Identify SNPs where the (potentially flipped) GWAS IUPAC code now matches the reference IUPAC code.
    matched<-ref_target[ref_target$IUPAC.x == ref_target$IUPAC.y,]
    # Combine the directly matched SNPs with the strand-corrected (flipped) SNPs.
    matched<-rbind(matched, flipped)

    # Flip the reference allele frequency (REF.FREQ) if GWAS A1 (A1.x) doesn't match reference A1 (A1.y).
    matched$REF.FREQ[matched$A1.x != matched$A1.y]<-1-matched$REF.FREQ[matched$A1.x != matched$A1.y]

    # Clean up columns after harmonization:
    # Keep the (potentially flipped) GWAS alleles as the final A1 and A2.
    matched$A1<-matched$A1.x
    matched$A1.y<-NULL
    matched$A1.x<-NULL
    matched$A2<-matched$A2.x
    matched$A2.y<-NULL
    matched$A2.x<-NULL
    # Remove IUPAC columns.
    matched$IUPAC.y<-NULL
    matched$IUPAC.x<-NULL
    # Use the reference BP (defaulting to GRCh37 here) as the final BP.
    # This assumes the reference panel primarily uses GRCh37 coordinates if merging by RSID.
    matched$BP<-matched$REF.BP_GRCh37
    # Remove all reference BP columns.
    matched$REF.BP_GRCh37<-NULL
    matched$REF.BP_GRCh38<-NULL
    matched$REF.BP_GRCh36<-NULL
    # Remove the original BP column from GWAS if it existed.
    matched$ORIGBP<-NULL
    # Use the reference chromosome identifier.
    matched$CHR<-matched$REF.CHR
    # Remove the prefixed reference CHR column.
    matched$REF.CHR<-NULL

    # Assign processed data for this chromosome merge (again, rbind seems redundant).
    GWAS_matched_chr<-rbind(matched)
    # Append the harmonized data for the current chromosome merge to the overall results.
    GWAS_matched<-rbind(GWAS_matched, GWAS_matched_chr)
  } # End of chromosome loop
} # End of if(is.na(target_build) & rsid_avail) block

# Check if harmonization actually produced results
if (!exists("GWAS_matched") || is.null(GWAS_matched) || nrow(GWAS_matched) == 0) {
    # Log an error if no variants could be matched to the reference
    sink(file = paste(opt$output,'.log',sep=''), append = T)
    cat('ERROR: No variants could be matched to the reference panel using either CHR/BP or RSID methods. Exiting.\n')
    if(is.na(target_build) & !rsid_avail & !chr_bp_avail) {
        cat('Reason: Neither sufficient CHR/BP information nor RSIDs were found in the input file.\n')
    } else if (is.na(target_build) & chr_bp_avail) {
         cat('Reason: Genome build could not be determined from CHR/BP, and RSID matching was either not possible or yielded no results.\n')
    } else {
         cat('Reason: Merging failed, possibly due to incorrect reference file paths, format mismatch, or no overlapping variants.\n')
    }
    sink()
    # Stop the script execution
    stop("Harmonization failed. Check log file for details.")
}


# Overwrite the original GWAS data table with the harmonized results.
GWAS<-GWAS_matched
# Clear the intermediate variable to save memory.
rm(GWAS_matched)

# Redirect output to log file (appending).
sink(file = paste(opt$output,'.log',sep=''), append = T)
# Log the number of variants remaining after the harmonization process.
cat('After matching variants to the reference,',dim(GWAS)[1],'variants remain.\n')
# Stop redirecting output.
sink()


#####
# Filter: Remove SNPs with INFO score below the specified threshold (--info).
#####

# Checks if an 'INFO' column exists in the harmonized GWAS data.
if(sum(names(GWAS) == 'INFO') == 1){
  # Filters the GWAS data, keeping only rows where INFO is greater than or equal to the threshold.
  GWAS<-GWAS[GWAS$INFO >= opt$info,]

  # Redirect output to log file (appending).
  sink(file = paste(opt$output,'.log',sep=''), append = T)
  # Log the number of variants remaining after the INFO score filtering.
  cat('After removal of SNPs with INFO < ',opt$info,', ',dim(GWAS)[1],' variants remain.\n', sep='')
  # Stop redirecting output.
  sink()
} else { # If 'INFO' column is not found
  # Redirect output to log file (appending).
  sink(file = paste(opt$output,'.log',sep=''), append = T)
  # Log a message indicating that INFO filtering was skipped.
  cat('INFO column is not present, skipping INFO score filtering.\n', sep='')
  # Stop redirecting output.
  sink()
}

#####
# Filter: Remove SNPs with reported MAF (Minor Allele Frequency) below the threshold (--maf).
# Assumes the frequency column is named 'FREQ' and represents the A1 allele frequency.
#####

# Checks if a 'FREQ' column (representing allele frequency from the input GWAS) exists.
if(sum(names(GWAS) == 'FREQ') == 1){
  # Filters the GWAS data, keeping rows where FREQ is >= maf AND <= (1-maf).
  # This ensures the minor allele frequency (min(FREQ, 1-FREQ)) is >= maf.
  GWAS<-GWAS[GWAS$FREQ >= opt$maf & GWAS$FREQ <= (1-opt$maf),]

  # Redirect output to log file (appending).
  sink(file = paste(opt$output,'.log',sep=''), append = T)
  # Log the number of variants remaining after filtering based on reported MAF.
  cat('After removal of SNPs with reported MAF < ',opt$maf,', ',dim(GWAS)[1],' variants remain.\n', sep='')
  # Stop redirecting output.
  sink()
} else { # If 'FREQ' column is not found
  # Redirect output to log file (appending).
  sink(file = paste(opt$output,'.log',sep=''), append = T)
  # Log a message indicating that reported MAF filtering was skipped.
  cat('Reported frequency column (FREQ) is not present, skipping reported MAF filtering.\n', sep='')
  # Stop redirecting output.
  sink()
}

#####
# Filter: Remove SNPs with reference panel MAF below the threshold (--maf).
# Assumes the reference frequency column added during harmonization is 'REF.FREQ'.
#####

# Checks if a 'REF.FREQ' column (reference allele frequency) exists.
if(sum(names(GWAS) == 'REF.FREQ') == 1){ # Corrected check from 'FREQ' to 'REF.FREQ' based on harmonization code
  # Filters the GWAS data, keeping rows where REF.FREQ is >= maf AND <= (1-maf).
  GWAS<-GWAS[GWAS$REF.FREQ >= opt$maf & GWAS$REF.FREQ <= (1-opt$maf),]

  # Redirect output to log file (appending).
  sink(file = paste(opt$output,'.log',sep=''), append = T)
  # Log the number of variants remaining after filtering based on reference MAF.
  cat('After removal of SNPs with reference MAF < ',opt$maf,', ',dim(GWAS)[1],' variants remain.\n', sep='') # Corrected log message
  # Stop redirecting output.
  sink()
} else { # If 'REF.FREQ' column is not found
  # Redirect output to log file (appending).
  sink(file = paste(opt$output,'.log',sep=''), append = T)
  # Log a message indicating that reference MAF filtering was skipped.
  cat('Reference frequency column (REF.FREQ) is not present, skipping reference MAF filtering.\n', sep='') # Corrected log message
  # Stop redirecting output.
  sink()
}

#####
# Filter: Remove SNPs with discordant MAF between reported ('FREQ') and reference ('REF.FREQ').
#####

# Checks if both 'FREQ' and 'REF.FREQ' columns exist.
if(sum(names(GWAS) == 'FREQ') == 1 & sum(names(GWAS) == 'REF.FREQ') == 1){ # Added check for REF.FREQ
  # Calculates the absolute difference between the reported frequency and the reference frequency.
  GWAS$diff<-abs(GWAS$FREQ-GWAS$REF.FREQ)

  # Create a scatter plot comparing reported vs reference frequencies for SNPs that EXCEED the difference threshold.
  # Save the plot as a PNG file.
  # Note: Uses bitmap() which might require specific graphics devices/packages (e.g., ghostscript) depending on the OS. Consider png() instead for wider compatibility.
  bitmap(paste0(opt$output,'.MAF_plot.png'), unit='px', res=300, width=1200, height=1200)
  # Plots REF.FREQ vs FREQ for discordant SNPs.
  plot(GWAS$REF.FREQ[GWAS$diff > opt$maf_diff],GWAS$FREQ[GWAS$diff > opt$maf_diff], xlim=c(0,1), ylim=c(0,1), xlab='Reference Allele Frequency', ylab='Sumstat Allele Frequency', main='Discordant Allele Frequencies') # Added title
  # Adds a diagonal line (y=x) for reference.
  abline(coef = c(0,1))
  # Closes the graphics device, saving the file.
  dev.off()

  # Filters the GWAS data, keeping only rows where the absolute MAF difference is LESS than the threshold.
  GWAS<-GWAS[GWAS$diff < opt$maf_diff,]
  # Removes the temporary 'diff' column.
  GWAS$diff<-NULL

  # Redirect output to log file (appending).
  sink(file = paste(opt$output,'.log',sep=''), append = T)
  # Log the number of variants remaining after removing SNPs with discordant MAFs. Note: Log message has '<' but code uses '<'. Assuming message should reflect the code.
  cat('After removal of SNPs with absolute MAF difference >= ',opt$maf_diff,', ',dim(GWAS)[1],' variants remain.\n', sep='') # Corrected log message based on code logic
  # Stop redirecting output.
  sink()
} else { # If either 'FREQ' or 'REF.FREQ' column is missing
  # Redirect output to log file (appending).
  sink(file = paste(opt$output,'.log',sep=''), append = T)
  # Log a message indicating that MAF discordance check was skipped.
  cat('Reported (FREQ) or Reference (REF.FREQ) MAF column is not present, so discordance cannot be determined.\n', sep='')
  # Stop redirecting output.
  sink()
}

#####
# Filter: Remove SNPs with out-of-bounds p-values (P > 1 or P <= 0).
#####

# Filters the GWAS data, keeping only rows where the P-value is strictly greater than 0 and less than or equal to 1.
GWAS<-GWAS[GWAS$P <= 1 & GWAS$P > 0,]

# Redirect output to log file (appending).
sink(file = paste(opt$output,'.log',sep=''), append = T)
# Log the number of variants remaining after removing SNPs with invalid P-values.
cat('After removal of SNPs with out-of-bound P values, ',dim(GWAS)[1],' variants remain.\n', sep='')
# Stop redirecting output.
sink()

#####
# Filter: Remove SNPs with duplicated rs numbers (SNP IDs).
#####

# Identifies SNP IDs that appear more than once in the 'SNP' column.
dups<-GWAS$SNP[duplicated(GWAS$SNP)]
# Filters the GWAS data, removing all rows where the SNP ID is in the list of duplicates.
GWAS<-GWAS[!(GWAS$SNP %in% dups),]

# Redirect output to log file (appending).
sink(file = paste(opt$output,'.log',sep=''), append = T)
# Log the number of variants remaining after removing duplicate SNPs.
cat('After removal of SNPs with duplicate IDs, ',dim(GWAS)[1],' variants remain.\n', sep='')
# Stop redirecting output.
sink()

#####
# Filter: Remove SNPs with sample size (N) < 3 Standard Deviations from the median N.
# This helps remove variants with unusually small or large sample sizes, which might indicate errors.
#####

# Checks if the 'N' column exists and has more than one unique value (filtering is meaningless otherwise).
if('N' %in% names(GWAS) && length(unique(na.omit(GWAS$N))) > 1){ # Added check for N column existence and made NA-robust
  # Calculates the standard deviation of the sample sizes (N), ignoring NAs.
  N_sd<-sd(GWAS$N, na.rm=TRUE)
  # Calculates the median sample size, ignoring NAs.
  N_median <- median(GWAS$N, na.rm=TRUE)
  # Defines the lower and upper bounds for acceptable N.
  N_lower_bound <- N_median - (3 * N_sd)
  N_upper_bound <- N_median + (3 * N_sd)
  # Filters the GWAS data, keeping only rows where N is within 3 SD of the median. Handles potential NAs in N column.
  GWAS<-GWAS[!is.na(GWAS$N) & GWAS$N < N_upper_bound & GWAS$N > N_lower_bound,]

  # Redirect output to log file (appending).
  sink(file = paste(opt$output,'.log',sep=''), append = T)
  # Log the number of variants remaining after filtering based on sample size deviation.
  cat('After removal of SNPs with N > ', N_upper_bound,' or < ', N_lower_bound,', ',dim(GWAS)[1],' variants remain.\n', sep='')
  # Stop redirecting output.
  sink()
} else { # If 'N' column is not present or has only one unique value
  # Redirect output to log file (appending).
  sink(file = paste(opt$output,'.log',sep=''), append = T)
  # Log a message indicating that N filtering was skipped.
  cat('N column is not present, is invariant, or contains only NAs. Skipping N filtering.\n', sep='') # Updated message
  # Stop redirecting output.
  sink()
}

#####
# Check for evidence of Genomic Control application and potentially recalculate P-values.
# Genomic control inflates standard errors (SE) to control for test statistic inflation.
# This check tries to reverse it if P-values seem inconsistent with BETA/OR and SE.
#####

# Check if Standard Error ('SE') column exists.
if(sum(names(GWAS) == 'SE') == 1){
  # Check if Odds Ratio ('OR') column exists.
  if(sum(names(GWAS) == 'OR') == 1){
    # Calculate BETA (log odds ratio) from OR. Ignore warnings for OR <= 0.
    suppressWarnings(GWAS$BETA<-log(GWAS$OR))
    # Calculate Z-score based on the calculated BETA and reported SE.
    GWAS$Z<-GWAS$BETA/GWAS$SE
    # Recalculate P-value (two-tailed) from the Z-score.
    GWAS$P_check<-2*pnorm(-abs(GWAS$Z))
    # Remove temporary Z and BETA columns.
    GWAS$Z<-NULL
    GWAS$BETA<-NULL # Keep BETA if originally present? Script removes it here.

    # Compare the mean of the original P-values with the mean of the recalculated P-values.
    # Use only non-NA values for comparison. Threshold is 0.01 difference in means.
    if(abs(mean(GWAS$P[!is.na(GWAS$P_check)], na.rm=TRUE) - mean(GWAS$P_check[!is.na(GWAS$P_check)], na.rm=TRUE)) > 0.01){
      # If the means differ substantially, replace the original P with the recalculated P.
      GWAS$P<-GWAS$P_check
      # Remove the temporary P_check column.
      GWAS$P_check<-NULL

      # Log that genomic control was likely detected and P-values were recomputed.
      sink(file = paste(opt$output,'.log',sep=''), append = T)
      cat('Potential genomic control detected based on OR/SE. P-value recomputed.\n', sep='')
      sink()
    } else { # If the means are similar
      # Log that genomic control was not detected.
      sink(file = paste(opt$output,'.log',sep=''), append = T)
      cat('Genomic control was not detected based on OR/SE consistency.\n', sep='')
      sink()
      # Remove the temporary P_check column.
      GWAS$P_check<-NULL
    }
  } else if(sum(names(GWAS) == 'BETA') == 1){ # Check if BETA column exists (if OR didn't)
    # Calculate Z-score directly from BETA and SE.
    GWAS$Z<-GWAS$BETA/GWAS$SE
    # Recalculate P-value (two-tailed) from the Z-score.
    GWAS$P_check<-2*pnorm(-abs(GWAS$Z))
    # Remove temporary Z column.
    GWAS$Z<-NULL

    # Compare the mean of the original P-values with the mean of the recalculated P-values.
    if(abs(mean(GWAS$P[!is.na(GWAS$P_check)], na.rm=TRUE) - mean(GWAS$P_check[!is.na(GWAS$P_check)], na.rm=TRUE)) > 0.01){
      # If means differ, replace original P with recalculated P.
      GWAS$P<-GWAS$P_check
      # Remove temporary P_check column.
      GWAS$P_check<-NULL

      # Log detection and recomputation.
      sink(file = paste(opt$output,'.log',sep=''), append = T)
      cat('Potential genomic control detected based on BETA/SE. P-value recomputed.\n', sep='')
      sink()
    } else { # If means are similar
      # Log non-detection.
      sink(file = paste(opt$output,'.log',sep=''), append = T)
      cat('Genomic control was not detected based on BETA/SE consistency.\n', sep='')
      sink()
      # Remove temporary P_check column.
      GWAS$P_check<-NULL
    }
  } else { # If neither OR nor BETA exists (but SE does)
    # Log that GC check isn't possible without effect size.
    sink(file = paste(opt$output,'.log',sep=''), append = T)
    cat('SE column present, but OR or BETA column missing. Cannot check for genomic control.\n', sep='')
    sink()
  }
} else { # If SE column is not present
  # Log that GC check isn't possible without SE.
  sink(file = paste(opt$output,'.log',sep=''), append = T)
  cat('SE column is not present, genomic control cannot be detected or corrected.\n', sep='')
  sink()
}

#####
# Insert SE (Standard Error) column if it's missing.
# Calculates SE from P-value and effect size (BETA or OR).
#####

# Checks if the 'SE' column is missing (count is 0).
if(sum(names(GWAS) == 'SE') == 0){
  # Check if 'BETA' column exists.
  if(sum(names(GWAS) == 'BETA') == 1){
    # Calculate Z-score from P-value (absolute value of inverse CDF of P/2).
    GWAS$Z<-abs(qnorm(GWAS$P/2))
    # Calculate SE as absolute value of BETA / Z. Handle potential division by zero if Z is 0 (P=1).
    GWAS$SE<-abs(GWAS$BETA/GWAS$Z)
    # Set SE to NA or Inf where Z was 0 or BETA was NA. Division by zero results in Inf. P=1 results in Z=0.
    GWAS$SE[is.infinite(GWAS$SE) | is.na(GWAS$SE)] <- NA
  } else if (sum(names(GWAS) == 'OR') == 1){ # Else, check if 'OR' column exists.
     # Calculate Z-score from P-value.
    GWAS$Z<-abs(qnorm(GWAS$P/2))
     # Calculate SE as absolute value of log(OR) / Z. Handles OR <= 0, NAs, division by zero.
    suppressWarnings(logOR <- log(GWAS$OR)) # Calculate log(OR), suppress warnings for OR<=0
    GWAS$SE<-abs(logOR / GWAS$Z)
    # Set SE to NA where calculation failed (OR<=0, P=1, NAs).
    GWAS$SE[is.infinite(GWAS$SE) | is.na(GWAS$SE) | is.nan(GWAS$SE)] <- NA # Add NaN check
  } else {
     # Log if SE cannot be calculated because effect size is missing
     sink(file = paste(opt$output,'.log',sep=''), append = T)
     cat('SE column missing, but cannot calculate it as neither BETA nor OR columns are present.\n', sep='')
     sink()
  }

  # If SE was successfully calculated in one of the branches above:
  if ('SE' %in% names(GWAS)) {
      # Remove the temporary Z-score column.
      GWAS$Z<-NULL

      # Log that the SE column was inserted.
      sink(file = paste(opt$output,'.log',sep=''), append = T)
      cat('SE column inserted based on available effect size (BETA/OR) and P-value.\n', sep='')
      sink()
  }
}

#####
# Filter: Remove SNPs with SE == 0 or SE is NA (often indicates issues).
#####

# Check if SE column exists before filtering
if ('SE' %in% names(GWAS)) {
    # Get number of rows before filtering
    n_before_se_filter <- nrow(GWAS)
    # Filters the GWAS data, keeping only rows where SE is not equal to 0 and is not NA.
    GWAS<-GWAS[!is.na(GWAS$SE) & GWAS$SE != 0,]
    # Get number of rows after filtering
    n_after_se_filter <- nrow(GWAS)

    # Log the result of the SE filtering.
    sink(file = paste(opt$output,'.log',sep=''), append = T)
    cat('After removal of SNPs with SE == 0 or SE == NA, ', n_after_se_filter,' variants remain (removed ', n_before_se_filter - n_after_se_filter, ').\n', sep='')
    sink()
} else {
    # Log if SE filtering couldn't be done.
    sink(file = paste(opt$output,'.log',sep=''), append = T)
    cat('SE column not present or could not be calculated, skipping removal of SNPs with SE == 0 or NA.\n', sep='')
    sink()
}


####
# Final check: Making sure A1,A2 columns are in uppercase.
# This might be required by downstream tools like LDPred2's snp_match() function.
####
# Convert A1 allele column to uppercase again (might be redundant if done earlier, but safe).
GWAS$A1 <- toupper(GWAS$A1)
# Convert A2 allele column to uppercase again.
GWAS$A2 <- toupper(GWAS$A2)

# Log this final conversion step.
sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Ensured A1 and A2 alleles are in uppercase for downstream compatibility.\n', sep='') # Improved message
sink()

#####
# Write out the cleaned and harmonized results.
#####

# Check if an output file with .gz extension already exists and remove it to prevent fwrite errors/appending.
if(file.exists(paste0(opt$output,'.gz'))){
  system(paste0(paste0('rm ',opt$output,'.gz')))
}
# Check if an output file without .gz extension already exists and remove it.
if(file.exists(paste0(opt$output))){
  system(paste0(paste0('rm ',opt$output)))
}

# Check the value of the --gz command-line option.
if(opt$gz == T){
  # If TRUE, write the final GWAS data table to a gzipped file using fwrite for speed. Tab-separated.
  fwrite(GWAS, paste0(opt$output,'.gz'), sep='\t')
} else {
  # If FALSE, write the final GWAS data table to a plain text file. Tab-separated.
  fwrite(GWAS, opt$output, sep='\t')
}

# Record the end time of the script execution.
end.time <- Sys.time()
# Calculate the total time taken for the script to run.
time.taken <- end.time - start.time
# Redirect output to the log file (appending).
sink(file = paste(opt$output,'.log',sep=''), append = T)
# Log the script's finish time.
cat('Analysis finished at',as.character(end.time),'\n')
# Log the total duration of the analysis, rounding to 2 decimal places.
cat('Analysis duration was',as.character(round(time.taken,2)),attr(time.taken, 'units'),'\n')
# Stop redirecting output - script finishes here.
sink()

# Implicitly exits Rscript.
```
