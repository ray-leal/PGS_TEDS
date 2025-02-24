#! /usr/bin/env Rscript
.libPaths()
library(optparse) 	# handles command-line arguments; getopt must be installed!
library(ggplot2) 	# used for visualisation - previously q-plot was used - ggplot2 = customisation but needs script updating to map variables to aesthetics (aes)
library(bigsnpr)	# handles the large genetic datasets (SNP arrays) - bigstatsr is required!
library(bigreadr)	# efficiently reads the big data text files, based on splitting "data.table::fread" command
library(data.table)	# (see above)

# --------------------------------------------------------------------------------------------
# REFERENCES:
# -------------------------------------------------------------------------------------------
# this scripts uses info from the following resources:
# recommended GWAS QC
# see 		https://github.com/privefl/paper-misspec/tree/main/code/prepare-sumstats
# see 		https://privefl.github.io/bigsnpr/articles/LDpred2.html
# see also  https://github.com/privefl/paper-infer/blob/main/code/example-with-provided-LD.R
# ew LD ref - hap-map + 

# --------------------------------------------------------------------------------------------
# DEFINE INPUT OPTIONS
# --------------------------------------------------------------------------------------------
# uses optparse package in R to allow passing arguments:
# NOTE: as of Feb 2025 filepath, updated to recovery folder for now; may need to revert if original file/folder location restored...

# -- sumstats: describes required columns based on trait type, 
# -- 	 geno: specifies that a .rds file format is needed 
# -- 	 type: indicates whether the trait is case/control (TRUE) or continuous (FALSE) - assumes case/control unless otherwise specified
# -- 	  out: specifies where the results will be stored
# -- 	  dir: specifies the directory containing the GWAS summary stats ***this has been specified for MDD in this test
# -- 	 misc: specifies the directory containing the LD reference data (hapmap3plus)
# -- 	cores: number of cores to use - 16 to 32 cores are best for PGS 

option_list = list(
  make_option(c("-s", "--sumstats"), type="character", default=NULL, 
              help=
                "Name of GWAS summary statistics.\n
                Note sumstats should have *at least* the following header:\n
                case/control traits: CHR BP A2 A1 NCAS NCON BETA SE \n
                continuous traits: CHR BP A2 A1 N BETA SE\n", metavar="character"),
  make_option(c("-g", "--geno"), type="character", 
              default="/cephfs/volumes/hpc_data_prj/teds/1971f5b8-1bb5-4a5f-8cb1-d27d54192a50/recovered/affy_oee_genotypes_final_sample_inclDZs/LDpred2geno/out/fullData_info75_maf005_geno02_mind02_HWE00001.rds", 
              help="(Path to) Genetic dataset in .rds format.\n 
              [default = %default]", metavar="character"),
  make_option(c("-t", "--type"), type="logical", 
              default = TRUE,
              help="Whether GWAS trait is case/control.\n 
              [default = %default]", metavar="logical"),
  make_option(c("-o", "--out"), type="character", 
              default="/cephfs/volumes/hpc_data_prj/teds/1971f5b8-1bb5-4a5f-8cb1-d27d54192a50/recovered/affy_oee_genotypes_final_sample_inclDZs/LDpred2geno/out/test/", 
              help="(Path to) Output directory.\n 
              [default = %default]", metavar="character"),
  make_option(c("-d", "--Sdir"), type="character", 
              default="/cephfs/volumes/hpc_data_prj/teds/1971f5b8-1bb5-4a5f-8cb1-d27d54192a50/recovered/affy_oee_genotypes_final_sample_inclDZs/SumStats_repository/mdd2025_no23andMe_eur_v3-49-24-11", 
              help="Sumstats directory.\n 
              [default = %default]", metavar="character"),
  make_option(c("-m", "--misc"), type="character", 
              default= "/cephfs/volumes/hpc_data_prj/teds/1971f5b8-1bb5-4a5f-8cb1-d27d54192a50/recovered/affy_oee_genotypes_final_sample_inclDZs/LDpred2geno/hapmap3plus/", 
              help="(Path to) Directory including LD reference info file(.rds) with ld scores, and LD matrices by chromosome.\n
              [default= %default]\n
              Note the directory should have  the following structure:\n 
              hapmap3plus
              ├── LDref
              │   ├── LD_with_blocks_chr1.rds
              │   ├── LD_with_blocks_chr2.rds
              │   ├── LD_with_blocks_chr3.rds
              etc...
              └── map_hm3_plus.rds", metavar="character"),
  make_option(c("-c", "--cores"), type="integer", 
              default= 32, 
              help="Number of cores. \n 
              [default = %default]", metavar="number")
); 

# create a parser using the options above, extracts the user provided values and stores values in 'opt'
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#print help
print_help(opt_parser)

#print options
print(opt)

out_path = opt$out
misc_path = opt$misc
NCORES = opt$cores
geno = opt$geno
sumstatDir = opt$Sdir

# --------------------------------------------------------------------------------------------
# Setup of logging file
# --------------------------------------------------------------------------------------------
# logs metadata e.g. start time, purpose of script, author etc.

start <- Sys.time() #start time

file_log <- paste0(opt$s,".log") 

file.create(file_log) #create log file

cat(" -------------------------------------------------","\n",
    " Generate LDpred2-auto and infinitesimal scores. ","\n",
    " bugs and questions: a.allegrini@ucl.ac.uk ","\n",
    "-------------------------------------------------","\n",
    " ","\n", sep = " ", file=file_log, append=TRUE)

# print help and  options
# redirects console output to the log log file set up above
sink(file = file_log, append = TRUE)
print_help(opt_parser)

cat("#####","SELECTED OPTIONS:","#####\n"," ", "\n", file=file_log, append=TRUE)
opt
cat("##########\n"," ", "\n", file=file_log, append=TRUE)
sink()

cat(paste0("Analyses started at ", start),"\n"," ","\n", file=file_log,append=TRUE)

cat("Reading: ", opt$s, " sumstats.", "\n"," ","\n", sep='',file=file_log,append=TRUE) 

# --------------------------------------------------------------------------------------------
# Read and process the summary statistics using fread2 command (reads text files)
# --------------------------------------------------------------------------------------------
# load sumstats and convert to LDpred header format + calculate effective sample size 
# reads the GWAS sumstats - looking for columns/headings called: 

	# CHR (chromosome) 
	# BP (base pair position) 
	# A1 (EFFECT allele)
	# A2 (NON-EFFECT allele) 
	# BETA (effect size) 
	# SE (std error) 
	# MAF (minor allele frequency)

# load in the sumstats table using fread2 function from the bigreadr package
sumstats <- bigreadr::fread2(input = paste0(sumstatDir,opt$s))
# write the information about the dimensions of the loaded table to the log file
cat("Loaded sumstats have: ", dim(sumstats)[1], " rows and ", dim(sumstats)[2], " columns.", "\n"," ","\n", sep="",file=file_log,append=TRUE) 
# write the column headers of the loaded table to the log file
cat("Sumstats header is: ", paste(names(sumstats), collapse = " ", sep = " "), "\n"," ","\n", sep="",file=file_log,append=TRUE) 

# --------------------------------------------------------------------------------------------
# basic QC 
# --------------------------------------------------------------------------------------------
# to add discard sample N whe .6 < .9 quant 
# also add per variant effective sample size if more than .5 discarded
# https://github.com/privefl/bigsnpr/issues/281
# All operations are logged to file_log to track what changes were made to the data 

cat("Starting sumstats QC...","\n"," ","\n", file=file_log,append=TRUE)  

# make sure MAF is actually MAF (i.e., max value is .5 or less) then filter out anything > .01
if("MAF" %in% colnames(sumstats)) {
  sumstats$MAF <- ifelse(sumstats$MAF <= .5, sumstats$MAF, (1-sumstats$MAF))
  cat("N = ", sum(sumstats$MAF < 0.01), "have been discarded because MAF < 0.01 ","\n", sep='',file=file_log,append=TRUE)
  sumstats <-  sumstats[sumstats$MAF >= 0.01 , ]
  
}else{
  cat("No MAF column provided.\n",sep='',file=file_log,append=TRUE) }

# If an INFO column exists (quality metric for imputed variants), filter on INFO = Removes SNPs with INFO < 0.6 (low-quality variants) 
if("INFO" %in% colnames(sumstats)) {
    cat("N = ", sum(sumstats$INFO < 0.6), "have been discarded because INFO < 0.6 ","\n", sep='',file=file_log,append=TRUE)
    sumstats <-  sumstats[sumstats$INFO >= 0.6,]
  }else{
  cat("No INFO column provided.\n",sep='',file=file_log,append=TRUE) }

# make sure all alleles are upper case
sumstats$A1 <- factor(toupper(sumstats$A1), c("A", "C", "G", "T"))
sumstats$A2 <- factor(toupper(sumstats$A2), c("A", "C", "G", "T"))

# effective sample size - distinquish between case control vs continuous trait
if(opt$type == T){ 
  
#check colummn names are correct 1 
# --------------------------------------------------------------------------------------------
if (!all(c("CHR","BP","A2","A1","NCAS","NCON","BETA","SE") %in% names(sumstats))) {
    stop(paste0("Sumstats header is not correct. Minium header required for case/control traits is:\n        ",
                paste(c("CHR","BP","A2","A1","NCAS","NCON","BETA","SE"), collapse=", ")))}
  
#select columns and rename (e.g. chr to CHR, pos to BP, a0 to A2, a1 to A1, Ncas to NCAS, Ncon to NCON, beta to BETA, beta_se to SE)
sumstats <- sumstats[,c("CHR","BP","A2","A1","NCAS","NCON","BETA","SE")]
names(sumstats) <- c("chr", "pos", "a0", "a1", "Ncas","Ncon","beta","beta_se")
 
# create and store Ncas and Ncon as variables
Ncas <- sumstats$Ncas 
Ncon <- sumstats$Ncon
  
sumstats$n_eff <- 4 / (1 / Ncas + 1 / Ncon)
  
}else{ #if not
  
# check colummn names are correct 2
# --------------------------------------------------------------------------------------------
  if (!all(c("CHR","BP","A2","A1","N","BETA","SE") %in% names(sumstats))) {
    stop(paste0("Sumstats header is not correct. Minium header required for continuous traits is:\n        ",
                paste(c("CHR","BP","A2","A1","N","BETA","SE"), collapse=", ")))}
  
#select columns and rename (e.g. chr to CHR, pos to BP, a0 to A2, a1 to A1, n_eff to N, beta to BETA, beta_se to SE)
  sumstats <- sumstats[,c("CHR","BP","A2","A1","N","BETA","SE")]
  names(sumstats) <- c("chr", "pos", "a0", "a1", "n_eff","beta","beta_se") 
}

#sometimes this is an issue (checks if chromosome number is stored as integer)
if(!is.integer(sumstats$chr)){
  sumstats$chr <- as.integer(sumstats$chr)
}

# --------------------------------------------------------------------------------------------
# Read in LD reference map (HapMap3+) 
# --------------------------------------------------------------------------------------------
# HapMap3+ is a combined dataset that includes haplotypes from the HapMap 3 project and the 1,000 Genomes Project. 
# The HapMap 3 project includes data from 1,301 samples across 11 populations, providing a diverse set of genetic variations. 
# The 1,000 Genomes Project adds even more genetic diversity, with data from multiple populations around the world. 
# info here: https://www.broadinstitute.org/medical-and-population-genetics/hapmap-3

# LOAD reference map.rds file: 
map_ldref <- readRDS(paste0(misc_path,"map_hm3_plus.rds")) 

# 			(snp_match) SNPs matching to Reference Panel 			 
# (return_flip_and_rev) track which SNPs need flipping or reversing
# returns in Tibble format (updated data frame format - do not change varname or type, indicate when variable is not present etc.)
(info_snp <- tibble::as_tibble(snp_match(sumstats, map_ldref, return_flip_and_rev = T)))

cat("There were: ", sum( vctrs::vec_duplicate_detect(sumstats[, c("chr","pos")]))," duplicated physical positions in GWAS data.\n" ," ","\n", sep='',file=file_log,append=TRUE)

# print results for number of SNPs matched, flipped, reversed, also appends to log file:
cat(" N = ", dim(info_snp)[1], "SNPs have been matched with reference data (i.e. HapMap3 + )\n",
    "N =  ",sum(info_snp$`_FLIP_`), "SNPs were flipped\n",
    "N =  ",sum(info_snp$`_REV_`), "were reversed.\n",
    "\n",file=file_log,append=TRUE)

# remove any missing values - drop NAs and make sure order is the same with SD file below
info_snp <- tidyr::drop_na(tibble::as_tibble(info_snp))

# calculate the chi-square stats for each SNP; ignores any missing values in the calculation (na.rm=T)
chi2 <- with(info_snp, (beta / beta_se)^2)
cat("GWAS chi^2 = ",mean(chi2,na.rm=T),".\n",sep='',file=file_log,append=TRUE)  

# --------------------------------------------------------------------------------------------
#SD REFERENCE
# --------------------------------------------------------------------------------------------

# Compute the reference standard deviation for SNP allele dosages (af_UKBB = allele frequency in UK BioBank)
sd_ldref <- with(info_snp, sqrt(2 * af_UKBB * (1 - af_UKBB)))

# process based on trait type (binary = T, continuous = F)
# print message to .log file if trait is binary. 
if(opt$type == T){  #if GWAS trait is binary, following block will run:
  cat("treating GWAS trait as binary.\n",sep='',file=file_log,append=TRUE)  
  sd_ss <- with(info_snp, 2 / sqrt(n_eff * beta_se^2 + beta^2)) #sumstats sd

}else{ #if it is continuous; the code inside this else block runs:

  cat("treating GWAS trait as continuous.\n",sep='',file=file_log,append=TRUE)  
  
  # consistent with e.g.: https://github.com/privefl/paper-misspec/blob/main/code/prepare-sumstats-bbj/height.R
  # although assumes sd(y) = 1 (beta std 1)
  # if not below sd(y) is reestimated 
  
  
# 99th percentile of sd_ss is computed and scaled to sqrt 0.5 (normalising the sd_ss acorss datasets)
  sd_ss = with(info_snp, 1 / sqrt(n_eff * beta_se^2 + beta^2))
  sd_ss = sd_ss / quantile(sd_ss, 0.99) * sqrt(0.5)
}

# Identify SNPs with unusual standard deviation values 
# e.g. under or over estimation vs ldref, sd less than 0.1 is too small/unreliable, sd less than 0.05 as potential rare alelle)   
is_bad <- sd_ss < (0.5 * sd_ldref) | sd_ss > (sd_ldref + 0.1) | sd_ss < 0.1 | sd_ldref < 0.05
cat(" ","\n", file=file_log,append=TRUE)

# --------------------------------------------------------------------------------------------
# Prepare data for plotting
# --------------------------------------------------------------------------------------------

# make (temp) out dir
dir_path <- paste0(out_path, opt$sumstats, "_", format(Sys.time(), "%d%B%Y"))
if (!dir.exists(dir_path)) {
  dir.create(dir_path, recursive = TRUE)
}
tmp <- paste0(out_path,opt$sumstats,'_',format(Sys.time(), '%d%B%Y'),'/')


# Create a proper data frame first
# --------------------------------------------------------------------------------------------
# sd_ldref (sd from alelle freq from LD ref), sd_ss (sd from GWAS summary stats), is_bad (vector marking SNPs as problematic)
plot_data <- data.frame(sd_ldref = sd_ldref, 
                       sd_ss = sd_ss, 
                       is_bad = is_bad)

# Set up ggplot object using plot_data
# --------------------------------------------------------------------------------------------
# Maps sd_ldref to the x-axis, sd_ss to the y-axis, and colors points by is_bad
# theme_bigstatsr = apply a predefined theme
# coord_equal = make sure bothe axes have same scale
#  scale_color_viridis_d = use a colourblind-friendly gradient (from viridis)for the points, reverse default direction
# geom_point(alpha = 0.5) = add points to the scatter plot with 50% transparency
# geom_abline(linetype = 2, color = "red") = add a diagonal reference line (dashed, red)
# labs = adds axis labels, color = "Removed?" = labels the legend to indicate flagged SNPs

p1 <- ggplot(data = plot_data, 
       aes(x = sd_ldref, y = sd_ss, color = is_bad)) +
  theme_bigstatsr() +
  coord_equal() +
  scale_color_viridis_d(direction = -1) +
  geom_point(alpha = 0.5) +  # Add this line to actually plot points
  geom_abline(linetype = 2, color = "red") +
  labs(x = "Standard deviations derived from allele frequencies of the LD reference",
       y = "Standard deviations derived from the summary statistics",
       color = "Removed?")

# save the plot to temp out folder
ggsave(paste0(tmp,"sd_", opt$s,".png"), plot = p1, width = 10, height = 7)

# Filter out bad SNPs from the dataset
df_beta <- info_snp[!is_bad, ] #remove bad SNPs

# Log the number of bad SNPs that were removed
cat(sum(is_bad, na.rm = T)," SNPs are bad.\n"," ","\n", sep='',file=file_log,append=TRUE)

# Log the number of SNPs remaining after QC
cat("After QC there are: ", dim(df_beta)[1], " SNPs.\n"," ","\n",sep='',file=file_log,append=TRUE)

# bit to be fixed/double checked:
# --------------------------------------------------------------------------------------------
# stop if more than 50% variants have discordant SD
# see also https://github.com/privefl/bigsnpr/issues/281

# count number of SNPs flagged as bad - ignoring the NA values; warn user if more than half of the variants have inconsistent standard deviations
if(sum(is_bad, na.rm = T) > (length(is_bad)*0.5)) {
  cat("WARNING: More than half the variants had a discordant SD. Imputing Neff.\n
      Double check your input sumstats: reference population, per-variant sample size, effective N.\n\n", file=file_log, append=TRUE) # see: https://github.com/privefl/bigsnpr/issues/281.

# for binary traits: Computes n_eff_imp (imputed effective sample size) using a formula from the referenced article. This equation adjusts for the variance explained by the SNP.
#see eq 4 and 5 here  https://www.sciencedirect.com/science/article/pii/S2666247722000525?via%3Dihub#sec3.2
  if(opt$type == T){
    info_snp$n_eff_imp <- (4 / sd_ldref^2 - info_snp$beta^2) / info_snp$beta_se^2
	#correcting difference in per allele effect size (if n_eff_imp is below 67% of the 90th percentile, it is set to NA)
    info_snp$n_eff_imp <- with(info_snp, ifelse(n_eff_imp >= .67 * quantile(n_eff_imp,.90), n_eff_imp, NA))  #doi: 10.1093/bioinformatics/btw613
	# Computes the median of imputed Neff while ignoring NA values = serves as a robust estimate of the sample size.
    info_snp$n_eff <- median(info_snp$n_eff_imp, na.rm=T)
	# Recalculate sd_ss using the imputed Neff
    sd_ss <- with(info_snp, 2 / sqrt(n_eff * beta_se^2 + beta^2))
    
  }else{

# for continuous traits: Uses a different formula for n_eff_imp in continuous traits. The denominator (sd_ldref^2 - beta^2) accounts for the expected variance in continuous outcomes.
	info_snp$n_eff_imp <- (1 / (sd_ldref^2 - info_snp$beta^2)) / info_snp$beta_se^2
	#correcting difference in per allele effect size (if n_eff_imp is below 67% of the 90th percentile, it is set to NA)
    info_snp$n_eff_imp <- with(info_snp, ifelse(n_eff_imp >= .67 * quantile(n_eff_imp,.90), n_eff_imp, NA))  
    # Computes the median of imputed Neff while ignoring NA values = serves as a robust estimate of the sample size.
	info_snp$n_eff <- median(info_snp$n_eff_imp, na.rm=T)
    # Uses the 1st percentile of 0.5 * (n_eff * beta_se^2 + beta^2) to estimate sd_y = Ensures consistent scaling of summary statistics.
    sd_y = with(info_snp, sqrt(quantile(0.5 * (n_eff * beta_se^2 + beta^2), 0.01))) # https://github.com/privefl/bigsnpr/issues/349
    # Recalculate sd_ss using the estimated sd_y
	sd_ss = with(info_snp, sd_y / sqrt(n_eff * beta_se^2 + beta^2))
	# Normalize sd_ss based on the 99th percentile
    sd_ss = sd_ss / quantile(sd_ss, 0.99) * sqrt(0.5)
  }
  
 # Counts how many SNPs were set to NA due to low sample size estimates. Logs this quality control (QC) step.
  cat("N = ",sum(is.na(info_snp$n_eff_imp)), " SNPs have been filtered to correct for low sample size.\n"," ","\n",sep='',file=file_log,append=TRUE)
  cat("new iputed median N = ", median(info_snp$n_eff_imp, na.rm=T) ,"\n", file=file_log,append=TRUE) 
  # flags SNPs under or over estimation vs ldref, sd less than 0.1 is too small/unreliable, sd less than 0.05 as potential rare alelle)   
  is_bad <- sd_ss < (0.5 * sd_ldref) | sd_ss > (sd_ldref + 0.1) | sd_ss < 0.1 | sd_ldref < 0.05
  cat(sum(is_bad, na.rm = T)," SNPs are bad.\n"," ","\n", sep='',file=file_log,append=TRUE)
  
# Create a proper data frame again
# --------------------------------------------------------------------------------------------

# sd_ldref: Standard deviation from the LD reference.
# sd_ss: Standard deviation from summary statistics.
# is_bad: Boolean flag indicating if an SNP is discordant.

plot_data_2 <- data.frame(sd_ldref = sd_ldref, 
                         sd_ss = sd_ss, 
                         is_bad = is_bad)

# Uses ggplot2 to create a scatter plot comparing sd_ldref and sd_ss.
# Colors points based on is_bad (whether an SNP was removed).
# geom_abline(linetype = 2, color = "red"): Adds a dashed red line representing equality (y = x), making it easy to spot outliers.
# theme_bigstatsr(): Applies a predefined theme.

p2 <- ggplot(data = plot_data_2, 
       aes(x = sd_ldref, y = sd_ss, color = is_bad)) +
  theme_bigstatsr() +
  coord_equal() +
  scale_color_viridis_d(direction = -1) +
  geom_point(alpha = 0.5) +  # Add this line to actually plot points
  geom_abline(linetype = 2, color = "red") +
  labs(x = "Standard deviations derived from allele frequencies of the LD reference",
       y = "Standard deviations derived from the summary statistics",
       color = "Removed?")

ggsave(paste0(tmp,"sd_", opt$s,".medianN.png"), plot = p2, width = 10, height = 7)

# Filters info_snp, keeping only the "good" SNPs (!is_bad).
  df_beta <- info_snp[!is_bad, ] #remove bad SNPs
  cat("After QC there are: ", dim(df_beta)[1], " SNPs.\n"," ","\n",sep='',file=file_log,append=TRUE)
}

# Load test data
# --------------------------------------------------------------------------------------------
obj.bigsnp <- snp_attach(paste0(geno))
# shortcut for geno test set
G <- obj.bigsnp$genotypes

# Extracts SNP metadata from the test dataset and creates a simplified map (map_test) with:
# chr: Chromosome number.
# pos: Physical position.
# a0, a1: Reference and alternative alleles.

map_test <- dplyr::transmute(obj.bigsnp$map,
                             chr = as.integer(chromosome), pos = physical.pos,
                             a0 = allele2, a1 = allele1)

#remove duplicates (if any) in test data 
dups <- vctrs::vec_duplicate_detect(map_test[, c("chr","pos")])
if (any(dups)) {
  cat("There are: ", sum( vctrs::vec_duplicate_detect(map_test[, c("chr","pos")])), " duplicated physical positions in test data.\n" ," ","\n", sep='',file=file_log,append=TRUE)
  cat("removing duplicated...\n"," ","\n",file=file_log,append=TRUE) 
  map_test <-  map_test[!dups, ] }

# match with variants in test data and prepare files for matrix multiplication later
cat("Matching with test data...\n"," ","\n",file=file_log,append=TRUE)

# Extracts the first four columns of df_beta (containing SNP metadata) and adds a placeholder beta = 1 for later matching.
map_pgs <- df_beta[1:4]; map_pgs$beta <- 1
# Uses snp_match() to match GWAS SNPs (map_pgs) with test dataset SNPs (map_test). return_flip_and_rev = T: Keeps track of strand flips and allele reversals.
map_pgs2 <- snp_match(map_pgs, map_test, return_flip_and_rev = T)

# match sumstats with reference data: log Total SNPs matched, Number of strand-flipped SNPs, Number of allele-reversed SNPs.
cat(" N = ", dim(map_pgs2)[1], "SNPs have been matched with test data\n",
    "N =  ",sum(map_pgs2$`_FLIP_`), "SNPs were flipped\n",
    "N =  ",sum(map_pgs2$`_REV_`), "were reversed.\n",
    "\n",file=file_log,append=TRUE)

# Keep Only SNPs That Exist in Test Data
#in_test <- vctrs::vec_in(df_beta[, c("chr", "pos")], map_test[, c("chr", "pos")]) #pgs (not pgs2) used
in_test <- vctrs::vec_in(df_beta[, c("chr", "pos")], map_pgs2[, c("chr", "pos")])

#cat("in_test vector = TRUE is: ", sum(in_test), " long.\n"," ","\n",sep='',file=file_log,append=TRUE)

# Match Test Data with GWAS Data
df_beta <- df_beta[in_test, ]
cat("There are: ", dim(df_beta)[1], " SNPs in common with test data.\n"," ","\n",sep='',file=file_log,append=TRUE)

# LDSC (Linkage Disequilibrium Score Regression)
# --------------------------------------------------------------------------------------------
# ld: Linkage disequilibrium scores.
# chi2 = (beta / beta_se)^2: Chi-square statistic per SNP.
# sample_size = n_eff: Effective sample size.
# ncores = NCORES: Number of CPU cores used.

cat("Running LDSC...\n"," ","\n",file=file_log,append=TRUE)

(ldsc <- with(df_beta, snp_ldsc(ld, ld_size = nrow(map_ldref),
                                chi2 = (beta / beta_se)^2,
                                sample_size = n_eff,
                                ncores = NCORES)))

# extracts the estimated heritabilty (h2)
h2_est <- ldsc[["h2"]]

# estimates heritability h2 using the LD Score Regression
# --------------------------------------------------------------------------------------------
cat(paste(names(ldsc),ldsc, sep ="=", collapse="; "),"\n"," ","\n", file=file_log, append=TRUE)
cat("Running sparse matrix...\n"," ","\n",file=file_log,append=TRUE)

# --------------------------------------------------------------------------------------------
# Check and remove existing .sbk file if it exists (added new)
# --------------------------------------------------------------------------------------------
sbk_file <- paste0(tmp, "corr_chr.sbk")
if (file.exists(sbk_file)) {
  file.remove(sbk_file)
  cat("Removed existing .sbk file\n", file=file_log, append=TRUE)
}

# --------------------------------------------------------------------------------------------
# LD Matrix to a Sparse matrix format
# --------------------------------------------------------------------------------------------
# LD matrix is a correlation matrix of SNPs and has millions of rows and columns.
# Most SNPs are only correlated with nearby SNPs due to local linkage disequilibrium. This results in a sparse matrix (i.e., most values are zero or near zero). 
# Sparse matrix: Store only the nonzero values, drastically reducing memory usage as needed for LDpred2

# loop through chromosomes 1 to 22
for (chr in 1:22) {
  
  cat(chr, ".. ", sep = "")
  
  ## find SNPs in df_beta
  ind.chr <- which(df_beta$chr == chr)
  ## find numeric IDs from df_beta
  ind.chr2 <- df_beta$`_NUM_ID_`[ind.chr]
  ## Matches SNP in df_beta (ind.chr2) with indices in LD reference (map_ldref).
  ind.chr3 <- match(ind.chr2, which(map_ldref$chr == chr))
 
# now load the LD correlation matrix for the current chromosome. Selects only the SNPs present in both df_beta and map_ldref.
  corr_chr <- readRDS(paste0(misc_path,"LDref/LD_with_blocks_chr", chr, ".rds"))[ind.chr3, ind.chr3]

# Convert the LD Matrix to a Sparse Format
# If chr == 1, creates a Sparse Filebacked Matrix (SFBM) to store LD correlations.
# For later chromosomes, appends their correlation matrices (corr_chr) to the existing SFBM.
  if (chr == 1) {
    corr <- as_SFBM(corr_chr, paste0(tmp,"corr_chr"), compact = TRUE)
  } else {
    corr$add_columns(corr_chr, nrow(corr))
  }
}

# Computes the size of the sparse LD matrix file in gigabytes (GB).
file.size(corr$sbk) / 1024^3 

cat("Running LDpred2-inf...\n"," ","\n",file=file_log,append=TRUE)

# Compute LDpred2-inf, a Bayesian polygenic risk score (PRS) model assuming an infinitesimal genetic architecture.
beta_inf <- snp_ldpred2_inf(corr, df_beta, h2 = h2_est)
# Adds chr, pos, a0, and a1 to the computed beta_inf values.
beta_inf2 <- cbind(df_beta[,1:4], beta_inf) #bind with chr:pos:a0:a1

#save weights
write.table(beta_inf2, paste0(tmp,opt$sumstats,"_beta_inf.txt"), col.names=T,row.names=F,quote=F)

# Computes polygenic scores (PRS) for test individuals:
	# G: Genotype matrix.
	# beta_inf: LDpred2-inf effect sizes.
	# map_pgs2[["_NUM_ID_"]]: Indices of matched SNPs.
	# ncores = NCORES: Parallel processing.

pred_inf <- big_prodVec(G, beta_inf * map_pgs2$beta,
                        ind.col = map_pgs2[["_NUM_ID_"]],
                        ncores = NCORES)

# Combines PRS predictions with test sample metadata:
final_pred_inf <- cbind(obj.bigsnp$fam,pred_inf)

#save final scores
write.table(final_pred_inf, paste0(tmp,opt$sumstats,"_pred_inf.txt"), col.names=T, row.names=F, quote = F)

cat("Running LDpred2-auto...\n"," ","\n", file=file_log,append=TRUE)

# --------------------------------------------------------------------------------------------
# LDpred2-auto
# --------------------------------------------------------------------------------------------

# Runs the LDpred2-auto model, a Bayesian approach that estimates polygenic risk scores (PRS).
# Uses the LD matrix (corr) and summary statistics (df_beta) as input.
# h2_init = h2_est: Uses the estimated heritability (h2_est) as an initial value.
# vec_p_init = seq_log(1e-4, 0.9, 30): Defines a log-spaced sequence of 30 values between 1e-4 and 0.9 for p, which represents the proportion of causal variants.
# allow_jump_sign = FALSE: Prevents large sign changes in effect sizes during estimation.
# shrink_corr = 0.95: Shrinks the correlation matrix to reduce overfitting.
# ncores = NCORES: Runs in parallel using multiple cores.


multi_auto <- snp_ldpred2_auto(corr, df_beta, h2_init = h2_est,
                               vec_p_init = seq_log(1e-4, 0.9, 30),
                               allow_jump_sign = FALSE, shrink_corr = 0.95,
                               ncores = NCORES)

#save weights
saveRDS(multi_auto, file = paste0(tmp,opt$sumstats,"_multi_beta_auto.rds"))

# Create data frames for p_est and h2_est from multi_auto results
auto_ref <- multi_auto[[1]]
df_p <- data.frame(iteration = seq_along(auto_ref$path_p_est), 
                   p_est = auto_ref$path_p_est)
df_h2 <- data.frame(iteration = seq_along(auto_ref$path_h2_est), 
                    h2_est = auto_ref$path_h2_est)

# Use the final estimates from the reference chain for the horizontal lines
p_est_value <- data.frame(yintercept = auto_ref$p_est)
h2_est_value <- data.frame(yintercept = auto_ref$h2_est)

# Check that chains converged
plot_grid(
  ggplot(df_p, aes(x = iteration, y = p_est)) +
    theme_bigstatsr() +
    geom_line() +  
    geom_hline(data = p_est_value, aes(yintercept = yintercept), col = "blue") +
    scale_y_log10() +
    labs(x = "Iteration", y = "p"),

  ggplot(df_h2, aes(x = iteration, y = h2_est)) +
    theme_bigstatsr() +
    geom_line() +  
    geom_hline(data = h2_est_value, aes(yintercept = yintercept), col = "blue") +
    labs(x = "Iteration", y = "h2"),

  ncol = 1, align = "hv"
)

ggsave(paste0(tmp,"auto_chains_", opt$s,".png"), width = 12, height = 10)

# Filter outlier betas and average remaining ones - new recommended way 
# where the range is greater than 95% of the highest quantile (indicating stable estimates).
(range <- sapply(multi_auto, function(auto) diff(range(auto$corr_est))))
(keep <- (range > (0.95 * quantile(range, 0.95))))
final_beta_auto <- rowMeans(sapply(multi_auto[keep], function(auto) auto$beta_est))

#save final betas (final estimated effect size) along with chromosome, position, and alleles.
final_beta_auto2 <- cbind(df_beta[,1:4], final_beta_auto) #bind with chr:pos:a0:a1
write.table(final_beta_auto2, paste0(tmp,opt$sumstats,"_final_beta_auto.txt"), col.names=T,row.names=F,quote=F)

# calculate score (matrix multiplication)
# Uses big_prodVec() for efficient matrix multiplication.
# Only uses SNPs present in both test and reference datasets (map_pgs2[["_NUM_ID_"]]).

final_pred_auto <- big_prodVec(G, final_beta_auto * map_pgs2$beta, 
                               ind.col = map_pgs2[["_NUM_ID_"]],
                               ncores = NCORES)
							   
# Merges the PRS with sample information (fam file).
final_pred_auto <- cbind(obj.bigsnp$fam,final_pred_auto)

# save final scores
write.table(final_pred_auto, paste0(tmp,opt$sumstats,"_pred_auto.txt"), col.names=T, row.names=F, quote = F)

#remove sbk
file.remove(paste0(tmp, "corr_chr.sbk"))

# Log Completion Time, mark end of analysis
end <- Sys.time()
cat(paste0("Analyses ended at ", end)," ","\n", file=file_log,append=TRUE)
cat(paste0("(Analyses took: ", round(as.numeric(difftime(end, start , units="mins")),digits=1) ," minutes)"),"\n"," ","\n", file=file_log,append=TRUE)
cat("###END###\n",file=file_log,append=TRUE)

#move logs and figures in out folder 
system(paste('mv ', file_log ,tmp, sep = " "))

