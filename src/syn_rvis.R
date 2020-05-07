################################################################################
# Calculate synRVIS scores                                                     #
# Author: Ryan Dhindsa                                                         #
################################################################################

# Imports ----------------------------------------------------------------------
library(data.table)
library(tidyverse)
library(lazyeval)
library(argparse)


# Globals ----------------------------------------------------------------------

FOURFOLD <- c("CTT", "CTC", "CTA", "CTG", 
              "GTT", "GTC", "GTA", "GTG", 
              "TCT", "TCC", "TCA", "TCG", 
              "CCT", "CCC", "CCA", "CCG", 
              "ACT", "ACC", "ACA", "ACG", 
              "GCT", "GCC", "GCA", "GCG", 
              "CGT", "CGC", "CGA", "CGG", 
              "GGT", "GGC", "GGA", "GGG")


# Functions --------------------------------------------------------------------
ReadSynGERP <- function(syn.gerp.path) {
  # Reads in synGERP file and calculates the percentile for each gene
  #
  # Args:
  #   syn.gerp.path: path to synGERP file
  #   
  # Returns:
  #   syn.gerp: a dataframe
  syn.gerp <- fread(syn.gerp.path)
  
  # Calculate percentile; we multiply by -1 because we want the
  #   lower percentile to indicate intolerance (positive GERP score ==
  #   higher conservation)
  syn.gerp$synGERP.centile <- PercentileRank(-1 * syn.gerp$synGERP)
  return(syn.gerp)
}

ReadVariants <- function(variants.path, csc.path, 
                         fourfold.only=F, fold.mafs=T) {
  # Read in tab-delimitted file containing variants
  # 
  # Args:
  #   variants.path: path to variant file
  #   csc.path: path to csc scores
  #   fourfold.only: whether to retain fourfold codons only
  #   fold.mafs: whether to fold the MAFs 
  #
  # Returns: 
  #   df: a dataframe
  
  df <- fread(variants.path)
  # remove multi-allelic sites
  df <- df[!grepl(",", df$ALT), ]
  df <- df[!grepl(",", df$REF), ]
  df$MAF <- as.numeric(df$MAF)
  
  # convert codons to upper case
  df$CODON_CHANGE <- toupper(df$CODON_CHANGE)
  df <- separate(data = df, col = CODON_CHANGE, 
                 into = c("REF_CODON", "ALT_CODON"), 
                 sep = "\\/")
  
  if(fold.mafs) {
    df$MAF <- ifelse(df$MAF > 0.5, 1 - df$MAF, df$MAF)
  }
  
  # Subset to fourfold degenerate sites (optional)
  if(fourfold.only) {
    df <- subset(df, df$REF_CODON %in% FOURFOLD & df$ALT_CODON %in% FOURFOLD)
  }
  
  # Annotate variants with CSC scores
  csc <- fread(csc.path)
  df <- AnnotateCSC(df, csc)
  return(df)
}


AnnotateCSC <- function(df, csc, score.use="293T_endo") {
  # Annotates dataframe with codon stability coefficient scores 
  #   (from Wu et al. supplement). 
  #
  # Args:
  #   df: dataframe containing variants
  #   csc: df containing codons and their CSC scores
  #   score.use: which CSC score you want to use (e.g. 293T_endo);
  #     must correspond to a column in the CSC file
  #
  # Returns:
  #   df: df with CSC scores appended
  csc <- subset(csc, select = c("codon", score.use))
  
  colnames(csc) <- c("CODON", "CSC")
  df <- AnnotateCodons(df, csc, "REF_CODON", "CODON", "REF")
  df <- AnnotateCodons(df, csc, "ALT_CODON", "CODON", "ALT")
  
  df$DELTA_CSC <- df$CSC_ALT - df$CSC_REF
  
  return(df)
}


AnnotateCodons <- function(df1, df2, by.1, by.2, suffix) {
  # Annotate codons in one df based on values in column of another df
  #
  # Args:
  #   df1: dataframe containing codons to annotate
  #   df2: dataframe containing codons and annotations
  #   by.1: which column to merge on in df1
  #   by.2: which column to merge on in df2
  #   suffix: suffix to add to column in annotation 
  #
  # Returns:
  #   annotated.df: dataframe with annotations added
  annotated.df <- merge(df1, df2, by.x = by.1, by.y = by.2, 
                        sort = F, suffixes = c(suffix), all.x = T)
  
  setnames(annotated.df, paste0(names(annotated.df), 
                                ifelse(names(annotated.df) %in% 
                                         setdiff(names(df2), names(df1)), 
                                       paste0("_", suffix), "")))
  return(annotated.df)
}


CodonPref <- function(df, csc.cutoff) {
  # Classifies codon changes as pref -> non-pref, non-pref --> pref, and 
  # Neutral
  #
  # Args:
  #   df: dataframe of SNVs
  #   csc.cutoff: the cutoff for CSC
  #
  # Returns:
  #   df: DF with annotated codon changes 
  
  # CSC pref
  df$CSC_CODON_PREF <- ifelse(df$CSC_REF > csc.cutoff & 
                                df$CSC_ALT < csc.cutoff, "O -> NO", 
                              ifelse(df$CSC_REF < csc.cutoff & 
                                       df$CSC_ALT > csc.cutoff, "NO -> O", 
                                     "Neutral"))
  df$CSC_CODON_PREF_CSQ <- paste(df$CSQ, df$CSC_CODON_PREF)
  
  return(df)
}


CalcXY <- function(df, maf.threshold) {
  # Takes in a dataframe of synonymous variants and calculates X (sum of
  # synonynmous variants per gene) and Y (# of syn variants with MAF >
  # maf.threshold)
  #
  # Args:
  #   df: a tab-delimitted dataframe containing synonymous variants
  #   maf.threshold: MAF threshold for Y
  #
  # Returns:
  #   xy: a dataframe with X and Y tallies
  xy <- df %>% 
    group_by(GENE) %>% 
    summarize(x = n(), 
              y = sum(CSC_CODON_PREF == "O -> NO" & MAF > maf.threshold))
  return(xy)
}


PercentileRank <- function(x) {
  # Calculates percentile of each observation in vector
  # Args: 
  #   x: a numeric vector
  #
  # Returns:
  #   none
  trunc(rank(x))/length(x)
}


CalcRVIS <- function(df) {
  # Creates linear model of X and Y tallies (from CalcXY()) and then
  # calculates studentized residuals
  #
  # Args:
  #   df: df containing X and Y tallies
  # 
  # Returns:
  #   df: a dataframe containing gene names and scores
  df <- df[complete.cases(df), ]
  
  lin.mod <- lm(formula = y ~ x, data = df)
  df$synRVIS <- rstudent(lin.mod)
  df$synRVIS.centile <- PercentileRank(df$synRVIS)
  
  return(df)
}


# Main -------
parser <- ArgumentParser()

parser$add_argument("-i", "--input",
                    required = TRUE,
                    help="Path to synonymous variants") 

parser$add_argument("-csc", "--csc-scores",
                    required = TRUE,
                    help="Path to CSC scores")

parser$add_argument("-maf", "--maf-cutoff",
                    required = TRUE,
                    help="MAF cutoff for defining common variants (Y-axis)")

parser$add_argument("-o", "--output",
                    required = TRUE,
                    help="Path to output directory")

args <- parser$parse_args()


syn.rvis <- ReadVariants(args$input, args$csc_scores, fold.mafs = T, 
                         fourfold.only = F) %>%
  CodonPref(csc.cutoff = 0) %>%
  CalcXY(maf.threshold = args$maf_cutoff) %>%
  CalcRVIS()

write.table(syn.rvis, args$output, quote = F, sep = "\t", row.names = F)
