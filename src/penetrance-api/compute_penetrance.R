#!/usr/bin/env Rscript

# Arjun Manrai
# NIH/NHLBI BDC
# Script to compute penetrance across BDC cohorts with dynamic G/P criteria

library(RSQLite)
suppressPackageStartupMessages(suppressWarnings(library(tidyverse)))
library(ggthemes)
library(optparse)

options(warn=-1)

# load helpers (incl. beta-binomial sampling function)
source('../utilities/sample_beta.R')

prefix <- "../../data/"

# define gene panel
# "incontrovertible list" from https://www.ahajournals.org/doi/10.1161/JAHA.119.015473
panel <- c('MYH7','TNNT2','MYBPC3','TPM1','MYL3','TNNI3','MYL2','ACTC1')

# option #2
#panel <- c("ABCC9", "ACTC1", "ACTN2", "ANKRD1", "BAG3", "CASQ2", "CAV3",
#  "CHRM2", "CRYAB", "CSRP3", "DES", "DMD", "DOLK", "DSC2", "DSG2",
#  "DSP", "DTNA", "EMD", "FHL2", "GATAD1", "GLA", "ILK", "JPH2",
#  "JUP", "LAMA4", "LAMP2", "LDB3", "LMNA", "CAVIN4", "MYBPC3",
#  "MYH6", "MYH7", "MYL2", "MYL3", "MYLK2", "MYOM1", "MYOZ2", "MYPN",
#  "NEBL", "NEXN", "PDLIM3", "PKP2", "PLN", "PRDM16", "PRKAG2",
#  "PTPN11", "RAF1", "RBM20", "RYR2", "SCN5A", "SGCD", "TAZ", "TCAP",
#  "TMEM43", "TNNC1", "TNNI3", "TNNT2", "TPM1", "TRDN", "TTN", "TTR", "VCL")

make.options <- function() {
  # make options
  option_list = list(
    make_option(c("-v", "--vus"), type="integer", default=0,
                help="Include variants of uncertain significance (1/0)?"),
    make_option(c("-l", "--likely_pathogenic"), type="integer", default=1,
                help="Include likely pathogenic variants (1/0)?"),
    make_option(c("-p", "--pathogenic"), type="integer", default=1,
                help="Include pathogenic variants (1/0)?"),
    make_option(c("--lvwt_min"), type="double", default=8,
                help="Minimum left ventricular wall thickness (in mm)"),
    make_option(c("--lvwt_max"), type="double", default=13,
                help="Maximum left ventricular wall thickness (in mm)"),
    make_option(c("--age_min"), type="integer", default=18,
                help="Minimum age (in years)"),
    make_option(c("--age_max"), type="integer", default=40,
                help="Maximum age (in years)"),
    make_option(c("-g","--gender"), type="character", default="all",
                help="Gender ('male','female','all')"),
    make_option(c("-r","--race"), type="character", default="all",
                help="Race"),
    make_option(c("-c","--cohorts"), type="character", default="jackson,cardia",
                help="Cohorts (as comma-separated list e.g. jackson,cardia,hvh)")
  )
  opt_parser = OptionParser(option_list=option_list)
  opt = parse_args(opt_parser)
  return(opt)
}

set.parameters <- function(opt) {
  # need to add in remaining options
  cohorts <- c("jackson","hvh")

  variants <- list(
    vus = opt$vus,
    likely_pathogenic = opt$likely_pathogenic,
    pathogenic = opt$pathogenic
  )

  phenos <- list(
    lvwt.min = opt$lvwt_min,
    lvwt.max = opt$lvwt_max,
    htn = 1 # 1 = include htn in controls; 0 = exlcude htn from controls
  )

  l <- list(cohorts=cohorts,variants=variants,phenos=phenos)

  return(l)
}

load.bdc <- function() {
  # connect to sql db
  db.location <- paste0(prefix, "prebuilt-dbs/bdc-test.db")
  conn <- dbConnect(RSQLite::SQLite(), db.location)
  df <- as_tibble(dbGetQuery(conn, "SELECT * FROM bdc"))
  dbDisconnect(conn)
  return(df)
}

load.variant.annotations <- function() {
  # connect to sql db
  db.location <- paste0(prefix, "prebuilt-dbs/variants.db")
  conn <- dbConnect(RSQLite::SQLite(), db.location)
  df <- as_tibble(dbGetQuery(conn, "SELECT * FROM clinvar"))
  dbDisconnect(conn)

  # heuristic [NEED TO UPDATE]
  rsID <- paste0('rs',sub(".*RS=","",df$Info))
  df <- as_tibble(cbind(rsID=rsID,Pathogenic=df$Pathogenic))
  df <- df %>%
    filter(nchar(rsID) < 20)

  df <- df %>%
    group_by(rsID) %>%
    slice(1) %>%
    ungroup()

  df$Likely_Pathogenic <- sample(c(0,1),nrow(df),replace=TRUE,prob=c(0.77,0.23))
  df$VUS <- sample(c(0,1),nrow(df),replace=TRUE,prob=c(0.6,0.4))

  return(df)
}

# compute penetrance dynamically with user-specified geno & pheno criteria
# use geno for variants considered
# use pheno criteria for CONTROLS [keep numerator terms fixed for simplicitly for now]
compute.pen <- function(df, var.annot, cohorts, variants, phenos) {
  # literature-based estimates
  caf.obs <- 917 # morita & fokstuen
  caf.sample.size <- 2912

  prev.obs <- 7 # maron
  prev.sample.size <- 4111

  controls <- df %>% #
    filter(cohort %in% cohorts) %>%
    filter(lvwt >= phenos$lvwt.min & lvwt <= phenos$lvwt.max)

  if (phenos$htn == 0) { # exclude htn from controls
    controls <- controls %>% filter(htn == 0)
  }

  # define which genetic variants to include
  rs <- c()
  if (variants$vus == 1) { # include vus
    rs.vus <- var.annot %>% filter(VUS == 1) %>% select(rsID)
    rs.vus <- intersect(rs.vus$rsID,colnames(controls))
    rs <- c(rs,rs.vus)
  }

  if (variants$likely_pathogenic == 1) { # include lp
    rs.lp <- var.annot %>% filter(Likely_Pathogenic == 1) %>% select(rsID)
    rs.lp <- intersect(rs.lp$rsID,colnames(controls))
    rs <- c(rs,rs.lp)
  }

  if (variants$pathogenic == 1) { # include p
    rs.p <- var.annot %>% filter(Pathogenic == 1) %>% select(rsID)
    rs.p <- intersect(rs.p$rsID,colnames(controls))
    rs <- c(rs,rs.p)
  }

  rs <- unique(rs)
  controls <- controls %>%
    select(c('id','cohort','lvwt','htn',rs))

  control.af <- rowSums(controls[,5:ncol(controls)])/(nrow(controls)*2)
  control.af <- 1-prod(1-control.af)
  control.n <- nrow(controls)

  #case.af <- rowSums(cases[,3:ncol(cases)])/(nrow(cases)*2) # if defining cases from cohorts not lit.
  #case.af <- 1-prod(1-case.af)
  #case.n <- nrow(cases)

  sample.case.af <- sample.beta(freq = caf.obs/caf.sample.size, n = caf.sample.size)
  case.af <- sample.case.af %>%
    quantile(c(0.05, 0.5, 0.95))
  sample.prev <- sample.beta(freq = prev.obs/prev.sample.size, n = prev.sample.size)
  prev <- sample.prev %>%
    quantile(c(0.05, 0.5, 0.95))
  sample.control.af <- sample.beta(freq = control.af, n = control.n)
  control.af <- sample.control.af %>%
    quantile(c(0.05, 0.5, 0.95))

  penetrance <- as_tibble(
    sample.case.af * sample.prev / (sample.control.af)
  )

  colnames(penetrance) <- 'Penetrance'

  return(penetrance)
}

plot.pen <- function(pen) {
  pen <- as_tibble(pen)
  p <- ggplot(data=pen,aes(x=Penetrance)) +
    geom_density() +
    theme_pander(nomargin = F) +
    labs(x = "Penetrance", y = "Density") +
    geom_vline(xintercept = as.double(quantile(pen$Penetrance,0.025)), linetype = "dashed", color = "red") +
    geom_vline(xintercept = as.double(quantile(pen$Penetrance,0.975)), linetype = "dashed", color = "red") +
    ggtitle("Penetrance Across Cohorts")
  return(p)
}

## main
df <- load.bdc()
df <- df[,1:25] # modify
var.annot <- load.variant.annotations()
opt <- make.options()

# set parameters from default/command line args
l <- set.parameters(opt)
cohorts <- l$cohorts
variants <- l$variants
phenos <- l$phenos

# verbose messages
#print(paste0('Computing penetrance with vus = ',variants$vus,"; lp = ",variants$likely_pathogenic,"; p = ",variants$pathogenic))
#print(paste0('...also using lvwt_min = ',phenos$lvwt.min," and lvwt_max = ",phenos$lvwt.max))

pen <- compute.pen(df, var.annot, cohorts, variants, phenos)

# save pen plot to file
p <- plot.pen(pen)

filename = paste0(paste(cohorts,collapse='_'),"_",paste(unlist(variants),collapse="_"),"_",paste(unlist(phenos),collapse="_"),'.pdf')
ggsave(filename=filename,p)

# print out pen median, low, high
low = as.double(round(quantile(pen$Penetrance,0.025),4))
map = as.double(round(quantile(pen$Penetrance,0.5),4))
high = as.double(round(quantile(pen$Penetrance,0.975),4))

print(paste(low,map,high,sep=';'))
