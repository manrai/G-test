# Arjun Manrai
# Draws on var-pheno empirical distributions to generate
#   multi-cohort data (as seen on BDC)

library(RSQLite)
library(tidyverse)
library(stringi)

prefix <- "../../data/"

## load data ##
# connect to sql db
db.location <- paste0(prefix, "prebuilt-dbs/variants.db")
conn <- dbConnect(RSQLite::SQLite(), db.location)

# get variants from subset of genes
#get.vars <- function() { # example with MYBPC3 and MYH7
#  res <- as_tibble(dbGetQuery(conn, "SELECT rsID,`Allele.Frequency`
#                              FROM variants
#                              WHERE (Chromosome = 11 AND
#                              Position >= 47352957 AND Position <= 47374253) OR
#                              (Chromosome = 14 AND
#                              Position >= 23881947 AND Position <= 23904870)"))
#  return(res)
#}

# get all available variants
get.vars <- function() { # example with MYBPC3 and MYH7
  res <- as_tibble(dbGetQuery(conn, "SELECT rsID,`Allele.Frequency`
                              FROM variants"))
  res <- res %>% 
    group_by(rsID) %>% 
    slice(1) %>%
    ungroup()

  return(res)
}

generate.genotypes <- function(vars) {
  vars <- vars %>%
    rowwise() %>%
    mutate(genotype = sample(c(0,1,2),
                             1,
                             prob=c(prob.0, prob.1, prob.2))) %>%
    ungroup()
  return(vars$genotype)
}

create.inds <- function(n, vars) {
  high.freq.vars <- vars %>%
    arrange(desc(Allele.Frequency)) %>%
    top_n(100)
  
  high.freq.vars <- high.freq.vars %>%
    mutate(prob.0 = (1-Allele.Frequency)^2,
           prob.1 = (1-Allele.Frequency)*Allele.Frequency,
           prob.2 = Allele.Frequency^2)
  
  rsIDs <- high.freq.vars$rsID
  ids <- stri_rand_strings(n,8,"[A-Za-z0-9]")
  cohort <- sample(c("jackson","cardia","hvh"),replace=TRUE,n)
  genotypes <- c()
  for (id in ids) {
    these.genotypes <- generate.genotypes(high.freq.vars)
    genotypes <- rbind(genotypes,these.genotypes)
  }
  genotypes <- as_tibble(genotypes)
  genotypes <- as_tibble(sapply(genotypes,as.integer))
  colnames(genotypes) <- rsIDs
  genotypes <- as_tibble(cbind(id=ids,cohort,genotypes))
}

generate.phenos <- function(genotypes,pheno.list = c('lvwt','htn')) {
  n <- nrow(genotypes)
  
  phenos <- c()
  for (pheno in pheno.list) {
    if (pheno == 'lvwt')  { # continuous
      phenos <- cbind(phenos,rnorm(n,10,2))
    }
    if (pheno == 'htn') { # binary
      phenos <- cbind(phenos,sample(c(0,1),n,replace=TRUE))
    }
  }
  
  df <- as_tibble(cbind(genotypes[,1:2],phenos,genotypes[,3:ncol(genotypes)]))
  colnames(df) <- c(colnames(genotypes[,1:2]),pheno.list,
                    colnames(genotypes[,3:ncol(genotypes)]))
  
  df
  
  return(df)
}

write.db <- function(df,remove.existing=1) {
  if (remove.existing) { # overwrite existing database
    file.remove("../../data/prebuilt-dbs/bdc-test.db")
  }
  conn <- dbConnect(RSQLite::SQLite(), "../../data/prebuilt-dbs/bdc-test.db")
  dbWriteTable(conn, "bdc", df)
}

## main
n <- 5000
vars <- get.vars()
genotypes <- create.inds(n,vars)
df <- generate.phenos(genotypes,pheno.list = c('lvwt','htn'))
write.db(df)

dbDisconnect(conn)
