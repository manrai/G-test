# Arjun Manrai
# Creates a regions file for bcftools from a sql db

library(RSQLite)
library(tidyverse)

prefix <- "../../data/"

is.pathogenic.simple <- function(info.col) { # simple example. replace in production version.
  #l <- strsplit(info.string,';')
  return(grepl("pathogenic",info.col))
}

# load merged clinvar variants
vcf.loc <- paste0(prefix,'/clinvar/clinvar_gnomad.vcf')
df <- read_tsv(vcf.loc,
               col_names = F)
colnames(df) <- c('Chromosome',
                  'Position',
                  'ID',
                  'Reference',
                  'Alternate',
                  'Qual',
                  'Filter',
                  'Info')
df$Chromosome <- as.integer(df$Chromosome)
df$Position <- as.integer(df$Position)
df$ID <- as.integer(df$ID)
df$Pathogenic <- is.pathogenic.simple(df$Info)
df %>% summary

# connect to sql db
db.location <- paste0(prefix, "prebuilt-dbs/variants.db")
conn <- dbConnect(RSQLite::SQLite(), db.location)

# write new sql table
dbWriteTable(conn, "clinvar", df)

# example queries
res.0 <- as_tibble(dbGetQuery(conn, "SELECT *
                              FROM clinvar
                              WHERE Chromosome = 11 AND
                              Position = 47369220"))
res.0

# join against gnomAD
res.1 <- as_tibble(dbGetQuery(conn, "SELECT v.Chromosome, v.Position, c.Pathogenic
                              FROM variants AS v
                              INNER JOIN clinvar AS c
                                ON (v.Chromosome = c.Chromosome AND
                                    v.Position = c.Position)
                              WHERE v.Chromosome = 11 AND
                              v.Position = 47369220"))
res.1

dbDisconnect(conn)




