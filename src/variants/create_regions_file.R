# Arjun Manrai
# Creates a regions file for bcftools from a sql db

library(RSQLite)
library(tidyverse)

prefix <- "../../data/"

# load SQL db
db.location <- paste0(prefix, "prebuilt-dbs/variants.db")
conn <- dbConnect(RSQLite::SQLite(), db.location)

res <- as_tibble(dbGetQuery(conn, "SELECT Chromosome, Position
                              FROM variants"))
res

dbDisconnect(conn)

out <- paste0(prefix,"clinvar/regions.tsv")
write_tsv(res, out, col_names = F)