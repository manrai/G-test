# Arjun Manrai
# Builds variants db

library(RSQLite)
library(tidyverse)

# load data into tibble
src.directory <- '../../data/gnomAD/'
filenames <- list.files(src.directory, pattern="*.csv", full.names=T)
df <- lapply(filenames, read.csv)
df <- as_tibble(bind_rows(df))

# store as SQLite db
remove_existing_db <- 1 # warning this will delete previous db!
if (remove_existing_db) {
  file.remove("../../data/prebuilt-dbs/variants.db")
}
conn <- dbConnect(RSQLite::SQLite(), "../../data/prebuilt-dbs/variants.db")
dbWriteTable(conn, "variants", df)
dbListTables(conn) # double check table was created

# example queries
res.0 <- as_tibble(dbGetQuery(conn, "SELECT COUNT(*)
                              FROM variants"))
res.0

res.1 <- as_tibble(dbGetQuery(conn, "SELECT Chromosome, Position, Reference, Alternate
                              FROM variants 
                              WHERE rsID = 'rs147315081'"))
res.1

res.2 <- as_tibble(dbGetQuery(conn, "SELECT rsID
                              FROM variants
                              WHERE Chromosome = 11 AND
                                Position = 47369220 AND
                                Reference = 'C' AND
                                Alternate = 'T'"))
res.2

# should remove periods from colnames
res.3 <- as_tibble(dbGetQuery(conn, "SELECT rsID, `Allele.Frequency`
                              FROM variants
                              WHERE Chromosome = 11 AND
                              Position = 47369220 AND
                              Reference = 'C' AND
                              Alternate = 'T'"))
res.3

res.4 <- as_tibble(dbGetQuery(conn, "SELECT rsID, `Allele.Frequency` as Freq
                              FROM variants
                              ORDER BY Freq DESC
                              LIMIT 25"))
res.4

dbDisconnect(conn)
