library(dplyr)
library(tidyr)
library(magrittr)

# Set directory relative to the path of THIS file (works in RStudio only)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))  
setwd("..") # navigate up to `inf-covid` dir

################################################################################
# Convert tract-level age-by-sex counts to PUMA-level age-by-sex counts.
################################################################################

E <-
    readRDS("data/pop_counts_tract_by_strata.RDS") %>% 
    # Keep only the final fields needed to construct the exposure/population matrix.
    select(c(PUMA5CE, stratum, pop)) %>% 
    group_by(PUMA5CE, stratum) %>% 
    summarize(pop = sum(pop), .groups="drop") %>% 
    # Convert from long format to one row per PUMA and one column per stratum.
    pivot_wider(names_from = stratum, values_from = pop) %>% 
    # Use PUMA code as the matrix row name rather than a data column.
    tibble::column_to_rownames("PUMA5CE") %>% 
    # Convert the wide data frame to a numeric matrix.
    as.matrix()

saveRDS(E, "data/E.RDS")