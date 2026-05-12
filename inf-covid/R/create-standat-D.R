library(dplyr)
library(tidyr)
library(magrittr)
library(stringr)

library(geosphere)
library(sf)

# Set directory relative to the path of THIS file (works in RStudio only)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))  
setwd("..") # navigate up to `inf-covid` dir

# Read the 2020 Michigan tract shapefile. The internal-point coordinate fields
# INTPTLON and INTPTLAT are used as tract centroids for distance calculations.
mi_tract <- st_read("data/mi_tracts_2020/mi_tract.shp")

# Load the Census tract-to-PUMA crosswalk for 2020
tract_to_puma_df <- read.csv("data/2020-Tract-To-PUMA.txt", stringsAsFactors = FALSE)

pop_counts_tract_by_strata <-
    readRDS("data/pop_counts_tract_by_strata.RDS") %>% 
    group_by(GEOID) %>% 
    summarize(pop = sum(pop), .groups="drop")
    
################################################################################
# Prepare tract-to-PUMA mapping for Michigan.
################################################################################

# Return a Michigan-only tract-to-PUMA crosswalk with GEOID components formatted
# as character strings using Census-standard zero padding.
get_mapping_df <- function() 
{
    tract_to_puma_df %>%
    subset(STATEFP==26) %>% 
    mutate(
        # State FIPS should be two digits; Michigan is 26.
        STATEFP = str_pad(as.character(STATEFP), width = 2, side = "left", pad = "0"),
        
        # County FIPS should be three digits.
        COUNTYFP = str_pad(as.character(COUNTYFP), width = 3, side = "left", pad = "0"),

        # Tract code should be six digits.
        TRACTCE = str_pad(as.character(TRACTCE), width = 6, side = "left", pad = "0"),
        
        # Keep PUMA code as character for the join; it is converted to integer later.
        PUMA5CE = as.character(PUMA5CE)
    ) %>% 
    return()
}

################################################################################
# Construct tract-to-tract distance data.
################################################################################

get_tract_dist_df <- function(mi_tract)
{
    #' Calculates the Haversine ("Great Circle") Distance between census tracts
    #' Returns:
    #'  - (data.frame) distance matrix in long format 
    
    # Extract tract internal-point longitude and latitude coordinates, coerce them
    # to numeric values, and compute all pairwise Haversine distances.
    # geosphere::distm() returns meters, so divide by 1000 to express distances
    # in kilometers. GEOIDs are used as row and column names to retain tract IDs.
    tract_dist_mat <- 
        as.matrix(mi_tract)[, c("INTPTLON", "INTPTLAT")] %>% 
        apply(2, as.numeric) %>%     
        distm(fun = distHaversine) %>% 
        magrittr::divide_by(1000) %>% 
        magrittr::set_rownames(mi_tract$GEOID) %>% 
        magrittr::set_colnames(mi_tract$GEOID)
    
    # Convert the square tract-by-tract distance matrix to long format with one
    # row per ordered tract pair (u, v). d_uv is the distance in kilometers.
    tract_dist_df <- 
        tract_dist_mat %>% 
        as.data.frame() %>%
        tibble::rownames_to_column("GEOID_u") %>%
        pivot_longer(-GEOID_u, names_to = "GEOID_v", values_to = "d_uv")
    
    return(tract_dist_df)
}

################################################################################
# Aggregate tract-pair distances to PUMA-pair population-weighted distances.
################################################################################

get_puma_dist_df <- function(tract_dist_df) 
{
    # Load tract-to-PUMA assignments from the NHGIS tract file, then join tract
    # population totals from the previously saved tract-level population counts.
    # The right join keeps the set of tracts appearing in pop_counts_tract.RDS.
    mi_tract_pop_counts <- 
        get_mapping_df() %>% 
        mutate(GEOID=paste0(STATEFP, COUNTYFP, TRACTCE), PUMA=PUMA5CE) %>%
        select(GEOID, PUMA) %>% 
        right_join(pop_counts_tract_by_strata, by="GEOID") %>% 
        mutate(pop = as.numeric(pop), PUMA=as.integer(PUMA))
    
    # Compute total population by PUMA. These totals are used below to normalize
    # the tract-pair weighted distance sums into population-weighted PUMA-pair
    # average distances.
    mi_PUMA_pop_counts <-
        mi_tract_pop_counts %>% 
        group_by(PUMA) %>% 
        summarize(pop = sum(pop)) %>% 
        ungroup()
    
    # Attach population and PUMA identifiers to each origin tract u and
    # destination tract v. Within each PUMA pair (i, j), sum e_u * e_v * d_uv
    # over all tract pairs, then divide by E_i * E_j to obtain the population-
    # weighted mean distance between residents of PUMA i and residents of PUMA j.
    puma_dist_df <- 
        tract_dist_df %>%
        left_join(rename(mi_tract_pop_counts, e_u = pop, PUMA_u = PUMA),
                  by = c("GEOID_u" = "GEOID")) %>%
        left_join(rename(mi_tract_pop_counts, e_v = pop, PUMA_v = PUMA),
                  by = c("GEOID_v" = "GEOID")) %>% 
        group_by(PUMA_u, PUMA_v) %>% 
        summarize(
            D_ij = sum(e_u * e_v * d_uv) 
        ) %>% 
        ungroup() %>% 
        left_join(mi_PUMA_pop_counts, by=c("PUMA_u"="PUMA")) %>% 
        left_join(mi_PUMA_pop_counts, by=c("PUMA_v"="PUMA"), suffix=c("_i","_j")) %>% 
        mutate(
            # Rename PUMA and population fields to match the final notation.
            PUMA_i = PUMA_u,
            PUMA_j = PUMA_v,
            E_i = pop_i,
            E_j = pop_j,
            # Normalize by the total number of ordered resident pairs across the
            # two PUMAs. For i == j, this includes within-tract terms with d_uv = 0.
            D_ij = D_ij / (E_i * E_j)
        ) %>% 
        select(PUMA_i, PUMA_j, E_i, E_j, D_ij)

    return(puma_dist_df)
}

################################################################################
# Build the PUMA-by-PUMA distance matrix, and save to disk
################################################################################

# Compute tract-pair distances, aggregate them to population-weighted PUMA-pair
# distances, and reshape the result into a square PUMA-by-PUMA matrix.
D <- 
    mi_tract %>% 
    get_tract_dist_df() %>% 
    get_puma_dist_df() %>% 
    select(c(PUMA_i, PUMA_j, D_ij)) %>%
    pivot_wider(
        names_from = PUMA_j,
        values_from = D_ij
    ) %>%
    tibble::column_to_rownames("PUMA_i") %>%
    as.matrix()     

saveRDS(D, "data/D.RDS")




