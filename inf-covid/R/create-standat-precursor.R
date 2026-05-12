library(dplyr)
library(tidyr)
library(magrittr)
library(stringr)

library(tidycensus)

# Set directory relative to the path of THIS file (works in RStudio only)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))  
setwd("..") # navigate up to `inf-covid` dir

# Register the Census API key from a local text file, stored as plain text
census_api_key(readLines("data/census-api-key.txt"))

# Load the Census tract-to-PUMA crosswalk for 2020
tract_to_puma_df <- read.csv("data/2020-Tract-To-PUMA.txt", stringsAsFactors = FALSE)

################################################################################
# Retrieve Census variable metadata and define age-by-sex variables of interest.
################################################################################

# Load metadata for the 2020 Decennial Census Demographic and Housing
# Characteristics file. This is useful for identifying the variable labels we need.
all_vars_df <- load_variables(2020, "dhc", cache = TRUE)
View(all_vars_df)

# P12 is the sex-by-age table. These variables select male and female age bins,
# excluding the table total and sex subtotals.
vars <- c("P12_003N", "P12_004N", "P12_005N", "P12_006N", "P12_007N", 
          "P12_008N", "P12_009N", "P12_010N", "P12_011N", "P12_012N", 
          "P12_013N", "P12_014N", "P12_015N", "P12_016N", "P12_017N", 
          "P12_018N", "P12_019N", "P12_020N", "P12_021N", "P12_022N", 
          "P12_023N", "P12_024N", "P12_025N", "P12_027N", "P12_028N", 
          "P12_029N", "P12_030N", "P12_031N", "P12_032N", "P12_033N", 
          "P12_034N", "P12_035N", "P12_036N", "P12_037N", "P12_038N", 
          "P12_039N", "P12_040N", "P12_041N", "P12_042N", "P12_043N", 
          "P12_044N", "P12_045N", "P12_046N", "P12_047N", "P12_048N", 
          "P12_049N")

# Download tract-level 2020 Decennial Census counts for Michigan for the selected
# sex-by-age variables. The output contains one row per tract-variable pair.
age20 <- get_decennial(geography = "tract", 
                       variables = vars, 
                       year = 2020,
                       sumfile = "dhc",
                       state = "MI", 
                       geometry = FALSE)

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
# # Aggregate to tract-by-stratum counts and save to disk
################################################################################

pop_counts_tract_by_strata <-
    age20 %>% 
    mutate(
        # Split the tract GEOID into components matching the crosswalk keys.
        STATEFP  = str_sub(GEOID, 1,  2),
        COUNTYFP = str_sub(GEOID, 3,  5),
        TRACTCE  = str_sub(GEOID, 6, 11)
    ) %>%
    left_join(
        # Attach each tract's 2020 PUMA assignment.
        get_mapping_df(), 
        by = c("STATEFP", "COUNTYFP", "TRACTCE")
    ) %>%
    mutate(
        # Extract the numeric suffix from variables such as P12_003N.
        # This suffix indexes the original Census age-sex cell.        
        age_sex = as.integer(substr(variable, 5, 7)),
        
        # Collapse the detailed P12 age categories into decade-scale age bins.
        # Male cells are 003--025 and female cells are 027--049; both sets are
        # mapped into the same age intervals.
        age = case_when(
            (age_sex %in% c( 3, 4,27,28)            ) ~ "[0,10)",
            (age_sex %in% c( 5, 6, 7,29,30,31)      ) ~ "[10,20)",
            (age_sex %in% c( 8, 9,10,11,32,33,34,35)) ~ "[20,30)",
            (age_sex %in% c(12,13,36,37)            ) ~ "[30,40)",
            (age_sex %in% c(14,15,38,39)            ) ~ "[40,50)",
            (age_sex %in% c(16,17,40,41)            ) ~ "[50,60)",
            (age_sex %in% c(18,19,20,21,42,43,44,45)) ~ "[60,70)",
            (age_sex %in% c(22,23,46,47)            ) ~ "[70,80)",
            (age_sex %in% c(24,25,48,49)            ) ~ "[80,Inf)"
        ),
        
        # Assign sex from the P12 variable suffix.
        sex = case_when(
            (age_sex %in%  3:25) ~ "M",
            (age_sex %in% 27:49) ~ "F"
        ),        
        
        # Combine sex and age bin into one population stratum label.
        stratum = as.factor(paste(sex, age, sep="-")),
        
        # Store the Census count as an integer population count.
        pop = as.integer(value),
        
        # Convert PUMA code to integer for the final output matrix row labels.
        PUMA5CE = as.integer(PUMA5CE),
        
        # Place newly created fields before the original variable column.
        .before = variable
    ) %>% 
    select(GEOID, TRACTCE, PUMA5CE, stratum, pop) %>%  
    group_by(GEOID, TRACTCE, PUMA5CE, stratum) %>% 
    summarize(pop = sum(pop), .groups="drop") %>% 
    arrange(GEOID, TRACTCE, PUMA5CE, stratum)


## Sanity check: should be {10,077,331}
sum(pop_counts_tract_by_strata$pop)

## Save to disk (will be used to construct D and E for the Stan model)
saveRDS(pop_counts_tract_by_strata, "data/pop_counts_tract_by_strata.RDS")

