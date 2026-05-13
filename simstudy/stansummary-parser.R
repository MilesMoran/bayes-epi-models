library(stringr)
library(tidyr)
library(dplyr)
library(magrittr)

read_stansummary <- function(path) {

    lines <- readLines(path)
    
    model.summary <- list()
    model.summary$name <- str_split_1(lines[1], ": ")[2]
    model.summary$n_chains <- str_match(lines[2], "([0-9]+) chains:")[,2] %>% as.integer()
    model.summary$iter     <- str_match(lines[2], "iter=([0-9]+)")[,2]    %>% as.integer()
    model.summary$warmup   <- str_match(lines[2], "warmup=([0-9]+)")[,2]  %>% as.integer()
    model.summary$thin     <- str_match(lines[2], "thin=([0-9]+)")[,2]    %>% as.integer()
    
    model.summary$time_warmup <- 
        str_match(lines[4], "Warmup took \\((.*?)\\) seconds")[,2] %>% 
        str_split(",", simplify=TRUE) %>% 
        as.integer()
    model.summary$time_sampling <- 
        str_match(lines[5], "Sampling took \\((.*?)\\) seconds")[,2] %>% 
        str_split(",", simplify=TRUE) %>% 
        as.integer()
    
    tbl_colnames <- str_split(lines[7], "[ |\t]+")[[1]]
    tbl_colnames[c(1,6,7,8)] <- c("Variable", "q5", "Median", "q95")
    
    start <- 9
    end <- length(lines)-4
    model.summary$tbl <- read.table(
        text = lines[start:end], 
        header     = FALSE, 
        na.strings = c("NA","na","nan"), 
        col.names  = tbl_colnames,
        colClasses = c("character", rep("numeric", 11))
    )
    
    return(model.summary)
}

