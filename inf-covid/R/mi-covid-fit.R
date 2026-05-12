library(dplyr)
library(tidyr)
library(magrittr)
library(stringr)
library(lubridate)

library(cmdstanr)
library(posterior)

# Set directory relative to the path of THIS file (works in RStudio only)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))  
setwd("..") # navigate up to `inf-covid` dir

#' At this time, exact reproduceability for this analysis is not an option, as
#' we are forbidden from sharing the original Michigan incidence data per our DUA. 
#' However, I detail here the schema required to "reconstruct" an artificial copy
#' of this data.frame object, in case you want to fit this model yourself. The 
#' `df` variable below should contain at least the following columns: 
#' - age_coarse (chr): 10-year age-group strings, one of "[0,10)", "[10,20)", etc.
#' - PUMACE20 (chr): PUMA areal unit UIDs as strings
#' - year_week (chr): year+week strings, formatted XXXX_YY, where XXXX denotes year and YY denotes ISO week
#' - n (int): incident counts (zeroes coded as 0 instead of NA, etc.)
df <- readRDS("")

################################################################################
### Helper Functions 
################################################################################

pad <- function(x) {
    if (nchar(x) == 1) 
        return(paste0("0",x))
    else
        return(x)
}

puma_pad <- function(x) {
    if (nchar(x) < 5)  {
        padder <- paste0(rep("0",5-nchar(x)), collapse = "")
        return(paste0(padder,x))
    }
    else
        return(x)
}

################################################################################
### Curate the data needed to run the Stan program
################################################################################

## Map the 10-year age groups of the original data into a "coarse" partition
## corresponding to I=6 age groups
age_coarse <- unique(df$age_coarse)
age_map <- data.frame(
    age_coarse = age_coarse,
    age_bin = c("[0,10)", "[10,20)", "[20,30)", "[30,40)", rep("[40,60)",2), rep("[60,Inf)",3))
)
df_agg <- df |> 
    left_join(age_map) |>
    subset(!(PUMACE20 %in% c("00100", "00200"))) |> 
    group_by(age_bin, PUMACE20, year_week) |>
    summarise(n = sum(n)) |> 
    mutate(
        year_week_pad = sapply(
            strsplit(year_week,"_"),
            \(x) paste0(x[1],"_",pad(x[2]))
        )
    )

## Convert the incident counts from the original `data.frame` object into a 
## single T×G×I array whose `dimnames` describe the subpopulation metadata
weeks <- unique(df_agg$year_week_pad)
strata <- unique(df_agg$age_bin)
pumas <- unique(df_agg$PUMACE20)
dy_mat <- array(NA_real_, c(length(weeks), length(pumas), length(strata)))
dy_mat[
    cbind(
        match(df_agg$year_week_pad, weeks),
        match(df_agg$PUMACE20, pumas),
        match(df_agg$age_bin, strata)
    )
] <- df_agg$n

## Subset to exclude the PUMAs of the Michigan upper-peninsula 
## Subset to include only the first T=45 time points
G_idx <- 3:68
D <- (readRDS("data/D.RDS")[G_idx, G_idx] / 10) |> log()
dy_mat <- dy_mat[1:45,,]

## Map the 10-year age groups of the population counts into a "coarse" partition
## corresponding to I=6 age groups
E_big <- readRDS("data/E.RDS")[G_idx,]
E_big <- E_big[,1:9] + E_big[,10:18]
E <- cbind(E_big[,1:4], rowSums(E_big[,5:6]), rowSums(E_big[,7:9]))
bins <- c(0,10,20,30,40,60,Inf)
colnames(E) <- paste0("[",head(bins,-1),",",tail(bins,-1),")")

dat <-  list("T"=dim(dy_mat)[1], "G"=dim(dy_mat)[2], "I"=dim(dy_mat)[3], "E"=E, "D"=D, "dy"=dy_mat,
         grainsize = 1)

################################################################################
### (Optional) Plot an aggregated view of the COVID-19 dataset
################################################################################

library(ggplot2)
library(ggmatplot)

plot_incidence <- 
    apply(dy_mat,c(1,3),sum)[1:45,] %>% 
    ggmatplot(plot_type = "l") + 
    theme_bw() + 
    scale_x_continuous(breaks=seq(1,45,5), labels=seq(11,55,5)) +
    scale_color_discrete(name="Age Group", labels=strata) +
    scale_linetype_discrete(name="Age Group", labels=strata) +
    labs(title = "Reported Incidence of COVID-19 in MI Lower Peninsula, March 2020 to Jan 2021", 
         y="", x="ISO Week")

ggsave("plots/plot-incidence-300dpi.png", plot_incidence, 
       scale=1, dpi=300, width=7.5, height=3.75, units="in")

################################################################################
### Fit the Stan model
################################################################################

init_params <- function() 
{
    list2env(dat, envir=environment()); ## unpack names of `dat` here

    # pre-determine the number of non-degenerate r[t,g,i]'s
    n_nonzero_y = 0;
    for(t in 2:T) {
    for(g in 1:G) {
    for(i in 1:I) {
        cumulative_incidence <- sum( dy[1:(t-1),g,i] );
        if(cumulative_incidence > 0) {
            n_nonzero_y = (n_nonzero_y+1);
        }
    }}}

    w_diag_prop <- 0.4;
    concentration <- 1.8;
    shape_diag <- I * w_diag_prop * concentration;
    shape_offd <- (1.0*I/(I-1)) * (1-w_diag_prop) * concentration;
    gamma_shapes <- shape_diag + shape_offd*((I-1):0)

    # initialize parameters by randomly-drawing from their priors (mostly).
    # The exception is C_comp, which is cumbersome to initialize according to
    # its prior, so we use an `rnorm` alternative that keeps the composition
    # close to uniform weights (which is OK for small I)
    init_list <- list(
        "log_R"  = rnorm(T-1, 0, 2),
        "theta"  = abs(rcauchy(1,0,5)), 
        "gamma"  = abs(rcauchy(1,0,1)),
        "k"      = abs(rcauchy(1,0,1e4)),
        "eta_log_tau" = rnorm(G-1,0,3),
        "eta_log_rho" = rnorm(G,0,1),
        "delta"  = abs(rnorm(G*I,0,1)),
        "u"      = rgamma(I, gamma_shapes, rate=2*I),
        "mean_log_rho" = rnorm(1,0,2),
        "sd_log_rho"   = abs(rnorm(1,0,1)),
        "log_r_raw"    = rnorm(n_nonzero_y,0,1),
        "C_comp" = rnorm(I*(I-1)/2, 0, 0.5)
    )
    return(init_list)
}

# Precompute the (seeded) initial values manually to ensure reproducibility
set.seed(1)
init_list <- replicate(
    n = 16,
    expr = init_params(),
    simplify = FALSE
)

mod <- cmdstan_model("stan/mi-covid-model.stan",
                      cpp_options = list(stan_threads = TRUE),
                      force_recompile=TRUE)

fit <- mod$sample(data = dat,
                  chains = 16,
                  parallel_chains = 8,
                  threads_per_chain = 1,
                  iter_sampling = 2000,
                  iter_warmup = 2000,
                  init = init_list,
                  refresh = 400,
                  step_size = 0.01,
                  max_treedepth = 12,
                  sig_figs = 9,
                  output_dir = "output",
                  output_basename = "priors-mi-covid-model",
                  seed=1) 

################################################################################
################################################################################

