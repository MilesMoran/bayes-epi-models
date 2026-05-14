library(dplyr)
library(tidyr)
library(magrittr)
library(stringr)
library(lubridate)

library(cmdstanr)
library(posterior)
library(hhh4contacts)

# Set directory relative to the path of THIS file (works in RStudio only)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))  
setwd("..") # navigate up to `inf-noro` dir

################################################################################
### Curate the data needed to run the Stan program
################################################################################

## Map the 5-year age groups of the original data into a "coarse" partition
## corresponding to I=6 age groups. Also subset the time-series to the 4-year 
##  (208-week) period spanning July 2011 to July 2015.
GROUPING <- c(1,2,2,4,4,2)
t_vec <- 27:234 
E <- aggregateCountsArray(pop2011, dim=2, grouping=GROUPING) 
strata <- dimnames(E)[[2]]

## Convert the incident counts from the original `surveillance::sts` object into 
## a single T×G×I array whose `dimnames` describe the subpopulation metadata
noroBEall <- noroBE(by="all", flatten=TRUE, agegroups=GROUPING, timeRange=c("2011-w27","2015-w26")) ## what's used in Meyer&Held

## Use the normalized contact matrix seen in Meyer and Held (2017) to define the
## fixed eigenvalues and eigenvectors needed in the "reduced" (Meyer-Held) model
W_strata_t <- contactmatrix("reciprocal", "all", GROUPING, normalize=TRUE)

## For the remaining model components, we match the original Meyer-Held model:
##  (1) Let the geog. "distance" matrix D be defined by adjacency order
##  (2) Let `log_e` denote the population fractions to be used as offsets
##  (3) Let `sin_omega_train` and `cos_omega_train` refer to the endemic
##      seasonality terms for the weeks under study
dat <- list(
    "y"     = array(data = attr(noroBEall, "observed"),
                    dim = c(208,12,6),
                    dimnames = list(week = dimnames(counts)$week[t_vec],
                                    district = dimnames(E)[[1]],
                                    agegroup = dimnames(E)[[2]])  ),
    "T"     = 208,
    "G"     = 12,
    "I"     = 6,
    "D"     = (neighbourhood(noroBEall)[1:12,1:12] + 1),
    "log_e" = log(E/sum(E)),
    "xmas"  = as.integer((t_vec %% 52) %in% c(0,1)),
    "sin_omega" = sin((t_vec-1)*(2*pi/52)),
    "cos_omega" = cos((t_vec-1)*(2*pi/52)),
    "EValC" = eigen(W_strata_t, symmetric = FALSE)$values,
    "EVecC" = eigen(W_strata_t, symmetric = FALSE)$vectors
)

rm(GROUPING, t_vec, E, noroBEall, W_strata_t)

################################################################################
### (Optional) Plot an aggregated view of the Berlin Norovirus dataset
################################################################################

library(ggplot2)
library(ggmatplot)

plot_incidence <-
    apply(dat$y, c(1,3), sum) %>% 
    ggmatplot(plot_type = "l") +
    theme_bw() + 
    scale_x_continuous(breaks=seq(1,208,26), 
                       labels=dimnames(dat$y)$week[seq(1,208,26)]) +
    scale_color_discrete(name="Age Group", labels=strata) + 
    scale_linetype_discrete(name="Age Group", labels=strata) + 
    labs(title = "Reported Incidence of Norovirus in Berlin, July 2011 to July 2015", 
         y="", x="ISO Week")
    
ggsave("plots/plot-incidence-300dpi.png", plot_incidence, 
       scale=1, dpi=300, width=7.5, height=3.75, units="in")

################################################################################
### Fit the Stan model(s)
#
# The `dat` list and the `init_params()` function are written to support fitting
# either the full or reduced model. When fitting the reduced model, parameters 
# only used by the full model (`theta`, `u`, `C_comp`, and `log_r_raw`) are 
# ignored. When fitting the full model, parameters only used by the reduced model 
# (`logpower`) and data only used by the reduced model (`EValC` and `EVecC`) are
# ignored. Switch between the fitting full and reduced models by changing the 
# filepath in `cmdstan_model(...)` to the other model's filepath.
################################################################################

init_params <- function() 
{
    list2env(dat, envir=environment()); ## unpack names of `dat` here
    
    # pre-determine the number of non-degenerate r[t,g,i]'s for the "full" model
    n_nonzero_y = sum(y[1:(T-1),,] > 0)
    
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
        "end_1"   = rnorm(1,3,2),
        "end_christmas"  = rnorm(  1, 0, 3),
        "end_strata_raw" = rnorm(I-1, 0, 3),
        "end_geogs_raw"  = rnorm(G-1, 0, 3),
        "end_sin" = rnorm(I, 0, 3),
        "end_cos" = rnorm(I, 0, 3),
        "ne_1"    = rnorm(1, 2, 5),
        "ne_strata_raw"  = rnorm(I-1, 0, 3),
        "ne_geogs_raw"   = rnorm(G-1, 0, 3),
        "ne_log_pop"     = rnorm(  1, 0, 2),
        "neweights_logd" = rnorm(  1, 0, 1),
        "theta"    = abs(rcauchy(1, 0, 1)), 
        "overdisp" = abs(rcauchy(1, 0, 1)),
        "log_r_raw"= rnorm(n_nonzero_y, 0, 1),
        "u"        = rgamma(I, gamma_shapes, rate=2*I),
        "C_comp"   = rnorm(I*(I-1)/2, 0, 0.5),
        "logpower" = rnorm(1, 0, 0.75)
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

mod <- cmdstan_model("stan/be-noro-model-full.stan", force_recompile=TRUE)

fit <- mod$sample(data = dat,
                  chains = 16,
                  parallel_chains = 16,
                  iter_sampling = 2000,
                  iter_warmup = 2000,
                  init = init_list,
                  refresh = 400,
                  step_size = 0.02,
                  max_treedepth = 12,
                  sig_figs = 9,
                  output_dir = "output",
                  output_basename = "be-noro-model-full",
                  seed = 1)

################################################################################
################################################################################
