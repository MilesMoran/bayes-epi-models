library(dplyr)
library(tidyr)
library(magrittr)
library(stringr)
library(lubridate)

library(cmdstanr)
library(posterior)
library(ggplot2)
library(knitr)
library(kableExtra)

# Set directory relative to the path of THIS file (works in RStudio only)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))  
setwd("..") # navigate up to `inf-covid` dir

################################################################################
### Load & summarize Stan draws
################################################################################

fnames_prior <- list.files("output", "^priors-mi-covid-model-\\d+\\.csv$", full.names=TRUE)
fnames_postr <- list.files("output", "^mi-covid-model-\\d+\\.csv$", full.names=TRUE)

## Load everything except the r[t,g,i]'s
variables <- c("W_strata", "log_alpha", "log_R", "theta", "gamma", "k", "log_tau", "u",
               "rho", "delta", "mean_log_rho", "sd_log_rho", "C_comp")
    
draws_prior <- read_cmdstan_csv(fnames_prior, variables)
draws_postr <- read_cmdstan_csv(fnames_postr, variables)

summ_prior <- summarize_draws(draws_prior$post_warmup_draws)
summ_postr <- summarize_draws(draws_postr$post_warmup_draws)

################################################################################
### Cursory inspections: make sure nothing is horribly wrong 
################################################################################

summ_postr %>% 
    arrange(desc(rhat)) %>% 
    select(variable, mean, median, rhat, ess_bulk, ess_tail) %>%  
    head(n=20L)

################################################################################
### Plot #1: Matrix plot of prior and posterior for `W_strata`
################################################################################

strata <- c("[0,10)", "[10,20)", "[20,30)", "[30,40)", "[40,60)", "[60,Inf)")

plot_wstrata <-
    summ_prior %>% 
    mutate(prior=TRUE) %>% 
    rbind(mutate(summ_postr, prior=FALSE)) %>% 
    subset(grepl("W_strata", variable)) %>% 
    select(variable, mean, q5, q95, prior) %>% 
    separate_wider_regex(
        variable, patterns=c("W_strata\\[", i = "\\d+", ",", i_ = "\\d+", "\\]")
    ) %>% 
    mutate(
        i = as.integer(i),
        i_ = as.integer(i_),
        Destination = strata[i],
        Origin = strata[i_],
        mean.val = mean,
        across(mean:q95, ~sprintf("%.2f", .x)),
        interval = paste0("(", q5, "-", q95, ")")
    ) %>% 
    ggplot(aes(x=Origin, y=Destination, fill=mean.val)) +
        geom_tile(color="white", linewidth=0.4) +
        geom_text(aes(label=mean), size=4, fontface="bold", vjust=-0.1) +
        geom_text(aes(label=interval), size=2.5, vjust=1.5) +
        facet_wrap(~rev(prior), labeller=\(x){ data.frame(c("Prior", "Posterior")) }) + 
        scale_y_discrete(limits=rev(strata)) +
        scale_fill_gradient2(name="Mean", low="blue", mid="white", high="red", midpoint=0) +
        coord_equal() +
        labs(x="", y="") +
        theme_bw(base_size = 12) +
        theme(
            panel.grid = element_blank(),
            axis.title.x = element_text(margin = margin(t = 8)),
            axis.title.y = element_text(margin = margin(r = 8)),
            strip.text = element_text(face = "bold"),
            legend.position = "none"
        )

ggsave("plots/plot-wstrata-300dpi.png", plot_wstrata,
       scale=1.5, dpi=300, width=7, height=3.5, units="in")

################################################################################
### Table #1: Prior and posterior summaries for a few primary estimands
################################################################################

### Posterior Medians ###

summ_prior %>% 
    mutate(prior=TRUE) %>% 
    rbind(mutate(summ_postr, prior=FALSE)) %>% 
    subset(variable %in% c("W_strata[1,1]", "W_strata[2,2]", "W_strata[1,2]", "k", "theta", "gamma")) %>% 
    select(variable, median, q5, q95, prior) %>% 
    mutate(
        across(median:q95, ~if_else(variable == "k", .x / 1e4, .x)),
        across(median:q95, ~format(.x, digits=1, nsmall=2, scientific=FALSE)),
        estm = paste0(median, " (", q5, "-", q95, ")")
    ) %>% 
    select(variable, estm, prior) %>% 
    pivot_wider(names_from=prior, values_from=estm, names_vary="slowest") %>% 
    slice(c(1,3,2,5,4,6)) %>% 
    mutate(variable = c("$w^{(I)}_{11}$", "$w^{(I)}_{22}$", "$w^{(I)}_{12}$",
                        "$\\gamma$", "$\\theta$", "$k/10^{4}$")) %>% 
    set_colnames(c("Param", "Prior", "Posterior")) %>% 
    kable(
        format = "latex",
        align = "crr",
        escape = FALSE
    ) %>%
    clipr::write_clip()

### Posterior Means ###

summ_prior %>% 
    mutate(prior=TRUE) %>% 
    rbind(mutate(summ_postr, prior=FALSE)) %>% 
    subset(variable %in% c("W_strata[1,1]", "W_strata[2,2]", "W_strata[1,2]", "k", "theta", "gamma")) %>% 
    select(variable, mean, q5, q95, prior) %>% 
    mutate(
        across(mean:q95, ~if_else(variable == "k", .x / 1e4, .x)),
        across(mean:q95, ~format(.x, digits=1, nsmall=2, scientific=FALSE)),
        estm = paste0(mean, " (", q5, "-", q95, ")")
    ) %>% 
    select(variable, estm, prior) %>% 
    pivot_wider(names_from=prior, values_from=estm, names_vary="slowest") %>% 
    slice(c(1,3,2,5,4,6)) %>% 
    mutate(variable = c("$w^{(I)}_{11}$", "$w^{(I)}_{22}$", "$w^{(I)}_{12}$",
                        "$\\gamma$", "$\\theta$", "$k/10^{4}$")) %>% 
    set_colnames(c("Param", "Prior", "Posterior")) %>% 
    kable(
        format = "latex",
        align = "crr",
        escape = FALSE
    ) %>%
    clipr::write_clip()
















