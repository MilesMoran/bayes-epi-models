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
setwd("..") # navigate up to `inf-noro` dir

################################################################################
### Load & summarize Stan draws
################################################################################

fnames_full_prior <- list.files("output", "^priors-be-noro-model-full-\\d+\\.csv$", full.names=TRUE)
fnames_full_postr <- list.files("output", "^be-noro-model-full-\\d+\\.csv$", full.names=TRUE)
fnames_redu_prior <- list.files("output", "^priors-be-noro-model-reduced-\\d+\\.csv$", full.names=TRUE)
fnames_redu_postr <- list.files("output", "^be-noro-model-reduced-\\d+\\.csv$", full.names=TRUE)

variables <- c("end_1", "end_christmas", "end_strata_raw", "end_geogs_raw", 
               "end_sin", "end_cos", "ne_log_pop", "ne_1", "ne_strata_raw", 
               "ne_geogs_raw", "neweights_logd", "overdisp", "W_strata")
variables_full <- c(variables, "theta", "C_comp", "u")
variables_redu <- c(variables, "logpower", "kappa")

draws_full_prior <- read_cmdstan_csv(fnames_full_prior, variables_full)
draws_full_postr <- read_cmdstan_csv(fnames_full_postr, variables_full)
draws_redu_prior <- read_cmdstan_csv(fnames_redu_prior, variables_redu)
draws_redu_postr <- read_cmdstan_csv(fnames_redu_postr, variables_redu)

summ_full_prior <- summarize_draws(draws_full_prior$post_warmup_draws)
summ_full_postr <- summarize_draws(draws_full_postr$post_warmup_draws)
summ_redu_prior <- summarize_draws(draws_redu_prior$post_warmup_draws)
summ_redu_postr <- summarize_draws(draws_redu_postr$post_warmup_draws)

summary_df <- rbind(mutate(summ_full_prior, model="full", prior=TRUE),
                    mutate(summ_full_postr, model="full", prior=FALSE),
                    mutate(summ_redu_prior, model="reduced", prior=TRUE),
                    mutate(summ_redu_postr, model="reduced", prior=FALSE))

################################################################################
### Cursory inspections: make sure nothing is horribly wrong 
################################################################################

summ_full_postr %>% 
    arrange(desc(rhat)) %>% 
    select(variable, mean, median, rhat, ess_bulk, ess_tail) %>%  
    head(n=10L)

summ_redu_postr %>% 
    arrange(desc(rhat)) %>% 
    select(variable, mean, median, rhat, ess_bulk, ess_tail) %>%  
    head(n=10L)


################################################################################
### Plot #1: Matrix plot of prior and posterior for `W_strata`
################################################################################

strata <- c("0–4", "5–14", "15–24", "25–44", "45–64", "65+") %>% 
            factor(., levels=.)

plot_wstrata <-
    summary_df %>% 
    subset(grepl("W_strata", variable)) %>% 
    select(variable, mean, q5, q95, model, prior) %>% 
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
        facet_wrap(
            model~rev(prior), 
            labeller=\(x){ 
                data.frame(c("Full Model Prior", "Full Model Posterior",
                             "Reduced Model Prior", "Reduced Model Posterior")) 
        }) + 
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
       scale=1.5, dpi=300, width=7, height=7, units="in")


################################################################################
### Table #1: Prior and posterior summaries for a few primary estimands
################################################################################

### Posterior Medians ###

summary_df %>% 
    subset(variable %in% c("theta", "overdisp")) %>% 
    select(model, prior, variable, median, q5, q95) %>% 
    mutate(
        across(median:q95, ~format(.x, digits=1, nsmall=2, scientific=FALSE)),
        estm = paste0(median, " (", q5, "-", q95, ")")
    ) %>% 
    select(model, prior, variable, estm) %>% 
    pivot_wider(names_from=c(model,prior), values_from=estm, names_vary="slowest") %>% 
    mutate(variable = c("$\\psi$", "$\\theta$")) %>% 
    set_colnames(c("Param", "Prior-F", "Posterior-F", "Prior-R", "Posterior-R")) %>% 
    kable(
        format = "latex",
        align = "crrrr",
        escape = FALSE
    ) %>% 
    clipr::write_clip()

### Posterior Means ###

summary_df %>% 
    subset(variable %in% c("W_strata[1,1]", "W_strata[2,2]", "W_strata[1,2]", "kappa")) %>% 
    select(model, prior, variable, median, q5, q95) %>% 
    mutate(
        across(median:q95, ~format(.x, digits=1, nsmall=2, scientific=FALSE)),
        estm = paste0(median, " (", q5, "-", q95, ")")
    ) %>% 
    select(model, prior, variable, estm) %>% 
    pivot_wider(names_from=c(model,prior), values_from=estm, names_vary="slowest") %>% 
    slice(c(1,3,2,4)) %>% 
    mutate(variable = c("$w^{(I)}_{11}$", "$w^{(I)}_{22}$", "$w^{(I)}_{12}$", "$\\kappa$")) %>% 
    set_colnames(c("Param", "Prior-F", "Posterior-F", "Prior-R", "Posterior-R")) %>% 
    kable(
        format = "latex",
        align = "crrrr",
        escape = FALSE
    ) %>% 
    clipr::write_clip()


