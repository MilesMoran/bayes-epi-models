library(hhh4contacts)
library(dplyr)

GROUPING <- c(1,2,2,4,4,2)

# noroBEall <- noroBE(by="all", flatten=TRUE, agegroups=GROUPING, timeRange=c("2011-w01","2016-w30")) ## the whole dataset
noroBEall <- noroBE(by="all", flatten=TRUE, agegroups=GROUPING, timeRange=c("2011-w27","2015-w26")) ## what's used in Meyer&Held
popBErbyg <- aggregateCountsArray(pop2011, dim=2, grouping=GROUPING)

DISTRICTS <- unique(stratum(noroBEall, 1))
G <- length(DISTRICTS)
STRATA <- unique(stratum(noroBEall, 2))
I <- length(STRATA)

Cgrouped <- contactmatrix(
    which = "reciprocal", # estimated by the Wallinga et al (2006) method
    type = "all",         # alternatively: "physical" only
    grouping = GROUPING   # age group specification
)
Cgrouped_norm <- (Cgrouped / rowSums(Cgrouped))

### MODEL FITS
## We estimate various hhh4() models with spatial power-law weights,
## population gravity and (power-adjusted) age-structured contact matrix
###

O <- neighbourhood(noroBEall)
all(O[1:12,1:12] == O[13:24,13:24])
neighbourhood(noroBEall) <- neighbourhood(noroBEall) + 1

DATAt <- list(t = epoch(noroBEall) - 1,
              christmas = 1*(epochInYear(noroBEall) %in% c(52, 1)))

## setup a model matrix with group indicators
MMG <- sapply(STRATA, 
              function (i) {
                index <- which(stratum(noroBEall, which = 2) == i)
                res <- col(noroBEall)
                res[] <- res %in% index
                res
              }, 
              simplify = FALSE, 
              USE.NAMES = TRUE)
str(MMG)

## setup model matrix with district indicators
MMR <- sapply(DISTRICTS, 
              function (g) {
                index <- which(stratum(noroBEall, which = 1) == g)
                res <- col(noroBEall)
                res[] <- res %in% index
                res
              }, 
              simplify = FALSE, 
              USE.NAMES = TRUE)
str(MMR)

## setup model matrix of group-specific seasonal terms
MMgS <- with(c(MMG, DATAt), unlist(lapply(X = STRATA,
                                          FUN = function (i) {
                                            gIndicator <- get(i)
                                            res <- list(gIndicator * sin(2 * pi * t/52),
                                                    gIndicator * cos(2 * pi * t/52))
                                            names(res) <- paste0(c("sin", "cos"), "(2 * pi * t/52).", i)
                                            res
                                          }), 
                                   recursive = FALSE, 
                                   use.names = TRUE)
             )
str(MMgS)

qSTRATA <- paste0("`", STRATA, "`")

## fit the endemic-only model
ma0 <- hhh4(noroBEall, 
            control = list(
                ## endemic formula: ~group + district + christmas + group:(sin+cos)
                end = list(f = reformulate(c(qSTRATA[-1], DISTRICTS[-1], "christmas",
                                             paste0("`", names(MMgS), "`")),
                                           intercept = TRUE),
                           offset = population(noroBEall) / rowSums(population(noroBEall))),
                # family = factor(stratum(noroBEall, which = 2)), # group-specific dispersion is used in Meyer-Held original paper
                family = "NegBin1", # single dispersion is what we use for parsimony
                data = c(MMG, MMR, DATAt, MMgS)
            ))

## fit the power-law model with the given contact matrix
ma_popPLC <- update(ma0,
                    ne = list(
                        ## epidemic formula: ~group + district + log(pop)
                        f = reformulate(c(qSTRATA[-1], DISTRICTS[-1], "log(pop)"),
                                        intercept = TRUE), 
                        weights = W_powerlaw(maxlag = 5, log = TRUE, normalize = FALSE,
                                             initial = c("logd" = log(2))),
                        scale = expandC(Cgrouped_norm, G),
                        normalize = TRUE),
                    data = list(pop = population(noroBEall)/rowSums(population(noroBEall)))
                    )

## AIC comparison
AIC(ma0, ma_popPLC)
cumsum(diff(AIC(ma0, ma_popPLC)$AIC)) |> t() |> t()

### fit power-adjusted contact matrix via profile likelihood
### CAVE: this takes a while (approx. 3 minutes)

ma_popPLCpower <- fitC(ma_popPLC, Cgrouped, normalize = TRUE, truncate = TRUE)

## AIC comparison
AIC(ma_popPLCpower, ma0)
diff(AIC(ma_popPLCpower, ma0)$AIC)

## model summary
summary(ma_popPLCpower, maxEV = TRUE, reparamPsi = TRUE,
        amplitudeShift = TRUE, idx2Exp = TRUE)
summary(ma_popPLCpower, maxEV = TRUE, reparamPsi = FALSE,
        amplitudeShift = FALSE, idx2Exp = c())



### Results from OG Meyer-Held Paper (unique dispersion for each age-group)
# overdisp.00-04                0.24±0.04 
# overdisp.05-14                1.98±0.49 
# overdisp.15-24                0.30±0.17 
# overdisp.25-44                0.03±0.04 
# overdisp.45-64                0.15±0.04 
# overdisp.65+                  0.40±0.03 

### Result from M-H model fit with single dispersion ψ
# overdisp 0.305±0.017 (0.288-0.322)

################################################################################ 
### Inspect age-group mixing weights estimates

plotC(Cgrouped)
powerC <- make_powerC(Cgrouped_norm, truncate = TRUE)
c(powerC(0.3243814)[1,1], powerC(0.4622776)[1,1], powerC(0.6001738)[1,1]) %>% format(digits=2, nsmall=2)
c(powerC(0.3243814)[2,2], powerC(0.4622776)[2,2], powerC(0.6001738)[2,2]) %>% format(digits=2, nsmall=2)
c(powerC(0.3243814)[2,1], powerC(0.4622776)[2,1], powerC(0.6001738)[2,1]) %>% format(digits=2, nsmall=2) # [2,1] because it's transposed

t(powerC(0.47)) # W_strata estimate




