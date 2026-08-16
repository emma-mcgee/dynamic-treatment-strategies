#############################################################################
#    Defining and Estimating Effects of Dynamic Treatment Strategies        #
#    Exercise 2 - Estimating Effects                                        #
#    ISPE Course - August 2026                                              #
#    Programmers: Emma McGee & Rienna Russo                                 #
#############################################################################

# Part 0: Data Setup --------------------------------------------------------

# Clear workspace
rm(list = ls())

# Load the following packages. If they do not exist on your computer,
# this code will automatically install them into your default library. 
if(!require(data.table)) { install.packages("data.table"); require(data.table)}
if(!require(speedglm)) { install.packages("speedglm"); require(speedglm)}
if(!require(tidyverse)){ install.packages("tidyverse"); require(tidyverse)}
if(!require(ggplot2)) { install.packages("ggplot2"); require(ggplot2)}


# Set working directory
setwd("/pathway_to_your_files")


# Load the data
bestmed <- read.csv("bestmed.csv", header=TRUE)
setDT(bestmed)


# Part 1: Inverse Probability Weighting --------------------------------------
## a) ------------------------------------------------------------------------
#   	Fit logistic regression models to estimate the denominator 
#     probabilities for the non-stabilized IP weights. 

# Weight model for probability of initiating medication at baseline
model_tx_b <- speedglm(strategy ~
                         # Baseline time-fixed covariates
                         Female + SDI  +
                         Age + I(Age^2) + Insurance +
                         
                         # Time-varying covariates
                         eGFR + I(eGFR^2) +
                         Obesity + hba1c,
                       data = bestmed[Time==0, ], 
                       family = binomial())

summary(model_tx_b)

# Weight model for probability of continuing medication for SGLT2i strategy 
model_tx_0 <- speedglm(sglt2i ~
                         # Baseline time-fixed covariates
                         Female + SDI  +
                         Age + I(Age^2) + Insurance +
                         
                         # Time-varying covariates
                         eGFR + I(eGFR^2) +
                         Obesity + hba1c +
                         
                         # Time 
                         Time + I(Time^2),
                       data = bestmed[excused==0 & strategy==0 & Time >0 & lag1_sglt2i==1 & lag1_glp1==0, ],
                       family = binomial())

summary(model_tx_0)


# Weight model for probability of continuing medication for GLP1 strategy
model_tx_1 <- speedglm(glp1 ~
                         # Baseline time-fixed covariates
                         Female + SDI  +
                         Age + I(Age^2) + Insurance +
                         
                         # Time-varying covariates
                         eGFR + I(eGFR^2) +
                         Obesity + hba1c +
                         
                         # Time 
                         Time + I(Time^2),
                       data = bestmed[excused==0 & strategy==1 & Time > 0 & lag1_glp1==1 & lag1_sglt2i==0,],
                       family = binomial())
summary(model_tx_1)


# Obtaining predicted values from each model to use in the denominators
# At baseline
bestmed[, p_tx_b := NA_real_]  
bestmed[Time == 0,
        p_tx_b := predict(model_tx_b, newdata = .SD, type = "response")]
# Over follow-up
bestmed[, `:=`(
  p_tx_0 = predict(model_tx_0, newdata = bestmed, type = "response"),
  p_tx_1 = predict(model_tx_1, newdata = bestmed, type = "response"))]
bestmed[Time == 0, `:=`(
  p_tx_0 = NA_real_,
  p_tx_1 = NA_real_
)]


# Calculating time-specific weight contributions in each time period
bestmed[, weight_t:=fcase(
  excused==1, 1,                                              # excused 
  strategy==1 & Time==0, 1/(p_tx_b),                          # initiated GLP1 strategy at baseline
  strategy==0 & Time==0, 1/(1-p_tx_b),                        # initiated SGLT2 strategy at baseline
  strategy==1 & Time >0 & glp1==1 & sglt2i==0 & excused==0, 1/(p_tx_1),    # following GLP1 strategy, not excused
  strategy==0 & Time > 0 & sglt2i==1 & glp1==0 & excused==0, 1/(p_tx_0)    # following SGLT2 strategy, not excused
)]
# weight contributions will be set to missing at other non-adherent times


## c) ------------------------------------------------------------------------
#   	Take the cumulative product of the time-specific probabilities estimated 
#     in (a) and construct time-varying IP weights for each eligible individual.

# Make sure rows are ordered from first to last time within ID
setorder(bestmed, id, Time)

# Take the cumulative product of the time-specific weight contributions 
# Beginning from the first time (start of follow-up), separately for each ID
bestmed[, w_a := cumprod(weight_t), by = "id"]


# Part 2: Censoring ----------------------------------------------------------
## b) ------------------------------------------------------------------------
#   	Implement censoring rules for each treatment strategy to create dataset   
#     in which all eligible individuals follow 'assigned' strategy. 

# Create variable to censor people
bestmed[, cens_time_dynamic:=min(ifelse((strategy==1 & glp1!=1 & excused==0 | 
                                         strategy==1 & glp1==1 & sglt2i==1 & excused==0 |
                                         strategy==0 & sglt2i!=1 & excused==0 | 
                                         strategy==0 & sglt2i==1 & glp1==1 & excused==0 ),
                                        Time, 999)), by="id"]

# Truncate follow-up once they stopped following their 'assigned' strategy 
bestmed_cens <- bestmed[Time<cens_time_dynamic,]

# Check censoring times
table_cen <- bestmed_cens %>%
  filter(Time == 0) %>%
  pull(cens_time_dynamic) %>%
  table()

print(table_cen)

# Check weights and censored dataset construction
View(bestmed_cens[,c("id","Time","cens_time_dynamic","strategy","excused",
                "glp1","sglt2i","p_tx_b","p_tx_1","p_tx_0","weight_t","w_a")])


## c) ------------------------------------------------------------------------
#   	Truncate IP weights at the 99th percentile.  

# Summarize weight distribution
summary(bestmed_cens$w_a)

# Truncate weights
bestmed_cens[, w_a:=ifelse(
  w_a>=quantile(w_a, p=0.99),
  quantile(w_a, p=0.99),
  w_a)]   

# Summarize weight distribution after truncating
summary(bestmed_cens$w_a)



# Part 3: Dynamic Marginal Structural Model ----------------------------------
## a) ------------------------------------------------------------------------
##    Fit an IP weighted pooled logistic regression model for the outcome of 
##    MACE
model_dynamic_msm <- speedglm(mace ~ Time*strategy + I(Time^2)*strategy,
                              family = binomial(),
                              data = bestmed_cens,
                              weights = w_a)
summary(model_dynamic_msm)

## b) ------------------------------------------------------------------------
##    Using the parameter estimates from the model in (a), estimate the 
##    marginal risk of MACE at each month of follow-up under each strategy.

# Create dataset of all timepoints (0–47) under each strategy (0,1) 
risk_results <- data.table(
  Time=rep(seq(0, 47), 2), 
  strategy=c(rep(0, 48), rep(1, 48)))

# Estimate discrete-time hazards under each strategy
risk_results[, hazard:= predict(model_dynamic_msm, type="response", newdata=risk_results)]

# Estimate survival from cumulative product of (1-hazard) for each strategy
risk_results[, survival:=cumprod(1-hazard), by=strategy]

# Estimate risks from survival probabilities
risk_results[, risk:=1-survival]

## c) ------------------------------------------------------------------------
##    Using the risks from (b), construct adjusted risk curves for the outcome
##    of MACE under each strategy.

# Prepare data
# Shift Time by +1 (outcomes appearing in interval k represent those that happened in interval k+1)
# Add Time=0, risk=0 for each strategy 
risk_plot <- risk_results %>% 
  mutate(Time=Time+1) %>% 
  add_row(Time=0, strategy=0, risk=0) %>%
  add_row(Time=0, strategy=1, risk=0) %>%
  arrange(strategy, Time)

# Plot risk curves
ggplot(risk_plot,
       aes(x = Time, y = risk,
           color    = factor(strategy),
           linetype = factor(strategy))) + 
  geom_line() + 
  xlab("Months") + 
  scale_x_continuous(breaks = seq(0, 48, 6)) +
  coord_cartesian(xlim = c(0, 48)) +
  ylab("Risk (%)") + 
  scale_y_continuous(limits = c(0, 0.075),
                     breaks = c(0, 0.025, 0.05, 0.075),
                     labels = c("0.0%", "2.5%", "5.0%", "7.5%")) + 
  theme_bw() +
  scale_color_manual(name   = "Strategy",
    values = c("1" = "#56B4E9", "0" = "#000080"), 
    labels = c("1" = "GLP-1RA", "0" = "SGLT2i")
  ) +
  scale_linetype_manual(
    name   = "Strategy",
    values = c("1" = "dashed", "0" = "solid"),
    labels = c("1" = "GLP-1RA", "0" = "SGLT2i")
  ) +
  theme(legend.position = "bottom")

## d) ------------------------------------------------------------------------
##    Using the risks from (b), estimate the 4-year risk difference 
##    and risk ratio for MACE

# 4-year risks, risk difference, and risk ratio
risk1 <- risk_results[Time==47 & strategy==1,]$risk
risk0 <- risk_results[Time==47 & strategy==0,]$risk
results <- c(risk1, risk0, risk1-risk0, risk1/risk0)

round(results, 4)


