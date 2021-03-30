
# Boiler plate stuff
# Start with empty environment 
rm(list =ls())
gc()

# Set working directory for project. Relative paths are in reference to here.
setwd("~/Desktop/NVASA_symposium")

# Options + Settings ------------------------------------------------------
options(scipen=999)

library(dplyr)
library(tidyr)
library(lubridate)
library(ggplot2)


## A brief background...

# ADM is a consulting company, which analyzes energy use and the effectiveness of
# energy saving programs for various U.S. utility companies. 
# 
# For reference, the average electricity consumption for a Nevada household is about 
# 31 kWh per day, or 935 kWh per month 
# (according to [electricitylocal.com](https://www.electricitylocal.com/states/nevada/)).

# >My data show that Northern Nevada in particular uses an average of `23 kWh` per day, 
# so about `700 kWh` per month, and `8,400 kWh` per year. 

  
  ## An example analysis
  
#   We will utilize what's called a **difference-in-differences** (DiD) approach-- which
# essentially describes the difference between the _control_ and _treatment_ groups'
# usage, in comparison to the data collected during the period prior _(pre)_ to treatment, 
# those collected after treatment _(post)_ or after intervention began. 

  
  
  ### First, we load in some tamed, slimmed down, and anonymous-ed monthly billing data:
# Load in billing data 

#' TODO replace with a public dataset for example
billing_data <- readRDS("./data/billing/model_data_slim.RDS")
head(billing_data)


# Get list of participants
participants <- billing_data %>% 
  distinct(ACCT, TREATMENT)

participants %>% group_by(TREATMENT) %>% count()


#Then, we **match** a control group

# Here, I have 5000 potential control group members. This means they are expected
# to have been set aside and not be targeted for any energy saving programs previous
# to that time. Usually, this would need to be verified, but I've already don't that...
# 
# Onward!
#   
#   The trouble is, that this particular set of controls were not chosen in relation
# to these treatments. Thus, in order to make sure we're comparing apples-to-apples,
# we choose control group members with similar average trends to specific treatment 
# group members before treatment is a factor.


create_seasons <- function(DATASET){
  # create average usage data set by season
  temp_data <- 
    DATASET %>% 
    mutate(SEASON = if_else(month(READING_END_DATE) %in% c(6:9), "Summer", "Winter")) %>%
    group_by(ACCT, SEASON) %>% 
    summarise(KWHD = mean(ADC, na.rm = TRUE),
              CDDD = mean(AVG_CDD, na.rm = TRUE),
              HDDD = mean(AVG_HDD, na.rm = TRUE),
              TREATMENT = first(TREATMENT))
    # data.table()
  
  gathered <- 
    temp_data %>% 
    pivot_longer(names_to = "TYPE", values_to = "VALUE", cols = KWHD:HDDD) # replaced gather()
  
  gathered$CAT <- paste0(gathered$SEASON, ".", gathered$TYPE)
  gathered$SEASON <- NULL  
  gathered$TYPE <- NULL  
  
  spread_data <- 
    gathered %>%
    pivot_wider(names_from = CAT, values_from = VALUE) #replaced spread
  
  data <- 
    spread_data %>% 
    mutate(KWHDN.Summer = Summer.KWHD/Summer.CDDD, # Normalize average usage data with weather variables
           KWHDN.Winter = Winter.KWHD/Winter.HDDD) 
  
  data <- 
    na.omit(data) %>% 
    filter(!(is.infinite(KWHDN.Summer) | is.infinite(KWHDN.Winter))) # Prevent NAs
  
    return(data)
}

# We will take our billing data, to filter for everything before treatment started.
# Then, prepare for **Matching** by differentiating seasonal usage--create normalized usage 
# terms for which to match usage.

# Matching based on usage only  ----------------------------------------
  # Get pre data
  pre_data <- filter(billing_data, POST == 0)
  
  # Prep for matching
  prepped <- create_seasons(pre_data) # lose a few here--for div/0 in the seasons calc
  
  head(prepped)
  

# We will utilize the R `MatchIt` package here, which has a myriad of options with
# ample statistical explanation for each.
# See the [MatchIt website](https://kosukeimai.github.io/MatchIt/index.html)
# and in particular, the "Reference" tab. 

# Load packages
library(MatchIt)
library(optmatch) # for optimal, vs nearest neighbor methods

# Match treatment and controls 
set.seed(1234) # seed matters for reproducibility due to "nearest neighbor" method

matched <- matchit(formula = TREATMENT ~ KWHDN_Summer + KWHDN_Winter, data = prepped, 
                   method = "full", distance = "glm", discard = "both")
                   # method = "nearest", distance = "glm", distance.options = list(),
                   # discard = "none")
# discard = "both" could result in better matches in general (as it will only keep the treatments for which we find good control matches), but no difference in results here.

# Built in function for formatting "longer" dataset
df_matched <- MatchIt::match.data(matched)[1:ncol(prepped)]
head(df_matched)

# count our new participants
participants_matched <-  df_matched %>% select(TREATMENT, ACCT) 
participants_matched %>% group_by(TREATMENT) %>% count()

  
#### Let's test those matches!
# Test Matches 
# set the names of our covariate variables
cov <- c("KWHDN_Summer","KWHDN_Winter") 
# 'apply' (loop over) the two season columns and perform a Student's t-test to
# compare usage between controls/treatments in our declared Winter/Summer months
t_test_seasons <- lapply(cov, function(X) {t.test(dt_matched[, X] ~ dt_matched$TREATMENT)})
# label the output of both tests for clarity
names(t_test_seasons) <- c("KWHDN_Summer","KWHDN_Winter")

t_test_seasons

# or more simply put:
t.test(df_matched$KWHDN_Summer ~ df_matched$TREATMENT)
t.test(df_matched$KWHDN_Winter ~ df_matched$TREATMENT)

  
  ### Regression Model
  
  # Here we will use the `lmer()` function from the `lme4` package. 
  # **Linear Mixed Effects Regression models** are built on the foundation of linear 
  # models, such that they are measuring the regression from the equation line:
  #   
  #   $AEC_{i,e}$ = $\beta_{1}CDD_{i,t} + \beta_{2}HDD_{i,t} + \beta_{3}Post_{i,t} + \beta_{4}Post_{i,t}*Treat_{i,t} + \beta_{5}Post_{i,t}*CDD_{i,t} + \beta_{6}Post_{i,t}*HDD_{i,t} + \alpha Customer_{i} +  E_{i,t}$
  #   
  #   These _interaction terms_ serve to incorporate both the _fixed effects_ and _random effects_ that
  # we need our model to consider. We need to account for for the confounders or casual terms, so that our results
  # are not just anecdotal evidence. (E.g., if you used less energy this year than last,
  #                                   is it only because the summer weather was more mild?)
  
  library(lme4)
  library(MuMIn) # for r squared

  #### But first we do aother statstical check:
  model_data <- inner_join(billing_data, participants_matched)
  
  # Get Pre Data -----------------------------------------------------------
  pre_data <- filter(model_data, POST == 0)
  
  # T-Test Pre ADC 
  t.test(pre_data[pre_data$TREATMENT == 0, ]$ADC, 
         pre_data[pre_data$TREATMENT == 1, ]$ADC, conf.level = 0.90)

  
  #### ... and a visual one!
  
  pre_data <- pre_data %>% mutate(DATE = format(as.Date(READING_END_DATE), "%Y-%m"))
  pre_data$TREATMENT <- as.factor(pre_data$TREATMENT)
  
  # Check plot
  pre_data %>%
    group_by(DATE, TREATMENT) %>%
    summarise(MEAN_ADC = mean(ADC, na.rm = TRUE)) %>%
    ggplot(., aes(x = as.factor(DATE), y = MEAN_ADC, colour = TREATMENT, linetype = TREATMENT)) +
    geom_line(aes(group = TREATMENT, size = TREATMENT)) +
    ylim(0, 40) +
    theme_bw() +
    scale_color_manual(name = "TREATMENT",
                       values = c('grey60', 'black'),
                       labels = c("Control", "Treatment")) + 
    scale_linetype_manual(name = "TREATMENT",
                          values = c("solid", "dashed"),
                          labels = c("Control", "Treatment")) +
    scale_size_manual(name = "TREATMENT",
                      values = c(2, 1),
                      labels = c("Control", "Treatment")) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.title = element_blank(),
          text = element_text(size = 14)) +
    labs(y = "Average Daily Consumption (kWh)",
         x = "") 

  
  #### It's time...
  mix_fit <- lmer(ADC ~ AVG_HDD + AVG_CDD + POST +
                    AVG_HDD:POST + AVG_CDD:POST +
                    TREATMENT:POST + 
                    (1 | ACCT), data = model_data)
  
  # Raw output from the Linear Mixed Regression Model
  mix_fit
  
  # Oddly, `summary` actually offers more statistical information. I also like to take
  # a look at the $r^2$ terms and the p-value.

  #results <- run_did_model(model_data)
  r2 <- r.squaredGLMM(mix_fit)
  pvalue <- drop1(mix_fit, test="Chisq")
  results <- summary(mix_fit)
  
  results
  r2 # conditional and marginal
  pvalue
  

  # Let's focus on this part:
  # 
  #                    Estimate  Std. Error    t value
  #   POST:TREATMENT 

# This is the interaction term representing the _post_ treatment period, for the 
# _treatment_ participants. Because the dataset has "dummy terms" for each of these
# variables, the coefficients only persist for `TREATMENT == 1` and `POST == 1`, 
# thus 0 for any other data. 
# 
# These results can be interpreted as follows:
# >Each participant in the treatment group saves about **0.52 kWh of energy per day**,
# for an annual total of about **190 kWh** per person. 
# 
# Recall that our annual average usage for NV was **8,400 kWh**, and **22 kWh** 
# per day per participant during the pre-treatment period for these folks in particular. 

# We consider these statistically significant results in part because the "t value",
# or Student's* test statistic (abs value 8.12) is greater than 1.645: the critical value for a 
#   two-sided t-test at the $\alpha$ = 0.10 level (90% confidence) for a large dataset 
#   (many degrees of freedom). 
  
  
  
# Also, as a "sanity check", we can check to make sure 0 isn't in the CI for post data:
post_data <- filter(model_data, POST == 1)

# T-Test Pre ADC 
t.test(post_data[post_data$TREATMENT == 0, ]$ADC, 
       post_data[post_data$TREATMENT == 1, ]$ADC, conf.level = 0.90)


#----


# #### In conclusion, we have sufficient evidence to _reject the null hypothesis_---
# that is, the assumption that the treatment and control groups are equivalent), and can 
# affirm there is a statistical difference in the electrical consumption between 
# the treatment and control groups. Indeed, the treatment group appears to use less
# energy than the control group and less than they themselves did prior to treatment.
# Thus, there is evidence that this energy saving intervention is a success.
# It should be noted that these savings are fairly minute in comparison to daily usage, however.... 
# (note average percent reduction)


