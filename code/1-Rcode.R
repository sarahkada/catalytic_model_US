# code from doi https://doi.org/10.7554/eLife.45474 continued and adapted by SK
rm(list=ls()); gc()
#pacman::p_unload('all')

# life expectancy average values in the US https://www.cdc.gov/nchs/fastats/life-expectancy.htm 
#### Fit models for specific countries
library(rstan)
library(coda)
library(ggplot2)
library(tidyverse)
library(loo)
#warning tidyr and stan have extract method so this is to avoid conflicts always use rstan::extract

options(mc.cores = parallel::detectCores()) # checks number of cores without having later to specify the cores argument
rstan_options(auto_write = TRUE) # extended packages to use stan
### Source function to prepare output
source("code/0-export-code.R")
iter = 5000
yearly_foi_sd = 0.2 # 20% upper bound for average FOI - decided after comparing 20%, 25% and 30% upper bounds (see notebooks)
life_exp <- 77; # life expectancy 
reporting_rate <- 0.1 # average reporting rate expected for all cases
# stan code:
stan_code <- "denguefoi_model.stan"

############################################################
### Catalytic model function and xtract output of interest -- PR data
load('data/cases_PR_2010_2019.Rdata') # dengue cases 
census_PR <- read.csv('data/census_PR_2010_2019.csv', header = TRUE) # population census
colnames(census_PR) <- gsub("X", "", colnames(census_PR)) # remove the X reas.csv adds to column header that have numbers
census_PR <- census_PR[,-1] # remove first column with age groups
j <- "PR"

## informed FOI prior
# calculate median age at secondary infection 
cases <- cases_PR_2010_2019[1:16,] %>% as_tibble() %>% 
  mutate(age_gp = age_gp) %>% 
  pivot_longer(!age_gp, names_to = "year", values_to = "cases") %>% 
  group_by(age_gp) %>%
  summarize(n = sum(cases))
mean_age <- cases %>% uncount(n) %>%
  pull(age_gp) %>%
  median() 

incidence <- calc_incidence(cases_PR_2010_2019, census_PR, reporting_rate)
foi_mean <- logitnorm::twCoefLogitnorm(foi_init(life_exp, mean_age, incidence), yearly_foi_sd)[1]
foi_sd <- logitnorm::twCoefLogitnorm(foi_init(life_exp, mean_age, incidence), yearly_foi_sd)[2] 

# run model
run <- run_catalytic_model(cases= cases_PR_2010_2019[1:16,], 
                           census = census_PR[1:16,],
                           max_age = 80,iter = iter,
                           foi_mean = foi_mean,
                           foi_sd = foi_sd)

model_output <- data_extract_multiyear(run$dat, run$output)
saveRDS(model_output,file= paste0("output/model_",j,".rds"))
