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
iter = 5000#50000
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
j <- "PR" # AS

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

# reporting parameters
reporting_mean = -2.2; reporting_sd= 0.7; # corresponds to 10% reporting median and 30% upper bound

# run model
run <- run_catalytic_model(cases= cases_PR_2010_2019[1:16,], 
                           census = census_PR[1:16,],
                           max_age = 80,iter = iter,
                           foi_mean = foi_mean,
                           foi_sd = foi_sd, 
                           reporting_mean = reporting_mean, 
                           reporting_sd = reporting_sd
                           )

model_output <- data_extract(run$dat, run$output)
saveRDS(model_output, file= paste0("output/model_",j,".rds"))

# some visualizations
j <- "PR" 
model_pr <- readRDS(file= paste0("output/model_",j,".rds"))

# plot seroprevalence
PR_ci <- model_pr$susc %>% as_tibble() %>% select(susc_ci_lr, susc_ci_ur, agegroup)

median_susceptibility <- model_pr$susc %>% filter(agegroup %in% c(1,2,3,4,5,7,9,11,13,15,17,19,21,23,25,30,35,40,45,50,55,60,65,70,75,80))

ggplot(median_susceptibility) +
  geom_line(aes(x = agegroup, y = susc_median, colour = "steelblue", alpha = .8), linetype = 1, size = .5) +
  geom_point(aes(x = agegroup, y = susc_median, alpha = .8),colour="white", size=2, shape=21) +
  geom_ribbon(data = PR_ci,
              aes(x =  agegroup, ymin = susc_ci_lr, ymax = susc_ci_ur), 
              inherit.aes = FALSE,
              fill = "steelblue",
              alpha = .4
  ) +
  scale_x_continuous(breaks = seq(0, 79, by =10)) +
  ylim(0, 1) +
  ggthemes::geom_rangeframe() +
  ylab("Proportion seropositive in PR") + xlab("Age") +
  ggthemes::theme_clean()  + theme(
    axis.line.x = element_line(
      colour = "black",
      size = 0.5,
      linetype = "solid"
    ),
    axis.line.y = element_line(
      colour = "black",
      size = 0.5,
      linetype = "solid"
    ),
    plot.background = element_blank(),
    legend.position = "top",
    legend.background = element_rect(color = NA),
    legend.text = element_text(size = 8)
  )

### Fit to reported cases 
prep_dat <- function(data){
  years <- data$data$T
  l_age <- data$data$lr_bound
  u_dage <- data$data$ur_bound
  indage <- (l_age+u_dage)/2-1
  # prepare cases
  cases <- as_tibble(data$data$secondary_cases) %>%
    mutate(., agegroup = data$data$lr_bound) # data frame with raw case counts
  cases <- gather(cases, names(cases)[1:years], key= year, value = cases)
  #gather posterior cases by year
  exp_reportedcases_median <- gather(data$exp_reportedcases[[1]], names(data$exp_reportedcases[[1]])[1:years], key=year, value = cases)
  # Low CI
  exp_reportedcases_ci_lower <-  data$exp_reportedcases[[2]] %>%
    mutate(., agegroup = data$data$lr_bound) 
  exp_reportedcases_ci_lower <- gather(exp_reportedcases_ci_lower, names(exp_reportedcases_ci_lower)[1:years], key=year, value = ci_lwr)
  # # Upper CI
  exp_reportedcases_ci_upper <- data$exp_reportedcases[[3]] %>%
    mutate(., agegroup = data$data$lr_bound) 
  exp_reportedcases_ci_upper <- gather(exp_reportedcases_ci_upper, names(exp_reportedcases_ci_upper)[1:years], key=year,value = ci_upr)
  
  #  list of lists
  exp_reportedcases <- exp_reportedcases_median %>% full_join(exp_reportedcases_ci_lower, by = c("agegroup", "year")) %>% full_join(exp_reportedcases_ci_upper, by = c("agegroup", "year"))
  
  return(list(cases = cases, exp_reportedcases =exp_reportedcases))
}

plot_fit_to_reported_cases <- function(data, plot_title = "Dengue case counts in", location = location, max_y_axis=max(cases_reported$cases)){
  
  expectedreportedcases <- prep_dat(data)$exp_reportedcases 
  cases_reported <- prep_dat(data)$cases
  
  p <- ggplot(data= expectedreportedcases, aes(y = cases, x= agegroup , group = factor(year))) +
    geom_line(color= "steelblue")+
    geom_ribbon(
      aes(x= agegroup, ymin = ci_lwr, ymax = ci_upr, group = factor(year)),
      fill = "steelblue",
      alpha = .5
    )  +
    geom_point(data = cases_reported, aes(x= agegroup, y = cases, group = factor(year), color = "Observed\ncases"),
               size=.7, shape=19,  alpha = .8) +
    facet_wrap(~year, ncol = 3) +
    #scale_fill_manual(values=c("palevioletred3","steelblue"))+
    scale_color_manual("",values=c('darkblue'),
                       labels=c("Observed\ncases"))+
    ylim(0,max_y_axis)+
    ylab(paste(plot_title, location)) + xlab("Age") +
    ggthemes::theme_clean()  + theme(axis.line.x = element_line(colour = "black", size = 0.5, 
                                                      linetype = "solid"), axis.line.y = element_line(colour = "black", 
                                                                                                      size = 0.5, linetype = "solid"), plot.background = element_blank(),
                           legend.background = element_rect(color = NA)) # this removes the black line surrounding
  return(p)
}
plot_fit_to_reported_cases(model_pr, location = "Puerto Rico")

## Figure 5. Long term FOI 
max_A <- 80; last_year_of_data <- 2019
prep_foi <- function(data1){
  # get years
  year <- seq(last_year_of_data, by=-1, length.out = max_A+data1$data$T)
  #date <- colnames(data1$data$secondary_cases)
  date <- paste0(year,"-01-01")
  
  # lambda 
  lambda <- data1$lambda %>% mutate(year = rev(date)) %>%
    mutate(location =data1$name)
  # 
  return(lambda)
}

foi <- prep_foi(model_pr)
ggplot(foi)+
  geom_line(aes(x = as.Date(year), y = lambda_median))+
  geom_ribbon(
    aes(x = as.Date(year), ymin = lambda_50_ci_lr, ymax = lambda_50_ci_ur),
    alpha = .6,
    fill= "steelblue"
  )  +
  geom_ribbon(
    aes(x = as.Date(year), ymin = lambda_ci_lr, ymax = lambda_ci_ur),
    alpha = .4, 
    fill = "steelblue"
  ) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 3))

## reporting rate plots
years <- colnames(model_pr$data$secondary_cases)
reporting <- model_pr$reportingrate %>% mutate(year = years)
reporting %>%
  ggplot(aes(y = reporting_rate_median, ymin = reporting_rate_ci_lr, ymax = reporting_rate_ci_ur))+
  geom_pointrange(aes(x = year, color = "steelblue"), alpha = .6)+
  ylab("Reporting rate") + xlab("Location") +
  ggthemes::theme_clean()  + theme(axis.text.x = element_text(angle = 45, hjust=1),axis.line.x = element_line(colour = "black", size = 0.5,
                                                                                                    linetype = "solid"), axis.line.y = element_line(colour = "black",
                                                                                                                                                    size = 0.5, linetype = "solid"), plot.background = element_blank(), legend.title=element_text(size=8), legend.position = "top", legend.background = element_rect(color = NA),legend.text=element_text(size=8)) # this removes the black

