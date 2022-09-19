### Extract output of interest
data_extract <- function(data, output){
  
  susceptibility <- rstan::extract(output, pars = "susc", inc_warmup=FALSE, permute=F)  
  susceptibilitytosecondary <- rstan::extract(output, pars="mono", inc_warmup=FALSE, permute=F)
  exp_reported_cases <- rstan::extract(output, pars = "exp_reported_cases", inc_warmup=FALSE, permute=F)
  posterior_exp_reported_cases <- exp_reported_cases
  lambda <- rstan::extract(output, pars = "lambda", inc_warmup=FALSE, permute=F)
  posterior_lambda <- lambda
  reportingrate <- rstan::extract(output, pars = "reporting_rate", inc_warmup=FALSE, permute=F)
  posterior_reportingrate <- reportingrate
  
  susc_median <- apply(susceptibility, 3, median) # median across chains and run
  susc_median <- matrix(1-susc_median, nrow= data$max_A, ncol= data$T) %>%  
    as_tibble() # arrange by year=col and row=nb f age groups and operation on 1-susceptibility
  colnames(susc_median) <- colnames(data$secondary_cases) # add year as colname
  # extract confidence intervals
  susc_ci <- apply(susceptibility, 3, function(x) quantile(x, c(0.025, 0.975))) # Quantiles across chains and run
  susc_ci_lower <- matrix(1-susc_ci[2,], nrow= data$max_A, ncol= data$T) %>%
    as_tibble() # arrange by year=col and row=nb f age groups
  susc_ci_upper <- matrix(1-susc_ci[1,], nrow= data$max_A, ncol= data$T) %>%  
    as_tibble() # arrange by year=col and row=nb f age groups
  colnames(susc_ci_lower) <- colnames(data$secondary_cases) # add year as colname
  colnames(susc_ci_upper) <- colnames(data$secondary_cases) # add year as colname
  
  # extract seropositivity for year 2010 to compare with Duncombe et al 2013 serostudy
  susc_2010 <- cbind(susc_median$`2010`,
                     susc_ci_lower$`2010`,
                     susc_ci_upper$`2010`,
                     age= c(1:data$max_A))
  colnames(susc_2010) <- c("susc_median", "susc_ci_lr", "susc_ci_ur", "agegroup")
  
  # keep last year of available data to plot susceptibility
  susc_median <- susc_median %>% select(last_col())
  susc_ci_lower <- susc_ci_lower %>% select(last_col())
  susc_ci_upper <- susc_ci_upper %>% select(last_col())
  
  susc <- cbind(susc_median,
                susc_ci_lower, 
                susc_ci_upper, 
                age= c(1:data$max_A))
  colnames(susc) <- c("susc_median", "susc_ci_lr", "susc_ci_ur", "agegroup")
  
  # extract median of susceptibility to secondary infection
  susctosecondary_median <- apply(susceptibilitytosecondary, 3, median) # median across chains and run
  susctosecondary_median <- matrix(susctosecondary_median, nrow= data$max_A, ncol= data$T) %>%
    as_tibble() # arrange by year=col and row=nb f age groups and operation on 1-susceptibility
  colnames(susctosecondary_median) <- colnames(data$secondary_cases) # add year as colname
  # extract confidence intervals
  susctosecondary_ci <- apply(susceptibilitytosecondary, 3, function(x) quantile(x, c(0.025, 0.975))) # Quantiles across chains and run
  susctosecondary_ci_lower <- matrix(susctosecondary_ci[1,], nrow= data$max_A, ncol= data$T) %>%
    as_tibble() # arrange by year=col and row=nb f age groups
  susctosecondary_ci_upper <- matrix(susctosecondary_ci[2,], nrow= data$max_A, ncol= data$T) %>%
    as_tibble() # arrange by year=col and row=nb f age groups
  colnames(susctosecondary_ci_lower) <- colnames(data$secondary_cases) # add year as colname
  colnames(susctosecondary_ci_upper) <- colnames(data$secondary_cases) # add year as colname
  
  # # keep last year of available data to plot susceptibility
  susctosecondary <- cbind(last(susctosecondary_median),
                           last(susctosecondary_ci_lower),
                           last(susctosecondary_ci_upper),
                           agegroup = c(1:data$max_A))
  colnames(susctosecondary) <- c("susctoscdry_median", "susctoscdry_ci_lr", "susctoscdry_ci_ur", "agegroup")
  
  # extract median of expected cases
  exp_reportedcases_median <- apply(exp_reported_cases, 3, median) # median across chains and run
  exp_reportedcases_median <- matrix(exp_reportedcases_median, nrow= data$AG, ncol= data$T) %>%  
    as_tibble() # arrange by year=col and row=nb f age groups and operation on 1-susceptibility
  colnames(exp_reportedcases_median) <- colnames(data$secondary_cases) # add year as colname
  exp_reportedcases_median <- mutate(exp_reportedcases_median, agegroup = data$lr_bound) # add a column with median age group
  # extract confidence intervals
  exp_reportedcases_ci <- apply(exp_reported_cases, 3, function(x) quantile(x, c(0.025, 0.975))) # Quantiles across chains and run
  exp_reportedcases_ci_lower <- matrix(exp_reportedcases_ci[1,], nrow= data$AG, ncol= data$T) %>%
    as_tibble() # arrange by year=col and row=nb f age groups
  exp_reportedcases_ci_upper <- matrix(exp_reportedcases_ci[2,], nrow= data$AG, ncol= data$T) %>%  
    as_tibble() # arrange by year=col and row=nb f age groups
  colnames(exp_reportedcases_ci_lower) <- colnames(data$secondary_cases) # add year as colname
  colnames(exp_reportedcases_ci_upper) <- colnames(data$secondary_cases) # add year as colname
  
  exp_reportedcases <- list(exp_reportedcases_median,exp_reportedcases_ci_lower, exp_reportedcases_ci_upper)
  
  # expectd incidence
  exp_inc <- rstan::extract(output, pars = "exp_inc", inc_warmup=FALSE, permute=F)
  exp_inc_median <- apply(exp_inc, 3, median) # median across chains and run
  exp_inc_median <- matrix(exp_inc_median, nrow= data$max_A, ncol= data$T) %>%  
    as_tibble() # arrange by year=col and row=nb f age groups and operation on 1-susceptibility
  colnames(exp_inc_median) <- colnames(data$secondary_cases) # add year as colname
  exp_inc_median <- mutate(exp_inc_median, agegroup = data$max_A) # add a column with median age group
  # extract confidence intervals
  exp_inc_ci <- apply(exp_inc, 3, function(x) quantile(x, c(0.025, 0.975))) # Quantiles across chains and run
  exp_inc_ci_lower <- matrix(exp_inc_ci[1,], nrow= data$max_A, ncol= data$T) %>%
    as_tibble() # arrange by year=col and row=nb f age groups
  exp_inc_ci_upper <- matrix(exp_inc_ci[2,], nrow= data$max_A, ncol= data$T) %>%  
    as_tibble() # arrange by year=col and row=nb f age groups
  colnames(exp_inc_ci_lower) <- colnames(data$secondary_cases) # add year as colname
  colnames(exp_inc_ci_upper) <- colnames(data$secondary_cases) # add year as colname
  
  exp_inc <- list(exp_inc_median,exp_inc_ci_lower, exp_inc_ci_upper)
  
  # extract reporting rate
  reporting_rate_median <- apply(reportingrate, 3, median)
  reporting_rate_ci <- apply(reportingrate, 3, function(x) quantile(x, c(0.025, 0.975))) # Quantiles across chains and run
  # reporting rate
  reporting_rate_ci_lower <- matrix(reporting_rate_ci[1,], nrow= data$T) %>%
    as_tibble() # arrange by year=col and row=nb f age groups
  reporting_rate_ci_upper <- matrix(reporting_rate_ci[2,], nrow= data$T) %>%  
    as_tibble() # arrange by year=col and row=nb f age groups
  #logit
  reporting_rate_logit <- rstan::extract(output, pars = "reporting_rate_logit", inc_warmup=FALSE, permute=F)
  posterior_reporting_rate_mean <- reporting_rate_logit
  reporting_rate_logit_median <- apply(reporting_rate_logit, 3, median)
  reporting_rate_logit_median <- matrix(reporting_rate_logit_median, nrow= (data$T)) %>%  
    as_tibble() # 
  reporting_rate_logit_ci <- apply(reporting_rate_logit, 3, function(x) quantile(x, c(0.025, 0.975))) # Quantiles across chains and run
  reporting_rate_logit_ci_lower <- matrix(reporting_rate_logit_ci[1,], nrow= data$T) %>%
    as_tibble() # arrange by year=col and row=nb f age groups
  reporting_rate_logit_ci_upper <- matrix(reporting_rate_logit_ci[2,], nrow= data$T) %>%  
    as_tibble() # arrange by year=col and row=nb f age groups
  reporting_rate_logit <- cbind(reporting_rate_logit_median,reporting_rate_logit_ci_lower, reporting_rate_logit_ci_upper)
  colnames(reporting_rate_logit) <- c("reporting_rate_logit_median", "reporting_rate_logit_ci_lr", "reporting_rate_logit_ci_ur")
  # hyperprior mean
  reporting_rate0 <- rstan::extract(output, pars = "reporting_rate0", inc_warmup=FALSE, permute=F)
  reporting_rate0_median <- apply(reporting_rate0, 3, median)
  reporting_rate0_ci <- apply(reporting_rate0, 3, function(x) quantile(x, c(0.025, 0.975))) 
  # extract confidence intervals
  reporting_rate0_ci_lower <- reporting_rate0_ci[1,]
  reporting_rate0_ci_upper <- reporting_rate0_ci[2,]
  reporting_rate0 <- cbind(reporting_rate0_median,reporting_rate0_ci_lower, reporting_rate0_ci_upper)
  colnames(reporting_rate0) <- c("reporting_rate0_median", "reporting_rate0_ci_lr", "reporting_rate0_ci_ur")
  
  reportingrate <- cbind(reporting_rate_median,reporting_rate_ci_lower, reporting_rate_ci_upper)
  colnames(reportingrate) <- c("reporting_rate_median", "reporting_rate_ci_lr", "reporting_rate_ci_ur")
  
  # extract FOI .
  lambda_median <- apply(lambda, 3, median)
  lambda_median <- matrix(lambda_median, nrow= (data$max_A+data$T)) %>%  
    as_tibble() # arrange by year=col and row=nb f age groups and operation on 1-susceptibility
  # extract confidence intervals
  lambda_ci <- apply(lambda, 3, function(x) quantile(x, c(0.025, 0.975))) # Quantiles across chains and run
  lambda_ci_lower <- matrix(lambda_ci[1,], nrow= (data$max_A+data$T)) %>% 
    as_tibble() # arrange by year=col and row=nb f age groups
  lambda_ci_upper <- matrix(lambda_ci[2,], nrow=(data$max_A+data$T)) %>% 
    as_tibble() # arrange by year=col and row=nb f age groups
  
  lambda_50_ci <- apply(lambda, 3, function(x) quantile(x, c(.25, .75))) # Quantiles across chains and run
  lambda_50_ci_lower <- matrix(lambda_50_ci[1,], nrow= (data$max_A+data$T)) %>% 
    as_tibble() # arrange by year=col and row=nb f age groups
  lambda_50_ci_upper <- matrix(lambda_50_ci[2,], nrow=(data$max_A+data$T)) %>% 
    as_tibble() # arrange by year=col and row=nb f age groups
  
  lambda <- cbind(lambda_median,lambda_50_ci_lower,lambda_50_ci_upper, lambda_ci_lower, lambda_ci_upper)
  colnames(lambda) <- c("lambda_median","lambda_50_ci_lr","lambda_50_ci_ur", "lambda_ci_lr", "lambda_ci_ur")
  
  # extract lambda priors
  # hyperprior mean
  lambda0_logit <- rstan::extract(output, pars = "lambda0_logit", inc_warmup=FALSE, permute=F)
  lambda0_logit_median <- apply(lambda0_logit, 3, median)
  lambda0_logit_ci <- apply(lambda0_logit, 3, function(x) quantile(x, c(0.025, 0.975))) 
  # extract confidence intervals
  lambda0_logit_ci_lower <- lambda0_logit_ci[1,]
  lambda0_logit_ci_upper <- lambda0_logit_ci[2,]
  lambda0_logit <- cbind(lambda0_logit_median,lambda0_logit_ci_lower, lambda0_logit_ci_upper)
  colnames(lambda0_logit) <- c("lambda0_logit_median", "lambda0_logit_ci_lr", "lambda0_logit_ci_ur")
  #logit
  lambda_logit <- rstan::extract(output, pars = "lambdaRE_logit", inc_warmup=FALSE, permute=F)
  posterior_lambda_logit <- lambda_logit
  lambda_logit_median <- apply(lambda_logit, 3, median)
  lambda_logit_median <- matrix(lambda_logit_median, nrow= (data$max_A+data$T)) %>%  
    as_tibble() # 
  lambda_logit_ci <- apply(lambda_logit, 3, function(x) quantile(x, c(0.025, 0.975))) # Quantiles across chains and run
  lambda_logit_ci_lower <- matrix(lambda_logit_ci[1,], nrow=(data$max_A+data$T)) %>%
    as_tibble() # arrange by year=col and row=nb f age groups
  lambda_logit_ci_upper <- matrix(lambda_logit_ci[2,], nrow= (data$max_A+data$T)) %>%  
    as_tibble() # arrange by year=col and row=nb f age groups
  lambda_logit <- cbind(lambda_logit_median,lambda_logit_ci_lower, lambda_logit_ci_upper)
  colnames(lambda_logit) <- c("lambda_logit_median", "lambda_logit_ci_lr", "lambda_logit_ci_ur")
  
  # extract cumul_FOI
  cum_lambda <- rstan::extract(output, pars = "cum_lambda", inc_warmup=FALSE, permute=F)
  cum_lambda_median <- apply(cum_lambda, 3, median)
  # extract confidence intervals
  cum_lambda_ci <- apply(cum_lambda, 3, function(x) quantile(x, c(0.025, 0.975))) # Quantiles across chains and run
  
  cum_lambda_median <- matrix(cum_lambda_median, nrow= data$max_A, ncol = data$T) %>%  
    as_tibble()
  colnames(cum_lambda_median) <- colnames(data$secondary_cases) # add year as colname
  cum_lambda_median <- mutate(cum_lambda_median, agegroup = c(1:data$max_A)) # add a column with median age group
  # extract confidence intervals
  cum_lambda_ci_lower <- matrix(cum_lambda_ci[1,], nrow=data$max_A, ncol = data$T) %>% 
    as_tibble() # arrange by year=col and row=nb f age groups
  cum_lambda_ci_upper <- matrix(cum_lambda_ci[2,], nrow=data$max_A, ncol = data$T) %>% 
    as_tibble() # arrange by year=col and row=nb f age groups
  
  colnames(cum_lambda_ci_lower) <- colnames(data$secondary_cases) # add year as colname
  colnames(cum_lambda_ci_upper) <- colnames(data$secondary_cases) # add year as colname
  cum_lambda <- list(cum_lambda_median, cum_lambda_ci_lower, cum_lambda_ci_upper)
  names(cum_lambda) <- c("cum_lambda_median", "cum_lambda_ci_lower", "cum_lambda_ci_upper")
  
  # extract disperion
  phi <- rstan::extract(output, pars = "phi", inc_warmup=FALSE, permute=F)
  phi_median <- apply(phi, 3, median)
  # extract confidence intervals
  phi_ci <- apply(phi, 3, function(x) quantile(x, c(0.025, 0.975))) # Quantiles across chains and run
  phi_ci_lower <- phi_ci[1]
  phi_ci_upper <- phi_ci[2]
  phi <- c(phi_median,phi_ci_lower,phi_ci_upper)
  names(phi) <- c("phi_median", "phi_ci_lower", "phi_ci_upper")
  
  # extract alpha (scale first age lambda)
  alpha <- rstan::extract(output, pars = "alpha", inc_warmup=FALSE, permute=F)
  posterior_alpha=alpha
  alpha_median <- apply(alpha, 3, median)
  # extract confidence intervals
  alpha_ci <- apply(alpha, 3, function(x) quantile(x, c(0.025, 0.975))) # Quantiles across chains and run
  alpha_ci_lower <- alpha_ci[1]
  alpha_ci_upper <- alpha_ci[2]
  alpha <- c(alpha_median,alpha_ci_lower, alpha_ci_upper)
  names(alpha) <- c("alpha_median", "alpha_ci_lower", "alpha_ci_upper")
  
  
  # model comparison 
  # Extract pointwise log-likelihood and compute LOO
  posterior_log_lik <- rstan::extract(output, pars = "log_lik", inc_warmup=FALSE, permute=F)
  loo2 <- loo(output, save_psis = TRUE)
  log_lik_1 <- extract_log_lik(output)
  waic_1 <- waic(log_lik_1)
  #loo2 <- loo::extract_log_lik(output, parameter_name = "log_lik")
  
  ## output
  list(data= data, susc=susc, 
       susc_2010=susc_2010,
       susctosecondary=susctosecondary, 
       exp_reportedcases=exp_reportedcases,
       exp_inc = exp_inc,
       lambda = lambda, reportingrate=reportingrate, 
       reporting_rate0=reporting_rate0,
       reporting_rate_logit=reporting_rate_logit,
       lambda0_logit=lambda0_logit,
       lambda_logit =lambda_logit,
       cum_lambda=cum_lambda, 
       phi = phi, alpha= alpha, posterior_log_lik= posterior_log_lik, loo2 = loo2, waic_1=waic_1,
       posterior_exp_reported_cases = posterior_exp_reported_cases,
       posterior_reportingrate = posterior_reportingrate,
       posterior_reporting_rate_mean=posterior_reporting_rate_mean,
       posterior_lambda=posterior_lambda,
       posterior_lambda_logit=posterior_lambda_logit,
       posterior_alpha = posterior_alpha)
}

############################################################
### Calculate mid-point of age group (needed to calculate likelihood)
roof_age <-
  c(0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85)
mid_age <- ceiling(roof_age[1:17] + diff(roof_age) / 2)

# Catalytic model function
run_catalytic_model <- function(cases, census, max_age, iter = 50000, chains = 4, thin= 10, AG, T, 
                                foi_mean, foi_sd, 
                                reporting_mean =  reporting_mean,
                                reporting_sd = reporting_sd, ...) {
  
  data_prov <-
    list(
      AG = nrow(cases), # A: the number of age classes - extracted from the number of age class in the case data
      max_A = max_age, # the max age
      T = ncol(cases), # T= the number of observation time points - years
      secondary_cases = cases, # the number of cases reported at each time point for each age group (until group 60 yo)
      pop = floor(census), #  total population at each time point for each age group
      lr_bound = seq(1,max_age, by=5),# age groups lower bounds
      ur_bound = seq(5,max_age+1, by=5), # age groups upper bounds
      foi_mean =  foi_mean,  # foi hyperparameter mean
      foi_sd = foi_sd, # foi SD (partitioned between hyperperameter and inter-annual variation)
      reporting_mean = reporting_mean,# reporting rate hyperparameter mean
      reporting_sd = reporting_sd) # reporting rate SD (partitioned between hyperperameter and inter-annual variation)
  
  ### Compile Stan code  
  comp2 <- stan(file=paste0("code/",stan_code), data=data_prov, chains = chains, ...) 
  
  ### Without starting values
  set.seed(.98342)
  #system.time(
  output <-
    stan(
      fit = comp2,
      data = data_prov,
      chains = chains,
      iter = iter,
      #warmup = warmup,
      thin = thin,
      control = list(adapt_delta = 0.9, max_treedepth = 12), ...
    )
  #)
  
  list(dat = data_prov, output = output)
}


##Set up prior for foi
## foi init
foi_init <- function(L, A, inc) {
  r0 <- 1 + L/A
  s <- 1 / r0
  inc / s
}

calc_incidence <- function(cases, pop, reporting_rate) {
  if (is.null(dim(pop))) {
    N <- sum(pop) # pop size
  } else {
    N <- sum(pop[, 1]) # pop size
  }
  total_cases <- sum(cases)
  nb_year <- ncol(cases) # years available of incidence data
  inc <- total_cases / nb_year / N / reporting_rate
  return(inc)
}

# for median age calculation
intervals <- seq(0, 80, by=5)
interval_medians <- (intervals[1:(length(intervals) - 1)] + intervals[2:length(intervals)])/2
age_gp <- interval_medians
