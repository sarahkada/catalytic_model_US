 //This code uses a negbin distribution of cases. Secondary cases only
// Reporting rate and FOI with hyperpriors
data {
  int <lower=0> AG; //the number of age classes
  int <lower=0> max_A; //the max age
  int <lower=0> T; //the number of observation time points
  int <lower=0>	secondary_cases[AG, T]; //the number of cases reported at each time point for each age group
  int <lower=0> pop[AG, T]; // the total population at each time point for each age group
  int <lower=0> lr_bound[AG]; // lower age groups bounds
  int <lower=0> ur_bound[AG]; // upper age groups bounds
  real foi_mean;// hyperparameter lambda mean prior
  real foi_sd;// hyperparameter lambda sd prior
  real reporting_mean;// hyperparameter reporting rate mean prior
  real reporting_sd;// hyperparameter reproting rate sd prior
}

parameters {
  real <lower=0> invroot_phi; // overdispersion 
  real reporting_rate0; // reporting rate hyperparameter
  real reporting_rate_logit[T]; // time varying reporting rate
  real lambda0_logit;
  real lambdaRE_logit[max_A + T]; // yearly FOI
  real <lower=0, upper=1> alpha; // proportion of first age group experiencing lambda[t-1]
}

transformed parameters {
  real <lower=0, upper=1> lambda[max_A + T]; // yearly FOI
  real <lower=0> cum_lambda[max_A, T];
  real reporting_rate[T]; // reporting rate
  real mono[max_A, T];
  real susc[max_A, T];
  real exp_inc[max_A, T];
  real exp_reported_cases[AG, T];
  real exp_inc_grouped[AG, T];
  real <lower=0> phi;
  
  phi = 1/(invroot_phi)^2;
  
  for (t in 1:T) {
    reporting_rate[t] = inv_logit(reporting_rate_logit[t]); // logit transform
  }
  
  for (t in 1:(max_A + T)) {
    lambda[t] = inv_logit(lambdaRE_logit[t]);
  }
  
  // Cumulative FOI by age group and year considered (T)
  cum_lambda[1, 1] = alpha * lambda[max_A];
  for (a in 2:max_A) {
    cum_lambda[a, 1] = sum(lambda[(max_A-a+1):max_A]);
  }
  for (t in 2:T) {
    cum_lambda[1, t] = alpha * lambda[max_A+t-1];//first age group born with no past exposure
    for (a in 2:max_A) {
      cum_lambda[a, t] = cum_lambda[a-1, t-1] + lambda[max_A+t-1];
    }
  }
  
  for (t in 1:T) {
    for (a in 1:max_A) {
      susc[a, t] = exp(-4 * cum_lambda[a, t]) ; // Eqn 1
      mono[a, t] = 4 * exp(-3 * cum_lambda[a, t]) * (1 - exp(-cum_lambda[a, t])) ; // Eqn 2
      exp_inc[a, t] = 3 * lambda[max_A+t] * mono[a, t];  //eqn3 : exp_inc : expected incidence
    }
  }
  
  for (t in 1:T) {
    for (a_gr in 1:AG) {
      exp_inc_grouped[a_gr, t] = mean(exp_inc[lr_bound[a_gr]:ur_bound[a_gr], t]);
      exp_reported_cases[a_gr, t] = exp_inc_grouped[a_gr, t] * reporting_rate[t] * pop[a_gr, t]; // eqn 4
    }
  }
  
}

model {
  // prior distribution
  
  alpha ~ beta(2, 1);
  
  reporting_rate0 ~ normal(reporting_mean, sqrt(reporting_sd^2 / 2));// 0,1
  for (t in 1:T) {
    reporting_rate_logit[t] ~ normal(reporting_rate0, sqrt(reporting_sd^2 / 2)); 
  }
  
  lambda0_logit ~ normal(foi_mean, sqrt(foi_sd ^2 / 2));
  // logitnorm::twCoefLogitnorm(0.05, 0.3, perc = 0.975)
  for (i in 1:(max_A+T)) {
    lambdaRE_logit[i] ~ normal(lambda0_logit, sqrt(foi_sd ^2 / 2)); 
  }
  //plot(density(boot::inv.logit(rnorm(1000, rnorm(1000,-5, 1),1))))
  invroot_phi ~ normal(0, 100);
  
  for (t in 1:T) {
    for (a_gr in 1:AG) {
      secondary_cases[a_gr, t] ~ neg_binomial_2(exp_reported_cases[a_gr, t], phi);
    }
  }
}

// Compute likelihood
generated quantities {
  matrix [AG, T] log_lik;
  for (a_gr in 1:AG) {
    for (t in 1:T) {
      log_lik[a_gr, t] = neg_binomial_2_lpmf(secondary_cases[a_gr, t] | exp_reported_cases[a_gr, t], phi);
    }
  }  
}











