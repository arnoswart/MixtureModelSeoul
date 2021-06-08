/*
 * Binary Mixture model for seoul rat datasets
 */
  
data {
    // External data to be supplied to the model
    int<lower=1>  n;
    int<lower=1>  n_sample_matrix;
    int<lower=2>  n_components;
    int<lower=0>  n_censored;
    int<lower=1>  n_study;
    
    real                                   logOD_corr[n];
    int<lower=1, upper=n_sample_matrix>    sample_matrix[n];
    int<lower=1, upper=n_study> study[n];
}

parameters {
  // Parameters to be estimated
  ordered[n_components]            mu_raw;
  vector<lower=0>[n_components]    sigma;
  vector<lower=0,upper=1>[n_study] p_study;
  real<lower=0, upper=0.4>         mu_serum;
}

transformed parameters {
  // Transformations of parameters
  vector[n_components] mu[n];
  real                 p[n];
  
  for( j in 1:n ){
    p[j] = p_study[study[j]];}
  
  for( j in 1:n ){
    mu[j] = mu_raw;
    if( sample_matrix[j]==2 ){
      mu[j] = mu[j] + mu_serum;
    }
  }
}

model { 
  // Priors
  mu_raw[1]   ~ normal( 0, 5 );
  mu_raw[2]   ~ normal( 0, 5 );
  mu_serum    ~ normal( 0, 5);
  sigma       ~ normal( 0, 5 );
  p_study     ~ beta( 1,1 );
  
  for( j in 1:n ){
    target +=log_mix(1-p[j], 
                     normal_lpdf(logOD_corr[j] | mu[j,1], sigma[1]),
                     normal_lpdf(logOD_corr[j] | mu[j,2], sigma[2]));}
  
  // Censored values for sample_matrix heart fluid  (=1) and wildrats (=3)
  // Means we can use mu_raw and p_sample_matrix
  target +=n_censored* log_mix(1- p_study[3], 
                               normal_lcdf(-2 | mu_raw[1], sigma[1]),
                               normal_lcdf(-2 | mu_raw[2], sigma[2]));
}
