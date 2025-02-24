// The input data is a vector 'y' of length 'N'.
data {
  int<lower=0> N;
  int<lower=0> M;
  vector[N] y;
  vector[N] x;
  array[N] int<lower=1,upper=M> group;
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  real alpha;
  real beta;
  real<lower=0> sigma_y;
  real<lower=0> sigma_alpha;
  real<lower=0> sigma_beta;
  real lognu_y;
  real lognu_alpha;
  vector[M] alpha_z;
  vector[M] beta_z;
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  vector[N] eta;
  for (i in 1:N) {
    eta[i] = Phi(alpha + sigma_alpha * alpha_z[group[i]] + (sigma_beta * beta_z[group[i]] + beta) * x[i]);
  }
  
  y ~ student_t(exp(lognu_y), eta, sigma_y);
  sigma_y ~ normal(0, 0.5);
  lognu_y ~ normal(0, 2.5);
  
  alpha ~ normal(-2, 2.5);
  sigma_alpha ~ normal(0, 2.5);
  lognu_alpha ~ normal(0, 5);
  alpha_z ~ student_t(exp(lognu_alpha), 0, 1);
  
  beta ~ normal(0, 2.5);
  sigma_beta ~ normal(0, 2.5);
  beta_z ~ std_normal();
}

