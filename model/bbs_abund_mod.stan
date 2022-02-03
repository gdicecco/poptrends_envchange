// Poisson model of BBS counts at each stateroute, groups by species

data {
  int<lower=1> J;            // number of spp
  int<lower=1> I;           // Number of stateroutes x species
  int<lower=1> N;           // total number of observations
  int<lower=1> Y[N];         // Count outcome: speciestotal   
  vector[N] year;     // Year
  int<lower=1, upper=I> index_group[N]; // Spp x stateroute ID
}

parameters {
    vector[I] beta;             // estimate slope
    vector[I] alpha;             // and intercept for each species/stateroute
}

transformed parameters{
  vector[N] lambda;  
  vector<lower=0>[N] mu;
  
    lambda = alpha[index_group] + beta[index_group] .* year;
    mu = exp(lambda);
}

model {
    Y ~ poisson(mu);
}


