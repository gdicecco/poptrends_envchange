// Poisson model of BBS counts at each stateroute, groups by species
// Model abundance trend estimates as linear function of trend in tmax, tmin, change in edge density at stateroute

data {
    int<lower=0> N;                     // number of obs
    int<lower=0> Nsp;                   // number of species
    int<lower=0> Nst;                   // number of species/stateroutes
    int<lower=1, upper=Nsp> sp[N];      // species ids for each obs
    int<lower=1, upper=Nst> cn_id[N];   // species/stateroute ids for each obs
    int<lower=1, upper=Nsp> cn_sp[Nst]; // species ids for each cn_id
    int<lower=1> y[N];                    // speciestotal counts at each year x stroute x spp
    vector[N] year;                     // time
    vector[Nst] tmax;                    // max temp
    vector[Nst] tmin;                    // min temp
    vector[Nst] edge;                   // change in edge density
    vector[Nsp] propFor;
    vector[Nsp] temp_range;
    vector[Nsp] Brange_Area_km2;
}

parameters {
    vector[Nst] beta;             // slope for each species/stateroute
    vector[Nst] alpha;             // intercept for each species/stateroute
    vector[Nsp] alpha2;           // intercept for each species in clim model
    vector[Nsp] beta2;           // species slope for tmax effect
    vector[Nsp] beta3;           // species slope for tmin effect
    vector[Nsp] beta4;           // species slope for edge dens effect
    vector[Nsp] kappa;
    vector[Nsp] kappa2;
    vector[Nsp] kappa3;
    vector[Nsp] kappa4;
    vector[Nsp] kappa5;
    vector[Nsp] kappa6;
    vector[Nsp] kappa7;
    vector[Nsp] kappa8;
    vector[Nsp] kappa9;
    vector[Nsp] kappa10;
    vector[Nsp] kappa11;
    vector[Nsp] kappa12;
    real<lower=0> sigma1;       // std dev of beta
    real<lower=0> sigma2;       // std dev of beta2
    real<lower=0> sigma3;       // std dev of beta3
    real<lower=0> sigma4;       // std dev of beta4
}

transformed parameters{
  vector[N] lambda;  
  vector[N] mu;
  vector[Nst] mu2;
  vector[Nsp] mu3;
  vector[Nsp] mu4;
  vector[Nsp] mu5;
  
  //Abundance indices
    lambda = alpha[cn_id] + beta[cn_id] .* year;
    mu = exp(lambda);
    
    //Predict abund indices based on env change
    mu2 = alpha2[cn_sp] + beta2[cn_sp] .* tmax + beta3[cn_sp] .* tmin + beta4[cn_sp] .* edge;
    
    //Predict env change effects based on traits
    mu3 = kappa + kappa2 .* propFor + kappa3 .* temp_range + kappa4 .* Brange_Area_km2;
    mu4 = kappa5 + kappa6 .* propFor + kappa7 .* temp_range + kappa8 .* Brange_Area_km2;
    mu5 = kappa9 + kappa10 .* propFor + kappa11 .* temp_range + kappa12 .* Brange_Area_km2;
    
}

model {
    y ~ poisson(mu);
    beta ~ normal(mu2, sigma1);
    beta2 ~ normal(mu3, sigma2);
    beta3 ~ normal(mu4, sigma3);
    beta4 ~ normal(mu5, sigma4);
    sigma1 ~ normal(0, 5); // stdev beta
    sigma2 ~ normal(0, 5); // stdev beta2
    sigma3 ~ normal(0, 5); // stdev beta3
    sigma4 ~ normal(0, 5); // stdev beta4
}


