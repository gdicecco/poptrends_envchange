## BRMS code

abund_mod <- bf(speciestotal ~ dpois(lambda),
                log(lambda) ~ 1 + year + first_yr + (1 + year + first_yr|aou),
                nl = T)

# abund_prior <- c(prior(normal(0, 10000), nlpar = "beta1"),
#                  prior(normal(0, 10000), nlpar = "beta2"))

abund_mod <- brm(formula = abund_mod,
                 family = zero_inflated_poisson("log"),
                 data = bbs_subset)

## Write this in R2jags -- tmaxtrend and abundtrend are param estimates, not measured variables so BRMS style doesn't work

# n <- number of stateroutes for a species

cat("model{
        for(i in 1:n){
        abund[i] ~ b1[is[i]] + b2[is[i]]*tmax[is[i]] + b3[is[i]]*tmin[is[i]] + b4[is[i]]*X[i, deltaED]
        
        }
        
        for(k in 1:nyears) {
        X[k, tmin] ~ b5[is[k]] + tmin*X[k, year]
        X[k, tmax] ~ b6[is[k]] + tmax*X[k, year]
        
        X[k, stoptotal] ~ dpois(lambda[k])
        lambda[k] ~ exp(b7[is[k]] + abund*X[k, year] + b8[is[k]]*X[k, firstyr])
        }
        
        sigma ~ dgamma(1, 1000)
        
        b3 ~ dnorm(0.0, 1000)
        }
        ", fill = T, file = "countsRE.txt"
)

## RStan code model
## Countijk ~ pois(lambdaijk)
## log(lambda ijk) = alpha + piijk*yeark + delta*FYOik



