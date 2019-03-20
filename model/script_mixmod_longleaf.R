#### Run mixed effects models on longleaf

library(nlme)

clim_hab_poptrend_mixedmod <- read.csv("clim_hab_poptrend_mixedmod.csv", stringsAsFactors = F)

randomslope_add <- lme(abundTrend ~ tmin + tmax + ppt + deltaED + deltaProp + Wintering_Slimit_general + Area_sensitivity, 
                       random = (~SPEC|SPEC), 
                       data = clim_hab_poptrend_mixedmod, method = "ML")
randomslope_inter <- lme(abundTrend ~ ppt + deltaProp + Wintering_Slimit_general + Area_sensitivity + tmin*deltaED + tmax*deltaED, 
                         random = (~SPEC|SPEC), 
                         data = clim_hab_poptrend_mixedmod, method = "ML")

save.image()
