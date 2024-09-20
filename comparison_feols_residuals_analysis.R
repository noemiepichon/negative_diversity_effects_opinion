## Comparison feols and residuals analysis
## NutNet dataset Dee et al. 2023
## PaNDiv datatset 2017-2023
## Eric Allan, Caterina Penone, Bernhard Schmid, Oscar Godoy, Noémie A. Pichon 


library(lmerTest)
library(MuMIn)
library(fixest)
library(ggplot2)
library(cowplot)
library(miceadds)
library(ggpubr)

rm(list=ls())

# setwd("")
bio_nutnet = read.csv("NutNetControlPlotData_derived.csv", header = T)
bio_pandiv = read.table("PaNDiv_biomass_2017to2023.txt", header = T)



#### NutNet data Dee et al 2023 ####


# remove NAs from the dataset
bio_nutnet = bio_nutnet[!is.na(bio_nutnet$live_mass),]
bio_nutnet = bio_nutnet[!is.na(bio_nutnet$rich),]


## Simple lmer model, positive effect of richness

MixedModSimple_Rich_NutNet <- lmer(log(live_mass) ~ log(rich) + (year) +
                        (1|newplotid) + (1|site_code), bio_nutnet, REML = F)
summary(MixedModSimple_Rich_NutNet)
r.squaredGLMM(MixedModSimple_Rich_NutNet)


## Lmer model, positive effect of richness

MixedMod_Rich_NutNet <- lmer(log(live_mass) ~ log(rich) + as.factor(country) + as.factor(habitat) + as.factor(year) + 
                        elevation + managed + burned + grazed + anthropogenic + 
                        TEMP_VAR_v2 + MIN_TEMP_v2 + MAX_TEMP_v2 + TEMP_WET_Q_v2 + TEMP_DRY_Q_v2 + TEMP_WARM_Q_v2 + TEMP_COLD_Q_v2 + 
                        pct_C + pct_N + ppm_P + ppm_K + ppm_Na + ppm_Mg + ppm_S + ppm_Na + ppm_Zn + ppm_Mn + ppm_Fe + ppm_Cu + ppm_B + 
                        pH + PercentSand + PercentSilt + PercentClay + 
                        (1|newplotid) + (1|site_code) , bio_nutnet, REML = F)
summary(MixedMod_Rich_NutNet)
r.squaredGLMM(MixedMod_Rich_NutNet)


## Feols model, negative effect of richness

bio_nutnet$site.by.yeardummy = paste(bio_nutnet$site_code, bio_nutnet$year, sep="_")

MainMod_Rich_NutNet <- feols(log(live_mass) ~ log(rich)  | newplotid + site.by.yeardummy, bio_nutnet) 
summary(MainMod_Rich_NutNet)


## Residuals analysis, same mean effect than feols, almost same standard error when clustering

bio_nutnet$plot_site = paste(bio_nutnet$plot, bio_nutnet$site, sep="_")
bio_nutnet$year_site = paste(bio_nutnet$year, bio_nutnet$site, sep="_")
biomodres = lm(log(live_mass) ~ plot_site + year_site, bio_nutnet)
bio_nutnet$bio_res = biomodres$residuals
richmodres = lm(log(rich) ~ plot_site + year_site, bio_nutnet)
bio_nutnet$rich_res = richmodres$residuals

ResMod_Rich_NutNet <- lm(bio_res ~ rich_res, bio_nutnet)

summary(ResMod_Rich_NutNet)
round(summary(ResMod_Rich_NutNet)$coefficients, digits = 10)

# try out clustered standard errors, same results with sandwich package
ResMod_Rich_NutNet_cl <- lm.cluster(bio_res ~ rich_res, 
                         cluster = 'newplotid',
                         # cluster = "site.by.yeardummy",
                         data = bio_nutnet)
round(summary(ResMod_Rich_NutNet_cl)[4], digits = 10)
#



#### PaNDiv biomass 2017 to 2023 ####


## Lmer model, positive effect of realised diversity

MixedMod_Rich_PaNDiv <- lmer((total.biomass) ~ (year) + log(realised.diversity) + nitrogen + fungicide + functional.composition +
                        (1|block) + (1|plot), bio_pandiv, REML = F) 
plot(MixedMod_Rich_PaNDiv)
summary(MixedMod_Rich_PaNDiv)
r.squaredGLMM(MixedMod_Rich_PaNDiv)


## Feols model, negative effect of realised diversity

bio_pandiv$block_plot = paste(bio_pandiv$block, bio_pandiv$plot, sep="_") 
bio_pandiv$block_year = paste(bio_pandiv$block, bio_pandiv$year, sep="_")

MainMod_Rich_PaNDiv <- feols((total.biomass) ~ log(realised.diversity) | block_plot + block_year, bio_pandiv) 

summary(MainMod_Rich_PaNDiv)
round(summary(MainMod_Rich_PaNDiv)$coefficients, digits = 10)


## Residuals analysis, same mean effect than feols, almost same standard error when clustering

biomodres = lm(total.biomass ~ block_plot + block_year, bio_pandiv)
bio_pandiv$bio_res = biomodres$residuals
richmodres = lm(log(realised.diversity) ~ block_plot + block_year, bio_pandiv)
bio_pandiv$rich_res = richmodres$residuals

ResMod_Rich_PaNDiv <- lm(bio_res ~ rich_res, bio_pandiv)

summary(ResMod_Rich_PaNDiv)
round(summary(ResMod_Rich_PaNDiv)$coefficients, digits = 10)

# try out clustered standard errors, same results with sandwich package
ResMod_Rich_PaNDiv_cl <- lm.cluster(bio_res ~ rich_res, 
                         cluster = 'block_plot',
                         data = bio_pandiv)
round(summary(ResMod_Rich_PaNDiv_cl)[4], digits = 10)
#




#### Plots PaNDiv ####

plot_div_bio = ggplot(bio_pandiv, aes(x = log(realised.diversity), y = total.biomass, color = as.factor(year)))+
  geom_point(size = 2, alpha = 0.3)+
  geom_smooth(method = "lm", se = F)+
  theme_cowplot()+
  scale_x_continuous(breaks = c(log(1), log(4), log(8), log(20)), labels = c(1,4,8,20))+
  theme(legend.title = element_blank())+
  labs(x = "Realised diversity", y = "Biomass (g)")+
  scale_colour_viridis_d()
plot_div_bio

plot_res = ggplot(bio_pandiv, aes(x = rich_res, y = bio_res))+#, color = as.factor(year)))+
  # geom_point(size = 2, alpha = 0.3)+
  geom_point(aes(x = rich_res, y = bio_res, color = as.factor(year)), size = 2, alpha = 0.3)+
  geom_smooth(method = "lm", se = F, color = "salmon")+
  theme_cowplot()+
  # scale_x_continuous(breaks = c(log(1), log(4), log(8), log(20)), labels = c(1,4,8,20))+
  theme(legend.title = element_blank())+
  labs(x = "Diversity residuals", y = "Biomass residuals")+
  scale_colour_viridis_d()
plot_res

ggarrange(plot_div_bio, plot_res,
          labels = c("a)", "b)"),
          ncol = 2, nrow = 1)
#