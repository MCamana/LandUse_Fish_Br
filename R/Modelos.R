##################################################################################
###This is the script to generate the effect of land use on two spatial scales,### 
###in addition to the history of land use on fish richness in Brazil##############
##################################################################################

##Required packages##
library(lme4)
library(lmerTest)
library(sjPlot) 
library(car)
library(MuMIn)
library(ggplot2)
library(effects)
library(visreg)
library(dplyr)
library(cowplot)
library(ggpubr)

###Data organization
getwd()
LUmont<-read.csv("R/DataOrganization/LUmont_corrected.csv", h=T, row.names = 1)
LU9<-read.table("R/DataOrganization/LUx9_corrected.txt", h=T, row.names = 1)
covars<-read.csv("R/DataOrganization/covars_corrected.csv", h=T, row.names = 1)
Time9<-read.table("R/DataOrganization/duration_all.txt", h=T)


all(LUmont$Name_fish==LU9$Name_fish)
all(Time9$Name_fish==LUmont$Name_fish)

predi.CLU<-LU9[,c("X9","Name_fish","Year_sampling","X0")]
colnames(predi.CLU)[c(1,4)]<-c("sub9","LUyear0.sub9")
head(predi.CLU)
predi.CLU$LUyear0.mont<-LUmont$X0

preditoras<-data.frame(predi.CLU,Time9[c("Duration_70","Duration_50","Duration_30")])

preditoras$NanosLU<-preditoras$Year_sampling-1985

cor(preditoras$Year_sampling,preditoras$NanosLU)

#List of fish occurrence
fishlist<-readRDS("R/DataOrganization/listaPeixes_corrected.RData")

#Calculating fish richness for each matrix
listarich<-lapply(fishlist,function(x){
  matsCommu<-x[,-c(1:7)]
  Rich<-rowSums(ifelse(matsCommu>0,1,0))
  cbind(x[,c(1:5)],Rich)
})

Fish.Rich<-do.call(rbind,listarich)
rownames(Fish.Rich)<-Fish.Rich$ID_Point
rownames(Fish.Rich) == rownames(LUmont)

all(rownames(Fish.Rich) == (preditoras$Name_fish))

###Joining environmental data with fish communities###
dados.mod<-data.frame(Fish.Rich,preditoras)

colnames(dados.mod)

#Scaling and centering data

dados.modSC<-scale(dados.mod[,c("Rich","LUyear0.mont","LUyear0.sub9","Duration_30",
                                "Duration_50","Duration_70")])

dados.modSC<-data.frame(as.factor(dados.mod[,"sub9"]),dados.modSC)
colnames(dados.modSC)[1]<-"sub9"
str(dados.modSC)

###Competing models### 

mod1.lmer<-lmer(Rich ~ LUyear0.mont * LUyear0.sub9 + (1 | sub9), data=dados.modSC)
mod2.lmer<-lmer(Rich ~ LUyear0.mont * Duration_30 + (1 | sub9), data=dados.modSC)
mod3.lmer<-lmer(Rich ~ LUyear0.mont * Duration_50 + (1 | sub9), data=dados.modSC)
mod4.lmer<-lmer(Rich ~ LUyear0.mont * Duration_70 + (1 | sub9), data=dados.modSC)
modNULL.lmer<-lmer(Rich ~ 1 + (1 | sub9), data=dados.modSC)

#Comparing models from two methods

anova(mod1.lmer,mod2.lmer,mod3.lmer,mod4.lmer,modNULL.lmer)

options(na.action = "na.fail")
resumodes<-model.sel(mod1.lmer,mod2.lmer,mod3.lmer,mod4.lmer,modNULL.lmer)
resumodes<-as.data.frame(resumodes)

#Individual effect for each model

tab_model(mod1.lmer, string.std="std.Beta")
tab_model(mod2.lmer, string.std="std.Beta")
tab_model(mod3.lmer, string.std="std.Beta")
tab_model(mod4.lmer, string.std="std.Beta")
tab_model(modNULL.lmer, string.std="std.Beta")

#Residual plots

plot(mod1.lmer)
plot(mod2.lmer)
plot(mod3.lmer)
plot(mod4.lmer)
plot(modNULL.lmer)

###Making the result plots for the significative models###

write.csv(dados.mod, "DataOrganization/dados_plot_bruto.csv")

dados.plot <- read.table("DataOrganization/dados_plot.csv", h = T, sep = ";")
colnames(dados.plot)

#Mod3, 1º plot
plot_mod3_dur <- ggplot(data = dados.plot, aes(x = LUyear0.mont, y = Rich)) +
  geom_point (aes(col = Duration_50_groups), size = 4, alpha = 0.6) +
  geom_smooth(method = "lm", aes(col = Duration_50_groups), alpha = 0.3, size = 2) +
  labs(x = "Current upland land use (%)", y = "Species richness") +
  theme_classic() +
  theme (axis.title.x = element_text(size = 18)) +
  theme (axis.text.y = element_text(size = 16)) +
  theme (axis.text.x = element_text(size = 16), axis.ticks.x = element_blank()) +
  theme (axis.title.y = element_text(size = 18)) +
  theme(legend.title = "element_text"(size = 14), legend.text = "element_text"(size = 12)) +
  scale_color_manual(name = "50% LU threshold", values=c("#018571", "#dfc27d"), 
              breaks=c("recent", "old"), labels=c("old" = "≥ 20 years", "recent" = "< 20 years"))

  plot_mod3_dur

#Mod3, 2º plot
plot_mod3_mont <- ggplot(data = dados.plot, aes(x = Duration_50, y = Rich)) +
  geom_point (aes(col = mont_groups), size = 4, alpha = 0.6) +
  geom_smooth(method = "lm", aes(col = mont_groups), alpha = 0.3, size = 2) +
  labs(x = "50% land use threshold (years)", y = "Species richness") +
  theme_classic() +
  theme (axis.title.x = element_text(size = 18)) +
  theme (axis.text.y = element_text(size = 16)) +
  theme (axis.text.x = element_text(size = 16), axis.ticks.x = element_blank()) +
  theme (axis.title.y = element_text(size = 18)) +
  theme(legend.title = "element_text"(size = 14), legend.text = "element_text"(size = 12)) +
  scale_color_manual(name = "Upland land use", values=c("#018571", "#dfc27d"),
              breaks=c("baixo", "alto"), labels=c("alto" = "High (≥  50%)", "baixo" = "Low (< 50%)"))

plot_mod3_mont

#Mod1, 1º plot
plot_mod1_sub9 <- ggplot(data = dados.plot, aes(x = LUyear0.mont, y = Rich)) +
  geom_point (aes(col = sub9_groups), size = 4, alpha = 0.6) +
  geom_smooth(method = "lm", aes(col = sub9_groups), alpha = 0.3, size = 2) +
  labs(x = "Current upland land use (%)", y = "Species richness") +
  theme_classic() +
  theme (axis.title.x = element_text(size = 18)) +
  theme (axis.text.y = element_text(size = 16)) +
  theme (axis.text.x = element_text(size = 16), axis.ticks.x = element_blank()) +
  theme (axis.title.y = element_text(size = 18)) +
  theme(legend.title = "element_text"(size = 14), legend.text = "element_text"(size = 12)) +
  scale_color_manual(name = "Regional land use", values=c("#018571", "#dfc27d"),
              breaks=c("baixo", "alto"), labels=c("alto" = "High (≥  50%)", "baixo" = "Low (< 50%)"))
  
plot_mod1_sub9

#Mod1, 2º plot
plot_mod1_mont <- ggplot(data = dados.plot, aes(x = LUyear0.sub9, y = Rich)) +
  geom_point (aes(col = mont_groups), size = 4, alpha = 0.6) +
  geom_smooth(method = "lm", aes(col = mont_groups), alpha = 0.3, size = 2) +
  labs(x = "Current regional land use (%)", y = "Species richness") +
  theme_classic() +
  theme (axis.title.x = element_text(size = 18)) +
  theme (axis.text.y = element_text(size = 16)) +
  theme (axis.text.x = element_text(size = 16), axis.ticks.x = element_blank()) +
  theme (axis.title.y = element_text(size = 18)) +
  theme(legend.title = "element_text"(size = 14), legend.text = "element_text"(size = 12)) +
  scale_color_manual(name = "Upland land use", values=c("#018571", "#dfc27d"), 
              breaks=c("baixo", "alto"), labels=c("alto" = "High (≥  50%)", "baixo" = "Low (< 50%)"))

  plot_mod1_mont

#Joining the plots in a single panel

rich_models_plots <- ggarrange(plot_mod3_dur, plot_mod3_mont, plot_mod1_sub9,
                               plot_mod1_mont, labels = c("A", "B", "C", "D"), label.x = 0.7, 
                               align = "hv", legend = "right")
rich_models_plots

#Saving
save_plot("output/rich_models_plots.png", rich_models_plots, dpi= 600, base_aspect_ratio = 3:1)

####Now, we will rerun the models using duration values for just 20 years of historical land use####
###Data organization
Time9_20y<-read.table("DataOrganization/duration_20y.txt", h=T)
nomes_time20 <- Time9_20y$Name_fish
predi.CLU_time20 <- subset(predi.CLU, Name_fish %in% nomes_time20)

all(predi.CLU_time20$Name_fish==Time9_20y$Name_fish)

colnames(Fish.Rich)
colnames(Fish.Rich)[colnames(Fish.Rich) == "ID_Point"] <- "Name_fish"

Fish.Rich_time20 <- subset(Fish.Rich, Name_fish %in% nomes_time20)
all(Fish.Rich_time20$Name_fish==Time9_20y$Name_fish)

preditoras_time20<-data.frame(predi.CLU_time20,Time9_20y[c("Duration_70_20","Duration_50_20","Duration_30_20")])

dados.mod_time20<-data.frame(Fish.Rich_time20,preditoras_time20)

colnames(dados.mod_time20)

#Scaling and centering data
dados.modSC_time20<-scale(dados.mod_time20[,c("Rich","LUyear0.mont","LUyear0.sub9","Duration_30_20",
                                              "Duration_50_20","Duration_70_20")])
dados.modSC_time20<-data.frame(as.factor(dados.mod_time20[,"sub9"]),dados.modSC_time20)
colnames(dados.modSC_time20)[1]<-"sub9"
str(dados.modSC_time20)

###Competing models###

mod1.lmer_time20<-lmer(Rich ~ LUyear0.mont * LUyear0.sub9 + (1 | sub9), data=dados.modSC_time20)
mod2.lmer_time20<-lmer(Rich ~ LUyear0.mont * Duration_30_20 + (1 | sub9), data=dados.modSC_time20)
mod3.lmer_time20<-lmer(Rich ~ LUyear0.mont * Duration_50_20 + (1 | sub9), data=dados.modSC_time20)
mod4.lmer_time20<-lmer(Rich ~ LUyear0.mont * Duration_70_20 + (1 | sub9), data=dados.modSC_time20)
modNULL.lmer_time20<-lmer(Rich ~ 1 + (1 | sub9), data=dados.modSC_time20)

#Comparing models from two methods

anova(mod1.lmer,mod2.lmer,mod3.lmer,mod4.lmer,modNULL.lmer)#

options(na.action = "na.fail")
resumodes_time20<-model.sel(mod1.lmer_time20,mod2.lmer_time20,mod3.lmer_time20,mod4.lmer_time20,modNULL.lmer_time20)
resumodes_time20<-as.data.frame(resumodes_time20)
resumodes_time20

#Individual effect for each model

tab_model(mod1.lmer_time20, string.std="std.Beta")
tab_model(mod2.lmer_time20, string.std="std.Beta")
tab_model(mod3.lmer_time20, string.std="std.Beta")
tab_model(mod4.lmer_time20, string.std="std.Beta")
tab_model(modNULL.lmer_time20, string.std="std.Beta")








