##################################################################################
###This is the script to generate the effect of land use on two spatial scales,### 
###in addition to the history of land use on fish beta diversity in Brazil########
##################################################################################

##Required packages##
library(lme4)
library(lmerTest)
library(betapart)
library(vegan)
library(car)
library(MuMIn)
library(sjPlot) 
library(car)
library(ggplot2)
library(cowplot)
library(ggpubr)

###Data organization
gc()
getwd()
list.files()
LUmont<-read.csv("DataOrganization/LUmont_corrected.csv", h=T, row.names = 1)
LU9<-read.table("DataOrganization/LUx9_corrected.txt", h=T, row.names = 1)
covars<-read.csv("DataOrganization/covars_corrected.csv", h=T, row.names = 1)
Time9<-read.table("DataOrganization/duration_all.txt", h=T)

#Jean - adicionei aqui
regional_up<-read.csv("DataOrganization/regional_upland.csv", sep = ";") #lendo o arquivo com a condição de uso das bacias


all(LUmont$Name_fish==LU9$Name_fish)
all(Time9$Name_fish==LUmont$Name_fish)

predi.CLU<-LU9[,c("X9","Name_fish","Year_sampling","X0")]
colnames(predi.CLU)[c(1,4)]<-c("sub9","LUyear0.sub9")
head(predi.CLU)
predi.CLU$LUyear0.mont<-LUmont$X0

preditoras<-data.frame(predi.CLU,Time9[c("Duration_70","Duration_50","Duration_30")])

preditoras$NanosLU<-preditoras$Year_sampling-1985

#List of fish occurrence
fishlist<-readRDS("DataOrganization/listaPeixes_corrected.RData")

# ###Calculando riqueza para cada matriz de peixes e montando uma lista nova
# listarich<-lapply(fishlist,function(x){
#   matsCommu<-x[,-c(1:7)]
#   Rich<-rowSums(ifelse(matsCommu>0,1,0))
#   cbind(x[,c(1:5)],Rich)
# })
# 
# Fish.Rich<-do.call(rbind,listarich)
# rownames(Fish.Rich)<-Fish.Rich$ID_Point
# rownames(Fish.Rich) == rownames(LUmont)
# 
# all(rownames(Fish.Rich) == (preditoras$Name_fish))

#Matrizes de peixes estão numa lista

fishlist<-readRDS("DataOrganization/listaPeixes_corrected.RData")
names (fishlist)

#Calculating fish beta diversity for each matrix 
listaBetas<-list()

for (j in 1:length(fishlist)) {
  
  dataSet.rod<-names(fishlist)[j]
  
  mat1<-fishlist[[dataSet.rod]]
  mat1<-droplevels(mat1)
  code.sel<-covars[,"Name_fish"] %in% mat1$ID_Point
  length(droplevels(mat1$ID_Point))
  cov.rod<-covars[code.sel,]
  nrow(cov.rod)
  
  code.rod<-cov.rod$X9
  unique(code.rod)
  #i<-1
  vet<-NA
  for (i in 1:length(unique(code.rod))) {
    code.rod1<-unique(code.rod)[i]
    sites.rod<-cov.rod[cov.rod[,"X9"]==code.rod1,]
    sites.rod<-sites.rod$Name_fish
    
    mat.rod<-mat1[mat1[,"ID_Point"] %in% sites.rod,]
    mat.rod
    
    ##Div. beta
    so.number<-mat.rod[,-c(1:7)]
    sem.zeros<-so.number[,colSums(so.number)>0]
    sem.zerosPA<-ifelse(sem.zeros>0,1,0)
    
    #Removendo sites sem peixes
    if (any(rowSums(sem.zerosPA)==0)) {
      print(dataSet.rod)
      print(mat.rod)
      print(rowSums(sem.zerosPA))
    }
    
    sem.zerosPA<-sem.zerosPA[rowSums(sem.zerosPA)>0,]
    
    ##DIVS BETA
    paresBETA<-beta.pair(sem.zerosPA)
    ###SIMP<-paresBETA$beta.sim
    SIMP<-paresBETA$beta.sor
    m.SIMP<-mean(SIMP)
    vet[i]<-m.SIMP
    names(vet)[i]<-code.rod1
    
  }
  
  vet1<-data.frame(vet)
  vet1$sub9<-rownames(vet1)
  vet1$dataset<-dataSet.rod
  
  listaBetas[[j]]<-vet1
  names(listaBetas)[j]<-dataSet.rod
}


#List with the beta diversity results
a<-do.call(rbind,listaBetas)
Turnover<-data.frame(a)
colnames(Turnover)[1]<-"BetaDivs"

######Dissimilarity between catchement land use inside sub 9

lista.Predi<-split(preditoras, preditoras$sub9)

predi.sub9<-lapply(lista.Predi, function(x){
  vec.rod<-x[1,c("sub9","Year_sampling","LUyear0.sub9","Duration_70","Duration_50","Duration_30")]
  dissLUmont<-mean(vegdist(x[,"LUyear0.mont"], method="euclidean"))
  vec.rod1<-data.frame(vec.rod,dissLUmont)
  })

predi.sub9<-do.call(rbind,predi.sub9)

nrow(predi.sub9)
nrow(Turnover)

#Ordering the lines

Turnover<-Turnover[order(as.numeric(Turnover$sub9)),]
Turnover$sub9 == predi.sub9$sub9

###Joining environmental data with fish communities###

###########################

#Jean - adicionei aqui

#merge data by the id of sub9 in both dataframes
rm(dados.mod)
dados.mod <- merge(Turnover, predi.sub9, by = "sub9", all = TRUE)
dados.mod <- merge(dados.mod, regional_up, by = "sub9", all = TRUE)

###########################

#Scaling and centering data

colnames(dados.mod)

dados.modSC<-scale(dados.mod[,c("BetaDivs","dissLUmont","LUyear0.sub9","Duration_30",
                                "Duration_50","Duration_70")])

dados.modSC<-data.frame(as.factor(dados.mod[,"dataset"]),dados.modSC)
colnames(dados.modSC)[1]<-"dataset"
str(dados.modSC)
colnames(dados.modSC)

####Competing models### 

mod1.lmer<-lmer(BetaDivs ~ LUyear0.sub9 + (1 | dataset), data=dados.modSC)
mod2.lmer<-lmer(BetaDivs ~ Duration_30 + (1 | dataset), data=dados.modSC)
mod3.lmer<-lmer(BetaDivs ~ Duration_50 + (1 | dataset), data=dados.modSC)
mod4.lmer<-lmer(BetaDivs ~ Duration_70 + (1 | dataset), data=dados.modSC)
modNULL.lmer<-lmer(BetaDivs ~ 1 + (1 | dataset), data=dados.modSC)

#Comparing models from two methods

anova(mod1.lmer,mod2.lmer,mod3.lmer,mod4.lmer,modNULL.lmer)

options(na.action = "na.fail")
resumodes<-model.sel(mod1.lmer,mod2.lmer,mod3.lmer,mod4.lmer,modNULL.lmer)
resumodes.beta<-as.data.frame(resumodes)
resumodes.beta

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

###########################

#Jean - adicionei aqui

colnames(dados.mod)

#Mod1
plot_beta_sub9 <- ggplot(data = dados.mod, aes(x = LUyear0.sub9, y = BetaDivs)) +
  geom_point (aes(col = status_all), size = 4, alpha = 0.8) +
  geom_smooth(method = "lm", alpha = 0.3, size = 2, color = "black") +
  labs(x = "Current regional land use (%)", y = "Beta diversity") +
  theme_classic() +
  theme (axis.title.x = element_text(size = 18)) +
  theme (axis.text.y = element_text(size = 16)) +
  theme (axis.text.x = element_text(size = 16), axis.ticks.x = element_blank()) +
  theme (axis.title.y = element_text(size = 18)) +
  theme(legend.title = "element_text"(size = 14), legend.text = "element_text"(size = 12)) +
  scale_color_manual(name = "50% LU threshold",
                     values=c("#01665e","#5ab4ac", "#d8b365", "#8c510a"), breaks=c("none", "only_upland", 
                    "only_regional","both"), labels=c("none"="None", "only_upland"="Upland", 
                                                      "only_regional"="Regional","both"="Both"))
plot_beta_sub9

#Mod3
plot_beta_dur50 <- ggplot(data = dados.mod, aes(x = Duration_50, y = BetaDivs)) +
  geom_point (aes(col = status_all), size = 4, alpha = 0.8) +
  geom_smooth(method = "lm", alpha = 0.3, size = 2, color = "black") +
  labs(x = " 50% land use threshold (years)", y = "Beta diversity") +
  theme_classic() +
  theme (axis.title.x = element_text(size = 18)) +
  theme (axis.text.y = element_text(size = 16)) +
  theme (axis.text.x = element_text(size = 16), axis.ticks.x = element_blank()) +
  theme (axis.title.y = element_text(size = 18)) +
  theme(legend.title = "element_text"(size = 14), legend.text = "element_text"(size = 12)) +
  scale_color_manual(name = "50% LU threshold",
                     values=c("#01665e","#5ab4ac", "#d8b365", "#8c510a"), breaks=c("none", 
                    "only_upland","only_regional","both"), labels=c("none"="None", "only_upland"="Upland", 
                                                                                                                     "only_regional"="Regional","both"="Both"))
plot_beta_dur50

#Joining the plots in a single panel

beta_models_plots <- ggarrange(plot_beta_sub9, plot_beta_dur50, 
                               labels = c("A", "B"), label.x = 0.9, align = "hv")
beta_models_plots


#Saving

save_plot("output/beta_models_plots.png", beta_models_plots, dpi= 600, base_aspect_ratio = 3:2)

####Now, we will rerun the models using duration values for just 20 years of historical land use####
###Data organization

Time9_20y<-read.table("DataOrganization/duration_20y.txt", h=T)
nomes_time20 <- Time9_20y$X9
predi.CLU_time20 <- subset(predi.CLU, sub9 %in% nomes_time20)
all(predi.CLU_time20$Name_fish==Time9_20y$X9)

colnames(Turnover)
colnames(Turnover)[colnames(Turnover) == "sub9"] <- "sub9"
colnames(Turnover)
 
Fish.beta_time20 <- subset(Turnover, sub9 %in% nomes_time20)

all(Fish.beta_time20$Name_fish==Time9_20y$Name_fish)
 
preditoras_time20<-data.frame(predi.CLU_time20,Time9_20y[c("Duration_70_20","Duration_50_20","Duration_30_20")])
colnames(preditoras_time20) 
dados.mod_time20<-data.frame(Fish.beta_time20,preditoras_time20)
# 
# colnames(dados.mod_time20)
# 
# dados.modSC_time20<-scale(dados.mod_time20[,c("Rich","LUyear0.mont","LUyear0.sub9","Duration_30_20",
#                                 "Duration_50_20","Duration_70_20")])
# 
# dados.modSC_time20<-data.frame(as.factor(dados.mod_time20[,"sub9"]),dados.modSC_time20)
# colnames(dados.modSC_time20)[1]<-"sub9"
# str(dados.modSC_time20)
# 
# mod1.lmer_time20<-lmer(Rich ~ LUyear0.mont * LUyear0.sub9 + (1 | sub9), data=dados.modSC_time20)
# mod2.lmer_time20<-lmer(Rich ~ LUyear0.mont * Duration_30_20 + (1 | sub9), data=dados.modSC_time20)
# mod3.lmer_time20<-lmer(Rich ~ LUyear0.mont * Duration_50_20 + (1 | sub9), data=dados.modSC_time20)
# mod4.lmer_time20<-lmer(Rich ~ LUyear0.mont * Duration_70_20 + (1 | sub9), data=dados.modSC_time20)
# modNULL.lmer_time20<-lmer(Rich ~ 1 + (1 | sub9), data=dados.modSC_time20)
# 
# options(na.action = "na.fail")
# resumodes_time20<-model.sel(mod1.lmer_time20,mod2.lmer_time20,mod3.lmer_time20,mod4.lmer_time20,modNULL.lmer_time20)
# resumodes_time20<-as.data.frame(resumodes_time20)
# write.csv2(resumodes_time20, "ResuMods_time20.csv")
# 
# tab_model(mod1.lmer_time20, string.std="std.Beta")
# tab_model(mod2.lmer_time20, string.std="std.Beta")
# tab_model(mod3.lmer_time20, string.std="std.Beta")
# tab_model(mod4.lmer_time20, string.std="std.Beta")
# tab_model(modNULL.lmer_time20, string.std="std.Beta")

