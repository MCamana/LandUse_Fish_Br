###Baselga beta div###
library(betapart)
#install.packages("betapart")

###Analises Turnover###
getwd()
##Lendo dados já organizados na pasta DataOrganization
LUmont<-read.csv("R/DataOrganization/LUmont_corrected.csv", h=T, row.names = 1)
LU9<-read.csv("R/DataOrganization/LUx9_corrected.csv", h=T, row.names = 1)
LU12<-read.csv("R/DataOrganization/LUx12_corrected.csv", h=T, row.names = 1)
covars<-read.csv("R/DataOrganization/covars_corrected.csv", h=T, row.names = 1)


nrow(LUmont)
nrow(LU9)

###Gerando as Divs betas para cada sub-bacia 9 dentro de cada dataset

#Matrizes de peixes estão numa lista
fishlist<-readRDS("R/DataOrganization/listaPeixes_corrected.RData")
names (fishlist)
#Dissimilaridade media de simpson entre pontos da sub9.
#Dissimilaridade de LU mont entre os mesmo pontos.
#LU total da sub9
listaBetas<-list()

for (j in 1:length(fishlist)) {
  
  dataSet.rod<-names(fishlist)[j]
  
  mat1<-fishlist[[dataSet.rod]]
  mat1<-droplevels(mat1)
  code.sel<-covars[,"Name_fish"] %in% mat1$ID_Point
  length(droplevels(mat1$ID_Point))
  cov.rod<-covars[code.sel,]
  nrow(cov.rod)
  
  #cov.rod$Name_fish == mat1$ID_Point
  code.rod<-cov.rod$X9
  unique(code.rod)
  
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
    sem.zerosPA<-ifelse(sem.zeros>=1,1,0)
    
    #Removendo sites sem peixes
    if (any(rowSums(sem.zerosPA)==0)) {
      print(dataSet.rod)
      print(mat.rod)
      print(rowSums(sem.zerosPA))
    }
    
    sem.zerosPA<-sem.zerosPA[rowSums(sem.zerosPA)>0,]
    
    ##DIVS BETA
    paresBETA<-beta.pair(sem.zerosPA)
    SIMP<-paresBETA$beta.sim
    m.SIMP<-mean(SIMP)
    vet[i]<-m.SIMP
    names(vet)[i]<-code.rod1
    
  }
  
  listaBetas[[j]]<-vet
  names(listaBetas)[j]<-dataSet.rod
}


##listaBetas contem os resultados das divs betas
a<-do.call(c,listaBetas)
a2<-data.frame(a)
colnames(a2)<-"BetaDivs"

###Removing site P0507 because it has 0 species sampled
LU9<-LU9[rownames(LU9)!="P0507",]
LUmont<-LUmont[rownames(LUmont)!="P0507",]

####Fazendo Mat de valores unicos de LU bacia 9#####
length(unique(LU9$X9))
mat.LU9.single<-matrix(NA,length(unique(LU9$X9)),ncol(LU9))
mat.LU9.single<-as.data.frame(mat.LU9.single)
colnames(mat.LU9.single)<-colnames(LU9)

mat.MONT.single<-matrix(NA,length(unique(LUmont$X9_sig)),ncol(LUmont))
mat.MONT.single<-as.data.frame(mat.MONT.single)
colnames(mat.MONT.single)<-colnames(LUmont)


for (i in 1:length(unique(LU9$X9))) {
  sub9.rod<-unique(LU9$X9)[i]
  LU9.rod<-LU9[LU9$X9==sub9.rod,]
  LU9.rod1<-LU9.rod[1,]
  mat.LU9.single[i,]<-LU9.rod1

  
  #Dissimilaridade media do uso do solo a montante dentro de subs9
  LUmont.rod<-LUmont[LUmont$X9_sig==sub9.rod,]
  dists.ano<-apply(LUmont.rod[,-c(1:6)],2,dist)
  
  if(nrow(LUmont.rod)>2){
    m.dists.ano<-colMeans(dists.ano,na.rm = TRUE)
    }else{m.dists.ano<-dists.ano} 
  
  mat.MONT.single[i,1:6]<-LUmont.rod[1,c(1:6)]
  mat.MONT.single[i,7:40]<-m.dists.ano

  }

###################################
# @@@@@@@@@@@ MODELOS @@@@@@@@@@@ #
###################################
mat.MONT.single
mat.LU9.single
a2

nrow(mat.MONT.single)
nrow(mat.LU9.single)
nrow(a2)

cbind(rownames(a2),mat.LU9.single$X9)
###################################
#Nao tem correlacao entre uso do solo na LU9 e diss do uso na LUmont
cor(mat.LU9.single$X11,mat.MONT.single$X11)

ncol(mat.LU9.single)
colnames(mat.LU9.single)
ncol(mat.MONT.single)
colnames(mat.MONT.single)


Coefs.divB.LU9<-matrix(NA,29,3)
Coefs.divB.LU9<-as.data.frame(Coefs.divB.LU9)
colnames(Coefs.divB.LU9)<-c("b_coef","p_value","N")

Coefs.divB.LUMont<-matrix(NA,29,3)
Coefs.divB.LUMont<-as.data.frame(Coefs.divB.LU9)
colnames(Coefs.divB.LUMont)<-c("b_coef","p_value","N")

#Fazer for da coluna 7 a 35 (anos X0 ao X28)
for (i in 7:35) {
LU9.ROD<-mat.LU9.single[,i]
diss.LUmont.ROD<-mat.MONT.single[,i]

datasets<-mat.LU9.single$Data

N<-sum((ifelse(LU9.ROD>0,1,0)),na.rm=TRUE)

mod.ano<-lm(scale(a2$BetaDivs)~scale(LU9.ROD)+scale(diss.LUmont.ROD))
#mod.ano<-lmer(scale(a2$BetaDivs)~scale(LU9.ROD)+scale(diss.LUmont.ROD)+(1|datasets))

teste.lm<-summary(mod.ano)

LU9.resu<-teste.lm$coefficients["scale(LU9.ROD)",][c("Estimate","Pr(>|t|)")]
LUmont.resu<-teste.lm$coefficients["scale(diss.LUmont.ROD)",][c("Estimate","Pr(>|t|)")]

Coefs.divB.LU9[i-6,]<-c(LU9.resu,N)
Coefs.divB.LUMont[i-6,]<-c(LUmont.resu,N)

}

summary (mod.ano)
cor(Coefs.divB.LU9$b_coef,Coefs.divB.LU9$N) 

library(lme4)
library(sjPlot)
library(lmerTest)
mrod2 <- lmer(scale(a2$BetaDivs)~scale(LU9.ROD)+scale(diss.LUmont.ROD) + (1 | datasets))
summary(mrod2)  

##############
#### PLOT ####

par(mfrow=c(1,1))
library(ggplot2)

# Adicionando uma nova coluna para representar os círculos
Coefs.divB.LU9$is_filled <- ifelse(Coefs.divB.LU9$p_value < 0.5, "<0.5", ">0.5")

# Gráfico de dispersão
ggplot(Coefs.divB.LU9, aes(x = seq_along(b_coef), y = b_coef, size = N, 
                            fill = is_filled)) +
  geom_point(shape = 21, alpha = 0.8) +
  scale_fill_manual(values = c("<0.5" = "black", ">0.5" = "white")) +
  labs(x = "Anos antes da amostragem", y = "Efeito do uso da terra sobre a B diversidade") +
  theme_minimal() 




