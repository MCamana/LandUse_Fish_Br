#Como LU dos diferentes anos desde 1985
#em 3 escalas espaciais diferentes estao relacioandos 
#com a riqueza de peixes coletados ao longo do Brasil
getwd()

##Lendo dados já organizados na pasta DataOrganization
# LUmont<-read.csv("LUmont.csv", h=T, row.names = 1)
# LU9<-read.csv("LUx9.csv", h=T, row.names = 1)
# LU12<-read.csv("LUx12.csv", h=T, row.names = 1)
# covars<-read.csv("covars.csv", h=T, row.names = 1)

LUmont<-read.csv("R/DataOrganization/LUmont_corrected.csv", h=T, row.names = 1)
LU9<-read.csv("R/DataOrganization/LUx9_corrected.csv", h=T, row.names = 1)
LU12<-read.csv("R/DataOrganization/LUx12_corrected.csv", h=T, row.names = 1)
covars<-read.csv("R/DataOrganization/covars_corrected.csv", h=T, row.names = 1)

##################
###BOXPLOTS LU####
##################
library(ggplot2)
library(tidyr)
library(dplyr)
library(cowplot)
#install.packages("cowplot")


##REMOVENDO ANOS NAO USADOS NAS BACIAS
LUmont <- LUmont[, !(names(LUmont) %in% c('X29','X30','X31','X32','X33'))]
LU9 <- LU9[, !(names(LU9) %in% c('X29','X30','X31','X32','X33'))]
LU12 <- LU12[, !(names(LU12) %in% c('X29','X30','X31','X32','X33'))]


#########################
##Organização dos dados##
#########################

#Matrizes de peixes estão numa lista
fishlist<-readRDS("DataOrganization/listaPeixes_corrected.RData")

###Calculando riqueza para cada matriz de peixes e montando uma lista nova
listarich<-lapply(fishlist,function(x){
  matsCommu<-x[,-c(1:7)]
  Rich<-rowSums(ifelse(matsCommu>0,1,0))
  cbind(x[,c(1:5)],Rich)
})

Fish.Rich<-do.call(rbind,listarich)
rownames(Fish.Rich)<-Fish.Rich$ID_Point
rownames(Fish.Rich) == rownames(LUmont)
all(rownames(Fish.Rich) == rownames(LUmont))

#Matriz com dados de riqueza
Fish.Rich

######################################################################
#############TESTES INICIAIS DO MODELO USANDO APENAS 1 LU#############
######################################################################
####MODELOS LINEAR MISTO ANINHADO####
####TESTANTO PARA O ANO X0
DFteste<-data.frame(Fish.Rich$Rich, covars, LUmont$X0)
str(DFteste)
colnames(DFteste)[1]<-"Rich"
colnames(DFteste)

#PREDITORAS FIXAS: uso da terra por ano LUmont, 
#Mont_temperaturamean, Mont_precmean, Geologia,topografiaCV, Area, Fresh_ecor  
#PREDITORAS ALEATÓRIAS ANINHADAS: (1 | 12)

library(lme4)

colnames(DFteste)
unique(DFteste$X9)
unique(DFteste$X12)

# ###TopografiaCV tem correlação de 0.65 com LU. Quando tira topografia LU fica significativo.
m0 <- lmer(Rich ~ LUmont.X0 + (1 | X9), data = DFteste,
           REML = FALSE, control = lmerControl(optimizer ="Nelder_Mead"))

m1 <- lmer(Rich ~ LUmont.X0 + (1 + LUmont.X0 | X9), data = DFteste,
           REML = FALSE, control = lmerControl(optimizer ="Nelder_Mead"))

m2 <- lmer(Rich ~ LUmont.X0 + (1 | X9) + (1 | Fresh_ecor), data = DFteste,
           REML = FALSE, control = lmerControl(optimizer ="Nelder_Mead"))

m3 <- lmer(Rich ~ LUmont.X0 + Mont_temperaturamean + Mont_precmean +
             topografiaCV + Area + (1 | X9), data = DFteste,
           REML = FALSE, control = lmerControl(optimizer ="Nelder_Mead"))

m4 <- lmer(Rich ~ LUmont.X0 + Mont_temperaturamean + Mont_precmean +
             topografiaCV + Area + (1 | X9) + (1 | Fresh_ecor), data = DFteste,
           REML = FALSE, control = lmerControl(optimizer ="Nelder_Mead"))

summary(m0)
summary(m1)
summary(m2)
summary(m3)
summary(m4)

anova(m0,m1,m2,m3,m4) #O modelos 4 é melhor

#install.packages("sjPlot")
library(sjPlot) #Para ver tabela de resultados
tab_model(m0, string.std="std.Beta")
tab_model(m1, string.std="std.Beta")
tab_model(m2, string.std="std.Beta")
tab_model(m3, string.std="std.Beta")
tab_model(m4, string.std="std.Beta")

##Verificando se os modelos estão inflados
library(car)
car::vif(m3) #o modelo não é inflado

#install.packages("lmerTest")
library(lmerTest) ### Carregar para pegar p valores

##############################
##### MODELOS DE RIQUEZA #####
##############################
#######FAZER O LOOPING########
##############################

##################
####MONTANTE######
##################
LUmat.mont<-LUmont
LUmat.12<-LU12
LUmat.9<-LU9

RIQmod0.mont<-matrix(NA,(ncol(LUmat.mont)-6),3)
colnames(RIQmod0.mont)<-c("Effect.size","pvalue","N")
rownames(RIQmod0.mont)<-colnames(LUmat.mont)[-c(1:6)]

RIQmod0.sub12<-matrix(NA,(ncol(LUmat.12)-6),3)
colnames(RIQmod0.sub12)<-c("Effect.size","pvalue","N")
rownames(RIQmod0.sub12)<- colnames(LUmat.12)[-c(1:6)]
 
RIQmod0.sub9<-matrix(NA,(ncol(LUmat.9)-6),3)
colnames(RIQmod0.sub9)<-c("Effect.size","pvalue","N")
rownames(RIQmod0.sub9)<- colnames(LUmat.9)[-c(1:6)]
 
listaAnos<-list()

#i<-7

#For comeca em 7 que e onde comecam os anos na matriz
for (i in 7:ncol(LUmat.mont)){
  LUano<-colnames(LUmat.mont[i])
  
  LUmont.rod<-LUmat.mont[,LUano]
  LU12.rod<-LUmat.12[,LUano]
  LU9.rod<-LUmat.9[,LUano]
  
  DFrod.mon<-data.frame(Fish.Rich$Rich, covars, LUmont.rod)
  DFrod.12<-data.frame(Fish.Rich$Rich, covars, LU12.rod)
  DFrod.9<-data.frame(Fish.Rich$Rich, covars, LU9.rod)
  
  #Filtrando-removendo dados de NA
  ptos.SNA.m<-na.exclude(DFrod.mon)
  ptos.SNA.12<-na.exclude(DFrod.12)
  ptos.SNA.9<-na.exclude(DFrod.9)
  
  #Removendo sub-bacias com menos que 5 pontos
  ptos<-table(ptos.SNA.m$X9)[table(ptos.SNA.m$X9)>=5]
  nomesSel<-names(ptos)
  sub.dados<-ptos.SNA.m[ptos.SNA.m$X9 %in% nomesSel,]
  DFrod.mon<-droplevels(sub.dados)

  ptos<-table(ptos.SNA.12$X9)[table(ptos.SNA.12$X9)>=5]
  nomesSel<-names(ptos)
  sub.dados<-ptos.SNA.12[ptos.SNA.12$X9 %in% nomesSel,]
  DFrod.12<-droplevels(sub.dados)

  ptos<-table(ptos.SNA.9$X9)[table(ptos.SNA.9$X9)>=5]
  nomesSel<-names(ptos)
  sub.dados<-ptos.SNA.9[ptos.SNA.9$X9 %in% nomesSel,]
  DFrod.9<-droplevels(sub.dados)


  #sub.dados<-DFrod
  mrod0.mon <- lmer(Fish.Rich.Rich ~ LUmont.rod + (1 | X9), data = DFrod.mon,
                REML = FALSE, control = lmerControl(optimizer ="Nelder_Mead"))
  #sub.dados<-DFrod
  mrod0.LU12 <- lmer(Fish.Rich.Rich ~ LU12.rod + (1 | X9), data = DFrod.12,
                REML = FALSE, control = lmerControl(optimizer ="Nelder_Mead"))
  #sub.dados<-DFrod
  mrod0.LU9 <- lmer(Fish.Rich.Rich ~ LU9.rod + (1 | X9), data = DFrod.9,
                REML = FALSE, control = lmerControl(optimizer ="Nelder_Mead"))
  
  resu.rod0.mon<-summary(mrod0.mon)
  BeP0<-resu.rod0.mon$coefficients["LUmont.rod",][c("Estimate","Pr(>|t|)")]
  
  resu.rod0.12<-summary(mrod0.LU12)
  BeP1<-resu.rod0.12$coefficients["LU12.rod",][c("Estimate","Pr(>|t|)")]
  
  resu.rod0.9<-summary(mrod0.LU9)
  BeP2<-resu.rod0.9$coefficients["LU9.rod",][c("Estimate","Pr(>|t|)")]
  
  RIQmod0.mont[(i-6),] <- c(BeP0,nrow(DFrod.mon))
  RIQmod0.sub12[(i-6),] <- c(BeP1,nrow(DFrod.12))
  RIQmod0.sub9[(i-6),] <- c(BeP2,nrow(DFrod.9))
  
  # RIQmod1.mont[(i-6),] <- c(BeP1,nrow(sub.dados))
  # RIQmod2.mont[(i-6),] <- c(BeP2,nrow(sub.dados))
  # RIQmod3.mont[(i-6),] <- c(BeP3,nrow(sub.dados))
  print(i)
    }

#Preditora
years<-0:28

###MODELO 0
#Objeto símbolos com p significativo e não significativo
simbolsA<-RIQmod0.mont[,"pvalue"]<0.05
preeA<-as.numeric(simbolsA)
preeA[preeA==0]<-21
preeA[preeA==1]<-19


x11()
par(mfrow=c(1,3))
plot(Effect.size~years, cex=RIQmod0.mont[,"N"]/50, data=RIQmod0.mont, 
     pch=preeA,col=1)

#x11()
plot(Effect.size~years, cex=RIQmod0.mont[,"N"]/50, data=RIQmod0.sub12, 
     pch=preeA,col=1)

#x11()
plot(Effect.size~years, cex=RIQmod0.mont[,"N"]/50, data=RIQmod0.sub9, 
     pch=preeA,col=1)


##################
####  LU 12 ######
##################
LUmat<-LU12

RIQmod0.LU12<-matrix(NA,32,3)
colnames(RIQmod0.LU12)<-c("Effect.size","pvalue","N")
rownames(RIQmod0.LU12)<-colnames(LUmat)[2:33]

RIQmod1.LU12<-matrix(NA,32,3)
colnames(RIQmod1.LU12)<-c("Effect.size","pvalue","N")
rownames(RIQmod1.LU12)<- colnames(LUmat)[2:33]

RIQmod2.LU12<-matrix(NA,32,3)
colnames(RIQmod2.LU12)<-c("Effect.size","pvalue","N")
rownames(RIQmod2.LU12)<- colnames(LUmat)[2:33]

RIQmod3.LU12<-matrix(NA,32,3)
colnames(RIQmod3.LU12)<-c("Effect.size","pvalue","N")
rownames(RIQmod3.LU12)<- colnames(LUmat)[2:33]

listaAnos<-list()

#i<-2

#Removendo os dois últimos anos do looping (32 e 33) pois só tem 12 observações
for (i in 2:(ncol(LUmat)-2)){
  LUano<-colnames(LUmat[i])
  LU<-LUmat[,LUano]
  DFrod<-data.frame(Fish.Rich$Rich, covars, LU)
  colnames(DFrod)[1]<-"Rich"
  colnames(DFrod)[12]<-"LUrod"
  colnames(DFrod)
  print(LUano)
  
  #Filtrando-removendo dados de NA
  ptos.SNA<-na.exclude(DFrod)
  #Removendo sub-bacias com menos que 5 pontos
  ptos<-table(ptos.SNA$X9)[table(ptos.SNA$X9)>=5]
  nomesSel<-names(ptos)
  sub.dados<-ptos.SNA[ptos.SNA$X9 %in% nomesSel,]
  sub.dados<-droplevels(sub.dados)
  table(sub.dados$X9)
  
  #sub.dados<-DFrod
  
  mrod0 <- lmer(Rich ~ LUrod + (1 | X9), data = sub.dados,
                REML = FALSE, control = lmerControl(optimizer ="Nelder_Mead"))
  mrod1 <- lmer(Rich ~ LUrod + (LUrod | X9), data = sub.dados,
                REML = FALSE, control = lmerControl(optimizer ="Nelder_Mead"))
  mrod2 <- lmer(Rich ~ LUrod + topografiaCV + Area + (1 | X9), data = sub.dados,
                REML = FALSE, control = lmerControl(optimizer ="Nelder_Mead"))
  mrod3 <- lmer(Rich ~ LUrod + topografiaCV + Area + (LUrod | X9), data = sub.dados,
                REML = FALSE, control = lmerControl(optimizer ="Nelder_Mead"))
  
  resu.rod0<-summary(mrod0)
  BeP0<-resu.rod0$coefficients["LUrod",][c("Estimate","Pr(>|t|)")]
  resu.rod1<-summary(mrod1)
  BeP1<-resu.rod1$coefficients["LUrod",][c("Estimate","Pr(>|t|)")]
  resu.rod2<-summary(mrod2)
  BeP2<-resu.rod2$coefficients["LUrod",][c("Estimate","Pr(>|t|)")]
  resu.rod3<-summary(mrod3)
  BeP3<-resu.rod3$coefficients["LUrod",][c("Estimate","Pr(>|t|)")]
  
  RIQmod0.LU12[(i-1),] <- c(BeP0,nrow(sub.dados))
  RIQmod1.LU12[(i-1),] <- c(BeP1,nrow(sub.dados))
  RIQmod2.LU12[(i-1),] <- c(BeP2,nrow(sub.dados))
  RIQmod3.LU12[(i-1),] <- c(BeP3,nrow(sub.dados))
}

years<-0:31

simbolsA<-RIQmod0.LU12[,"pvalue"]<0.05
preeA<-as.numeric(simbolsA)
preeA[preeA==0]<-21
preeA[preeA==1]<-19

plot(Effect.size~years, cex=RIQmod0.LU12[,"N"]/50, data=RIQmod0.LU12, 
     pch=preeA,col=1)

simbolsB<-RIQmod1.LU12[,"pvalue"]<0.05
preeB<-as.numeric(simbolsB)
preeB[preeB==0]<-21
preeB[preeB==1]<-19

plot(Effect.size~years, cex=RIQmod1.LU12[,"N"]/50, data=RIQmod1.LU12, 
     pch = preeB)

simbolsC<-RIQmod2.LU12[,"pvalue"]<0.05
preeC<-as.numeric(simbolsC)
preeC[preeC==0]<-21
preeC[preeC==1]<-19

plot(Effect.size~years, cex=RIQmod2.LU12[,"N"]/50, data=RIQmod2.LU12, 
     pch = preeC)

simbolsD<-RIQmod3.LU12[,"pvalue"]<0.05
preeD<-as.numeric(simbolsD)
preeD[preeD==0]<-21
preeD[preeD==1]<-19

plot(Effect.size~years, cex=RIQmod3.LU12[,"N"]/50, data=RIQmod3.LU12, 
     pch = preeD)


##################
####  LU 9  ######
##################
LUmat<-LU9

RIQmod0.LU9<-matrix(NA,32,3)
colnames(RIQmod0.LU9)<-c("Effect.size","pvalue","N")
rownames(RIQmod0.LU9)<-colnames(LUmat)[2:33]

RIQmod1.LU9<-matrix(NA,32,3)
colnames(RIQmod1.LU9)<-c("Effect.size","pvalue","N")
rownames(RIQmod1.LU9)<- colnames(LUmat)[2:33]

RIQmod2.LU9<-matrix(NA,32,3)
colnames(RIQmod2.LU9)<-c("Effect.size","pvalue","N")
rownames(RIQmod2.LU9)<- colnames(LUmat)[2:33]

RIQmod3.LU9<-matrix(NA,32,3)
colnames(RIQmod3.LU9)<-c("Effect.size","pvalue","N")
rownames(RIQmod3.LU9)<- colnames(LUmat)[2:33]

listaAnos<-list()

#i<-2

#Removendo os dois últimos anos do looping (32 e 33) pois só tem 12 observações
for (i in 2:(ncol(LUmat)-2)){
  LUano<-colnames(LUmat[i])
  LU<-LUmat[,LUano]
  DFrod<-data.frame(Fish.Rich$Rich, covars, LU)
  colnames(DFrod)[1]<-"Rich"
  colnames(DFrod)[12]<-"LUrod"
  colnames(DFrod)
  print(LUano)
  
  #Filtrando-removendo dados de NA
  ptos.SNA<-na.exclude(DFrod)
  #Removendo sub-bacias com menos que 5 pontos
  ptos<-table(ptos.SNA$X9)[table(ptos.SNA$X9)>=5]
  nomesSel<-names(ptos)
  sub.dados<-ptos.SNA[ptos.SNA$X9 %in% nomesSel,]
  sub.dados<-droplevels(sub.dados)
  table(sub.dados$X9)
  
  #sub.dados<-DFrod
  
  mrod0 <- lmer(Rich ~ LUrod + (1 | X9), data = sub.dados,
                REML = FALSE, control = lmerControl(optimizer ="Nelder_Mead"))
  mrod1 <- lmer(Rich ~ LUrod + (LUrod | X9), data = sub.dados,
                REML = FALSE, control = lmerControl(optimizer ="Nelder_Mead"))
  mrod2 <- lmer(Rich ~ LUrod + topografiaCV + Area + (1 | X9), data = sub.dados,
                REML = FALSE, control = lmerControl(optimizer ="Nelder_Mead"))
  mrod3 <- lmer(Rich ~ LUrod + topografiaCV + Area + (LUrod | X9), data = sub.dados,
                REML = FALSE, control = lmerControl(optimizer ="Nelder_Mead"))
  
  resu.rod0<-summary(mrod0)
  BeP0<-resu.rod0$coefficients["LUrod",][c("Estimate","Pr(>|t|)")]
  resu.rod1<-summary(mrod1)
  BeP1<-resu.rod1$coefficients["LUrod",][c("Estimate","Pr(>|t|)")]
  resu.rod2<-summary(mrod2)
  BeP2<-resu.rod2$coefficients["LUrod",][c("Estimate","Pr(>|t|)")]
  resu.rod3<-summary(mrod3)
  BeP3<-resu.rod3$coefficients["LUrod",][c("Estimate","Pr(>|t|)")]
  
  RIQmod0.LU9[(i-1),] <- c(BeP0,nrow(sub.dados))
  RIQmod1.LU9[(i-1),] <- c(BeP1,nrow(sub.dados))
  RIQmod2.LU9[(i-1),] <- c(BeP2,nrow(sub.dados))
  RIQmod3.LU9[(i-1),] <- c(BeP3,nrow(sub.dados))
  
  }

years<-0:31

simbolsA<-RIQmod0.LU9[,"pvalue"]<0.05
preeA<-as.numeric(simbolsA)
preeA[preeA==0]<-21
preeA[preeA==1]<-19

plot(Effect.size~years, cex=RIQmod0.LU9[,"N"]/50, data=RIQmod0.LU9, 
     pch=preeA,col=1)

simbolsB<-RIQmod1.LU9[,"pvalue"]<0.05
preeB<-as.numeric(simbolsB)
preeB[preeB==0]<-21
preeB[preeB==1]<-19

plot(Effect.size~years, cex=RIQmod1.LU9[,"N"]/50, data=RIQmod1.LU9, 
     pch = preeB)

simbolsC<-RIQmod2.LU9[,"pvalue"]<0.05
preeC<-as.numeric(simbolsC)
preeC[preeC==0]<-21
preeC[preeC==1]<-19

plot(Effect.size~years, cex=RIQmod2.LU9[,"N"]/50, data=RIQmod2.LU9, 
     pch = preeC)

simbolsD<-RIQmod3.LU9[,"pvalue"]<0.05
preeD<-as.numeric(simbolsD)
preeD[preeD==0]<-21
preeD[preeD==1]<-19

plot(Effect.size~years, cex=RIQmod3.LU9[,"N"]/50, data=RIQmod3.LU9, 
     pch = preeD)


############PLOTS COMPOSTOS############
x11(height=11,width=9)

png(filename="fig1.png", width = 8, height = 10, units = "in", 
    bg = "white", res=600)
par(mfrow=c(4,3), oma=c(2, 1.3, 1, 0.1), mar=c(1, 2, 1, 0.1), mgp=c(1.5, 0.3, 0), cex=0.8, las=0, tcl=-0.3)

#MOD0
simbolsD<-RIQmod0.mont[,"pvalue"]<0.05
preeD<-as.numeric(simbolsD)
preeD[preeD==0]<-21
preeD[preeD==1]<-19
plot(Effect.size~years, ylab="", xlab="", cex=RIQmod0.mont[,"N"]/50, data=RIQmod0.mont, 
     pch = preeD)
mtext("Effect size", las=0, cex=1.2, 2, 1.8, outer=F)
mtext("Upland", las=0, cex=1.5, 3, 0.5, outer=F, font = 2)

simbolsD<-RIQmod0.LU12[,"pvalue"]<0.05
preeD<-as.numeric(simbolsD)
preeD[preeD==0]<-21
preeD[preeD==1]<-19
plot(Effect.size~years, xlab="", ylab="",cex=RIQmod0.LU12[,"N"]/50, data=RIQmod0.LU12, 
     pch = preeD)
mtext("SB-12", las=0, cex=1.5, 3, 0.5, outer=F, font = 2)

simbolsD<-RIQmod0.LU9[,"pvalue"]<0.05
preeD<-as.numeric(simbolsD)
preeD[preeD==0]<-21
preeD[preeD==1]<-19
plot(Effect.size~years, xlab="", ylab="", cex=RIQmod0.LU9[,"N"]/50, data=RIQmod0.LU9, 
     pch = preeD)
mtext("SB-9", las=0, cex=1.5, 3, 0.5, outer=F, font = 2)

#MOD2
simbolsD<-RIQmod2.mont[,"pvalue"]<0.05
preeD<-as.numeric(simbolsD)
preeD[preeD==0]<-21
preeD[preeD==1]<-19
plot(Effect.size~years, ylab="", xlab="", cex=RIQmod2.mont[,"N"]/50, data=RIQmod2.mont, 
     pch = preeD)
mtext("Effect size", las=0, cex=1.2, 2, 1.8, outer=F)

simbolsD<-RIQmod2.LU12[,"pvalue"]<0.05
preeD<-as.numeric(simbolsD)
preeD[preeD==0]<-21
preeD[preeD==1]<-19
plot(Effect.size~years, xlab="", ylab="",cex=RIQmod2.LU12[,"N"]/50, data=RIQmod2.LU12, 
     pch = preeD)

simbolsD<-RIQmod2.LU9[,"pvalue"]<0.05
preeD<-as.numeric(simbolsD)
preeD[preeD==0]<-21
preeD[preeD==1]<-19
plot(Effect.size~years, xlab="", ylab="", cex=RIQmod2.LU9[,"N"]/50, data=RIQmod2.LU9, 
     pch = preeD)
#MOD1
simbolsD<-RIQmod1.mont[,"pvalue"]<0.05
preeD<-as.numeric(simbolsD)
preeD[preeD==0]<-21
preeD[preeD==1]<-19
plot(Effect.size~years, xlab="", ylab="", cex=RIQmod1.mont[,"N"]/50, data=RIQmod1.mont, 
     pch = preeD)
mtext("Effect size", las=0, cex=1.2, 2, 1.8, outer=F)

simbolsD<-RIQmod1.LU12[,"pvalue"]<0.05
preeD<-as.numeric(simbolsD)
preeD[preeD==0]<-21
preeD[preeD==1]<-19
plot(Effect.size~years, xlab="", ylab="", cex=RIQmod1.LU12[,"N"]/50, data=RIQmod1.LU12, 
     pch = preeD)

simbolsD<-RIQmod1.LU9[,"pvalue"]<0.05
preeD<-as.numeric(simbolsD)
preeD[preeD==0]<-21
preeD[preeD==1]<-19
plot(Effect.size~years, xlab="", ylab="", cex=RIQmod1.LU9[,"N"]/50, data=RIQmod1.LU9, 
     pch = preeD)

#MOD3
simbolsD<-RIQmod3.mont[,"pvalue"]<0.05
preeD<-as.numeric(simbolsD)
preeD[preeD==0]<-21
preeD[preeD==1]<-19
plot(Effect.size~years, ylab="", xlab="", cex=RIQmod3.mont[,"N"]/50, data=RIQmod3.mont, 
     pch = preeD)
mtext("Years", las=0, cex=1.2, 1, 1.8, outer=F)
mtext("Effect size", las=0, cex=1.2, 2, 1.8, outer=F)

simbolsD<-RIQmod3.LU12[,"pvalue"]<0.05
preeD<-as.numeric(simbolsD)
preeD[preeD==0]<-21
preeD[preeD==1]<-19
plot(Effect.size~years, ylab="", xlab="", cex=RIQmod3.LU12[,"N"]/50, data=RIQmod3.LU12, 
     pch = preeD)
mtext("Years", las=0, cex=1.2, 1, 1.8, outer=F)

simbolsD<-RIQmod3.LU9[,"pvalue"]<0.05
preeD<-as.numeric(simbolsD)
preeD[preeD==0]<-21
preeD[preeD==1]<-19
plot(Effect.size~years, ylab="", xlab="", cex=RIQmod3.LU9[,"N"]/50, data=RIQmod3.LU9, 
     pch = preeD)
mtext("Years", las=0, cex=1.2, 1, 1.8, outer=F)

dev.off()


# x11(height=5,width=6)
# par(mfrow=c(1,1), mar=c(4, 3.5, .1, 0.1), cex=1, las=0, tcl=-0.3)
# vecpoints<-ifelse(RIQmod[,2]<0.05,19,1)
# plot(RIQmod[,1]~c(0:33), xlab = "", ylab="Effect size (richness)", xaxt="n", pch=vecpoints, mgp=c(2.2, 0.5, 0),tcl=-0.3,lwd=2, cex=1.2, cex.lab=1.2)
# #abline(lm(RIQmod[,1]~c(1985:2013)))
# summary(lm(RIQmod[,1]~c(0:33)))
# axis(side = 1, at = c(0:33), labels = seq(0,33,1) ,las=2, mgp=c(2.5, 0.5, 0),tcl=-0.3)
# mtext("Years", las=0, cex=1.2, 1, 2.8, outer=F)
# dev.copy2pdf(device = x11, file="FIGS/EffectRIQ_anos_47.pdf")








#########PLOTS CAMANA
dados_long_mont <- gather(LUmont, key = "Variable", value = "Value", -Name_fish) # Transformar os dados no formato longo para o ggplot
dados_long_mont$Variable <- factor(dados_long_mont$Variable, levels = colnames(LUmont)[-1]) # Definir a ordem das variáveis para o eixo X

# Criar o gráfico boxplot usando ggplot2
plot_mont <- ggplot(dados_long_mont, aes(x = Variable, y = Value)) +
  geom_boxplot(fill = "chartreuse4") +
  labs(y = "Native Land Cover (%)",x = "", title = "A")+
  scale_x_discrete(labels = as.character(seq(0, 31))) +
  theme_classic() +
  theme (axis.title.x = element_text(size = 16)) +
  theme (axis.text.y = element_text(size = 14)) +
  theme (axis.text.x = element_text(size = 12, angle = 0), axis.ticks.x = element_blank()) +
  theme (axis.title.y = element_text(size = 16)) +
  theme(plot.title = element_text(hjust = 1, vjust = 1))

##12##
dados_long_12 <- gather(LU12, key = "Variable", value = "Value", -Name_fish) # Transformar os dados no formato longo para o ggplot
dados_long_12$Variable <- factor(dados_long_12$Variable, levels = colnames(LUmont)[-1]) # Definir a ordem das variáveis para o eixo X

# Criar o gráfico boxplot usando ggplot2
plot_12 <- ggplot(dados_long_12, aes(x = Variable, y = Value)) +
  geom_boxplot(fill = "chartreuse4") +
  labs(y = "Native Land Cover (%)",x = "", title = "B")+
  scale_x_discrete(labels = as.character(seq(0, 31))) +
  theme_classic() +
  theme (axis.title.x = element_text(size = 16)) +
  theme (axis.text.y = element_text(size = 14)) +
  theme (axis.text.x = element_text(size = 12, angle = 0), axis.ticks.x = element_blank()) +
  theme (axis.title.y = element_text(size = 16)) +
  theme(plot.title = element_text(hjust = 1, vjust = 1))

##9##
dados_long_9 <- gather(LU9, key = "Variable", value = "Value", -Name_fish) # Transformar os dados no formato longo para o ggplot
dados_long_9$Variable <- factor(dados_long_9$Variable, levels = colnames(LUmont)[-1]) # Definir a ordem das variáveis para o eixo X

# Criar o gráfico boxplot usando ggplot2
plot_9 <-ggplot(dados_long_9, aes(x = Variable, y = Value)) +
  geom_boxplot(fill = "chartreuse4") +
  labs(y ="Native Land Cover (%)", x = "Years Before Sampling", y = "", title = "C")+
  scale_x_discrete(labels = as.character(seq(0, 31))) +
  theme_classic() +
  theme (axis.title.x = element_text(size = 16)) +
  theme (axis.text.y = element_text(size = 14)) +
  theme (axis.text.x = element_text(size = 12, angle = 0), axis.ticks.x = element_blank()) +
  theme (axis.title.y = element_text(size = 16)) +
  theme(plot.title = element_text(hjust = 1, vjust = 1))

##Exportando
plot_escalas <- plot_grid(plot_mont, plot_12,plot_9, ncol = 1, 
                          labels = "", rel_heights = c(100, 100, 100))
x11()
plot_escalas

ggsave("output/combined_plot.png", plot_escalas, width = 12, height = 15, units = "in")

