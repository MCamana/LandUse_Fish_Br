library(lme4)
library(lmerTest)

LUmont<-read.csv("DataOrganization/LUmont_corrected.csv", h=T, row.names = 1)
LU9<-read.csv("DataOrganization/LUx9_corrected.csv", h=T, row.names = 1)
LU12<-read.csv("DataOrganization/LUx12_corrected.csv", h=T, row.names = 1)
covars<-read.csv("DataOrganization/covars_corrected.csv", h=T, row.names = 1)

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

##REMOVENDO ANOS NAO USADOS NAS BACIAS
LUmont <- LUmont[, !(names(LUmont) %in% c('X29','X30','X31','X32','X33'))]
LU9 <- LU9[, !(names(LU9) %in% c('X29','X30','X31','X32','X33'))]
LU12 <- LU12[, !(names(LU12) %in% c('X29','X30','X31','X32','X33'))]

#######################
########MODELS MONT####
LUmat<-LUmont

RIQmod<-matrix(NA,(ncol(LUmat)-6),3)
colnames(RIQmod)<-c("Effect.size","pvalue","N")
rownames(RIQmod)<-colnames(LUmat)[-c(1:6)]

#listaAnos<-list()

#i<-7

#For comeca em 7 que e onde comecam os anos na matriz
for (i in 7:ncol(LUmat)){
  LUano<-colnames(LUmat[i])
  LUrod<-LUmat[,LUano]
  DFrod<-data.frame(Fish.Rich$Rich, covars, LUrod)
  
   # #Filtrando-removendo dados de NA
   # ptos.SNA<-na.exclude(DFrod)
   # #Removendo sub-bacias com menos que 5 pontos
   # ptos<-table(ptos.SNA$X9)[table(ptos.SNA$X9)>=5]
   # nomesSel<-names(ptos)
   # sub.dados<-ptos.SNA[ptos.SNA$X9 %in% nomesSel,]
   # DFrod<-droplevels(sub.dados)
  
  mrod <- lmer(Fish.Rich.Rich ~ LUrod + (1 | X9), data = DFrod,
               REML = FALSE, control = lmerControl(optimizer ="Nelder_Mead"))
  resu.rod<-summary(mrod)
  
  BeP<-resu.rod$coefficients["LUrod",][c("Estimate","Pr(>|t|)")]
  
  RIQmod[(i-6),] <- c(BeP,nrow(DFrod))
  
  print(i)
}


#Preditora
years<-0:28

###MODELO 0
#Objeto símbolos com p significativo e não significativo
simbolsA<-RIQmod[,"pvalue"]<0.05
preeA<-as.numeric(simbolsA)
preeA[preeA==0]<-21
preeA[preeA==1]<-19

preeA.mont<-preeA

# x11()
# par(mfrow=c(1,2))
plot(Effect.size~years, cex=2, data=RIQmod, pch=preeA, col=1)

RIQmod.mont<-RIQmod

######################
######## MODELS 12 ###
LUmat<-LU12

RIQmod<-matrix(NA,(ncol(LUmat)-6),3)
colnames(RIQmod)<-c("Effect.size","pvalue","N")
rownames(RIQmod)<-colnames(LUmat)[-c(1:6)]

#listaAnos<-list()

#For comeca em 7 que e onde comecam os anos na matriz
for (i in 7:ncol(LUmat)){
  LUano<-colnames(LUmat[i])
  LUrod<-LUmat[,LUano]
  DFrod<-data.frame(Fish.Rich$Rich, covars, LUrod)
  
   # #Filtrando-removendo dados de NA
   # ptos.SNA<-na.exclude(DFrod)
   # #Removendo sub-bacias com menos que 5 pontos
   # ptos<-table(ptos.SNA$X9)[table(ptos.SNA$X9)>=5]
   # nomesSel<-names(ptos)
   # sub.dados<-ptos.SNA[ptos.SNA$X9 %in% nomesSel,]
   # DFrod<-droplevels(sub.dados)
  
  mrod <- lmer(Fish.Rich.Rich ~ LUrod + (1 | X9), data = DFrod,
               REML = FALSE, control = lmerControl(optimizer ="Nelder_Mead"))
  resu.rod<-summary(mrod)
  
  BeP<-resu.rod$coefficients["LUrod",][c("Estimate","Pr(>|t|)")]
  
  RIQmod[(i-6),] <- c(BeP,nrow(DFrod))
  
  print(i)
}

#Preditora
years<-0:28

###MODELO 0
#Objeto símbolos com p significativo e não significativo
simbolsA<-RIQmod[,"pvalue"]<0.05
preeA<-as.numeric(simbolsA)
preeA[preeA==0]<-21
preeA[preeA==1]<-19

preeA.12<-preeA

#x11()
#par(mfrow=c(1,3))
plot(Effect.size~years, cex=2, data=RIQmod, pch=preeA.12, col=1)

RIQmod12<-RIQmod

######################
######## MODELS 9 ####
LUmat<-LU9

RIQmod<-matrix(NA,(ncol(LUmat)-6),3)
colnames(RIQmod)<-c("Effect.size","pvalue","N")
rownames(RIQmod)<-colnames(LUmat)[-c(1:6)]

#listaAnos<-list()

#i<-7

#For comeca em 7 que e onde comecam os anos na matriz
for (i in 7:ncol(LUmat)){
  LUano<-colnames(LUmat[i])
  LUrod<-LUmat[,LUano]
  DFrod<-data.frame(Fish.Rich$Rich, covars, LUrod)
  
   # #Filtrando-removendo dados de NA
   # ptos.SNA<-na.exclude(DFrod)
   # #Removendo sub-bacias com menos que 5 pontos
   # ptos<-table(ptos.SNA$X9)[table(ptos.SNA$X9)>=5]
   # nomesSel<-names(ptos)
   # sub.dados<-ptos.SNA[ptos.SNA$X9 %in% nomesSel,]
   # DFrod<-droplevels(sub.dados)
  
  mrod <- lmer(Fish.Rich.Rich ~ LUrod + (1 | X9), data = DFrod,
               REML = FALSE, control = lmerControl(optimizer ="Nelder_Mead"))
  resu.rod<-summary(mrod)
  
  BeP<-resu.rod$coefficients["LUrod",][c("Estimate","Pr(>|t|)")]
  
  RIQmod[(i-6),] <- c(BeP,nrow(DFrod))
  
  print(i)
}

#Preditora
years<-0:28

###MODELO 0
#Objeto símbolos com p significativo e não significativo
simbolsA<-RIQmod[,"pvalue"]<0.05
preeA<-as.numeric(simbolsA)
preeA[preeA==0]<-21
preeA[preeA==1]<-19

preeA.9<-preeA

#x11()
#par(mfrow=c(1,3))
plot(Effect.size~years, cex=2, data=RIQmod, pch=preeA.9, col=1)

RIQmod9<-RIQmod

#################
##### PLOTS #####

###Plot para os 366 pontos
#getwd()
png(filename="fig5mais.png", width = 10, height = 4, units = "in", 
    bg = "white", res=600)

par(mfrow=c(1,3), oma=c(2, 1.3, 0.1, 0.1), mar=c(1, 2, 0.1, 0.1), mgp=c(1.5, 0.3, 0), cex=0.8, las=0, tcl=-0.3)


#par(mfrow=c(1,3), cex.lab = 1.5, mar=c(4,4,1,1))

plot(Effect.size ~ years, cex = RIQmod.mont[,"N"]/100, data = RIQmod.mont, pch = preeA.mont, ylim = c(-0.11, -0.03),
     xlab = "", ylab ="") 

plot(Effect.size ~ years, cex = RIQmod12[,"N"]/100, data = RIQmod12, pch = preeA.12, ylim = c(-0.11, -0.03),
     xlab = "", ylab = "") 

plot(Effect.size ~ years, cex = RIQmod9[,"N"]/100, data = RIQmod9, pch = preeA.9, ylim = c(-0.11, -0.03),
     xlab = "", ylab = "")

mtext("Years before the sampling", side = 1, line = 1, outer = T)
mtext("Effect Size (R²)", side = 2, line = 0, outer = T)
dev.off()


# ###Plot para bacias 9 com mais de 5 pontos de coleta 
# getwd()
# png(filename="output/fig5mais.png", width = 12, height = 5, units = "in", 
#     bg = "white", res=600)
# 
# par(mfrow=c(1,3))
# par(cex.lab = 1.5)
# 
# plot(Effect.size ~ years, cex = RIQmod.mont[,"N"]/50, data = RIQmod.mont, pch = 19, ylim = c(-0.11, -0.03),
#      xlab = "", ylab = "Effect Size (R²)")
# mtext("Years before the sampling", side = 1, line = 2.5)
# 
# plot(Effect.size ~ years, cex = RIQmod12[,"N"]/50, data = RIQmod12, pch = 19, ylim = c(-0.11, -0.03),
#      xlab = "", ylab = "")
# mtext("Years before the sampling", side = 1, line = 2.5)
# 
# plot(Effect.size ~ years, cex = RIQmod9[,"N"]/50, data = RIQmod9, pch = 19, ylim = c(-0.11, -0.03),
#      xlab = "", ylab = "")
# mtext("Years before the sampling", side = 1, line = 2.5)
# 
# dev.off()
