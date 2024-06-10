####Este e o script para avaliar como as coberturas dos diferentes anos desde 1985
#### em tr?s escalas espaciais diferentes estao relacioandos 
##### com a Riqueza de peixes coletados ao longo do Brasil ####


###Esse é um script adaptado do script do mestrado####

####Libraries####
library(here)
# library(iNEXT)
# library(ape)
# library(vegan)


here()

dat_env_mont <- read.table (here::here("data","processed","historic_use_mont.txt"), header=T) ##Buscamos os dados brutos na pasta
str(dat_env_mont)
head(dat_env_mont)

#Códigos usados para obter os sites dos dados de peixes
codes<-dat_env_mont$Name_fish

#Carregando abiot para pegar os códigos originais das planilhas de peixes
abiot<-read.csv (here::here("data","raw","Abiotic.csv"), header=T)

abiot.codes<-abiot[abiot$ID_Point %in% codes,]

length(codes)
length(abiot.codes$ID_Point)
#Conferindo se todos os codes de uso do solo estão nos abiot
all(sort(abiot.codes$ID_Point) == sort(codes))


head(abiot.codes)
abiot.codesSEL<-abiot.codes[,c("ID_Point","ID_Name","Name","Latitude", "Longitude")]

####CARREGANDO MATRIZES DE PEIXES
#Assign folder object to the path where Biodata is. See example below.
Biofolder <- here::here("data","raw","fishes_sites")
fileslist <- list.files(path = Biofolder, pattern = "*.csv")
datalist <- lapply(fileslist, function(x)read.csv(paste(Biofolder,"/", x, sep=''), stringsAsFactors = T))

#Split list of csv files within biodata folder
nameslista<-gsub('.csv', '', fileslist)

names(datalist)<-nameslista


####Selecionando dados de peixes com base nos códicos em abiot.codesSEL
abiot.codesSEL

datalist$Fish_Amazonia_MT

names(datalist)





####################
#######RIQUEZA######
####################




####Avaliando o nível de completness das comunidades
##out2 <- iNEXT(fishes, datatype="abundance")
##ggiNEXT(out2, type=2, se=TRUE)
#ggiNEXT(out2, type=1, se=TRUE)
##out2$DataInfo #Usamos ponto de corte de nível de sample coverage (completness SC) = 0.95 (95%)
##NomeRemover<-out2$DataInfo[out2$DataInfo$SC<0.94,]$site
##NomeRemover<-as.character(NomeRemover)

###Filtrando todos pontos menos 1 que nao tem amostragem completa
##C.dat_env<-dat_env[dat_env$Name_fish!=NomeRemover,]
##C.fishes<-fishes[rownames(fishes)!=NomeRemover,]


#################################################
### Avaliando a correlaçao espacial de Moran ####
#################################################
coords<-read.csv("ptos_coordenadas.csv", h=T, sep=",")
#C.coords<-coords[coords$Codigo_Ca!=NomeRemover,]

plot(coords$X, coords$Y, cex=0.5, pch=19)
text(coords$X, coords$Y, labels=coords$Codigo_Ca)

all(C.coords$Codigo_Ca==C.dat19$Pontos)
all(C.coords$Codigo_Ca==rownames(C.peixes))

geo<-cbind(C.coords$X, C.coords$Y)
rownames(geo)<-C.coords$Codigo_Ca
samples.dist <- as.matrix( dist(geo) )
samples.dist.inv <- 1/samples.dist
diag(samples.dist.inv) <- 0

#Exemplo com a riqueza
riqTOT<-specnumber(fishes)
mod0 <- lm (fishes$Richness~dat_env$X0)
summary(mod0)
Moran.I(residuals(mod0), samples.dist.inv ,alternative="greater") #Riqueza está espacialmente autocorrelacionada

#Testando com composicao
#peixeslog<-log(C.peixes+1)
#PCOAfish<-capscale(peixeslog~1, distance="bray")
#summary(PCOAfish)
#Eixos<-scores(PCOAfish)$sites

Eixos
#Tem autocorrelacao espacial para composicao
modpcoa <- lm (Eixos[,1]~X1993+bacia+log(area), data=C.dat19)
summary(modpcoa)
Moran.I(residuals(modpcoa), samples.dist.inv ,alternative="greater") #Riqueza está espacialmente autocorrelacionada

matMor<-matrix(NA,nrow (Eixos),2)
colnames(matMor)<-c("MoranI","pvalue")
rownames(matMor)<-rownames(Eixos)

Eixos1<-Eixos
dat1<-C.dat19
geo1<-geo
peixes1<-C.peixes

vecEX<-c("CGR1903","CGR1904","QUA1403","SPS1404","SAT1403","CSR1904","SPS1401","CCA1301","BAG1502","CCE1303","BAG1503","CCE1301")
peixes1<-peixes1[!rownames(peixes1) %in% vecEX,]
dat1<-dat1[!dat1$Pontos %in% vecEX,]
geo1<-geo[!rownames(geo) %in% vecEX,]
riqTOT1<-specnumber(peixes1)

peixeslog1<-log(peixes1+1)
PCOAfish<-capscale(peixeslog1~1, distance="bray")
summary(PCOAfish)
Eixos1<-scores(PCOAfish)$sites

samples.dist1 <- as.matrix( dist(geo1) )
samples.dist.inv1 <- 1/samples.dist1
diag(samples.dist.inv1) <- 0

modpcoa1 <- lm (Eixos1[,1]~X1993+bacia+log(area), data=dat1)
summary(modpcoa1)
Moran.I(residuals(modpcoa1), samples.dist.inv1 ,alternative="greater") #Riqueza está espacialmente autocorrelacionada

modriq <- lm (riqTOT1~X1993+bacia+log(area), data=dat1)
Moran.I(residuals(modriq), samples.dist.inv1 ,alternative="greater") #Riqueza está espacialmente autocorrelacionada

plot(geo1[,1],geo1[,2], pch=19, cex=0.5)
text(geo1[,1],geo1[,2], label=rownames(geo1), cex = 0.8)

matMor<-matrix(NA,nrow (Eixos1),2)
colnames(matMor)<-c("MoranI","pvalue")
rownames(matMor)<-rownames(Eixos1)


###ABAIXO JA FOI RODADO PARA DECIDIR A RETIRADA DE PONTOS AUTOCORRELACIONADOS
# for (i in 1:nrow(Eixos1)) {
#   Rod<-Eixos1[,1]
#   Rod1<-Rod[-i]
#   Rod1env<-dat1[-i,]
# 
#   samples.dist <- as.matrix( dist(geo1[-i,]) )
#   samples.dist.inv <- 1/samples.dist
#   diag(samples.dist.inv) <- 0
#   
#   modpcoa <- lm (Rod1~X1993+bacia+log(area), data= Rod1env)
#   MorR<-Moran.I(residuals(modpcoa), samples.dist.inv ,alternative="greater") #Riqueza está espacialmente autocorrelacionada
#   matMor[i,]<-c(MorR$observed,MorR$p.value)
#     }
# matMor<-as.data.frame(matMor)
# teste<-matMor[order(matMor$MoranI),]
# 
# plot(teste$MoranI, 1:nrow(teste))
# text(teste$MoranI, 1:nrow(teste), label=rownames(teste), cex = 0.8)


#####FINAL FICAMOS COM peixes1 e dat1

#Riqueza total
riqTOT<-specnumber(peixes1)
####################################################

plot (riqTOT~dat1$area)
abline(lm(riqTOT~(dat1$area)))


#Avaliar normalidade dos dados de Riqueza
#install.packages("fitdistrplus")
# require(fitdistrplus)
# normal=fitdist(riqTOT,"norm")
# lnormal=fitdist(riqTOT,"lnorm")
# par(mfrow=c(1,2), mar=c(4,4,2,2))
# cdfcomp(list(normal,lnormal),horizontals=F, lwd=2,addlegend=T,legendtext=c("Normal","LNormal"))
# qqcomp(list(normal,lnormal),addlegend=T,legendtext=c("Normal","Lnormal"))###a Lognormal foi a melhor
# par(mfrow=c(1,1), cex.lab=1, cex.axis=1)

## Conferindo se os nomes das linhas (sitios amostrais) sao iguais
rownames(dat_env_mont)<-dat_env_mont$Name_fish
dat_env_mont<-dat_env_mont[,!colnames(dat_env_mont) %in% "Name_fish"]
rownames(dat_env_mont)==rownames(fishes)
str(dat_env_mont)

## Agora, criamos modelos de cada gravar o coef. de inclinação padronizado (tamanho do efeito) so da cobertura vegetal isolada
##Este primeiro modelos de teste será para a escala de montante
#RIQUEZA TOTAL
RIQmod_mont<-matrix(NA,34,2)
colnames(RIQmod)
colnames(RIQmod)<-c("Effect.size","pvalue")
rownames(RIQmod)<- colnames(dat_env_mont)[1:34]

for (i in 1:34){
  mod0 <- lm (fishes$Richness~dat_env_mont[,i])
  RIQmod[i,1] <- coef(mod0)[2]
  pval<-anova(mod0)[5][1,]
  RIQmod[i,2] <- pval
  }

RIQmod

x11(height=5,width=6)
par(mfrow=c(1,1), mar=c(4, 3.5, .1, 0.1), cex=1, las=0, tcl=-0.3)
vecpoints<-ifelse(RIQmod[,2]<0.05,19,1)
plot(RIQmod[,1]~c(0:33), xlab = "", ylab="Effect size (richness)", xaxt="n", pch=vecpoints, mgp=c(2.2, 0.5, 0),tcl=-0.3,lwd=2, cex=1.2, cex.lab=1.2)
#abline(lm(RIQmod[,1]~c(1985:2013)))
summary(lm(RIQmod[,1]~c(0:33)))
axis(side = 1, at = c(0:33), labels = seq(0,33,1) ,las=2, mgp=c(2.5, 0.5, 0),tcl=-0.3)
mtext("Years", las=0, cex=1.2, 1, 2.8, outer=F)
dev.copy2pdf(device = x11, file="FIGS/EffectRIQ_anos_47.pdf")

### Grafico para mat suplementar 1993 e 2013 ####
degrass1993<-100-dat1$X1993
degrass2013<-100-dat1$X2013

x11(height=5,width=6)
par(mfrow=c(1,1), mar=c(3.2, 3.2, .1, 0.1), cex=1, las=0, tcl=-0.3)
plot(riqTOT~degrass1993, xlab='Native vegetation conversion (%)',ylab = 'Species richness', mgp=c(2, 0.5, 0),tcl=-0.3, lwd=1, bg="orange", pch=21, cex=1.5, cex.lab=1.2)
points(riqTOT~degrass2013, xlab='Native vegetation conversion (%)',ylab = 'Species richness',  mgp=c(2, 0.5, 0),tcl=-0.3, lwd=1, bg="blue", pch=21, cex=1.5, cex.lab=1.2)
abline (lm(riqTOT~degrass1993),lwd=3, col="orange")
abline (lm(riqTOT~degrass2013),lwd=3, col="blue", lty = 2)
legend("topright", legend=c("1993", "2013"),col=c("black", "black"), cex=1, pt.bg=c("orange","blue"), pch=c(21,21),pt.cex=1.2, bty="o")
dev.copy2pdf(device = x11, file="FIGS/RiqXCob_47.pdf")


#@@@@@@@@@@@#
#    RDA    #
#@@@@@@@@@@@#
#Matriz de peixes SEM os pontos de 2019 
peixeslog<-log(peixes1+1)

##LOOP para obter R2 da RDA
COMPmod<-matrix(NA, 29,2)
colnames(COMPmod)
colnames(COMPmod)<-c("Effect.size","pvalue")
rownames(COMPmod)<- colnames(dat1)[1:29]

for (i in 1:29){
  mod0 <- rda(peixeslog~dat1[,i]+area+Condition(bacia), data=dat1)
  COMPmod[i,1] <- RsquareAdj(mod0)$r.squared
  pval<-anova(mod0, by="terms")$`Pr(>F)`[1]
  COMPmod[i,2] <- pval
  names(COMPmod)[i] <- colnames(dat19)[i]
}

COMPmod<-as.data.frame(COMPmod)
COMPmod
mod2001 <- rda(peixeslog~dat1$X2001+area+Condition(bacia), data=dat1)
anova(mod1993, by="terms")

mod2013 <- rda(peixeslog~dat1$X2013+area+Condition(bacia), data=dat1)
anova(mod2013, by="terms")

x11(height=5,width=6)
par(mfrow=c(1,1), mar=c(4, 3.5, .1, 0.1), cex=1, las=0, tcl=-0.3)
vecpoints1<-ifelse(COMPmod[,2]<0.05,19,1)
plot(COMPmod[,1]~c(1985:2013), xlab = "", ylab="Effect size (composition)", xaxt="n", mgp=c(2.2, 0.5, 0),tcl=-0.3,lwd=2, cex=1.2, cex.lab=1.2, pch=vecpoints1)
axis(side = 1, at = c(1985:2013), labels = seq(1985,2013,1) ,las=2, mgp=c(2.5, 0.5, 0),tcl=-0.3)
mtext("Years", las=0, cex=1.2, 1, 2.8, outer=F)
summary(lm(COMPmod[,1]~c(1985:2013)))
#abline(lm(COMPmod[,1]~c(1985:2013)))
#plot(lm(COMPmod[,1]~c(1985:2013)))
dev.copy2pdf(device = x11, file="FIGS/RDA_anos.pdf")


####### PLOT COMPOSTO RIQ E COMP ########
x11(height=10,width=6)
par(mfrow=c(2,1), mar=c(4, 3.5, .1, 0.1), cex=1, las=0, tcl=-0.3)
vecpoints<-ifelse(RIQmod[,2]<0.05,19,1)
plot(RIQmod[,1], xlab = "", ylab="Effect size (richness)", xaxt="n", mgp=c(2.2, 0.5, 0),tcl=-0.3,lwd=2, cex=1.2, cex.lab=1.4, ylim=c(0.15,0.56), pch=vecpoints)
axis(side = 1, at = seq(1,29,1), labels = seq(1985,2013,1) ,las=2, mgp=c(2.5, 0.5, 0),tcl=-0.3)
legend("topright", legend=c("A"),cex=1.5, bty="n", text.font=2)

vecpoints1<-ifelse(COMPmod[,2]<0.05,19,1)
plot(COMPmod[,1], xlab = "", ylab="Effect size (composition)", xaxt="n", mgp=c(2.2, 0.5, 0),tcl=-0.3,lwd=2, cex=1.2, cex.lab=1.4, ylim=c(0.05,0.072), pch=vecpoints1)
axis(side = 1, at = seq(1,29,1), labels = seq(1985,2013,1) ,las=2, mgp=c(2.5, 0.5, 0),tcl=-0.3)
mtext("Years", las=0, cex=1.4, 1, 3, outer=F)
legend("topright", legend=c("B"),cex=1.5, bty="n",text.font=2)
dev.copy2pdf(device = x11, file="FIGS/Plot_composto.pdf")

png("FIGS/Plot_composto.png", width = 6, height =10, units = 'in', res = 600)
par(mfrow=c(2,1), mar=c(4, 3.5, .1, 0.1), cex=1, las=0, tcl=-0.3)
vecpoints<-ifelse(RIQmod[,2]<0.05,19,1)
plot(RIQmod[,1], xlab = "", ylab="Effect size (richness)", xaxt="n", mgp=c(2.2, 0.5, 0),tcl=-0.3,lwd=2, cex=1.2, cex.lab=1.4, ylim=c(0.15,0.56), pch=vecpoints)
axis(side = 1, at = seq(1,29,1), labels = seq(1985,2013,1) ,las=2, mgp=c(2.5, 0.5, 0),tcl=-0.3)
legend("topright", legend=c("A"),cex=1.5, bty="n", text.font=2)

vecpoints1<-ifelse(COMPmod[,2]<0.05,19,1)
plot(COMPmod[,1], xlab = "", ylab="Effect size (composition)", xaxt="n", mgp=c(2.2, 0.5, 0),tcl=-0.3,lwd=2, cex=1.2, cex.lab=1.4, ylim=c(0.05,0.072), pch=vecpoints1)
axis(side = 1, at = seq(1,29,1), labels = seq(1985,2013,1) ,las=2, mgp=c(2.5, 0.5, 0),tcl=-0.3)
mtext("Years", las=0, cex=1.4, 1, 3, outer=F)
legend("topright", legend=c("B"),cex=1.5, bty="n",text.font=2)
dev.off()





