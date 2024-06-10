####Organização dos dados para gerar as matrizes para analisar####
####Libraries####
#install.packages("here")
library(here)
here()

###LAND USE MONTANTE
dataLUmont <- read.table (here::here("data","processed","historic_use_mont.txt"), header=T) ##Buscamos os dados brutos na pasta
str(dataLUmont)
head(dataLUmont)
rownames(dataLUmont)<-dataLUmont$Name_fish

##Ordenamos as linhas de acordo com os nomes em ordem alfabética dos códigos
dataLUmont<-dataLUmont[sort(rownames(dataLUmont)),]


###LAND USE ESCALA 12
dataLU12 <- read.table (here::here("data","processed","historic_use_12.txt"), header=T)
rownames(dataLU12)<-dataLU12$Name_fish
rownames(dataLU12)
nrow(dataLUmont)
nrow(dataLU12)

##Ordenamos as linhas de LU12
dataLU12<-dataLU12[sort(rownames(dataLU12)),]

##Conferir se é igual a montante
rownames(dataLU12) == rownames(dataLUmont)
all(rownames(dataLU12) == rownames(dataLUmont))


###LAND USE ESCALA 9
dataLU9 <- read.table (here::here("data","processed","historic_use_9.txt"), header=T)
rownames(dataLU9)<-dataLU9$Name_fish
rownames(dataLU9)

nrow(dataLUmont)
nrow(dataLU12)
nrow(dataLU9)

##Ordenamos as linhas de LU9
dataLU9<-dataLU9[sort(rownames(dataLU9)),]

##Conferir se é igual a montante e igual ao 12
rownames(dataLU9) == rownames(dataLUmont)
rownames(dataLU9) == rownames(dataLU12)
all(rownames(dataLU9) == rownames(dataLUmont))
all(rownames(dataLU9) == rownames(dataLU12))


#Códigos usados para obter os sites dos dados de peixes
codes<-dataLUmont$Name_fish
rownames(dataLUmont)==codes

#Carregando abiot para pegar os códigos originais das planilhas de peixes
abiot<-read.csv (here::here("data","raw","Abiotic.csv"), header=T)

#Filtrando só os dados ambientais codes
abiot.codes<-abiot[abiot$ID_Point %in% codes,]

length(codes)
length(abiot.codes$ID_Point)
#Conferindo se todos os codes de uso do solo estão nos abiot
#Alterando a ordem para mesma de codes
rownames(abiot.codes)<-abiot.codes$ID_Point
abiot.codes<-abiot.codes[sort(abiot.codes$ID_Point),]

#Tudo na mesma ordem
all(rownames(abiot.codes) == codes)

#Selecionando dados abióticos importantes
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
#i=1
####Selecionando dados de peixes com base nos códicos em abiot.codesSEL
fishlist<-list()
for (i in 1:length(unique(abiot.codesSEL$ID_Name))){
  nomeROD <- unique(abiot.codesSEL$ID_Name)[i]
  matROD <- datalist[[nomeROD]]
  abiROD <- abiot.codes[abiot.codesSEL$ID_Name==nomeROD,]
  
  print(i)
  print(nomeROD)
  
  matROD.sel1 <- matROD[matROD$ID_Point %in% abiROD$ID_Point,]
  
  nrow(matROD.sel1)
  nrow(abiROD)
  
  print(all(matROD.sel1$Name == abiROD$Name))

  fishlist[[i]] <- matROD.sel1
  names(fishlist)[i] <- nomeROD
    
  }

names(fishlist)

#Conferindo se os códigos batem
teste<-lapply(fishlist,function(x){x[,c(1:3)]})

mat.teste<-do.call(rbind, teste)
nrow(mat.teste)
length(codes)

all(mat.teste$ID_Point == codes)

##############################

####################
#######PEIXES######
####################
###REMOVENDO ESPÉCIES SEM OCORRÊNCIA
fishlist2<-lapply(fishlist, function(x){
rod<-x[,-c(1:7)]
sele<-colSums(rod)>0
rodsele<-rod[,sele]
cbind(x[,c(1:7)],rodsele)
})

##CONFERINDO
lapply(fishlist2,function(x)colSums(x[,-c(1:7)]))

###fishlist2 tem todas matrizes de peixes filtradas com somentes os pontos de codes e removendo spp com zero ocorrências.


###Carregando as covariáveis
covars<-read.table (here::here("data","processed","covars.txt"), h=T)
here::here()
###Rodando as análise para RIQUEZA
head(covars)

#Conferir a ordenação das linhas nas matrizes
dataLUmont$Name_fish
covars$Name_fish

rownames(covars)<-covars$Name_fish
covars<-covars[sort(rownames(covars)),]
rownames(covars) == rownames(dataLUmont)
all(rownames(covars) == rownames(dataLUmont))
rownames(covars) == mat.teste$ID_Point



cbind(rownames(covars),rownames(dataLUmont))

nrow(covars)
nrow(dataLUmont)
table(covars$Name_fish)


#Salvando todos os dados ordenados
write.csv(dataLUmont, "LUmont_corrected.csv")
write.csv(dataLU12, "LUx12_corrected.csv")
write.csv(dataLU9, "LUx9_corrected.csv")
write.csv(covars, "covars_corrected.csv")
saveRDS (fishlist2,file="listaPeixes_corrected.RData")



