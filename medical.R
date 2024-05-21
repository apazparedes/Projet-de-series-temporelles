library(zoo)
library(tseries)
library(fUnitRoots)
library(lubridate)
library(aTSA)
library(astsa)
library(forecast)
library(ggplot2)

#Fabrication d'instruments et de fournitures à usage médical et dentaire

#Lecture de la série et transformation de la base en un objet ZOO :
fichier3<-"~/Projet-de-series-temporelles/Data/valeurs_mensuelles3.csv"
data3 <- read.csv2(fichier3)
zoo_data3 <- as.yearmon(data3$lib, "%Y-%m")
xm3<- zoo(as.numeric(data3$indice), order.by = zoo_data3)

#La série choisie représente l'indice de production de la fabrication d'instruments et de fournitures à usage médical et dentaire.
plot(xm3, main="IPI de la fabrication d'instruments et de fournitures à usage 
     médical et dentaire",
     xlab="Années",
     ylab="")
acf(xm3)

#Au vu de la représentation graphique de la série ainsi que de sa fonction d'autocorrélation, on peut supposer une tendance linéaire qui serait positive.
#Nous pouvons intuiter ce résultat en régressant la série sur le temps : 
summary(lm(xm3 ~ time(xm3)))

#Le coéfficient associé à la tendance linéaire est bien positif et sembl significatif. Cependant on ne peut pas le confirmer car les résidus sont potentiellement corrélés ce qui rendrait le test non valide si c'était le cas.


#Pour éliminer la tendance linéaire de notre série nous la différencions à l'ordre 1:
#Différenciation de la série à l'odre 1:
diff3<- diff(xm3,1)

k<-plot(cbind(xm3, diff3), 
     main = "IPI de la fabrication d'instruments et de fournitures à usage 
     médical et dentaire",
     xlab = "Années", ylab = c(xm3="Série",
                              diff3="Série différenciée"),  
     names.arg = years,       
     col.axis = "black")   
#La série en différence première semble relativement stable autour d'une constante nulle. Elle pourrait être stationnaire.


#On réalise un test de racines unitaires de "Dickey-Fuller". Cependant pour interpréter le test on veut d'abord vérifier que les résidus du modèle sont non corrélés.
adf <- adfTest(xm3, lag=0, type="ct") #
adf


Qtests <- function(series, k, fitdf=0) {
  pvals <- apply(matrix(1:k), 1, FUN=function(l) {
    pval <- if (l<=fitdf) NA else Box.test(series, lag=l, type="Ljung-Box", fitdf=fitdf)$p.value
    return(c("lag"=l,"pval"=pval))
  })
  return(t(pvals))
}

Qtests(adf@test$lm$residuals, 24, fitdf = length(adf@test$lm$coefficients))
#

#L’absence d’autocorrélation des résidus est rejetée à chaque fois à partir d'un retard de 4. Le test ADF avec aucun retard n’est donc pas valide. 
#Nous devons donc ajouter des retards de ∆Xt jusqu'à ce que les résidus ne soient plus autocorrélés.

series <- xm3; kmax <- 24; adftype="ct"
adfTest_valid <- function(series, kmax, adftype){
  k <- 0
  noautocorr <- 0
  while (noautocorr==0){
    cat(paste0("ADF with ",k," lags: residuals OK? "))
    adf <- adfTest(series, lags=k, type=adftype)
    pvals <- Qtests(adf@test$lm$residuals, 24, 
                    fitdf = length(adf@test$lm$coefficients))[,2]
    if (sum(pvals<0.05,na.rm=T)==0) {
      noautocorr <- 1; cat("OK \n")
    } else cat("nope \n")
    k <- k+1
  }
  return(adf)
}

adf <- adfTest_valid(xm3,24,adftype="ct")

#Il a fallu considérer 5 retards pour supprimer l'autocorrélation des résidus.
#
Qtests(adf@test$lm$residuals, 24, fitdf = length(adf@test$lm$coefficients))
#
adf
#L'hypothèse nulle de racine unitaire n'est pas rejetée. Le modèle a donc au moins une tendance linéaire d'ordre 1.
#Régressons désormais la série différentiée sur le temps. 
summary(lm(diff3 ~ time(diff3)))

#Il n'y a bien ni constante ni tendance significative. Effectuons donc le test ADF dans le cas sans constante ni
#tendance, en vérifiant l’absence autocorrélation des résidus.
#
adf <- adfTest_valid(diff3,24,"nc")
#

Qtests(adf@test$lm$residuals, 24, fitdf = length(adf@test$lm$coefficients))
#L'hypothèse de non-corrélation des résidus n'est jamais rejetée à un niveau significatif. On peut donc réaliser le test ADF en considérant 4 retards.
adf
#
#Le test rejette bien la racine unitaire. La série différenciée est bien stationnaire.
par(mfrow=c(1,2))
pacf(diff3,24,
     main="PACF  de la série différentiée",
     ylab="",
     xlab="Retards");
acf(diff3,24,main="ACF de la série différentiée",
    ylab="",
    xlab="Retards") 


pmax=3;qmax=2

# Q5 ####

# Q5 ####
pqs <- expand.grid(0:pmax,0:qmax) #combinaisons possibles de p<=p* et q<=q*
mat <- matrix(NA, nrow=pmax+1, ncol=qmax+1)
rownames(mat) <- paste0("p=",0:pmax) #renomme les lignes
colnames(mat) <- paste0("q=",0:qmax) #renomme les colonnes
AICs <- mat #matrice ou assigner les AIC
BICs <- mat #matrice ou assigner les BIC
for (row in 1:dim(pqs)[1]){
  p <- pqs[row,1]
  q <- pqs[row,2]
  estim <- try(arima(diff3,c(p,0,q), include.mean=F)) #tente d'estimer l'ARIMA
  AICs[p+1,q+1] <- if (class(estim)=="try-error") NA else estim$aic
  BICs[p+1,q+1] <- if (class(estim)=="try-error") NA else BIC(estim)
}
AICs
BICs
AICs==min(AICs)
#
BICs==min(BICs)
#
 
arima011 <- arima(xm3,c(0,1,1),include.mean=F)
arima212
arima011

Qtests(arima011$residuals,24,fitdf=3)

adj_r2 <- function(model){
  ss_res <- sum(model$residuals^2)
  ss_tot <- sum(diff[-c(1:max(p,q))]^2)
  p <- model$arma[1]
  q <- model$arma[2]
  n <- model$nobs-max(p,q)
  adj_r2 <- 1-(ss_res/(n-p-q-1))/(ss_tot/(n-1))
  return(adj_r2)
}

adj_r2(arima011)

arimabis <- sarima(xm3,p=0,d=1,q=1,details=FALSE)
as.Date(1990-01)
xm4<-ts(as.numeric(data3$indice),start=c(1990,02),end=c(2019,13),frequency=12)



forebis<-sarima.for(xdata=xm4,n.ahead=24,p=0,d=1,q=1,plot=TRUE,plot.all = FALSE,
                    xlab="Années",
                    main="Prédictions à horizon de deux années")

forebismonth <-sarima.for(xdata=xm4,n.ahead=2,p=0,d=1,q=1,plot=TRUE,plot.all = FALSE,
                    xlab="Années",main="Prédictions à horizon de deux mois")
                  


foreplot<- plot(y=c(xm4,forebis$pred),x=time)


time<-period(num=390,units="month")

y=ym(9001)
ybis<-y+months(0:11)



for (i in seq(0,29)) {
  yter<-ybis+years(i)
  print(yter)}

y<-ym(2001)+months(0:23)
time<-as.Date(time(xm3))
time<-c(time,y)

for (i in seq(0,29)) {
  print(i)}

plot(arima011$residuals,xlab="Années",ylab="",main="Résidus du modèle ARIMA(0,1,1)");
plot(fore,xlab="Années",ylab="",main="Forecast du modèle ARIMA(0,1,1)")