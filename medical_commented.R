

#==============================================================================#
##### Packages and library ####
#==============================================================================#

library(zoo)
library(tseries)
library(fUnitRoots)
library(lubridate)
library(aTSA)
library(astsa)
library(forecast)
library(ggplot2)

#==============================================================================#
##### Lecture de la série et transformation de la base en un objet ZOO : ####
#==============================================================================#

fichier<-"~/Projet-de-series-temporelles/Data/valeurs_mensuelles3.csv"
data <- read.csv2(fichier)
zoo_data <- as.yearmon(data$lib, "%Y-%m")
xm<- zoo(as.numeric(data$indice), order.by = zoo_data)

#==============================================================================#
##### Partie 1: ####
#==============================================================================#
#==============================================================================#
##### Question 2:  ####
#==============================================================================#

plot(xm, main="IPI de la fabrication d'instruments et de fournitures à usage 
     médical et dentaire",
     xlab="Années",
     ylab="")
acf(xm)

#==============================================================================#
"Au vu de la représentation graphique de la série ainsi que de sa fonction
d'autocorrélation, on peut supposer une tendance linéaire qui serait positive.
Nous pouvons intuiter ce résultat en régressant la série sur le temps :"
#==============================================================================#
summary(lm(xm ~ time(xm)))

#==============================================================================#
"Le coefficient associé à la tendance linéaire est bien positif et semble 
significatif. Cependant on ne peut pas le confirmer car les résidus sont
potentiellement corrélés ce qui rendrait le test non valide si c'était le cas.
Pour éliminer la tendance linéaire de notre série nous la différencions à 
l'ordre 1:"
#==============================================================================#

#==============================================================================#
##### Différenciation de la série à l'ordre 1 #####
#==============================================================================#

diff<- diff(xm,1)

#==============================================================================#
##### Question 3:  ####
#==============================================================================#
#==============================================================================#
##### Représentation graphique de la série et de la série différenciée #####
#==============================================================================#

k<-plot(cbind(xm, diff), 
        main = "IPI de la fabrication d'instruments et de fournitures à usage 
     médical et dentaire",
        xlab = "Années", ylab = c(xm="Série",
                                  diff="Série différenciée"),  
        names.arg = years,       
        col.axis = "black")   

#==============================================================================#
"La série en différence première semble relativement stable autour d'une 
constante nulle. Elle pourrait être stationnaire."

"On réalise un test de racines unitaires de 'Dickey-Fuller'. Cependant pour 
interpréter le test on veut d'abord vérifier que les résidus du modèle sont 
non corrélés."
#==============================================================================#

#==============================================================================#
#####Premier test ADF #####
#==============================================================================#

adf <- adfTest(xm, lag=0, type="ct") #
adf
#==============================================================================#
##### Attention le résultat de ce test n'est pas encore interprétable ! #####
#==============================================================================#

#==============================================================================#
##### Premier Q-test pour vérifier l'autocorrélation des résidus #####
#==============================================================================#

#Création de la fonction test
Qtests <- function(series, k, fitdf=0) {
  pvals <- apply(matrix(1:k), 1, FUN=function(l) {
    pval <- if (l<=fitdf) NA else Box.test(series, lag=l, type="Ljung-Box", 
                                           fitdf=fitdf)$p.value
    return(c("lag"=l,"pval"=pval))
  })
  return(t(pvals))
}

#On exécute le test
Qtests(adf@test$lm$residuals, 24, fitdf = length(adf@test$lm$coefficients))

#==============================================================================#
"#L’absence d’autocorrélation des résidus est rejetée à chaque fois à partir 
d'un retard de 4 sur la série des résidus. Le test ADF avec aucun retard n’est 
donc pas valide. Nous devons donc ajouter des retards de ∆Xt jusqu'à ce que les 
résidus ne soient plus autocorrélés."
#==============================================================================#

#==============================================================================#
#### Test pour savoir à partir de combien de retard le test ADF est valide ####
#==============================================================================#

series <- xm; kmax <- 24; adftype="ct"
#==============================================================================#
"Création de la fonction test : la fonction s'arrête lorsque on a le retard pour
lequel les résidus sont non corrélés"
#==============================================================================#
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

#On exécute le test
adf <- adfTest_valid(xm,24,adftype="ct")

#==============================================================================#
"Il a fallu considérer 5 retards pour supprimer l'autocorrélation des résidus."
#==============================================================================#

#==============================================================================#
####Q-test et test ADF en prenant en compte le nombre de retard nécessaire####
#==============================================================================#

Qtests(adf@test$lm$residuals, 24, fitdf = length(adf@test$lm$coefficients))
#
adf

#==============================================================================#
"L'hypothèse nulle de racine unitaire n'est pas rejetée. e modèle a donc au 
moins une tendance linéaire d'ordre 1. Régressons désormais la série 
différentiée sur le temps."
#==============================================================================#

#==============================================================================#
####Régression linéaire de la série différentiée sur le temps####
#==============================================================================#
summary(lm(diff ~ time(diff)))

#==============================================================================#
"Il n'y a bien ni constante ni tendance significative. Effectuons donc le test 
ADF dans le cas sans constante ni tendance, en vérifiant l’absence 
d'autocorrélation des résidus."
#==============================================================================#


#==============================================================================#
####Test ADF sans constante ni tendance####
#==============================================================================#
#
adf <- adfTest_valid(diff,24,"nc")
#

#==============================================================================#
####Q test sur série différenciée et ADF sans tendance####
#==============================================================================#
Qtests(adf@test$lm$residuals, 24, fitdf = length(adf@test$lm$coefficients))
#==============================================================================#
"L'hypothèse de non-corrélation des résidus n'est jamais rejetée à un niveau 
significatif. On peut donc réaliser le test ADF en considérant 4 retards."
#==============================================================================#

#==============================================================================#
##### Test ADF en considérant 4 retards #####
#==============================================================================#

adf
#
#==============================================================================#
"Le test rejette bien la racine unitaire. La série différenciée est bien 
stationnaire."
#==============================================================================#

#==============================================================================#
##### Partie 2 : #####
#==============================================================================#

#==============================================================================#
##### Question 4 : #####
#==============================================================================#

#==============================================================================#
##### Étape 1 : Identification #####
#==============================================================================#
par(mfrow=c(1,2))
pacf(diff,24,
     main="PACF  de la série différentiée",
     ylab="",
     xlab="Retards");
acf(diff,24,main="ACF de la série différentiée",
    ylab="",
    xlab="Retards") 

#==============================================================================#
"Les fonctions d’autocorrélations et d’autocorrélations partielles sont 
significatives à l’ordre 2 et 3 respectivement. On choisit donc qmax=2 et pmax=3 
pour notre modèle ARMA(p,q). Tous les modèles possibles sont des modèles 
ARIMA(p, 1, q) avec p ≤ pmax et q ≤ qmax."
#==============================================================================#

pmax=3;qmax=2

#==============================================================================#
##### Étape 2 : Estimation #####
#==============================================================================#

arima011 <- sarima(xm,p=0,d=1,q=1,details=FALSE)
arima111 <- sarima(xm,p=1,d=1,q=1,details=FALSE)
arima012 <- sarima(xm,p=0,d=1,q=2,details=FALSE)

#==============================================================================#
##### Étape 3 : Validité ? #####
#==============================================================================#

Qtests(arima111$fit$residuals,24,fitdf=3)
Qtests(arima012$fit$residuals,24,fitdf=3)
Qtests(arima011$fit$residuals,24,fitdf=3)

#==============================================================================#
"L'absence d'autocorrélation n'est jamais rejetée pour les trois modèles. Nous 
devons maintenant choisir le meilleur"
#==============================================================================#


#==============================================================================#
##### Étape 4 : Sélection  #####
#==============================================================================#

pqs <- expand.grid(0:pmax,0:qmax) #combinaisons possibles de p<=p* et q<=q*
mat <- matrix(NA, nrow=pmax+1, ncol=qmax+1)
rownames(mat) <- paste0("p=",0:pmax) #renomme les lignes
colnames(mat) <- paste0("q=",0:qmax) #renomme les colonnes
AICs <- mat #matrice ou assigner les AIC
BICs <- mat #matrice ou assigner les BIC
for (row in 1:dim(pqs)[1]){
  p <- pqs[row,1]
  q <- pqs[row,2]
  estim <- try(arima(diff,c(p,0,q), include.mean=F)) #tente d'estimer l'ARIMA
  AICs[p+1,q+1] <- if (class(estim)=="try-error") NA else estim$aic
  BICs[p+1,q+1] <- if (class(estim)=="try-error") NA else BIC(estim)
}
AICs
BICs
AICs==min(AICs)
#
BICs==min(BICs)
#

#==============================================================================#
"Le modèle MA1 est celui qui minimise les critères AICs et BICs. Nous chosissons 
donc celui-là."
#==============================================================================#


#==============================================================================#
##### Partie 3 :  #####
#==============================================================================#

as.Date(1990-01)
IPI<-ts(as.numeric(data$indice),start=c(1990,02),end=c(2019,13),frequency=12)
foreyear<-sarima.for(xdata=IPI,n.ahead=24,p=0,d=1,q=1,plot=TRUE,
                     plot.all = FALSE,
                    xlab="Années",
                    main="Prédictions à horizon de deux années")

foremonth <-sarima.for(xdata=IPI,n.ahead=2,p=0,d=1,q=1,plot=TRUE,
                          plot.all = FALSE,
                          xlab="Années",
                          main="Prédictions à horizon de deux mois")

