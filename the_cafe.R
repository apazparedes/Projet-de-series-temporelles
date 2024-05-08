library(zoo)
library(tseries)
library(fUnitRoots)
library(lubridate)
library(aTSA)
library(forecast)
#Transformation du thé et du café

fichier2<-"~/2A ENSAE/S2/série tempo linéaires/Projet-de-series-temporelles/Data/valeurs_mensuelles2.csv"
data2 <- read.csv2(fichier2)
zoo_data2 <- as.yearmon(data2$lib, "%Y-%m")
xm2<- zoo(as.numeric(data2$indice), order.by = zoo_data2)


diff2<- diff(xm2,1)
k<-plot(cbind(xm2,diff2))


summary(lm(xm2 ~ time(xm2)))

adf <- adfTest(xm2, lag=0, type="ct") #
adf
#
Qtests <- function(series, k, fitdf=0) {
  pvals <- apply(matrix(1:k), 1, FUN=function(l) {
    pval <- if (l<=fitdf) NA else Box.test(series, lag=l, type="Ljung-Box", fitdf=fitdf)$p.value
    return(c("lag"=l,"pval"=pval))
  })
  return(t(pvals))
}

Qtests(adf@test$lm$residuals, 24, fitdf = length(adf@test$lm$coefficients))
#
#
series <- xm2; kmax <- 24; adftype="ct"
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

adf <- adfTest_valid(xm2,24,adftype="ct")
#
Qtests(adf@test$lm$residuals, 24, fitdf = length(adf@test$lm$coefficients))
#
adf
#
summary(lm(diff2 ~ time(diff2)))
#
adf <- adfTest_valid(diff2,24,"nc")
#
Qtests(adf@test$lm$residuals, 24, fitdf = length(adf@test$lm$coefficients))
adf
#

par(mfrow=c(1,2))
pacf(diff2,24);acf(diff2,24) 


pmax=3;qmax=2

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
arima212 <- arima(xm,c(2,1,2),include.mean=F)
arima212


Qtests(arima212$residuals,24,fitdf=3)

adj_r2 <- function(model){
  ss_res <- sum(model$residuals^2)
  ss_tot <- sum(diff[-c(1:max(p,q))]^2)
  p <- model$arma[1]
  q <- model$arma[2]
  n <- model$nobs-max(p,q)
  adj_r2 <- 1-(ss_res/(n-p-q-1))/(ss_tot/(n-1))
  return(adj_r2)
}

adj_r2(arima212)
# 
dev.off()
plot(arima212$residuals)
axis(side=1)