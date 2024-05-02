library(zoo)
library(tseries)
library(fUnitRoots)

data=read.csv2("~/2A ENSAE/S2/série tempo linéaires/Projet-de-series-temporelles/Data/valeurs_mensuelles.csv")

date= as.yearmon(data$lib, "%Y-%m")

data$indice=as.numeric(data$indice)
serie = zoo(data$indice, order.by = date)

plot(serie, ylim=c(77, 107), 
     xlab="Année", ylab="Indice (base 100 en 2021)", 
     main="IPI - Fabrication de cacao, chocolat et confiseries")

acf(serie)

diff_serie= diff(serie,1)

acf(diff_serie)

plot.ts(diff_serie)


summary(lm(serie~time(serie)))

adf <- adfTest(serie, lag=0, type="nc") #
adf
#
Qtests <- function(series, k, fitdf=0) {
  pvals <- apply(matrix(1:k), 1, FUN=function(l) {
    pval <- if (l<=fitdf) NA else Box.test(series, lag=l, type="Ljung-Box", fitdf=fitdf)$p.value
    return(c("lag"=l,"pval"=pval))
  })
  return(t(pvals))
}
?adfTest

str(adf)
Qtests(adf@test$lm$residuals, 24, fitdf = length(adf@test$lm$coefficients))
#
#
series <- serie ; kmax <- 24; adftype="ct"
adfTest_valid <- function(series, kmax, adftype){
  k <- 0
  noautocorr <- 0
  while (noautocorr==0){
    cat(paste0("ADF with ",k," lags: residuals OK? "))
    adf <- adfTest(series, lags=k, type=adftype)
    pvals <- Qtests(adf@test$lm$residuals, 24, fitdf = length(adf@test$lm$coefficients))[,2]
    if (sum(pvals<0.05,na.rm=T)==0) {
      noautocorr <- 1; cat("OK \n")
    } else cat("nope \n")
    k <- k+1
  }
  return(adf)
}
adf <- adfTest_valid(serie,24,adftype="ct")

# 3 lag pour que les résidus soient des BB

#
Qtests(adf@test$lm$residuals, 24, fitdf = length(adf@test$lm$coefficients))
#
adf
#
summary(lm(diff_serie ~ time(diff_serie)))
#
adf <- adfTest_valid(diff_serie,24,"nc")
#
Qtests(adf@test$lm$residuals, 24, fitdf = length(adf@test$lm$coefficients))
adf

# 14 lag pour avoir des résidus non corrélés 