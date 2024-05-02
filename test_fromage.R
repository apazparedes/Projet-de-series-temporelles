library(zoo)
library(tseries)
library(fUnitRoots)
#data= read.csv("C:\\Users\\elero\\OneDrive\\Documents\\2A ENSAE\\S2\\série tempo linéaires\\Projet-de-series-temporelles\\Data\\valeurs.csv")

date= as.yearmon(data$lib, "%Y-%m")

data$indice=as.numeric(data$indice)
serie = zoo(data$indice)


serie_ord <- zoo(data$indice, order.by=date)

serie_short = serie_ord[time(serie_ord) <= as.yearmon("2019-12")]


plot(serie_short, ylim=c(77, 107), 
     xlab="Année", ylab="Indice (base 100 en 2021)", 
     main="IPI - Exploitation de laiteries et fabrication de fromage")


acf(serie_short)
pacf(serie_short)

diff_serie= diff(serie_short,1)

summary(lm( serie_short ~ time(serie_short)))

diff_serie2= diff(serie_short,12)


adf=adfTest(diff_serie, lag=0, type = "ct")


Qtests <- function(series, k, fitdf=0) {
  pvals <- apply(matrix(1:k), 1, FUN=function(l) {
    pval <- if (l<=fitdf) NA else Box.test(series, lag=l, type="Ljung-Box", fitdf=fitdf)$p.value
    return(c("lag"=l,"pval"=pval))
  })
  return(t(pvals))
}

Qtests(adf@test$lm$residuals,24,length(adf@test$lm$coefficients))


acf(diff_serie)
pacf(diff_serie)
acf(diff_serie2)

diff_serie3= diff(diff_serie,12)
acf(diff_serie3)

diff_serie4= diff(diff_serie,1)
acf(diff_serie4)

plot(diff_serie)

pp.test(diff_serie)
adf.test(diff_serie)

adfTest(diff_serie, lag=0, type = "ct")

adfTest(diff_serie, lag=5, type = "ct")



#en utilisant le corr Q4 TD5

adf <- adfTest(serie_short, lag=0, type="ct") #
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
series <- serie_short; kmax <- 24; adftype="ct"
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
adf <- adfTest_valid(serie_short,24,adftype="ct")

#21 lag pour que les résidus soient des BB

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

#autre manière de transformer pr avoir qqch de stationnaire ? 




#test enlever deterministic trend 

trend= lm(serie_short~time(serie_short))
detrend=residuals(trend)
plot.ts(detrend)

acf(detrend)

diff_detrend=diff(detrend)

acf(diff_detrend)

diffbis_detrend= diff(diff_detrend)

acf(diffbis_detrend)

derniere=diff(log(serie_short))

plot.ts(derniere)

acf(derniere)

summary(lm(derniere ~ time(derniere)))
#
adf <- adfTest_valid(derniere,24,"nc")
#
Qtests(adf@test$lm$residuals, 24, fitdf = length(adf@test$lm$coefficients))
adf
