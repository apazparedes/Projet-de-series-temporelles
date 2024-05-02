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


Qtests <- function(series, k, fitdf=0) {
  pvals <- apply(matrix(1:k), 1, FUN=function(l) {
    pval <- if (l<=fitdf) NA else Box.test(series, lag=l, type="Ljung-Box", fitdf=fitdf)$p.value
    return(c("lag"=l,"pval"=pval))
  })
  return(t(pvals))
}


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
