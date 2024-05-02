
fichier<-"~/2A ENSAE/S2/série tempo linéaires/Projet-de-series-temporelles/Data/valeurs_mensuelles.csv"
data <- read.csv2(fichier)

colnames(data) <- c("lib", "indice","codes")
data<-subset(data,codes!="")
