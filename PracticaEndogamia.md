# Práctica Endogamia
R Core Development Team
Packages: “adegenet”, “pegas” y “hierfstat”


Primero debemos instalar las librerías que vamos a utilizar.
```
install.packages (“adegenet”)
install.packages (“pegas”)
install.packages (“hierfstat”)
```
Si la librería ya está instalada. Únicamente necesitamos llamarla:
```
library(adegenet)
library (pegas)
library (hierfstat)
```
Empezaremos con análisis básicos de genética de poblaciones para determinar la presencia de desequilibrio de Hardy-Weinberg en las poblaciones. Utilizaremos el set de datos proporcionado para el curso y la libreria “adegenet”: Micros_Genepop_Mexicanus.gen
```
a <- read.genepop (“Micros_Genepop_Mexicanus.gen”, ncode = 3)
a
```
```
toto <- hw.test(a, res="matrix")
dim(toto)
colnames(toto)
idx <- which(toto<0.001,TRUE)
idx
toto <- hw.test(a, res="full")
mapply(function(i,j) toto[[i]][[j]], idx[,2], idx[,1], SIMPLIFY=FALSE)
```
El desequilibrio de Hardy-Weinberg indica que alguna fuerza evolutiva está actuando en nuestra población.

## Para estimar el coeficiente de endogamia por población
```
b <- basic.stats(a)
Ho.locus <- b$Ho
Hs.locus <- b$Hs
Fis.LA <- 1-(mean(Ho.locus[,1])/mean(Hs.locus[,1]))
Fis.LS <- 1-(mean(Ho.locus[,2])/mean(Hs.locus[,2]))
Fis.AT <- 1-(mean(Ho.locus[,3])/mean(Hs.locus[,3]))
``` 
Escribimos los promedios por población en un sólo objeto
```
Fis <- rbind(Fis.LA, Fis.LS, Fis.AT)
```

## Determinaremos si los individuos presentes en nuestra población presentan endogamia.
Necesitamos conocer el valor de N para el set de datos
```
a
```
// 58 individuals; 10 loci; 67 alleles; size: 33.4 Kb

Pero como son diploides, entonces necesitamos saber el valor de 2N

2*58

## Vamos a estimar la endogamia promedio para los perritos de las praderas Mexicanos
```
temp <- inbreeding(a, N=116)
class(temp)
head(names(temp))
head (temp [[1]],20)
```
temp es una lista de valores muestreados a partir de la distribución de la verosimilitud de cada individuo; Obtenemos los valores promedio para todos los individuos utilizando sapply:
```
Fbar <- sapply (temp, mean)
```
Obtenemos un histograma
```
hist(Fbar, col="firebrick", main="Endogamia promedio en los perritos llaneros Mexicanos")
```
Observamos que algunos individuos tienen valores de endogamia altos (>0.6)
```
which (Fbar>0.6)
F <- inbreeding (a, res.type="function") [which(Fbar>0.6)]
F
```
El resultado obtenido de F es un tanto críptico, pero podemos observarlo a partir de un gráfico.
```
plot(F$LA20, main=paste("Endogamia en el individuo", names(F)), xlab="Endogamia (F)", ylab="Densidad de la Probabilidad")
```
## Podemos crear un histograma de la F individual para cada población.

Primero creamos objetos con el subset de individuos de cada población.
```
LA_F <- Fbar[c(1:24)]
LS_F <- Fbar[c(25:42)]
AT_F <- Fbar[c(43:58)]
```
Hacemos un histograma con cada uno de los subsets.
```
par(mfrow=c(3,1))
hist (LA_F, col = "darkorchid1", xlab="Individual inbreeding coefficient (F)", ylab= "Number of individuals", main = "Los Angeles (LA)", ylim=c(0, 6), xlim = c(0,0.7), border = "white", breaks = 14)
abline (v=mean(LA_F), lwd=3, lty=2)
abline (v=median(LA_F), lwd=3, lty=1)
hist (LS_F, col = "chocolate4", xlab="Individual inbreeding coefficient (F)", ylab= "Number of individuals", main = "La Soledad (LS)", ylim=c(0, 9), xlim = c(0,0.7), border = "white")
abline (v=mean(LS_F), lwd=3, lty=2)
abline (v=median(LS_F), lwd=3, lty=1)
hist (AT_F, col = "darkgoldenrod1", xlab="Individual inbreeding coefficient (F)", ylab= "Number of individuals", main = "Artesillas (AT)", ylim=c(0, 4), xlim = c(0,0.7), border = "white", breaks = 14)
abline (v=mean(AT_F), lwd=3, lty=2)
abline (v=median(AT_F), lwd=3, lty=1)
```

