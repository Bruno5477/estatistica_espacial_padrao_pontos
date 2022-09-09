library(raster)
library(geoR)
library(readxl)
library(rgdal)
setwd("~/Etatistica/2021.2/ESTATÍSTICA ESPACIAL I/Desafio 2/bases")
# Lendo a base de dados
base <- read_excel("Temp_MG.xlsx")
dim(base)
names(base)


# Lendo o SHAPE


mapaMG <- readOGR("shape_minas/Mg_region.shp",layer="Mg_region")

# Mapa com pontos
par(mar=c(4,4,0,0),cex=1.4)
plot(mapaMG,xlab="longitude",ylab="latitude")
points(base$long,base$lat,col=2,pch=21,bg=5)
axis(1)
axis(2)


# Transformando em geodata (para usar as funcoes do pacote geoR)
dados <- as.geodata(base, coords.col = 4:3, data.col = 2)

plot(dados)

par(mfrow=c(1,1))
points(dados,cex.min=1, cex.max=5, pt.sizes="quintiles",
       col=terrain.colors(5))

#Obtenção do variograma
variograma <- variog(dados)
plot(variograma)

#Definição da família e busca por valores iniciais
eyefit(variograma)

#Estimando o efeito pepita:
mv <- likfit(dados, ini=c(15.42,0.0001159668), trend = trend.spatial(~coords, dados), cov.model="gaussian",
             fix.nugget = TRUE, nugget=5)
summary(mv)



plot(variograma)
lines(mv,lty=1,lwd=2,col="blue")

### Interpolação ponderada pelo inverso da distância

library(gstat)
library(sp)

base.idw <- as.data.frame(cbind(dados$coords,dados$data))
coordinates(base.idw) <- c("long", "lat")

summary(dados$coords)

### Krigagem

long.pred = seq(-50.63,-40.25,l=65)
lat.pred = seq(-23.04,-14.41,l=65)
grid.pred <- expand.grid(long.pred,lat.pred)


# krigagem
krigagem <- krige.conv(dados, loc=grid.pred, 
                       krige=krige.control(type.krige="SK",
                                           trend.d = trend.spatial(~coords, dados),
                                           trend.l = trend.spatial(~Var1+Var2, grid.pred), 
                                           beta=mv$beta, 
                                           nugget=mv$nugget, 
                                           cov.pars=mv$cov.pars))

hist(krigagem$predict)
# media
par(mfrow=c(1,1),mar=c(4,4,0.5,0.5))
image(krigagem, loc=grid.pred,
      col=terrain.colors(256))
plot(mapaMG,xlab="longitude",ylab="latitude",add=TRUE)


# variancia
par(mfrow=c(1,1),mar=c(4,4,0.5,0.5))
image(krigagem, loc=grid.pred, coords=dados$coords, 
      values=krigagem$krige.var, col=terrain.colors(256))
plot(mapaMG,xlab="longitude",ylab="latitude",add=TRUE)


##### Validação cruzada

valid <- xvalid(dados, model=mv)
plot(valid, dados)

#EQM de predição
EQM.pred <- mean(valid$error^2)

