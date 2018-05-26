library(readxl)
ETANOL2 <- read_excel("ETANOL2.xlsx")
ETANOL2<-as.data.frame(ETANOL2)

#cargando los paquetes que vamos a utilizar
require(gamlss)
require(MASS)
require(corrplot)
require(splines)
require(plotly)
require(car)
require(scatterplot3d)
require(leaps)
require(persp3d)
#=================================================================
##ANALISIS DESCRIPTIVO DE LOS DATOS
#=================================================================

summary(ETANOL2)#estadisticas de resumen
str(ETANOL2)
ecor<-cor(ETANOL2)#para observar la correlacion de los datos
corrplot(ecor, type="upper", order="hclust", tl.col="black", tl.srt=45)

plot(density(ETANOL2$concentracion))

#realizando el gráfico pairs
# Mejorando el gráfico
panel.reg <- function (x, y)
{
  points(x, y, pch=20)
  abline(lm(y ~ x), lwd=2, col='dodgerblue2')
}
# Funci?n para crear el histograma
panel.hist <- function(x, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col="dodgerblue2", ...)
}
# Funcion para obtener la correlacion
panel.cor <- function(x, y, digits=2, prefix="", cex.cor)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- cor(x, y)
  txt <- format(c(r, 0.123456789), digits=digits)[1]
  txt <- paste(prefix, txt, sep="")
  if(missing(cex.cor)) cex <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex * abs(r))
}

pairs(ETANOL2,
      upper.panel = panel.reg,
      diag.panel = panel.hist,
      lower.panel = panel.cor)


###gráfico de disperción en 3d
pdf(file="dispercion.pdf", height=8, width=8)
scatterplot3d(x=ETANOL2$velocidad, y=ETANOL2$Temperatura, z=ETANOL2$concentracion, pch=16, cex.lab=1,
              highlight.3d=TRUE, type="h", main="Disperción de los datos",
              xlab='Velocidad de agitación (rpm)',
              ylab='Temperatura °C',
              zlab='Etanol (g/L)',box=T, cex.symbols = 2)
dev.off()

#=============================================================
####TOMANDO LA CONCENTRACIÓN COMO VARIABLE RESPUESTA
#============================================================

##modelos lm

mod01<-lm(concentracion ~ bs(Temperatura)+ velocidad, data=ETANOL2)
summary(mod01)
yg01<-mod01$fitted.values #extraer los valores ajustados
mse01 <- mean((ETANOL2$concentracion - yg01)^2)#
cor(yg01, ETANOL2$concentracion) #


mod0<-lm(concentracion ~ Temperatura, data=ETANOL2)#malo
summary(mod0) #0.492
anova(mod0)
plot(mod0)##gráficos de residuales
plot(x=ETANOL2$Temperatura,y=ETANOL2$concentracion, type="p", pch=19)
abline(mod0, col="red")
yg0<-mod0$fitted.values #extraer los valores ajustados
mse0 <- mean((ETANOL2$concentracion - yg0)^2)#4.139223
cor(yg0, ETANOL2$concentracion) #0.3770

#---este modelo es el mas malo
mod1<-lm(concentracion ~ velocidad, data=ETANOL2)
summary(mod1) #r2a=malisimo (negativo) -0.06
plot(x=ETANOL2$velocidad,y=ETANOL2$concentracion, type="p", pch=19)
abline(mod1, col="red")
yg1<-mod1$fitted.values
mse1 <- mean((ETANOL2$concentracion - yg1)^2)#6762.026
plot(mod1)
cor(yg1, ETANOL2$concentracion)

##---
mod2<-lm(concentracion ~ Temperatura+velocidad,data = ETANOL2)
summary(mod2) #r2a=0.4839
anova(mod2)
yg2<-mod2$fitted.values #extraer los valores ajustados
mse2 <- mean((ETANOL2$concentracion - yg2)^2)#3.7414
plot(mod2)
cor(yg2, ETANOL2$concentracion)


#modelos polinomiales  
#=================================================================================
#####MODELO ELEGIDO
#============================================================================

mod3<-lm(concentracion ~ Temperatura+velocidad+ I(Temperatura^2)+I(velocidad^2), data=ETANOL2)
summary(mod3)  #r2a=0.9506
anova(mod3)
yg3<-mod3$fitted.values #extraer los valores ajustados
mse3 <- mean((ETANOL2$concentracion - yg3)^2)#0.2686
cor(yg3, ETANOL2$concentracion)
plot(mod3)
#worm plot
pdf(file="wp.pdf", height=8, width=8)
wp(mod3)
dev.off()


Temp<- seq(from=41.4, to=86.6, length.out=30)
veloc<- seq(from=358.6, to=641.4, length.out=30)
Rend <- function(temp, veloc) {
  res <- coef(mod3) * c(1, temp, veloc, temp^2, veloc^2)
  sum(res)
}

Rend <- Vectorize(Rend)
Etanol  <- outer(Temp, veloc, Rend)
pdf(file="superficie.pdf", height=8, width=8)
persp(x=Temp, y=veloc, z=Etanol,zlab = "Etanol (g/L)", ylab="",
      xlab = "Temperatura", theta=20, phi=25, 
      ticktype = "detailed", col='rosybrown1')
dev.off()

tempe<-ETANOL2$Temperatura
velo<-ETANOL2$velocidad

minus_rend <- function(x) {
  tempe <- x[1]
  velo <- x[2]
  new.data <- data.frame(tempe=c(1, tempe), velo=c(1, velo))
  -predict(mod3, new.data)[2]
}

inicio <- c(0, 0)  # valores iniciales para la busqueda
names(inicio) <- c('Temperatura', 'velocidad') # Colocando nombres

res <- nlminb(start=inicio, objective=minus_rend,
              lower=c(41.4, 358.6), # minimos de las variables
              upper=c(86.6, 641.4), # maximos de las variables
              control=list(trace=1))

res # Para ver todo el contenido de res
res$par  # Valores optimos
-res$objective  # Valor del objetivo







mod4<-lm(concentracion ~ Temperatura+velocidad+ I(Temperatura^2), data=ETANOL2)
summary(mod4)  #r2a=0.89
anova(mod4)
yg4<-mod4$fitted.values #extraer los valores ajustados
mse4 <- mean((ETANOL2$concentracion - yg4)^2)#0.6890
cor(yg4, ETANOL2$concentracion)
min(ETANOL2$Temperatura)
max(ETANOL2$Temperatura)


mod5<-lm(concentracion~poly(Temperatura, degree = 3)+poly(velocidad, degree=3), data=ETANOL2)
summary(mod5) #r2a=0.9898
anova(mod5)
yg5<-mod5$fitted.values #extraer los valores ajustados
mse5 <- mean((ETANOL2$concentracion - yg5)^2)#0.0368
cor(yg5, ETANOL2$concentracion)

Temp<- seq(from=41.4, to=86.6, length.out=30)
veloc<- seq(from=358.6, to=641.4, length.out=30)
conc <- function(temp, veloc) {
  res <- coef(mod5) * c(1, temp, veloc, temp^2, veloc^2, temp^3,veloc^3)
  sum(res)
}

conc <- Vectorize(conc)
Etanol<- outer(Temp, veloc, conc )

persp(x=Temp, y=veloc, z=Etanol,
      theta=50, phi=30, ticktype = "detailed", col='salmon1')







mod6<-lm(concentracion~poly(Temperatura, degree =3 )+velocidad,data=ETANOL2)
summary(mod6)#0.8864
yg6<-mod6$fitted.values #extraer los valores ajustados
mse6 <- mean((ETANOL2$concentracion - yg6)^2)#3.741
cor(yg6, ETANOL2$concentracion)


#modelo con interaccion
mod7<-lm(concentracion ~ Temperatura+velocidad+ I(Temperatura^2)+I(velocidad^2)+Temperatura*velocidad, data=ETANOL2)
summary(mod7)#0.9407
anova(mod7)
yg7<-mod7$fitted.values #extraer los valores ajustados
mse7 <- mean((ETANOL2$concentracion - yg7)^2)#0.2685
cor(yg7, ETANOL2$concentracion)

#modelos gamlss
d<-fitDist(ETANOL2$concentracion, type = "realplus", k=2)
d$fits #WEI3(2) GA(2) LOGNO(2) IG(2) IGAMMA(2) EXP(2) PARETO(2)  

mod8<- gamlss(ETANOL2$concentracion~Temperatura +velocidad,
              data = ETANOL2, sigma.formula =~1,
              mu.formula=~ Temperatura+velocidad+I(Temperatura)^2,
              k=2, family=GA)
summary(mod8)#45.856
mu8 <- predict(object=mod8, newdata=ETANOL2, what='mu', type='response')
cor(ETANOL2$concentracion, mu8)
mse8 <- mean((ETANOL2$concentracion - mu8)^2)



mod9<-gamlss(ETANOL2$concentracion ~Temperatura+velocidad,
            data=ETANOL2, family=LOGNO, mu.formula=~Temperatura +
              velocidad, sigma.formula = ~Temperatura+velocidad)
summary(mod9)#45.885  
mu9 <- predict(object=mod9, newdata=ETANOL2, what='mu', type='response')
cor(ETANOL2$concentracion, mu9)
mse9 <- mean((ETANOL2$concentracion - mu9)^2)
wp(mod9)

mod10<-gamlss(ETANOL2$concentracion ~ Temperatura+velocidad,
             data=ETANOL2,family =IG, mu.formula=~Temperatura +
               velocidad, sigma.formula = ~Temperatura+velocidad,
             k=2)
summary(mod10)#44.7558
mu10 <- predict(object=mod10, newdata=ETANOL2, what='mu', type='response')
cor(ETANOL2$concentracion, mu10)
mse10 <- mean((ETANOL2$concentracion - mu10)^2)
wp(mod10)

mod11<-gamlss(ETANOL2$concentracion ~Temperatura+velocidad, data=ETANOL2,
              family=EXP, mu.formula=~Temperatura +velocidad+I(Temperatura))
summary(mod11)

mu11 <- predict(object=mod11, newdata=ETANOL2, what='mu', type='response')
cor(ETANOL2$concentracion, mu11)
mse11 <- mean((ETANOL2$concentracion - mu11)^2)
wp(mod11)

              
#=============================================================
####TOMANDO EL RENDIMIENTO COMO VARIABLE RESPUESTA
#============================================================

##grafico de dispercion de los datos

p <- plot_ly(data=ETANOL2, x=~velocidad, y=~Temperatura, z=~rendimiento) %>%
  add_markers() %>%
  layout(scene=list(xaxis=list(title='vel de agitacion'),
                    yaxis=list(title='Temperatura'),
                    zaxis=list(title='Rendimiento')))
p

##modelos lm
mod0<-lm(rendimiento ~ Temperatura, data=ETANOL2)
summary(mod0)
anova(mod0)
plot(mod0)
yg0<-mod0$fitted.values #extraer los valores ajustados
mse0 <- mean((ETANOL2$rendimiento - yg0)^2)#calcualando el mse

mod1<-lm(rendimiento ~ velocidad, data=ETANOL2)
summary(mod1) #r2a=malisimo (negativo)
yg1<-mod1$fitted.values
mse1 <- mean((ETANOL2$rendimiento - yg1)^2)
plot(mod1)

mod2<-lm(rendimiento ~ Temperatura+velocidad,data = ETANOL2)
summary(mod2) #r2a=0.4482
anova(mod2)
yg2<-mod2$fitted.values #extraer los valores ajustados
mse2 <- mean((ETANOL2$rendimiento - yg2)^2)#calcualando el mse***
plot(mod2)

#modelos polinomiales  *********
mod3<-lm(rendimiento ~ Temperatura+velocidad+ I(Temperatura^2)+I(velocidad^2), data=ETANOL2)
summary(mod3)  #r2a=0.9288
anova(mod3)
yg3<-mod3$fitted.values #extraer los valores ajustados
mse3 <- mean((ETANOL2$rendimiento - yg3)^2)#23.837
plot(mod3)

mod4<-lm(rendimiento ~ Temperatura+velocidad+ I(Temperatura^2), data=ETANOL2)
summary(mod4)  #r2a=0.86
anova(mod4)
yg4<-mod4$fitted.values #extraer los valores ajustados
mse4 <- mean((ETANOL2$rendimiento - yg4)^2)#52.006
plot(mod4)


mod5<-lm(rendimiento~poly(Temperatura, degree = 3)+poly(velocidad, degree=3), data=ETANOL2)
summary(mod5) #r2a=0.9813
anova(mod5)
yg5<-mod5$fitted.values #extraer los valores ajustados
mse5 <- mean((ETANOL2$rendimiento - yg5)^2)#4.16
plot(mod5)

mod6<-lm(rendimiento~poly(Temperatura)+velocidad,data=ETANOL2)
summary(mod6)#0.448
yg6<-mod6$fitted.values #extraer los valores ajustados
mse6 <- mean((ETANOL2$rendimiento - yg6)^2)#246.5192
plot(mod6)

#modelo con interaccion
mod7<-lm(rendimiento ~ Temperatura+velocidad+ I(Temperatura^2)+I(velocidad^2)+Temperatura*velocidad, data=ETANOL2)
summary(mod7)#0.9158
anova(mod7)
yg7<-mod7$fitted.values #extraer los valores ajustados
mse7 <- mean((ETANOL2$rendimiento - yg7)^2)#23.51153
plot(mod7)

#modelos gamlss
d<-fitDist(ETANOL2$rendimiento, type = "realplus", k=log(length(ETANOL2$rendimiento)))
d$fits #BCPEo (4) 94.494,GG (3) 97.864, GGB2(4) 100.26, WEI3 (2) 104.24

mod8<-gamlss(rendimiento~., family= GG, data=ETANOL2, k=2)
summary(mod8) #AIC 63.85794
wp(mod8)
yg8<-mod8$f #extraer los valores ajustados
mse8 <- mean((ETANOL2$rendimiento - yg5)^2)#4.16
plot(mod8)


mod9<-gamlss(rendimiento~., family= BCPEo(mu.link = "log", nu.link = "identity",
                                          sigma.link = "log", tau.link = "log"),
             data=ETANOL2, k=2)
summary(mod9)
wp(mod9)
plot(mod9)

mod8<-gamlss(rendimiento~., family= WEI3, data=ETANOL2, k=2)
summary(mod8)
