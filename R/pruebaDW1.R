#' Función `pruebaDW1()`
#'
#' Realiza test Durbin-Watson para autocorrelación de orden 1 en un MRLM. Usa la función `durbinWatsonTest()` de la librería `car`.
#'
#' @param modelo Un objeto tipo `lm` en el que se guardó el ajuste de un MRLM.
#' @return un objeto `data.frame` con las siguientes variables:
#' * `rho(1) estimado` estimación de la autocorrelación de orden 1 `rho(1)` usando los residuales del modelo.
#' * `Estadistico D-W` El valor del estadístico de la prueba.
#' * `VP H1: rho(1)>0` El valor P para el test `H0: rho(1)=0`, vs. `H1:rho(1)>0`.
#' * `VP H1: rho(1)<0` El valor P para el test `H0: rho(1)=0`, vs. `H1:rho(1)<0`.
#' @export
#' @examples
#' #Simulando serie con modelo de tendencia lineal y errores RB normales
#' t <- 1:300
#' yt <- 3500+200*t+rnorm(300,sd=4000)
#' yt <- ts(yt,freq=12,start=c(1980,1))
#' #Gráfico de la serie simulada
#' plot(yt)
#' #Ajuste del modelo de tendencia lineal y errores RB normales
#' mrl <- lm(yt~t)
#' summary(mrl)
#' #Aplicando función para test Durbin-Watson
#' pruebaDW1(modelo=mrl)
#' #Simulando serie con modelo de tendencia lineal y errores correlacionados
#' zt <- 3500+200*t+arima.sim(n=300,list(ar=0.7),sd=4000)
#' zt <- ts(zt,freq=12,start=c(1980,1))
#' #Gráfico de la serie simulada
#' plot(zt)
#' #Ajuste del modelo de tendencia lineal y errores RB normales
#' mrl2 <- lm(zt~t)
#' summary(mrl2)
#' #Aplicando función para test Durbin-Watson
#' pruebaDW1(modelo=mrl2)
pruebaDW1 <- function(modelo){
dwneg <- car::durbinWatsonTest(modelo,max.lag=1,method="normal",alternative="negative")
dwpos <- car::durbinWatsonTest(modelo,max.lag=1,method="normal",alternative="positive")
res <- data.frame(dwneg$r,dwneg$dw,dwpos$p,dwneg$p,row.names="Resultados")
names(res) <- c("rho(1) estimado","Estadistico D-W","VP H1: rho(1)>0","VP H1: rho(1)<0")
res
}
