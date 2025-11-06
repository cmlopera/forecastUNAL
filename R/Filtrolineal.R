#' Función `Filtrolineal()`
#'
#' Ajusta y pronostica series no estacionales a través de filtros lineales basados en medias móviles simples unilaterales o bilaterales, usando la estrategia circular.
#' @param serie Un objeto tipo `ts` con los valores de la serie a ajustar y pronosticar.
#' @param pesos.wi Vector de pesos de las medias móviles para la estimación del nivel. Para medias móviles simples unilaterales debe ser de longitud `m+1`, para `i=0,1,...,m`, y para bilaterales debe ser de longitud `2m+1`, para `i=-m,...,m`.
#' @param tipo El tipo de medias móviles para la estimación del nivel, `tipo=1` para medias móviles simples unilaterales, `tipo=2` para medias móviles simples bilaterales.
#' @param h Número de períodos futuros a pronosticar, por defecto es 1.
#' @param level El nivel de confianza para los intervalos de pronóstico, por defecto es 0.95.
#'
#' @return Una lista con los siguientes componentes:
#' * `nivel` Objeto tipo `ts` de las medias móviles calculadas con la estrategia circular.
#' * `fitted` Objeto tipo `ts` con los valores ajustados de la serie.
#' * `residuals` Objeto tipo `ts` con los valores de la serie de tiempo de los residuos del ajuste.
#' * `forecast` Objeto tipo `mts` con los valores de los pronósticos puntuales y por Intervalos para los `h` períodos.
#' * `RMSE` Estimacion de la raíz cuadrada del MSE del ajuste.
#' @export
#' @examples
#' #Serie de la tasa de variación trimestral del PIB in UK 1990:Q1 a 2004:Q1
#' datos <- ts(GDPchange,freq=4,start=c(1990,1))
#' #gráfico de la serie
#' plot(datos,type="b",pch=19,lwd=3)
#'
#' m <- 3 #tamaño muestra de validación cruzada
#' n <- length(datos)-m #tamaño muestra de ajuste
#' t <- 1:n
#' yt <- ts(datos[t],freq=4,start=c(1990,1))
#' tnuevo <- (n+1):length(datos)
#' ytnuevo <- ts(datos[tnuevo],freq=4,start=c(2003,3)) #pronósticos inician en Q3-2003
#'
#' #Ajuste y pronósticos con Filtro de Henderson de 2m+1=7 términos
#' w7 <- c(-0.059,0.059,0.294,0.413,0.294,0.059,-0.059) #Vector de pesos wi
#' FH7 <- Filtrolineal(serie=yt,pesos.wi=w7,tipo=2,h=m) #Obtención de ajustes y pronósticos
#' ythatFH7 <- fitted(FH7);ythatFH7 #Valores ajustados
#'
#' #Serie y su ajuste
#' plot(datos,lwd=2)
#' lines(ythatFH7,col=2,lwd=2)
#' legend("bottomright",legend=c("Obs.","FH7"),col=c(1,2),lty=1,cex=2,lwd=2)
#'
#' #RMSE de ajuste
#' RMSEFH7 <- FH7$RMSE
#' RMSEFH7
#'
#' #pronósticos puntuales y por IPs del 95%
#' pronFH7 <- FH7$forecast
#' pronFH7
#' #Precisión pronósticos puntuales
#' forecast::accuracy(FH7$forecast[,1],ytnuevo)
#' #Precisión pronósticos por IPs
#' precisintervals(real=ytnuevo,LIP=FH7$forecast[,2],LSP=FH7$forecast[,3])
#'
#' #Gráficos de residuos de ajuste
#' diagresTSModel(modelo=FH7,RMSE=RMSEFH7,all=FALSE,single=TRUE)
Filtrolineal <- function(serie,pesos.wi,tipo,h=1,level=0.95){
nivel <- stats::filter(serie,filter=pesos.wi,method="convolution",sides=tipo,circular=TRUE)
s <- stats::frequency(serie)
inicio <- stats::start(nivel)
ythat <- stats::ts(append(NA,nivel[-length(nivel)]),frequency=s,start=inicio)
residuos <- serie-ythat
RMSE <- forecast::accuracy(ythat,serie)[2]
tpron <- (length(serie)+1):(length(serie)+h)
if(stats::end(serie)[2]<s){
start.pron <- c(stats::end(serie)[1],stats::end(serie)[2]+1)
}
if(stats::end(serie)[2]==s){
start.pron <- c(stats::end(serie)[1]+1,1)
}
ytpron <- stats::ts(rep(nivel[length(nivel)],h),frequency=s,start=start.pron)
df <- length(serie)-2
lwr <- ytpron-stats::qt((1-level)/2,df=df,lower.tail=FALSE)*RMSE*sqrt(1+sum(pesos.wi)^2)
upr <- ytpron+stats::qt((1-level)/2,df=df,lower.tail=FALSE)*RMSE*sqrt(1+sum(pesos.wi)^2)
tablapron <- stats::ts(data.frame("Point Forecast"=ytpron,lwr,upr),frequency=s,start=start.pron)
result <- list(nivel=nivel,fitted=ythat,residuals=residuos,forecast=tablapron,RMSE=RMSE)
result
}

