#' Función `DescompLoessv02()`
#'
#' Estima y predice serie estacional aditiva y multiplicativa combinando filtro de la descomposición clásica y la regresión loess.
#' @param serie.ajuste Un objeto `ts` de serie de tiempo univariada con los valores de la serie a ajustar.
#' @param tipo.descomp Tipo de descomposición `"additive"` ó `"multiplicative"`, en el último caso es bajo modelo parcialmente multiplicativo.
#' @param grado Grado del polinomio local a ser ajustado sobre la serie previamente desestacionalizada con la estacionalidad estimada con el filtro de la descomposición clásica de la función R `decompose()`, bajo el modelo aditivo o multiplicativo, según valor del argumento `tipo.descomp`. Su valor por defecto es 1.
#' @param criterio Tipo de criterio de información `"aicc"` ó `"gcv"`, a usar para la elección del parámetro de suavizamiento en la rutina de loess óptimo mediante la función R `loess.as()` de la librería `fANCOVA`. Por defecto es `"aicc"`.
#' @param h Número de predicciones a calcular después del ajuste. Por defecto es 1.
#' @param level Nivel de confianza de los intervalos de predicción. Por defecto 0.95.
#' @return Una lista con los siguientes objetos:
#' * `deltasi` Vector con las estimaciones de los efectos estacionales.
#' * `alfa.optim` Valor del parámetro de suavizamiento loess óptimo.
#' * `nep` El número de parámetros equivalentes loess.
#' * `p` El total de parámetros aproximados en el ajuste de tendencia y estacionalidad.
#' * `St` Objeto `ts` de la estimación de la componente estacional.
#' * `Tt` Objeto `ts` de la estimación de la tendencia loess.
#' * `ytd` La serie desestacionalizada.
#' * `fitted` Objeto `ts` de la serie estimada combinando aditiva o multiplicativamente, según el tipo de descomposición, las estimaciones de la tendencia y de la estacionalidad.
#' * `forecast` Objeto `ts` multivariado con las predicciones puntuales y por intervalo para los `h` períodos pedidos.
#' * `residuals` Objeto `ts` con los residuos de ajuste obtenidos como la diferencia entre los valores observados y los ajustados para la serie.
#' * `MSE` Estimación de la varianza del error de ajuste.
#' @export
#' @examples
#' #Serie de tiempo con inicio en enero de 1965
#' ventascompanyX
#'
#' m <- 6 #tamaño muestra de validación cruzada
#' n <- length(ventascompanyX)-m #tamaño muestra de ajuste
#' t <- 1:n #índice de tiempo en periodos de ajuste
#'
#' #Serie de tiempo en periodos de ajuste
#' yt <- ts(ventascompanyX[t],frequency=12,start=c(1965,1))
#'
#' #Ajuste por descomposición multiplicativa y tendencia loess lineal
#' #se usa criterio AICC para el parámetro de suavizamiento loess
#' ajusteDLL1 <- DescompLoessv02(serie.ajuste=yt,tipo.descomp="multiplicative",
#'                               grado=1,criterio="aicc",h=m,level=0.95)
#'
#' #Gráfico de la serie y su ajuste
#' plot(ventascompanyX)
#' lines(fitted(ajusteDLL1),col=2,lwd=2)
#' legend("topleft",legend=c("Serie","Ajuste"),col=c(1,2),lty=1)
#'
#' #Pronósticos puntuales y por IPs en los períodos de validación cruzada
#' ajusteDLL1$forecast
#'
#' #Estimación de sigma
#' sigma <- sqrt(ajusteDLL1$MSE)
#'
#' #residuos del ajuste
#' diagresTSModel(modelo=ajusteDLL1,RMSE=sigma,all=FALSE,single=FALSE)
DescompLoessv02 <- function(serie.ajuste,tipo.descomp,grado=1,criterio,h=1,level=0.95){
s <- stats::tsp(serie.ajuste)[3]
tajuste <- 1:length(serie.ajuste)
tpron <- (length(serie.ajuste)+1):(length(serie.ajuste)+h)
if(stats::end(serie.ajuste)[2]<s){
start.pron <- c(stats::end(serie.ajuste)[1],stats::end(serie.ajuste)[2]+1)
}
if(stats::end(serie.ajuste)[2]==s){
start.pron <- c(stats::end(serie.ajuste)[1]+1,1)
}

descom <- stats::decompose(serie.ajuste,type=tipo.descomp)
St <- descom$seasonal
estacionini <- stats::cycle(serie.ajuste)[1]
if(estacionini==1){
deltas_i <- descom$figure
}
if(estacionini!=1){
j <- estacionini;deltas_i <- c(descom$figure[(s-j+2):s],descom$figure[1:(s-j+1)])
}
efectos <- deltas_i
if(s==12){
names(efectos) <- c("s1","s2","s3","s4","s5","s6","s7","s8","s9","s10","s11","s12")
}
if(s==4){
names(efectos) <- c("s1","s2","s3","s4")
}
efectos <- data.frame(efectos)
names(efectos) <- ""
cat("Efectos estacionales estimados")
print(efectos)

indi <- forecast::seasonaldummy(serie.ajuste,h)
indi2 <- cbind(indi,1-apply(indi,1,sum))
Stnuevo <- stats::ts(indi2%*%deltas_i,frequency=s,start=start.pron)

if(tipo.descomp=="additive"){
ytd <- serie.ajuste-St
}
if(tipo.descomp=="multiplicative"){
ytd <- serie.ajuste/St
}
ajusteLoess1 <- fANCOVA::loess.as(tajuste,ytd,degree=grado,criterion=criterio,family="gaussian",plot=F)
cat("\n")
cat("Resumen loess sobre serie desestacionalizada:")
cat("\n")
print(summary(ajusteLoess1))
alfa.optim1 <- ajusteLoess1$pars$span
nep1 <- round(ajusteLoess1$enp)
p1 <- nep1+s-1 #numero aproximado de parametros en ajuste Modelo
Tt1 <- stats::ts(stats::fitted(ajusteLoess1),frequency=s,start=stats::start(serie.ajuste))
loess <- stats::predict(stats::loess(ytd~tajuste,span=alfa.optim1,degree=grado,control=stats::loess.control(surface="direct")),data.frame(tajuste=tpron),se=TRUE)
Ttnuevo1 <- stats::ts(loess$fit,freq=s,start=start.pron)
if(tipo.descomp=="additive"){
ythat1 <- St+Tt1
ytpron1 <- Ttnuevo1+Stnuevo
et1 <- serie.ajuste-ythat1
df1 <- length(serie.ajuste)-nep1+s-1
MSE1 <- sum(et1^2)/df1
varpred <- MSE1+(loess$se.fit)^2
LIP <- ytpron1-stats::qt((1-level)/2,df=df1,lower.tail=F)*sqrt(varpred)
LSP <- ytpron1+stats::qt((1-level)/2,df=df1,lower.tail=F)*sqrt(varpred)
}
if(tipo.descomp=="multiplicative"){
ythat1 <- St*Tt1
ytpron1 <- Ttnuevo1*Stnuevo
et1 <- serie.ajuste-ythat1
df1 <- length(serie.ajuste)-nep1+s-1
MSE1 <- sum(et1^2)/df1
varpred <- MSE1+(Stnuevo^2)*(loess$se.fit)^2
LIP <- ytpron1-stats::qt((1-level)/2,df=df1,lower.tail=F)*sqrt(varpred)
LSP <- ytpron1+stats::qt((1-level)/2,df=df1,lower.tail=F)*sqrt(varpred)
}
tablapron1 <- stats::ts(cbind(forecast=ytpron1,lower=LIP,upper=LSP),freq=s,start=start.pron)
result <- list(deltasi=deltas_i,alfa.optim=alfa.optim1,nep=nep1,p=p1,St=St,Tt=Tt1,ytd=ytd,fitted=ythat1,forecast=tablapron1,residuals=et1,MSE=MSE1)
result
}
