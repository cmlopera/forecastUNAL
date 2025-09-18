#' Función `LoessOptimo()`
#'
#' Para ajuste y pronóstico de una serie no estacional mediante loess óptimo según citerios de información.
#' @param serie Un objeto tipo `ts` de serie de tiempo univariada con los valores de la serie a ajustar.
#' @param grado Grado del polinomio local a ser ajustado sobre la serie. `grado=1` es para loess lineal y `grado=2` es para loess cuadrático. Su valor por defecto es 1.
#' @param criterio Tipo de criterio de información `"aicc"` ó `"gcv"`, a usar para la elección del parámetro de suavizamiento en la rutina de loess óptimo mediante la función R `loess.as()` de la librería `fANCOVA`. Por defecto es `"aicc"`.
#' @param h Número de predicciones a calcular después del ajuste. Por defecto es 1.
#' @param level Nivel de confianza de los intervalos de predicción. Por defecto 0.95.
#' @returns Una lista con los siguientes objetos:
#' * `alfa.optim` valor del parámetro de suavizamiento loess óptimo.
#' * `nep` El número de parámetros equivalentes loess.
#' * `MSE` Estimación de la varianza del error de ajuste.
#' * `fitted` Objeto `ts` de la serie estimada.
#' * `residuals` Objeto `ts` con los residuos de ajuste obtenidos como la diferencia entre los valores observados y los ajustados para la serie.
#' * `forecast` Objeto `ts` multivariado con las predicciones puntuales y por intervalo para los ha periodos pedidos.
#' @export
#' @examples
#' serie <- ts(DatosChippy$serie,freq=12,start=c(1990,1))
#' curva <- ts(DatosChippy$Chippy,freq=12,start=c(1990,1))
#'
#' #Gráfica de la serie y tendencia verdadera
#' plot(serie,lwd=1)
#' lines(curva,col=2,lwd=2)
#' legend("topright",legend=c(expression(sin(cos(t)*exp(-t/2)))),col=2,lty=1,lwd=2)
#'
#' #Ajuste con validación cruzada con últimos 4 datos
#' h <- 4
#' n <- length(serie)-h
#' t <- 1:n #tiempos de ajuste
#' #valores a ajustar
#' yt <- ts(serie[t],freq=frequency(serie),start=start(serie))
#' tnuevo <- (n+1):length(serie) #tiempos de predicciones
#' #Fecha inicio predicciones
#' startpred=end(ts(serie[1:(n+1)],freq=frequency(serie),start=start(serie)))
#' #Valores observados en tiempos de predicciones
#' ytnuevo <- ts(serie[tnuevo],freq=frequency(serie),start=startpred)
#'
#' ajusteLoess <- LoessOptimo(serie=yt,grado=2,criterio="aicc",h=h,level=0.95)
#' ythat <- fitted(ajusteLoess) #serie ajustada
#'
#' #Serie y su ajuste
#' plot(serie)
#' lines(curva,col=2,lwd=2) #Tendencia verdadera
#' lines(ythat,col="blue1",lwd=2)
#' legend("topright",legend=c("Tendencia verdadera","Ajuste"),
#' col=c(2,"blue1"),lty=1,lwd=2)
#'
#' #Número aproximado de parámetros
#' p <- round(ajusteLoess$nep)
#' p
#' #AIC y BIC aproximados version exp(C*n(p))
#' Criterios <- critinfresid(residuales=residuals(ajusteLoess),npar=p)
#' Criterios
#'
#' #Pronósticos con IPs del 95% de confianza
#' pronosticos <- ajusteLoess$forecast
#' pronosticos
#'
#' #MSE aproximado del ajuste
#' MSE <- ajusteLoess$MSE
#'
#' #Gráfico serie de residuos de ajuste
#' plot(residuals(ajusteLoess))
#' abline(h=c(-2*sqrt(MSE),0,2*sqrt(MSE)),col=2)
LoessOptimo <- function(serie,grado=1,criterio,h=1,level=0.95){
  ta <- 1:length(serie)
  s <- stats::frequency(serie)
  n <- length(serie)
  ajuste <- fANCOVA::loess.as(ta,serie,degree=grado,criterion=criterio,family="gaussian",plot=F)
  print(summary(ajuste))
  alpha <- ajuste$pars$span
  ythat <- stats::ts(stats::fitted(ajuste),freq=s,start=stats::start(serie))
  MSE <- sum(stats::residuals(ajuste)^2)/(n-round(ajuste$enp))
  residuos <- stats::ts(stats::residuals(ajuste),freq=s,start=stats::start(serie))
  tpron <- (length(serie)+1):(length(serie)+h)
  if(stats::end(serie)[2]<s){
    start.pron <- c(stats::end(serie)[1],stats::end(serie)[2]+1)
  }
  if(stats::end(serie)[2]==s){
    start.pron <- c(stats::end(serie)[1]+1,1)
  }
  aux1p <- stats::loess(serie~ta,degree=grado,span=alpha,family="gaussian",control=stats::loess.control(surface="direct"))
  aux2p <- stats::predict(aux1p,data.frame(ta=tpron),se=T)
  ytpron <- aux2p$fit
  LI <- ytpron-stats::qt((1-level)/2,df=aux2p$df,lower.tail=F)*aux2p$se.fit
  LP <- ytpron+stats::qt((1-level)/2,df=aux2p$df,lower.tail=F)*aux2p$se.fit
  tablapron <- stats::ts(data.frame("Point forecast"=ytpron,lwr=LI,upr=LP),freq=s,start=start.pron)
  resul <- list(alfa.optim=alpha,nep=ajuste$enp,MSE=MSE,fitted=ythat,residuals=residuos,forecast=tablapron)
  resul
}
