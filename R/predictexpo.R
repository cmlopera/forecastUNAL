#' Función `predictexpo()`
#'
#' Para la predicción puntual y por IPs sobre modelos de regresión exponencial ajustados con la función `regexponencialv02()`. También calcula I.Cs para la respuesta media en lugar de IPs.
#'
#' @param modelo Objeto tipo `nls()` generado por la función `regexponencialv02()`.
#' @param new.data `data.frame` con los valores de los predictores en los tiempos de predicción. El nombre de las variables debe ser igual al usado en el ajuste del modelo.
#' @param level Un valor en (0,1) para el nivel de confianza de los intervalos. Por defecto es 0.95.
#' @param interval Una cadena de caracteres correspondiendo a uno de los siguientes casos: `"none"`, `"confidence"`, `"prediction"`. Por defecto es `"none"`, es decir, calcula solo estimaciones o pronósticos puntuales. Para intervalos de confianza en la respuesta media use `"confidence"` y para intervalos de predicción use `"prediction"`.
#'
#' @return Un objeto tipo `data.frame` con tres variables:
#' * `forecast` La estimación o pronóstico puntual de la respuesta.
#' * `lower`, `upper` Los límites inferior y superior de los intervalos pedidos.
#' @export
#' @examples
#' #Simulando serie aditiva con modelo de tendencia exponencial lineal y errores correlacionados
#' n <- 250
#' t <- 1:n
#' error <- arima.sim(n=n,list(ar=0.8),sd=108)
#' Yt <- ts(exp(6.5+0.0105*t)+error,freq=1,start=1)
#' plot.ts(Yt)
#' lines(t,exp(6.5+0.0105*t),type="l",col="red2")
#'
#' #Matriz de predictores
#' X <- Mipoly(tiempo=t,grado=1)
#' #Ajuste del modelo exponencial lineal
#' modelo <- regexponencialv02(respuesta=Yt,data=X)
#' summary(modelo)
#' #cinco tiempos de predicción
#' Xnuevo=Mipoly(tiempo=(n+1):(n+5),grado=1)
#' #Predicciones puntuales y por IPs del 95% de confianza
#' predictexpo(modelo=modelo,new.data=Xnuevo,interval="prediction")
predictexpo <- function(modelo,new.data,level=0.95,interval = c("none", "confidence", "prediction")){
stopifnot("`modelo` debe ser un objeto clase `nls` " =class(modelo)=="nls",
          "`interval` debe ser uno de los siguientes argumentos `none`, `confidence`, `prediction`"=interval%in%c("none", "confidence", "prediction"),
          "`level` debe ser valor en (0,1)"=level>0 & level<1)
#aprovechando normalidad del vector de parametros estimados y considerando g(beta)=exp(t(x0[i])%*%beta) con x0[i]
#el vector de predictores en i esima prediccion, incluye al 1 como primera entrada
x0 <- cbind(1,as.matrix(new.data)) #matriz con filas siendo los puntos de prediccion o valores de variables explicatorias en la prediccion
mse <- summary(modelo)$sigma^2 #estimacion varianza del error del modelo de regresion exponencial
yhat <- stats::predict(modelo,newdata=new.data)
varbetas <- as.matrix(stats::vcov(modelo)) #matriz de varianzas covarianzas del vector de parametros estimados modelo exponencial
Sigmapred <- (yhat^2)*diag(x0%*%varbetas%*%t(x0)) #aproximacion por metodo delta de la varianza de la respuesta predicha o estimada
if(interval=="none"){
res <- stats::predict(modelo,newdata=new.data)
}
if(interval=="confidence"){
#ICs para la respuesta media con metodo delta en la varianza de la estimacion puntual de la respuesta media
LI <- yhat-stats::qt((1-level)/2,df=summary(modelo)$df[2],lower.tail=F)*sqrt(Sigmapred)
LS <- yhat+stats::qt((1-level)/2,df=summary(modelo)$df[2],lower.tail=F)*sqrt(Sigmapred)
res <- data.frame(estimate=yhat,lower=LI,upper=LS)
}
if(interval=="prediction"){
#IPs con metodo delta en la varianza de la prediccion puntual de la respuesta
LI <- yhat-stats::qt((1-level)/2,df=summary(modelo)$df[2],lower.tail=F)*sqrt(mse+Sigmapred)
LS <- yhat+stats::qt((1-level)/2,df=summary(modelo)$df[2],lower.tail=F)*sqrt(mse+Sigmapred)
res <- data.frame(forecast=yhat,lower=LI,upper=LS)
}
res
}

