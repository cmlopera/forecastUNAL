#' Función `critinfresid()`
#'
#' Calcula criterios de información como en Diebold(2001).
#'
#' @param residuales Un vector de diferencias entre las observaciones y su ajuste, en la escala original de los datos.
#' @param npar Un valor entero igual al número de parámetros del modelo ajustado.
#'
#' @return Un vector de longitud tres con `p` (el número de parámetros) `AIC` (criterio Akaike) y `BIC` (criterio Bayesiano)
#' @export
#'
#' @examples
#' #Simulando serie aditiva con modelo de tendencia lineal y errores iid normales
#' t <- 1:300
#' yt <- 3500+200*t+rnorm(300,sd=4000)
#' yt <- ts(yt,freq=12,start=c(1980,1))
#' #Ajuste modelo de tendencia lineal
#' mrl <- lm(yt~t)
#' summary(mrl)
#' ythat <- ts(fitted(mrl),freq=12,start=c(1980,1))
#' #Gráfico serie observada y su ajuste
#' plot(yt)
#' lines(ythat,col=2)
#' #Calculando AIC y BIC sobre modelo ajustado
#' p <- length(coef(mrl)) #Número de parámetros
#' critinfresid(residuales=residuals(mrl),npar=p)
#'
#' #Simulando serie multiplicativa con modelo log lineal y errores correlacionados de media cero
#' t <- 1:540
#' logYt <- 2.5+0.0105*t+arima.sim(n=length(t),list(ar=0.6),sd=0.27)
#' Yt <- ts(exp(logYt),freq=1,start=1)
#' #Ajuste modelo log lineal con errores iid normales
#' modelo <- lm(log(Yt)~t)
#' summary(modelo)
#' #Valores estimados de la serie en escala original de los datos
#' Ythat <- ts(exp(fitted(modelo))*exp(summary(modelo)$sigma^2/2),freq=1,start=start(Yt))
#' #Gráfico serie observada y su ajuste
#' plot(Yt)
#' lines(Ythat,col=2)
#' legend("topleft",legend=c("Serie","Ajuste"),col=c(1,2),lty=1)
#' #Calculando AIC y BIC sobre modelo ajustado
#' seudores <- Yt-Ythat #Seudo residuos para calcular criterios de información
#' p <- length(coef(modelo)) #Número de parámetros
#' critinfresid(residuales=seudores,npar=p)
critinfresid <- function(residuales,npar) {
  AIC <- exp(log(mean(residuales^2))+2*npar/length(residuales))
  BIC <- exp(log(mean(residuales^2))+npar*log(length(residuales))/length(residuales))
  crit <- list(p=npar,AIC=AIC,BIC=BIC)
  unlist(crit)
}
