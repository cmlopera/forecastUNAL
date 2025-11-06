#' Función `BPLBtest()`
#'
#' Realiza pruebas Ljung-Box y Box-Pierce para `m=6, 12, 18,...,[maxlag/6]`,  con  `maxlag` el máximo orden de rezagos. Usa la función `Box.test()` de la librería `stats`.
#'
#' @param serie Vector de valores del proceso sobre el que se quiere evaluar incorrelación con tests Ljung-Box o Box-Pierce.
#' @param maxlag Entero indicando el máximo orden de rezago temporal hasta el cual se prueba la significancia de la función de autocorrelación.
#' @param type Cadena de caracteres para indicar el tipo de prueba, `"Ljung"` (por defecto) para las pruebas Ljung-Box y `"Box"` para las pruebas Box-Pierce.
#' @return Un objeto `data.frame` con tres variables:
#' * `Xsquared` los valores de los estadísticos de prueba para `m=6, 12,...,[maxlag/6]`.
#' * `df` los grados de libertad del estadístico de prueba, corresponden a los valores de `m=6, 12,...,[maxlag/6]`.
#' * `pvalue` los valores P para las pruebas con `m=6, 12,...,[maxlag/6]`.
#' @export
#' @examples
#' #Simulando serie de un proceso estacionario de media cero con variables correlacionadas
#' zt <- ts(arima.sim(n=300,list(ar=0.7),sd=4000),freq=1,start=1)
#' #Aplicando pruebas Ljung-Box con máximo orden de correlación de 36
#' BPLBtest(serie=zt,maxlag=36)
#' #Aplicando pruebas Box-Pierce con máximo orden de correlación de 36
#' BPLBtest(serie=zt,maxlag=36,type="Box")
#' #Simulando serie de un proceso ruido blanco gaussiano
#' zt <- ts(rnorm(n=300,sd=4000),freq=1,start=1)
#' #Aplicando pruebas Ljung-Box con máximo orden de correlación de 36
#' BPLBtest(serie=zt,maxlag=36)
#' #Aplicando pruebas Box-Pierce con máximo orden de correlación de 36
#' BPLBtest(serie=zt,maxlag=36,type="Box")
BPLBtest <- function(serie,maxlag,type=c("Ljung","Box")){
type <- match.arg(type)
aux <- floor(maxlag/6)
Xsquared <- c(rep(NA,aux))
df=c(rep(NA,aux))
pvalue <- c(rep(NA,aux))
for(i in 1:aux){
test <- stats::Box.test(serie,lag=(6*i),type=type)
Xsquared[i] <- test[[1]]
df[i] <- test[[2]]
pvalue[i] <- test[[3]]
}
lag <- 6*c(1:aux)
teste <- as.data.frame(cbind(Xsquared,df,pvalue))
rownames(teste) <- lag
teste
}
