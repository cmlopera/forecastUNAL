#' `Función interpdeltas()`
#'
#' Extrae las estimaciones y construye gráfico de los efectos estacionales en modelos de regresion global estacionales que usan indicadoras con nivel de referencia el `s` ésimo período calendario. Particularmente en modelos polinomiales estacionales con error RB y con error ARMA, modelos log polinomiales estacionales con error RB y con error ARMA y modelos exponenciales polinomiales estacionales con error RB. En los modelos polinomiales estacionales bien sea con error RB o ARMA, extrae los valores estimados de los parámetros asociados a las variables indicadoras y crea la gráfica de estos valores vs. periodo calendario, asignando valor de cero en el último período. En los modelos log polinomiales estacionales con error RB o ARMA y en los modelos exponenciales polinomiales estacionales con error RB, extrae los valores estimados de los parámetros asociados a las variables indicadoras, pero es necesario informar que estos modelos son multiplicativos para que sean exponenciadas estas estimaciones, las multiplique por 100% y las grafique contra el período calendario, asignando el valor de 100% al último periodo del año.
#'
#' @param modelo El nombre del objeto R con el modelo ajustado. Objetos tipo `lm()`, `nls()` o `Arima()`.
#' @param ordenp,ordenq,ordenP,ordenQ Argumentos de valor entero, para especificar los valores de los ordenes `p`, `q`, `P`, `Q` de modelos `ARMA(p,q)(P,Q)[s]`. Por defecto estan fijos en 0, es decir, asumiendo error RB. Si el error es un `AR(p)` el usuario debe especificar el valor de `p`; si el error es `MA(q)` debe especificar el valor de `q`; si es `ARMA(p,q)` debe especificar `p` y `q`, etc.
#' @param gradopoly El grado del polinomio global. Siempre debe ser especificado.
#' @param aditivo Argumento logico `(TRUE, FALSE)`, por defecto es `TRUE`, este valor debe usarse en modelos polinomiales estacionales con indicadoras con error RB o error ARMA. En modelos log polinomiales estacionales con error RB o error ARMA, o modelos exponenciales polinomiales estacionales con error RB, este argumento debe especificarse igual a `FALSE`.
#' @param plotit Argumento logico `(TRUE, FALSE)`, por defecto es `TRUE` indicando que se desea la grafica de los efectos estacionales para su interpretación en la escala de la serie, según si el modelo es aditivo o multiplicativo.
#'
#' @return Con `plotit=TRUE`, despliega la gráfica correspondiente y en la consola R el vector de parámetros estacionales estimados (modelos aditivos) o de los valores exponenciados de estas estimaciones y multiplicados por 100% (modelos multiplicativos). Además, crea un objeto lista con los siguientes componentes:
#' * `periodo` Un vector de valores enteros de 1 a `s` (la longitud deñ año calendario).
#' * `deltasi` Un vector con las estimaciones de los parámetros asociados a las indicadoras y adiciona el valor de cero al final. Este vector aparece cuando `aditivo=TRUE`.
#' * `expdeltasi100` Un vector con los valores exponenciados de las estimaciones de los parametros asociados a las indicadoras, multiplicados por 100 y adiciona el valor de 100 al final. Este vector aparece cuando `aditivo=FALSE`.
#' @export
#' @examples
#' library(forecast)
#' #Ej. 1 Modelo polinomial estacional con indicadoras (nivel ref=12)
#' datos1 <- salario
#' n1 <- length(datos1)-12
#' t1 <- 1:n1
#' yt1 <- ts(datos1[t1],freq=12,start=start(datos1))
#' poli1 <- Mipoly(tiempo=t1,grado=4)
#' mes <- seasonaldummy(yt1)
#' X1 <- data.frame(poli1,mes)
#' modelo1 <- lm(yt1~.,data=X1)
#' summary(modelo1)
#' #Extrayendo y graficando las estimaciones de los coeficientes de regresión delta_i
#' interpdeltas(modelo1,gradopoly=4,aditivo=TRUE,plotit=TRUE)
#'
#' #Ej. 2 log cúbico estacional con indicadoras (nivel ref=Q4), error ARMA(1,2)x(0,1)[4]
#' datos2 <- ts(CEMENTO,fre=4,start=c(1956,1))
#' n2 <- length(datos2)-4
#' t2 <- 1:n2
#' yt2 <- ts(datos2[t2],freq=4,start=start(datos2))
#' poli2 <- Mipoly(tiempo=t2,grado=3)
#' trimestre <- seasonaldummy(yt2)
#' X2 <- data.frame(poli2,trimestre)
#' modelo2 <- Arima(log(yt2),order=c(1,0,2),
#' seasonal=list(order=c(0,0,1)),xreg=as.matrix(X2),method="ML")
#' k2 <- length(coef(modelo2)[coef(modelo2)!=0]);k2 #número de parámetros
#' df2 <- n2-k2 #grados de libertad
#' lmtest::coeftest(modelo2,df=df2)
#' interpdeltas(modelo=modelo2,ordenp=1,ordenq=2,ordenP=0,ordenQ=1,
#'              gradopoly=3,aditivo=FALSE,plotit=TRUE)
#'
#' #Ej. 3, exponencial cúbico estacional con indicadoras (nivel ref=Q4)
#' modelo3 <- regexponencialv02(respuesta=yt2,data=X2)
#' summary(modelo3)
#' interpdeltas(modelo=modelo3,gradopoly=3,aditivo=FALSE,plotit=TRUE)
interpdeltas <- function(modelo,ordenp=0L,ordenq=0L,ordenP=0L,ordenQ=0L,gradopoly,aditivo=TRUE,plotit=TRUE){
clasemod <- class(modelo)
stopifnot("`modelo` debe ser un objeto clase `lm`,`nls`, `forecast_ARIMA`, `ARIMA` o `Arima` " =clasemod%in%c("lm","nls","forecast_ARIMA","ARIMA","Arima"))
aux <- ordenp+ordenq+ordenP+ordenQ+gradopoly+1
if(aditivo==TRUE){
deltas <- c(stats::coef(modelo)[-c(1:aux)],0)
names(deltas)=paste0("delta",1:length(deltas))
print(deltas[1:(length(deltas)-1)])
if(plotit==TRUE){
plot(1:length(deltas),deltas,type="b",pch=1,lwd=3,xlab="periodo calendario",ylab=expression(widehat(delta)[i]),xaxt="n")
  graphics::axis(1,at=1:length(deltas),1:length(deltas))
}
periodo <- 1:length(deltas)
res <- list(periodo=periodo,deltasi=deltas)
res
}
else{
if(aditivo==FALSE){
expdeltas100 <- exp(c(stats::coef(modelo)[-c(1:aux)],0))*100
names(expdeltas100) <- paste0("exp","(delta",1:length(expdeltas100),")*100%")
print(expdeltas100[1:(length(expdeltas100)-1)])
if(plotit==TRUE){
plot(1:length(expdeltas100),expdeltas100,type="b",pch=1,lwd=3,xlab="periodo calendario",ylab=expression(exp(widehat(delta)[i])*100),xaxt="n")
graphics::axis(1,at=1:length(expdeltas100),labels=1:length(expdeltas100))
}
periodo <- 1:length(expdeltas100)
res <- list(periodo=periodo,expdeltasi100=expdeltas100)
res
}
}
}
