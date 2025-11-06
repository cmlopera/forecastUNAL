#' Funcion `estimarecursiva()`
#'
#' Realiza la estimación recursiva de un modelo de regresión lineal, obteniendo además los residuos recursivos, la gráfica `CUSUMt-Recursivo` y el test `CUSUMt-Recursivo` mediante las funciones `recresid()`, `efp()` y `sctest()`, respectivamente, de la librería `strucchange`.
#'
#' @param respuesta Un vector numérico o serie de tiempo con los valores de la respuesta del modelo de regresión lineal múltiple.
#' @param dataX Un `data.frame` cuyas columnas son los predictores del modelo de regresión lineal múltiple, el número de filas debe ser igual a la longitud de la variable `respuesta`.
#' @param min.n Un escalar para indicar el mínimo tamaño de muestra con el que debe iniciar la estimación recursiva, el cual no puede ser inferior a `ncol(dataX) + 2`.
#' @param nrow1,ncol1 Valores enteros que definen el número de `filas, columnas` en el diseño de la ventana gráfica para representar las gráficas de los residuales recursivos y del estadístico `CUSUMt-Recursivo`. Si ambos valores son iguales a 1, cada gráfica sale en una ventana por separado.
#' @param nrow2,ncol2 Valores enteros que definen el número de `filas, columnas` en el diseño de la ventana gráfica para representar las gráficas de las estimaciones recursivas de los parámetros cuando el argumento `plot.recursive = TRUE`. Cuando el número de gráficas a ser representadas  supera al número de celdas definido en el diseño de la ventana de salida, entonces se generarán las ventanas gráficas adicionales que sean necesarias con el mismo diseño.
#' @param plot.recursive Argumento lógico, sus valores posibles son `TRUE` (por defecto) y `FALSE`, para indicar si se desean o no las gráficas de las estimaciones recursivas de los parámetros del modelo de regresión lineal múltiple.
#'
#' @return La función muestra en la consola `R` la estimación del modelo global sobre el total de observaciones leídas en el vector `respuesta` y los resultados del test `CUSUMt-recursivo`. Genera también el gráfico de los residuos recursivos y el gráfico del estadístico `CUSUMt-recursivo`. Además, produce un objeto tipo lista con los siguientes componentes:
#' * `n` Una matriz de una sola columna con sus filas siendo los valores de tamaño de muestra desde `min.n` hasta `length(respuesta)`.
#' * `estimacion_recursiva` Un arreglo de matrices, donde cada matriz tiene tres columnas: la estimación de los parámetros y los límites inferior y superior de confianza del 95% de cada parámetro. El número de matrices es igual al total de tamaños de muestras entre `min.n` y `length(respuesta)`.
#' * `ajusteglobal` Un objeto tipo `lm()` del ajuste de la `respuesta` vs. los predictores en `dataX`.
#' * `resid_recursivos` Un vector de residuos recursivos.
#' * `test_CUSUM` Los resultados (estadístico de prueba y valor-P) del test `CUSUMt-recursivo`.
#' @export
#' @examples
#' library(forecast)
#' datos <- ts(CEMENTO,frequency=4,start=c(1956,1)) # serie CEMENTO con todos los datos
#' #Variables para modelo log cuadrático estacional con indicadoras, nivel referencia trimestre 4
#' t <- 1:length(datos) #Definiendo índice de tiempo para todos los datos
#' poli <- Mipoly(tiempo=t,grado=2)
#' trimestre <- seasonaldummy(datos) #Creando variables indicadoras, nivel de  referencia trimestre 4
#' X <- data.frame(poli,trimestre) #Matriz de predictores
#' resultados <- estimarecursiva(respuesta=log(datos),dataX=X,min.n=24,
#'               nrow1=2,ncol1=1,nrow2=2,ncol2=2,plot.recursive=TRUE)
estimarecursiva <- function(respuesta,dataX,min.n,nrow1=2,ncol1=1,nrow2=2,ncol2=2,plot.recursive=TRUE){
names.vars <- names(dataX)
names.param <- c("Intercept", paste0("param_", names.vars))
p <- ncol(dataX)+1
m <- matrix(min.n:length(respuesta),ncol=1)
max.n <- length(respuesta)
posix <- (min.n+max.n)/2
ajusteglobal <- stats::lm(respuesta~.,data=dataX)
rr <- strucchange::recresid(ajusteglobal)
RRCUSUM <- strucchange::efp(respuesta~.,data=dataX, type = "Rec-CUSUM")
bound <- as.numeric(strucchange::boundary(RRCUSUM, alpha = 0.05))
datares <- data.frame(tiempo=seq((p+1):max.n),recursive_res=rr)
datarescusum <- data.frame(tiempo=as.numeric(stats::time(RRCUSUM$process)),EPF=as.numeric(RRCUSUM$process),LCL=-bound,UCL=bound)
cat("Ajuste Global")
cat("\n")
print(summary(ajusteglobal))
estim <- function(n){
datan <- data.frame(dataX[1:n,])
names(datan) <- names.vars
modelo <- stats::lm(respuesta[1:n]~.,data=datan)
resul <- cbind(stats::coef(modelo),stats::confint(modelo))
}
test <- strucchange::sctest(respuesta~.,data=dataX)
cat("\n")
cat("Resultados Test Cusum Recursivo")
cat("\n")
print(test)

g1 <- ggplot2::ggplot(datares,ggplot2::aes(x=tiempo,y=recursive_res))+ggplot2::geom_line()+ggplot2::geom_hline(yintercept=0,color='red')
g2 <- ggplot2::ggplot(datarescusum,ggplot2::aes(x=tiempo,y=EPF))+ggplot2::geom_line()+
  ggplot2::geom_line(ggplot2::aes(y=LCL,color='red'))+
  ggplot2::geom_line(ggplot2::aes(y=UCL,color='red'))+
  ggplot2::geom_hline(yintercept=0,color='black')+
  ggplot2::theme(legend.position = "none")+
  ggplot2::ylab("Empirical fluctuation process")

print(gridExtra::marrangeGrob(list(g1,g2),nrow = nrow1,ncol= ncol1))

b <- array(apply(m,1,estim),dim=c(p,3,nrow(m)),dimnames=list(names.param,c("Estimate","LCL","UCL"),paste0("n=",m)))
if(plot.recursive){
myplots <- list()
for(i in 1:p){
pari <- t(b[i,,])
href <- stats::coef(ajusteglobal)[i]
grafi <- eval(substitute(ggplot2::ggplot(pari,ggplot2::aes(x=min.n:max.n))+
                           ggplot2::geom_line(ggplot2::aes(y=Estimate,color='EstimRec'))+
                           ggplot2::geom_line(ggplot2::aes(y=LCL,color='ICs'))+
                           ggplot2::geom_line(ggplot2::aes(y=UCL,color='ICs'))+
                           ggplot2::geom_hline(yintercept=href,color='red')+
                           ggplot2::annotate("text", x=posix, y=href,size = 2.3,label= "Estimacion global")+
                           ggplot2::scale_color_manual(name="",
                     breaks=c('EstimRec','ICs'),
                     values=c('EstimRec'='black','ICs'='blue'))+
                       ggplot2::theme(legend.position = "top")+
                       ggplot2::xlab("n")+
                       ggplot2::ylab(names.param[i]),list(i = i)))
myplots[[i]] <- grafi
}

print(gridExtra::marrangeGrob(myplots, nrow = nrow2, ncol = ncol2))
}
result <- list(n=m,estimacion_recursiva=b,ajusteglobal=ajusteglobal,resid_recursivos=rr,test_CUSUM=test)
}
