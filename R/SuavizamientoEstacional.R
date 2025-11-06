#' Función `SuavizamientoEstacional()`
#'
#' Sobre un objeto `ts` realiza ajuste y predicción por suavizamiento exponencial Holt-Winters aditivo o multiplicativo. El suavizamiento puede ser óptimo o con valores predefinidos para los parámetros de suavizamiento. Usa la función R `HoltWinters()` pero presenta en orden de los periodos del año las últimas estimaciones suavizadas de los efectos estacionales.
#'
#' @param serie Un objeto tipo `ts` de frecuencia `s`.
#' @param alpha Parámetro de suavizamiento del nivel, debe ser un valor entre 0 y 1. Por defecto es `NULL`, en este caso, será hallado minimizando el SSE del ajuste.
#' @param beta Parámetro de suavizamiento de la pendiente, debe ser un valor entre 0 y 1. Por defecto es `NULL`, en este caso, será hallado minimizando el SSE del ajuste.
#' @param gamma Parámetro de suavizamiento de los efectos estacionales, debe ser un valor entre 0 y 1. Por defecto es `NULL`, en este caso, será hallado minimizando el SSE del ajuste.
#' @param seasonal Cadena de caracteres para indicar el tipo de descomposición, por defecto es `"additive"` para el caso aditivo. Para el caso multiplicativo se define como `"multiplicative"`.
#' @param h Número de pronósticos que se realizarán inmediatamente después del ajuste. Por defecto es 1.
#' @param optim.start Vector con componentes de nombres `alpha`, `beta` y `gamma`, de los valores iniciales de los parámetros de suavizamiento pasados al optimizador.
#' @param optim.control Lista con parámetros de control adicionales pasados a `optim`.
#' @return La función produce un objeto tipo lista con los siguientes componentes:
#' * `coefficients` Un `data.frame` con una sola columna de longitud `s+2`, con los valores finales suavizados del nivel, la pendiente y los efectos estacionales (estos ultimos en el orden `i=1, 2, ...,s`).
#' * `fitted` Objeto `ts` de los valores ajustados, con la misma frecuencia de `serie` pero iniciando en el tiempo `s+1`.
#' * `residuals` Objeto `ts` de los residuos de ajuste con la misma frecuencia de `serie`.
#' * `forecast` Objeto serie de tiempo multivariada (`mts`), con los pronósticos puntuales (columna 1), y límites de predicción del 95% de confianza (columnas 2 y 3), por defecto, para `h=1` períodos después del ajuste.
#' * `MSE` El MSE aproximado del ajuste.
#' @export
#' @examples
#' #Serie estacional trimestral multiplicativa de N=155 observaciones, inicio 1956-Q1
#' CEMENTO
#' datos <- ts(CEMENTO,freq=4,start=c(1956,1))
#' plot(datos)
#' t <- 1:151
#' yt <- ts(CEMENTO[t],freq=4,start=c(1956,1)) #serie de primeros 151 valores
#' #Fecha inicio predicciones
#' startpred <- end(ts(CEMENTO[1:152],freq=4,start=c(1956,1)))
#' #serie últimos cuatro valores
#' ytf <- ts(CEMENTO[152:length(CEMENTO)],freq=4,start=startpred)
#'
#' #Ajuste por SEHW multiplicativo óptimo con n=151 y predicción con h=4
#' modelo <- SuavizamientoEstacional(yt,seasonal="multiplicative",h=4)
#' ythat <- fitted(modelo) #serie de valores ajustados
#'
#' #Gráfica del ajuste
#' plot(datos)
#' lines(ythat,col=2)
#' legend("topleft",legend=c("Original","Ajustada SEHW"),col=c(1,2),lty=1)
#'
#' #Calculando de forma aproximada AIC y BIC usando exp(C*n(p))
#' numpar <- frequency(datos)+1 #Aprox. del número de parámetros
#' Criterios <- critinfresid(residuales=residuals(modelo),npar=numpar)
#' Criterios
#'
#' #Gráficos de residuos
#' MSE <- modelo$MSE; MSE #MSE aproximado del ajuste total del Suavizamiento
#' diagresTSModel(modelo=modelo,RMSE=sqrt(MSE),all=FALSE,single=FALSE)
#'
#' #Predicciones puntuales y por IPs del 95% de confianza
#' predicciones <- modelo$forecast
#' predicciones
#' #Precisión de las predicciones
#' forecast::accuracy(predicciones[,1],ytf)
#' precisintervals(real=ytf,LIP=predicciones[,2],LSP=predicciones[,3])
SuavizamientoEstacional <- function(serie,alpha=NULL,beta=NULL,gamma=NULL,seasonal="additive",
                                    h=1,optim.start=c(alpha = 0.3, beta = 0.1, gamma = 0.1),
                                    optim.control = list()){
suaviza <- stats::HoltWinters(serie,alpha=alpha,beta=beta,gamma=gamma,seasonal=seasonal,
                              optim.start=optim.start,optim.control=optim.control)
predicc <- stats::predict(suaviza,n.ahead=h,prediction.interval=TRUE)[,c(1,3,2)]
ythat <- stats::fitted(suaviza)[,1]
res <- stats::residuals(suaviza)
s <- stats::frequency(suaviza$x)
df <- length(serie)-2*s-((s-1)+2)
MSE <- suaviza$SSE/df

if(stats::end(suaviza$x)[2]<s){
estacionini <- stats::end(suaviza$x)[2]+1
}
if(stats::end(suaviza$x)[2]==s){
estacionini <- 1
}

if(estacionini==1){
efectossuav <- suaviza$coef[-c(1,2)]
}
if(estacionini!=1){
j <- estacionini
efectossuav <- c(suaviza$coef[-c(1,2)][(s-j+2):s],suaviza$coef[-c(1,2)][1:(s-j+1)])
}
if(s==12){
names(efectossuav) <- c("s1","s2","s3","s4","s5","s6","s7","s8","s9","s10","s11","s12")
}
if(s==4){
names(efectossuav) <- c("s1","s2","s3","s4")
}
pars <- data.frame(alpha=suaviza$alpha,beta=suaviza$beta,gamma=suaviza$gamma,row.names ="")
names(pars) <- c("alpha","beta","gamma")
coefi <- data.frame(c(suaviza$coef[1:2],efectossuav))
names(coefi) <- ""
cat("\n")
cat("Call:")
cat("\n")
print(suaviza$call)
cat("\n")
cat("Smoothing parameters:")
print(t(pars))
cat("\n")
cat("Coefficients:")
print(coefi)
result <- list(coefficients=coefi,fitted=ythat,residuals=res,forecast=predicc,MSE=MSE)
result
}
