#' Función `regexpoErrorARMAv03()`
#'
#' Ajusta y pronostica de forma aproximada modelos de regresión exponenciales como el considerado en la función de usuario `regexponencialv02()`, pero con error ARMA, así: primero estima los parámetros de la regresión exponencial por mínimos cuadrados no lineales, luego ajusta a los residuos de este modelo un ARMA estacionario de media cero, mediante la funcion `Arima()` de la libreria `forecast`. Los valores estimados de la serie son aproximados sumando las estimaciones de la regresión no lineal y del modelo ARMA. El anterior procedimiento también es aplicado en la construcción de pronósticos puntuales. Intervalos de pronóstico son construidos aproximando el error estándar del error de predicción como la raíz cuadrada de la suma de la varianza estimada del error de predicción sobre el modelo ARMA y la varianza estimada (por el método delta) de la estructura exponencial evaluada en los tiempos de predicción.
#'
#' @param respuesta Un objeto `ts` univariado con los datos a ser ajustados.
#' @param data,newdata Objetos tipo `data.frame`, cuyas columnas corresponden a los valores de las variables predictoras del modelo en el ajuste y en el pronóstico, respectivamente.
#' @param level Nivel de confianza para intervalos de predicción, dado como un valor en `(0,1)`, por defecto es 0.95.
#' @param control Una lista opcional de configuraciones de control para la función `nls()`. Consulte `nls.control()` para conocer los nombres de los parámetros de control configurables y sus efectos.
#' @param order La especificación de la parte no estacional del modelo ARIMA. Los tres componentes `(p,d,q)` son los órdenes AR, de diferenciación y el orden MA, respectivamente.
#' @param seasonal La especificación de la parte estacional del modelo ARIMA, más el período (que por defecto es `frequency(respuesta)`). Este argumento debe ser dado como una lista, con el orden de los componentes y el período, pero una especificación de solo un vector numérico de longitud 3 se convertirá en una lista adecuada con las especificaciones de los órdenes estacionales `(P,D,Q)` de la parte AR, de la diferenciación y del orden MA, de periodo `s=frequency(respuesta)`.
#' @param fixed Un vector numérico opcional de la misma longitud que la suma de los órdenes de la ecuación del ARMA estacionario de media cero a ajustar sobre los residuos del ajuste exponencial. Ver en la ayuda de la función `arima()` los detalles de este argumento.
#' @param method Método de ajuste, a saber, máxima verosimilitud `("ML")`, Máxima verosimilitud combinada con mínimos cuadrados condicionales `("CSS-ML")` o mínimos cuadrados condicionales `("CSS")`. El valor predeterminado (a menos que existan valores faltantes) es `"CSS-ML"`, en el cual se usan mínimos cuadrados condicionales para encontrar valores iniciales y luego aplica máxima verosimilitud.
#' @param optim.method El valor pasado al argumento `‘method’` para `‘optim’` en la funcion `Arima()`. Por defecto es `"BFGS"`.
#' @param optim.control Una lista de parámetros de control para `‘optim’`, usados en la función `Arima()`.
#'
#' @return La función produce una lista con los siguientes componentes:
#' * `coefficients` Una matriz con la tabla de parámetros estimados, sus errores estándar, el estadístico `T0` y valor P asociado. Tenga en cuenta que no hay una estimación conjunta de los parámetros de regresión de la funcion exponencial y de los parámetros del modelo ARMA del error estructural, de modo que los errores estándar provienen del ajuste separado de las dos estructuras, pero los valores p son calculados bajo una distribución t-student cuyos grados de libertad son `df=n - (total de parámetros)`. Esta tabla puede ser extraida con la función `coef()` aplicada al objeto R donde se guarde la estimación del modelo.
#' * `fitted` Objeto tipo `"ts"` univariado con los valores ajustados de la respuesta. Estos valores pueden ser extraidos mediante la función `fitted()` aplicada    sobre el objeto R donde se guarde la estimación del modelo.
#' * `residuals` Objeto tipo `"ts"` univariado con los residuos del ajuste de la respuesta. Estos valores pueden ser extraidos mediante la función `residuals()` aplicada sobre el objeto R donde se guarde la estimación del modelo.
#' * `sigma2` Estimación de la varianza de las innovaciones del modelo ARMA definido para el error estructural del modelo de regresión exponencial. Puede ser extraido como `nombre_objeto$sigma2`, donde `'nombre_objeto'` es el nombre del objeto R donde se guarda la estimación del modelo.
#' * `forecast` Objeto tipo `"mts"`, `"ts"`, `"matrix"`, `"array"`, con los pronósticos puntuales y por intervalos de la respuesta para `h=nrow(newdata)` períodos después del ajuste.
#' @export
#' @examples
#' library(forecast)
#' #Serie estacional trimestral multiplicativa
#' #de N=155 observaciones, inicio 1956-Q1
#' datos <- ts(CEMENTO,freq=4,start=c(1956,1))
#' n <- length(datos)-4
#' t <- 1:n
#' #Variables para ajuste y pronóstico modelo exponencial
#' #cúbico estacional con indicadoras (nivel ref=Q4) y error AR(6)
#' #serie de primeros n valores
#' yt <- ts(CEMENTO[t],freq=4,start=start(datos))
#' poli <- Mipoly(tiempo=t,grado=3)
#' trimestre <- seasonaldummy(yt)
#' #data.frame de los predictores en el ajuste
#' X <- data.frame(poli,trimestre)
#' tnuevo <- (n+1):length(datos)
#' #Fecha inicio predicciones
#' startpred <- end(ts(datos[1:(n+1)],freq=4,start=c(1956,1)))
#' #serie de los últimos cuatro valores
#' ytf <- ts(datos[tnuevo],freq=4,start=startpred)
#' polinuevo <- Mipoly(tiempo=tnuevo,grado=3)
#' trimestrenuevo <- seasonaldummy(yt,h=4)
#' #data.frame de los predictores  en la predicción
#' Xnuevo <- data.frame(polinuevo,trimestrenuevo)
#'
#' #Ajuste y pronóstico del modelo
#' modelo <- regexpoErrorARMAv03(respuesta=yt,data=X,newdata=Xnuevo,
#'            order=c(6,0,0),method="ML")
#' coef(modelo)
#'
#' #Criterios de información
#' k <- modelo$p;k
#' critinfresid(residuales=residuals(modelo),npar=k)
#'
#' #Serie ajustada
#' ythat <- fitted(modelo)
#' plot(datos)
#' lines(ythat,col=2)
#' legend("topleft",legend=c("Observaciones","ajuste"),lty=1,col=1:2,lwd=2)
#'
#' #Pronósticos y su precisión
#' predmod <- modelo$forecast; predmod
#' ytpron <- predmod[,1]
#' accuracy(ytpron,ytf)
#' precisintervals(real=ytf,LIP=predmod[,2],LSP=predmod[,3])
#'
#' #Gráficos y resultados para validar supuestos sobre errores de ajuste
#' diagresTSModel(modelo=modelo,RMSE=sqrt(modelo$sigma2),all=TRUE,single=FALSE)
regexpoErrorARMAv03 <- function(respuesta,data,newdata,
                                level=0.95,
                                control=stats::nls.control(),
                                order= c(0L, 0L, 0L),
                                seasonal=list(order = c(0L, 0L, 0L),period=NA),
                                fixed=NULL,method=c("CSS-ML", "ML", "CSS"),
                                optim.method="BFGS",
                                optim.control = list()){
  yt <- respuesta
  model.aux <- stats::lm(log(yt)~.,data)
  names.vars <- c(1,names(data))
  names.param <- c("Intercept",paste0("param_",names(data)))
  miformula <- stats::as.formula(paste("yt~",paste(paste("exp(",paste(names.param,names.vars,sep="*",collapse="+"),sep=""),")",sep="")))
  coef0 <- as.list(stats::coef(model.aux));names(coef0)=as.list(names.param)

    modelo <- stats::nls(miformula,start=coef0,control=control,data=data)
  serie.et <- stats::ts(stats::residuals(modelo),freq=stats::frequency(yt),
                        start=stats::start(yt))
  ciclos <- forecast::Arima(serie.et,order=order,seasonal=seasonal,
                            include.mean=F,
                            fixed=fixed,
                            method=method,
                            optim.method=optim.method,
                            optim.control=optim.control)

  p1 <- length(stats::coef(modelo)[stats::coef(modelo)!=0])
  p2 <- p1+length(stats::coef(ciclos)[stats::coef(ciclos)!=0])
  df <- length(yt)-p2
  tabla <- rbind(lmtest::coeftest(ciclos,df=df),lmtest::coeftest(modelo,df=df))[,1:4]
  hatexpo <- stats::ts(stats::fitted(modelo),freq=stats::frequency(yt),
                       start=stats::start(yt))
  hatciclos <- stats::fitted(ciclos)
  ythat <- hatexpo + hatciclos
  sigma2 <- ciclos$sigma2
  s <- stats::tsp(yt)[3]
  tajuste <- 1:length(yt)
  h <- nrow(newdata)
  tpron <- (length(yt)+1):(length(yt)+h)
  if(stats::end(yt)[2]<s){
    start.pron <- c(stats::end(yt)[1],stats::end(yt)[2]+1)
  }
  if(stats::end(yt)[2]==s){
    start.pron <- c(stats::end(yt)[1]+1,1)
  }
  resid <- stats::residuals(ciclos)

  pronexpo <- stats::ts(stats::predict(modelo,newdata=newdata),
                        freq=stats::frequency(yt),start=start.pron)

  pronciclos <- stats::ts(forecast::forecast(ciclos,h=h)$mean,
                          freq=stats::frequency(yt),start=start.pron)

  ytpron <- pronexpo + pronciclos

  #IPs aproximados
  #el vector de predictores en i esima prediccion, incluye al 1 como primera entrada
  #matriz con valores de variables explicatorias en la prediccion
  x0 <- cbind(1,as.matrix(newdata))
  #estimacion varianza del error del modelo de regresion exponencial
  mse <- summary(modelo)$sigma^2
  #vector prediccion estructura exponencial
  predexp <- stats::predict(modelo,newdata=newdata)

  #matriz de varianzas covarianzas parametros estimados modelo exponencial
  varbetas <- as.matrix(stats::vcov(modelo))
  #aproximacion varianza estructura exponencial estimada
  Sigmapred <- (predexp^2)*diag(x0%*%varbetas%*%t(x0))
  se.res <- stats::predict(ciclos,n.ahead=h)$se
  Sigmatot <- Sigmapred+(se.res)^2
  LI <- ytpron-stats::qt((1-level)/2,df=df,lower.tail=F)*sqrt(Sigmatot)
  LS <- ytpron+stats::qt((1-level)/2,df=df,lower.tail=F)*sqrt(Sigmatot)
  tablapron <- stats::ts(data.frame(forecast=ytpron,lower=LI,upper=LS),
                         freq=stats::frequency(yt),start=start.pron)

  result <- list(coefficients=tabla,fitted=ythat,residuals=resid,sigma2=sigma2,
                 p=p2,forecast=tablapron)
  result
}
