#' Función `regexponencialv02()`
#'
#' Ajusta modelos de regresión de la forma `g(x)=exp(t(x)%*%beta)`, con `t(x)` la transpuesta de un vector de predictores `x` y `beta` el vector de parámetros.
#'
#' @param respuesta El vector de valores de la variable respuesta.
#' @param data Un objeto `data.frame` con los predictores a usar.
#' @param control Una lista opcional de configuraciones de control. Consulte `nls.control()` para los valores de control configurables y su efecto.
#'
#' @return Son los mismos predeterminados para los objetos tipo `nls` generados por la función R `nls()`.
#' @export
#'
#' @examples
#' #Simulando serie aditiva con modelo de tendencia exponencial lineal y errores iid normales
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
#'
#' #Gráficos de residuos de ajuste
#' diagresTSModel(modelo=modelo,RMSE=summary(modelo)$sigma,all=FALSE,single=FALSE)
regexponencialv02 <- function(respuesta,data,control=stats::nls.control()){
yt <- respuesta
model.aux <- stats::lm(log(yt)~.,data)
names.vars <- c(1,names(data))
names.param <- c("Intercept",paste0("param_",names(data)))
miformula <- stats::as.formula(paste("yt~",paste(paste("exp(",paste(names.param,names.vars,sep="*",collapse="+"),sep=""),")",sep="")))
coef0 <- as.list(stats::coef(model.aux));names(coef0) <- as.list(names.param)
modelo <- stats::nls(miformula,start=coef0,control=control,data=data)
modelo
}
