#' Función `diagresTSModel()`
#'
#' Genera gráficos de residuos sobre modelos de series de tiempo. Adicionalmente, puede generar la ACF, PACF muestrales, el gráfico de probabilidad normal, test Shapiro-Wilks y Tests Ljung-Box.
#'
#' @param modelo Un objeto R con el ajuste de un modelo de serie de tiempo. Puede ser de clase `lm`, `nls`, `Arima`, `HoltWinters`, entre otros manejados en la asignatura.
#' @param RMSE Un valor numérico correspondiente a una estimación de la desviación estándar del error de ajuste. Se requiere para las líneas horizontales pasando por `-2*RMSE` y `+2*RMSE.`
#' @param all Argumento lógico. Por defecto es FALSE indicando que solo los gráficos de residuos vs. tiempo y vs. valores ajustados serán producidos. Si es `TRUE`, exhibe además los gráficos de ACF y PACF muestrales, de probabilidad normal, test Shapiro-Wilks y Tests Ljung-Box para `m=6, 12,...`, hasta máximo múltiplo de 6 menor o igual a `lagmax`.
#' @param single Un valor lógico. Por defecto es `TRUE` indicando que cada gráfica debe ser presentada en ventana independiente. Si se ajusta a `FALSE` y `all=FALSE`, los dos gráficos de residuos son presentados en una sola ventana en un arreglo de `2x1`, y si `all=TRUE`, los cinco gáficos que se generan son presentados en una sola ventana, en un arreglo de `3x2`.
#' @param lagmax Valor numérico indicando el máximo orden de rezago a considerar en gráficos de la ACF y PACF muestrales y en las pruebas Ljung-Box. Por defecto su valor es 36.
#' @return La función no genera ningún tipo de objeto R sobre el cual luego pueda obtenerse algún valor o resultado, solo las gráficas pedidas y resultados de pruebas exhibidas en la consola R.
#' @export
#' @examples
#' library(forecast)
#' library(lmtest)
#' library(car)
#' datos <- ts(CEMENTO,fre=4,start=c(1956,1))
#' n <- length(datos)-4
#' t <- 1:n
#' yt <- ts(datos[t],frequency=4,start=c(1956,1))
#' poli <- Mipoly(tiempo=t,grado=3)
#' trimestre <- seasonaldummy(yt)
#' X <- data.frame(poli,trimestre) #La matriz de los predictores en el ajuste
#'
#' #Modelo log cúbico estacional con indicadoras (nivel de ref.=Q4)
#' modelo1 <- lm(log(yt)~.,data=X)
#' summary(modelo1)
#' diagresTSModel(modelo=modelo1,RMSE=summary(modelo1)$sigma)
#' diagresTSModel(modelo=modelo1,RMSE=summary(modelo1)$sigma,all=TRUE,single=TRUE)
#' diagresTSModel(modelo=modelo1,RMSE=summary(modelo1)$sigma,all=TRUE,single=FALSE)
#' diagresTSModel(modelo=modelo1,RMSE=summary(modelo1)$sigma,all=FALSE,single=FALSE)
#'
#' #Modelo exponencial cúbico estacional con indicadoras (nivel de ref.=Q4)
#' modelo2b <- regexponencialv02(respuesta=yt,data=X)
#' summary(modelo2b)
#' diagresTSModel(modelo=modelo2b,RMSE=summary(modelo2b)$sigma,all=TRUE,single=FALSE)
#'
#' #Modelo log cúbico estacional con indicadoras (nivel de ref.=Q4) y error AR(6)
#' modelo2 <- Arima(log(yt),order=c(6,0,0),xreg=as.matrix(X),method="ML")
#' df2 <- n-length(coef(modelo2)[coef(modelo2)!=0])
#' coeftest(modelo2,df=df2)
#' diagresTSModel(modelo=modelo2,RMSE=sqrt(modelo2$sigma2),all=TRUE,single=FALSE)
#'
#' #Suavizamiento Holt-Winters multiplicativo con predicción de h=4 trimestres
#' modelo3 <- SuavizamientoEstacional(yt,seasonal="multiplicative",h=4)
#' diagresTSModel(modelo=modelo3,RMSE=sqrt(modelo3$MSE),all=TRUE,single=FALSE)
#'
#' #Combinación filtros de descomposición multiplicativa y
#' #Loess cuadrático óptimo con span según criterio AICC
#' modelo4 <- DescompLoessv02(serie.ajuste=yt,tipo.descomp="multiplicative",
#'            grado=2,criterio="aicc",h=4,level=0.95)
#' diagresTSModel(modelo=modelo4,RMSE=sqrt(modelo4$MSE),all=TRUE,single=FALSE)
diagresTSModel <- function(modelo,RMSE,all=FALSE,single=TRUE,lagmax=36){
  acflimits <- function(x){
    stopifnot(x$type=="correlation")
    ci.type <- "ma"
    ci <- 0.95
    C0 <- stats::qnorm((1+ci)/2)/sqrt(x$n.used)
    maxlag <- max(x$lag)
    minlag <- min(x$lag)
    if(minlag==0){
      aux <- C0*sqrt(1+2*cumsum((x$acf[2:(maxlag)])^2))
      limites <- c(C0,aux)
    }

    if(minlag>0){
      aux <- C0*sqrt(1+2*cumsum((x$acf[1:(maxlag-1)])^2))
      limites <- c(C0,aux)
    }
    res <- list(lag=x$lag,acf=x$acf,clim=limites,citype=ci.type)
    res
  }


  resdata <- data.frame(time=1:length(stats::residuals(modelo)),
                        residuals=as.numeric(stats::residuals(modelo)),
                        fitted=as.numeric(stats::fitted(modelo)))

  g1 <- ggplot2::ggplot(resdata,ggplot2::aes(x=time,y=residuals))+
    ggplot2::geom_line(,na.rm =TRUE)+
    ggplot2::geom_hline(yintercept=-2*RMSE,color=4)+
    ggplot2::geom_hline(yintercept=2*RMSE,color=4)+
    ggplot2::geom_hline(yintercept=0,color=4)+
    ggplot2::ggtitle("residuos vs. tiempo")+
    ggplot2::theme(plot.title=ggplot2::element_text(size=10))+
    ggplot2::theme(axis.title=ggplot2::element_text(size=10))

  g2 <- ggplot2::ggplot(resdata,ggplot2::aes(x=fitted,y=residuals))+
    ggplot2::geom_point(,na.rm =TRUE)+
    ggplot2::geom_hline(yintercept=-2*RMSE,color=4)+
    ggplot2::geom_hline(yintercept=2*RMSE,color=4)+
    ggplot2::geom_hline(yintercept=0,color=4)+
    ggplot2::ggtitle("residuos vs. ajustados")+
    ggplot2::theme(plot.title=ggplot2::element_text(size=10))+
    ggplot2::theme(axis.title=ggplot2::element_text(size=10))

  res2 <- acflimits(stats::acf(as.numeric(stats::residuals(modelo)),lag.max=lagmax,na.action=na.pass,plot=F))

  if(min(res2$lag)==0){
    dataACF <- data.frame(lag=res2$lag[-1],acf=res2$acf[-1],clim=res2$clim)
    g3<-ggplot2::ggplot(dataACF,ggplot2::aes(x=lag,xend=lag,y=0,yend=acf))+
      ggplot2::scale_x_continuous(breaks = seq(5, max(res2$lag), by = 5))+
      ggplot2::theme(panel.grid.minor.x = ggplot2::element_blank())+
      ggplot2::geom_segment(,na.rm =TRUE)+
      ggplot2::geom_line(ggplot2::aes(x=lag,y=clim,color=4),linetype = 2)+
      ggplot2::geom_line(ggplot2::aes(x=lag,y=-clim,color=4),linetype = 2)+
      ggplot2::geom_hline(yintercept=0,color='black')+
      ggplot2::theme(legend.position = "none")+
      ggplot2::ylab("ACF")+
      ggplot2::ggtitle("ACF residuos")+
      ggplot2::theme(plot.title=ggplot2::element_text(size=10))+
      ggplot2::theme(axis.title=ggplot2::element_text(size=10))
  }

  if(min(res2$lag)>0){
    dataACF <- data.frame(lag=res2$lag,acf=res2$acf,clim=res2$clim)
    g3<-ggplot2::ggplot(dataACF,ggplot2::aes(x=lag,xend=lag,y=0,yend=acf))+
      ggplot2::scale_x_continuous(breaks = seq(5, max(res2$lag), by = 5))+
      ggplot2::theme(panel.grid.minor.x = ggplot2::element_blank())+
      ggplot2::geom_segment(,na.rm =TRUE)+
      ggplot2::geom_line(ggplot2::aes(x=lag,y=clim,color=4),linetype = 2)+
      ggplot2::geom_line(ggplot2::aes(x=lag,y=-clim,color=4),linetype = 2)+
      ggplot2::geom_hline(yintercept=0,color='black')+
      ggplot2::theme(legend.position = "none")+
      ggplot2::ylab("ACF")+
      ggplot2::ggtitle("ACF residuos")+
      ggplot2::theme(plot.title=ggplot2::element_text(size=10))+
      ggplot2::theme(axis.title=ggplot2::element_text(size=10))
  }

  g4 <- forecast::ggPacf(as.numeric(stats::residuals(modelo)),lag.max=lagmax)+
    ggplot2::ggtitle("PACF residuos")+
    ggplot2::theme(plot.title=ggplot2::element_text(size=10))+
    ggplot2::theme(axis.title=ggplot2::element_text(size=10))

  g5<-ggplot2::ggplot(resdata, ggplot2::aes(sample = residuals)) +
    ggplot2::stat_qq(,na.rm = TRUE) +
    ggplot2::stat_qq_line(,na.rm = TRUE) +
    ggplot2::ggtitle("Normal Probability Plot")+
    ggplot2::theme(plot.title=ggplot2::element_text(size=10))+
    ggplot2::stat_qq_line(color=4,na.rm = TRUE)+
    ggplot2::xlab("Theorical Quantiles")+
    ggplot2::ylab("Sample Quantiles")+
    ggplot2::theme(axis.title=ggplot2::element_text(size=10))


  if(all==TRUE && single==FALSE){
    cat("Resultados Tests Ljung-Box")
    cat("\n")
    cat("\n")
    print(BPLBtest(stats::residuals(modelo),maxlag=lagmax))
    cat("\n")
    cat("Resultados Test de Normalidad Shapiro-Wilk")
    cat("\n")
    print(stats::shapiro.test(stats::residuals(modelo)))

    print(gridExtra::marrangeGrob(list(g1,g3,g5,g2,g4),nrow = 3,ncol= 2))
  }
  if(all==TRUE && single==TRUE){
    cat("Resultados Tests Ljung-Box")
    cat("\n")
    cat("\n")
    print(BPLBtest(stats::residuals(modelo),maxlag=lagmax))
    cat("\n")
    cat("Resultados Test de Normalidad Shapiro-Wilk")
    cat("\n")
    print(stats::shapiro.test(stats::residuals(modelo)))

    print(gridExtra::marrangeGrob(list(g1,g2,g3,g4,g5),nrow = 1,ncol= 1))
  }

  if(all==FALSE && single==FALSE){
    print(gridExtra::marrangeGrob(list(g1,g2),nrow = 2,ncol= 1))
  }

  if(all==FALSE && single==TRUE){
    print(gridExtra::marrangeGrob(list(g1,g2),nrow = 1,ncol= 1))
  }
}

