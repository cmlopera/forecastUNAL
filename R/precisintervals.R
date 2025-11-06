#' Función `precisintervals()`
#'
#' Evalúa precisión de intervalos de predicción (IPs).
#'
#' @param real Vector numérico con los valores observados en los periodos de pronóstico.
#' @param LIP Vector numérico con los valores de los límites inferiores de los IPs.
#' @param LSP Vector numérico con los valores de los límites superiores de los IPs.
#' @param alpha Escalar entre 0 y 1 correspondiendo al nivel `1-alpha` de confianza.
#'
#' @return Un vector con `AmplitudProm` la amplitud promedio, `Cobertura(%)` la cobertura (porcentaje) y `ScoreProm` el Score promedio de los IPs.
#' @export
#'
#' @examples
#' Real <- c(181958.0,185303.0,183429.0)
#' Pronost <- data.frame(rbind(c(175546.9,163886.4,188037.2),
#' c(175938.6,164246.2,188463.5),c(176323.1,164599.1,188882.2)))
#' names(Pronost) <- c("Pred","LimInf","LimSup")
#' Pronost
#' precisintervals(real=Real,LIP=Pronost[,2],LSP=Pronost[,3])
precisintervals <- function(real,LIP,LSP,alpha=0.05){
  a <- LSP-LIP
  am <- mean(a)
  Indi <- ifelse(real>=LIP & real<=LSP,1,0)
  p <- mean(Indi)*100
  Scalpha <- (LSP-LIP)+(2/alpha)*(LIP-real)*ifelse(real<LIP,1,0)+(2/alpha)*(real-LSP)*ifelse(real>LSP,1,0)
  Scalphamean <- mean(Scalpha)
  res <- list(AmplitudProm=am,"Cobertura(%)"=p,ScoreProm=Scalphamean)
  unlist(res)
}
