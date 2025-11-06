#' Función `Mipoly()`
#'
#' Genera las variables de los polinomios para modelos de regresión con funciones de tendencia polinomiales. Se basa en la función `poly()` de `stats`.
#' @param tiempo Un vector con los valores del indice de tiempo t.
#' @param grado El grado p del polinomio deseado, debe ser un valor entero positivo a partir de 1.
#' @return Un `data.frame` cuyas variables son `t` y sus potencias denominadas `t2`, `t3`, `...`, `tp`, según el grado polinomial p.
#' @export
#' @examples
#' t <- 1:20 #índice de tiempo primeros 20 valores
#' poli2 <- Mipoly(tiempo=t,grado=2) #para predictores con polinomio de grado 2
#' poli2
#' poli3 <- Mipoly(tiempo=t,grado=3) #para predictores con polinomio de grado 3
#' poli3
Mipoly <- function(tiempo,grado){
 if(grado>1){
 poli <- data.frame(stats::poly(tiempo,degree=grado,raw=TRUE))
 names(poli) <- c("t",paste0("t",2:ncol(poli)))
}
 if(grado==1){
 poli <- data.frame(t=tiempo)
}
 poli
}
