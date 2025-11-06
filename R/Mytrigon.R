#' Función `Mytrigon()`
#'
#' Genera las variables trigonométricas `sin(2*pi*Fj*tiempo)` y `cos(2*pi*Fj*tiempo)` en frecuencias Fj.
#' @param tiempo Vector de valores enteros positivos con el índice de tiempo para los cuales se debe calcular las funciones trigonométricas.
#' @param Frecuencias Vector con las frecuencias `Fj` correspondientes a las ondas sinusoidales armónicas a considerar.
#' @param indicej Vector de valores enteros indicando el sub índice `j` de las funciones trigonométricas para cada frecuencia en `Frecuencias` de manera que a `Frecuencias[i]` le asigna el indice `indicej[i]`.
#' @return Un `data.frame` con las variables `sinFj`, `cosFj`, siendo `j` el valor asignado en `indicej`. Para la frecuencia 1/2 solo crea la componente coseno, que corresponde a `cos(pi*tiempo)`. Los pares de variables `sinj`, `cosj` aparecen en orden de menor a mayor frecuencia pero conservando el sufijo `j` asignado en `indicej`.
#' @export
#' @examples
#' t <- 1:20 #índice de tiempo primeros 20 valores
#' #Trigonométricas en Fj=j/12 para j=1 a 6 y F7=0.35
#' trigono <- Mytrigon(tiempo=t,Frecuencias=c(c(1:6)/12,0.35),indicej=c(1,2,3,4,5,6,7))
#' trigono #Variables sinF7, cosF7 corresponden a la frecuencia F7=0.35
Mytrigon <- function(tiempo,Frecuencias,indicej){
F <- sort(Frecuencias)
if(max(F)!=0.5){
K <- length(F)
trig <- matrix(,ncol=2*K,nrow=length(tiempo))
for(i in 1:K){
trig[,2*i-1]=sin(2*F[i]*pi*tiempo)
trig[,2*i]=cos(2*F[i]*pi*tiempo)
}
trig <- data.frame(trig)
J <- indicej[order(Frecuencias)]
names(trig) <- paste0(rep(c("sinF","cosF"),times=K),rep(J,each=2))
}
if(max(F)==0.5){
K <- length(F)-1
trig <- matrix(,ncol=2*K+1,nrow=length(tiempo))
for(i in 1:K){
trig[,2*i-1] <- sin(2*F[i]*pi*tiempo)
trig[,2*i] <- cos(2*F[i]*pi*tiempo)
}
trig[,2*K+1] <- cos(2*F[K+1]*pi*tiempo)
trig <- data.frame(trig)
J <- indicej[order(Frecuencias)]
names(trig) <- c(paste0(rep(c("sinF","cosF"),times=K),rep(J[1:K],each=2)),paste0("cosF",J[K+1]))
}
trig
}

