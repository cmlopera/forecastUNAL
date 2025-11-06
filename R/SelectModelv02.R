#' Función `SelectModelv02()`
#'
#' Usando criterios de información selecciona el mejor modelo AR(p) completo. Está basado en la Función SelectModel() de la librería FitAR (ya no disponible en Cran-R) de A.I. McLeod and Y. Zhang.
#'
#' @param z Una serie de tiempo.
#' @param lag.max Un entero indicando el máximo orden autorregresivo para el proceso `AR(p)`.
#' @param Criterion El criterio de información a usar para la selección del orden `p`. Por defecto es `"BIC"` (Bayesian Information Criterion). Otras opciones son `"AIC"` (Akaike Information Criterion), `"UBIC"` (Uncertainty-Penalized Bayesian Information Criterion), `"EBIC"` (Extended BIC), y `"GIC"` (Generalizad Information Criterion).
#' @param Best El número final de modelos a ser seleccionados.
#' @param Candidates El número de inicialmente seleccionados usando el criterio aproximado.
#' @param t Parámetro de ajuste en los criterios`"EBIC"`, `"GIC"`.
#'
#' @return La función crean un objeto tipo `"matrix"`, `"array"` de `Best` filas por 3 columnas, con los siguientes valores en sus columnas: `p` que corresponde a los valores del orden autorregresivo para los `Best` primeros casos de mínimo valor del criterio usado; los valores del criterio de información correspondiente, el exacto: `Criterion-Exact` y el aproximado: `Criterion-Approx`, donde `Criterion` corresponde al acrónimo del criterio de información usado.
#' @export
#' @examples
#' #Simulando serie de un proceso AR(1) estacionario de media cero
#' zt <- ts(arima.sim(n=300,list(ar=0.7),sd=4000),freq=1,start=1)
#' SelectModelv02(z=zt,Criterion="AIC")
#' SelectModelv02(z=zt,Criterion="BIC")
#'
#' #Simulando serie de un proceso AR(6) estacionario de media cero
#' xt <- arima.sim(n=300,list(ar=c(0.65,0.32,-0.27,0.25,0,-0.23)),sd=0.025)
#' SelectModelv02(z=xt,Criterion="AIC")
#' SelectModelv02(z=xt,Criterion="BIC")
#' SelectModelv02(z=xt,Criterion="EBIC")
#' SelectModelv02(z=xt,Criterion="UBIC")
#' SelectModelv02(z=xt,Criterion="GIC")
SelectModelv02 <- function(z, lag.max=15, Criterion="default", Best=3, Candidates=5, t="default"){
  ARModel="AR"
  stopifnot(length(z)>0, length(z)>lag.max, lag.max>1, Best>0, Candidates>0)
  is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
  stopifnot(is.wholenumber(lag.max))
  `ARToPacf` <- function(phi){ # in DetAR, GetFitARz, SelectModel
    phik=phi
    L=length(phi)
    if(L==0) return(0)
    pi=numeric(L)
    for (k in 1:L){
      LL=L+1-k
      a <- phik[LL]
      pi[L+1-k] <- a
      phikp1 <- phik[-LL]
      if(is.na(a) || abs(a)==1)
        stop("transformation is not defined, partial correlation = 1")
      phik <- (phikp1+a*rev(phikp1))/(1-a^2)
    }
    pi
  }
  `AR1Est` <- function(z, MeanValue=0){ # in GetFitARz
    stopifnot(length(z)>1)
    n <- length(z)
    m <- MeanValue
    mz <- rep(m,n)
    a <- sum((z-mz)^2)
    b <- sum((z[-1]-mz[-1])*(z[-n]-mz[-n]))
    c <- sum((z[c(-1,-n)]-mz[c(-1,-n)])^2)
    i <- complex(1,0,1)
    x <- ((-16)*b^3+18*a*b*c+24*b^3*n-27*a*b*c*n-9*b*c^2*n-12*b^3*n^2+9*a*b*c*n^2+27*b*c^2*n^2+2*b^3*n^3-18*b*c^2*n^3)
    y <- (-b^2*(-2+n)^2+3*c*(-1+n)*(-a-c*n))
    f <- complex(1,x^2+4*y^3,0)
    z <- (x+sqrt(f))^(1/3)
    g <- x^2+4*y^3
    z1 <- (x+(-g)^(1/2)*i)^(1/3)
    part1 <- (n-2)*b/(3*c*(n-1))
    part2 <- (1-sqrt(3)*i)*y/(3*2^(2/3)*c*(n-1)*z)
    part3 <- (1+sqrt(3)*i)*z/(6*2^(1/3)*c*(n-1))
    Re(part1+part2-part3)
  }
  `ChampernowneD` <- function(z, p, MeanZero=FALSE){ # in GetFitARz, LoglikelihoodAR
    n<-length(z)
    if(MeanZero) y<-z
    else y<-z-mean(z)
    x0<-x<-y
    for (i in 1:p)
      x<-c(x,c(rep(0,i),y[1:(n-i)]))
    x<-matrix(x, nrow=p+1, ncol=length(x0), byrow=TRUE)
    C<-c(x%*%x0)
    A<-stats::toeplitz(C)
    E<-matrix(0, nrow=p+1, ncol=p+1)
    for (j in 1:p)
      for (i in 1:j){
        E[i+1,j+1] <- E[i,j]+y[i]*y[j]+y[n+1-i]*y[n+1-j]
      }
    for (j in 1:(p+1))
      for (i in 1:j)
        E[j,i]=E[i,j]
    A-E
  }
  `DetAR` <- function(phi){ # in FastLoglikelihoodAR, LoglikelihoodAR
    z<-ARToPacf(phi)
    1/prod((1-z^2)^(1:length(phi)))
  }
  `FastLoglikelihoodAR` <- function(phi, n, CD){ # in GetFitARz
    phis<-c(1,-phi)
    LL <- -log(DetAR(phi))/2-(n/2)*log(sum(crossprod(phis,CD)*phis)/n)
    if (!is.finite(LL)){
      LL<--1e35 }
    LL
  }
  `LoglikelihoodAR` <- function(phi, z, MeanValue=0){ # in GetFitARz
    if(length(phi)==0) phi=0
    phis<-c(1,-phi)
    y<-z-MeanValue
    n<-length(z)
    -log(DetAR(phi))/2 - (n/2)*log(sum(crossprod(phis,ChampernowneD(y,length(phis)-1,MeanZero=TRUE))*phis)/n)
  }
  `PacfToAR` <- function(zeta){ # in GetFitARz
    L=length(zeta)
    if (L==0) return(numeric(0))
    if (L==1) return(zeta)
    phik=zeta[1]
    for (k in 2:L){
      phikm1=phik
      phik=c(phikm1-zeta[k]*rev(phikm1),zeta[k])
    }
    phik
  }
  `GetFitARz` <- function(z, pvec, MeanValue=0, ...){ # in GetFitAR, SelectModel
    stopifnot(length(z)>0, length(z)>=2*max(pvec), length(pvec)>0,all(pvec>=0))
    pVec<-pvec[pvec>0]
    y <- z-MeanValue
    n<-length(y)
    if (length(pVec)==0 || length(pvec)==0)
      return(list(loglikelihood=-(n/2)*log(sum(y^2)/n),zetaHat=NULL,phiHat=NULL,convergence=0,algorithm="cubic"))
    PMAX <- max(pVec)
    if (PMAX == 1){
      phiHat <- AR1Est(y)
      LogL <- LoglikelihoodAR(phiHat,z)
      return(list(loglikelihood=LogL, zetaHat=phiHat, phiHat=phiHat, convergence=0))
    }
    PEFF<-length(pVec)
    CD<-ChampernowneD(y,PMAX,MeanZero=TRUE)
    xpar<-numeric(PMAX)
    EntropyAR<-function(x){
      if (max(abs(x))>0.999)
        out<-1e35
      else {
        xpar[pVec]<-x
        out<- -FastLoglikelihoodAR(PacfToAR(xpar),n,CD)
      }
      out
    }
    xinit<-ARToPacf(stats::ar.burg(y, aic=FALSE, order.max=PMAX, demean=FALSE)$ar)[pVec]
    #Sometimes there are problems with "L-BFGS-B" -- it frequently tests the endpoints which is bad news due
    #to numerical problems such as ARToPacf(PacfToAR(rep(0.99,20))) is not correct!
    #So it is better to use "BFGS" with a penalty instead.
    #ans<-optim(xinit,EntropyAR,method="L-BFGS-B", lower=rep(-0.9999,PEFF), upper=rep(0.9999,PEFF),control=list(trace=6),...)
    ans<-stats::optim(xinit,EntropyAR,method="BFGS", control="trace", ...)
    alg<-1
    if(ans$convergence != 0) {
      alg<-2
      warning("Convergence problem. convergence=", ans$convergence)
      warning("Trying Nelder-Mead algorithm ...")
      ans<-stats::optim(xinit,EntropyAR,method="Nelder-Mead", ...)
      if(ans$convergence != 0)
        warning("Still convergence problem, convergence= ", ans$convergence)
    }
    zetaHat<-ans$par
    zetas<-numeric(PMAX)
    zetas[pVec]<-zetaHat
    list(loglikelihood=-ans$value, zetaHat=ans$par, phiHat=PacfToAR(zetas),convergence=ans$convergence, algorithm=c("BFGS","Nelder-Mead")[alg],pvec=pvec)
  }
  BestCandidates<-Candidates
  IsValidCriterionQ <- Criterion %in% c("default", "AIC", "BIC", "UBIC", "EBIC", "BICq", "GIC")
  if (!IsValidCriterionQ)
    stop("Criterion = ", Criterion, " not known.")
  if (Best > BestCandidates)
    BestCandidates<-Best
  SubsetQ <- FALSE
  method<-Criterion
  if (Criterion == "default")
    method <- "BIC"
  if (!SubsetQ && Criterion=="UBIC")
    method <- "BIC"
  #set tuning parameter
  P<-0.01
  Q<-0.25
  G<-1
  if (method=="EBIC"  && t!="default")  G <- t
  if (method=="QBIC"  && t!="default")  Q <- t
  if (method=="GIC"   && t!="default")  P <- t
  if (P>=0.25 || P<=0)
    stop("error: GIC tuning parameter invalid")
  if (Q<=0 || Q>=1)
    stop("error: BICq tuning parameter invalid")
  #approximate likelihood "AR"
  zta<-ARToPacf(stats::ar.burg(z,aic=FALSE,order.max=lag.max)$ar)
  n<-length(z)
  LagRange<-1:lag.max
  if (method=="UBIC"){
    mColNames<-list(c("p", "UBIC-Exact", "UBIC-Approx"))
    PENALTY1 <- log(n) + 2*lchoose(lag.max, 1)
    penalty<-log(n)*(1+LagRange)+2*lchoose(lag.max, LagRange)
  }
  if (method=="EBIC"){
    mColNames<-list(c("p", "EBIC-Exact", "EBIC-Approx"))
    PENALTY1 <- log(n) + 2*G*lchoose(lag.max, 1)
    penalty<-log(n)*(1+LagRange)+2*G*lchoose(lag.max, LagRange)
  }
  if (method=="BICq"){
    mColNames<-list(c("p", "BICq-Exact", "BICq-Approx"))
    PENALTY1 <- log(n) - 2*log(Q/(1-Q))
    penalty<-(1+LagRange)*PENALTY1
  }
  if (method=="BIC"){
    PENALTY1 <- log(n)
    mColNames<-list(c("p", "BIC-Exact", "BIC-Approx"))
    penalty<-(1+LagRange)*PENALTY1
  }
  if (method=="AIC"){
    PENALTY1 <- 2
    mColNames<-list(c("p", "AIC-Exact", "AIC-Approx"))
    penalty<-(1+LagRange)*PENALTY1
  }
  if (method=="GIC"){
    mColNames<-list(c("p", "GIC-Exact", "GIC-Approx"))
    PENALTY1 <- stats::qchisq(p=(1+sqrt(1-4*P))/2, df=1)
    penalty<-(1+LagRange)*PENALTY1
  }
  LagsEntering<-1:lag.max
  LLapprox <- n*log(cumprod(1-zta[LagsEntering]^2))
  AnIC <- LLapprox + penalty
  #
  IndCandidates<-order(AnIC)[1:BestCandidates]
  AnICexact<-numeric(BestCandidates+1)
  AnICexact<-numeric(BestCandidates)
  AnICApprox<-numeric(BestCandidates)
  for (i in 1:BestCandidates){
    p<-LagsEntering[IndCandidates[i]]-1
    AnICApprox[i]<-AnIC[p+1]
    ans<-GetFitARz(z-mean(z), 0:p)
    LL<-ans$loglikelihood
    #mean included in all models, k=p+1
    if (method=="AIC")
      penalty<-2*(p+1)
    if (method=="BIC")
      penalty<-log(n)*(p+1)
    if (method=="BICq")
      penalty<-log(n)*(1+p)-2*((p+1)*log(Q/(1-Q)))
    if (method=="EBIC") #is equivalent to BIC really!
      penalty<-log(n) + 2*G*lchoose(p+1, 1)
    if (method=="UBIC")
      penalty<-log(n) + 2*lchoose(p+1, 1)
    if (method=="GIC")
      penalty<-(1+p)*PENALTY1
    if (method=="BIC")
      penalty<-log(n)*(p+1)
    AnICexact[i]<- -2*LL+penalty
  }
  m<-c(LagsEntering[IndCandidates]-1,AnICexact,AnICApprox)
  m<-matrix(m,ncol=3)
  if (Best==1) {
    m<-m[order(AnICexact),drop=FALSE]
    m<-m[1:Best,drop=FALSE]
  }
  else {
    m<-m[order(AnICexact),]
    m<-m[1:Best,]
  }
  if (Best > 1)
    dimnames(m)<-c(list(1:Best), mColNames)
  if (Best > 1)
    m
  else
    if (is.list(m))
      m[[1]]$p
  else
    as.vector(m[1])
}


