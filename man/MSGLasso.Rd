\name{MSGLasso}
\alias{MSGLasso}

\title{A function to fit the Multivariate Sparse Group Lasso with an arbitrary group structure (MSGLasso)}

\description{A function to fit the Multivariate Sparse Group Lasso with an arbitrary group structure using the mixed coordinate descent algorithm.}

\usage{
MSGLasso(X.m, Y.m, grp.WTs, Pen.L, Pen.G, PQ.grps, GR.grps, grp_Norm0,
 lam1, lam.G, Beta0 = NULL)
}

\arguments{
  \item{X.m}{numeric predictor matrix (n by p): columns correspond to predictor variables and rows correspond to samples. Missing values are not allowed.}
  \item{Y.m}{numeric predictor matrix (n by q): columns correspond to response variables and rows correspond to samples. Missing values are not allowed.}
  \item{grp.WTs}{user specified adaptive group weighting matrix of g by r, for putting different penalization levels on different groups. Missing values are not allowed.}
  \item{Pen.L}{user specified single-entry level penalization indictor mateix of p by q. 1 for being penalized and 0 for not. Missing values are not allowed.}
  \item{Pen.G}{user specified group level penalization indictor mateix of g by r. 1 for being penalized and 0 for not. Missing values are not allowed.}
  \item{PQ.grps}{the group attributing matrix of (p+q) by (gmax+1), where gmax is max number of different groups a single variable belongs to. Each row corresponds to a (predictor or response) varaible, and starts with group indexes the variable belongs to and followed by 999.}
  \item{GR.grps}{the variable attributing matrix of (g+r)*(cmax+1), where cmax is max number of variables a single group contains. Each row corresponds to a (predictor or response) group, and starts with variable indexes the group contains to and followed by 999.}
  \item{grp_Norm0}{a numeric matrix (g by r) containing starting L2 group norm values. Should be calculated from the Beta starting value matrix Beta0.}
  \item{lam1}{lasso panelty parameter scaler.}
  \item{lam.G}{group penalty parameter matrix (g by r).}
  \item{Beta0}{numeric matrix (p by q) containing starting beta values. By defalt use a zero matrix.}
}
\details{Uses the mixed coordinate descent algorithm for fitting the multivariate sparse group lasso in a multivariate-response-multiple-predictor linear regression setting, with an arbitrary group structure on the regression coefficient matrix (Li, Nan and Zhu 2014).}

\value{
    A list with five components:
    \item{Beta }{the estimated regression coefficient matrix (p by q).}
    \item{grpNorm }{the L2 group norm matrix (g by r) of the estimated regression coefficient matrix.}
    \item{E }{residual matrix (n by q).}
    \item{rss.v }{a vector of length q recording the resisual sum square for each of the q responses.}
    \item{rss }{a scaler of overall residual sum of square.}
    \item{iter }{a positive interger recording the number of iterations till convergence.}
}

\references{
    Y. Li, B. Nan and J. Zhu (2015)
    Multivariate sparse group lasso for the multivariate multiple linear regression with an arbitrary group structure. Biometrics. DOI: 10.1111/biom.12292
}

\author{Yanming Li, Bin Nan, Ji Zhu}

\examples{

#####################################################
# Simulate data
#####################################################

set.seed(sample(1:100,1))
G.arr <- c(0,20,20,20,20,20,20,20,20,20,20)

data("Beta.m")

######## generate data set for model fitting

simDataGen<-function(N, Beta, rho, s, G.arr, seed=1){

P<-nrow(Beta)
Q<-ncol(Beta)
gsum<-0
X.m<-NULL

set.seed(seed)

Sig<-matrix(0,P,P)
jstart <-1

for(g in 1:length(G.arr)-1){
X.m<-cbind(X.m, matrix(rnorm(N*G.arr[g+1]),N,G.arr[g+1], byrow=TRUE))

for(i in 2:P){ for(j in jstart: (i-1)){

    Sig[i,j]<-rho^(abs(i-j))

    Sig[j,i]<-Sig[i,j]

}}
jstart <- jstart + G.arr[g+1]
}


diag(Sig)<-1
R<-chol(Sig)

X.m<-X.m\%*\%R

SVsum <-0

for (q in 1:Q){SVsum <-SVsum+var(X.m \%*\% Beta[,q])}
sdr =sqrt(s*SVsum/Q)

E.m <- matrix(rnorm(N*Q,0,sdr),N, Q, byrow=TRUE)

Y.m<-X.m\%*\%Beta+E.m

return(list(X=X.m, Y=Y.m, E=E.m))
}

N <-150

rho=0.5; 
s=4;

Data <- simDataGen(N, Beta.m,rho, s, G.arr, seed=sample(1:100,1))
X.m<-Data$X
Y.m<-Data$Y


############################################################
## fit model for one set of (lam1, lam.G) using example data
############################################################

P <- dim(Beta.m)[1]
Q <- dim(Beta.m)[2]
G <- 10
R <- 10

gmax <- 1
cmax <- 20
GarrStarts <- c(0,20,40,60,80,100,120,140,160,180)
GarrEnds <- c(19,39,59,79,99,119,139,159,179,199)
RarrStarts <- c(0,20,40,60,80,100,120,140,160,180)
RarrEnds <- c(19,39,59,79,99,119,139,159,179,199)

tmp <- FindingPQGrps(P, Q, G, R, gmax, GarrStarts, GarrEnds, RarrStarts, RarrEnds)
PQgrps <- tmp$PQgrps

tmp1 <- Cal_grpWTs(P, Q, G, R, gmax, PQgrps)
grpWTs <- tmp1$grpWTs

tmp2 <- FindingGRGrps(P, Q, G, R, cmax, GarrStarts, GarrEnds, RarrStarts, RarrEnds)
GRgrps <- tmp2$GRgrps

Pen_L <- matrix(rep(1,P*Q),P,Q, byrow=TRUE)
Pen_G <- matrix(rep(1,G*R),G,R, byrow=TRUE)
grp_Norm0 <- matrix(rep(0, G*R), nrow=G, byrow=TRUE)

MSGLassolam1 <- 1.6
MSGLassolamG <- 0.26
MSGLassolamG.m <- matrix(rep(MSGLassolamG, G*R),G,R,byrow=TRUE)
 
system.time(try <-MSGLasso(X.m, Y.m, grpWTs, Pen_L, Pen_G, PQgrps, GRgrps, grp_Norm0, 
   MSGLassolam1, MSGLassolamG.m, Beta0=NULL))

\dontrun{
################################################
## visulizing model fitting results
################################################

########visulizing selection effect using heatmaps

MYplotBW <- function(Beta){
colorNP <- ceiling(abs(max(Beta)))+2
ColorValueP <- colorRampPalette(c("gray50", "black"))(colorNP)
colorNN <- ceiling(abs(min(Beta)))+2
ColorValueN <- colorRampPalette(c("gray50", "white"))(colorNN)
P <-nrow(Beta)
Q <-ncol(Beta)
Xlim <- c(0,2*(P+1))
Ylim <- c(0,2*(Q+1))
plot(0, type="n", xlab="", ylab="", xlim=Xlim, ylim=Ylim, cex.lab=1.0, 
  bty="n", axes=FALSE)
for (p in 1:P){
for (q in 1:Q){
        k0 <- Beta[p,q]
if(k0==0){
      rect(2*(P-p+1)-1,2*(Q-q+1)-1, 2*(P-p+1)+1, 2*(Q-q+1)+1, col="white", border=NA)
}
if(k0>0){
      k <- ceiling(k0)+1
      if(k>2) {k <- k+1}
          rect(2*(P-p+1)-1,2*(Q-q+1)-1, 2*(P-p+1)+1, 2*(Q-q+1)+1, 
               col="black", border=NA)
}
if(k0<0){
      k <- ceiling(abs(k0))+1
      if(k>2) {k <- k+1}
          rect(2*(P-p+1)-1,2*(Q-q+1)-1, 2*(P-p+1)+1, 2*(Q-q+1)+1, 
               col="black", border=NA)
}
}
}

rect(1,1,2*P, 2*Q, lty=2)
}

MYplotBW(try$Beta)

rect(1,1,40,40, lty=2)
rect(41,41,80,80, lty=2)
rect(81,81,120,120, lty=2)
rect(121,121,160,160, lty=2)
rect(161,161,200,200, lty=2)
rect(201,201,240,240, lty=2)
rect(241,241,280,280, lty=2)
rect(281,281,320,320, lty=2)
rect(361,1,400,400, lty=2)

######## visulizing the true Beta matrix

#X11()

MYplotBW(Beta.m)

rect(1,1,40,40, lty=2)
rect(41,41,80,80, lty=2)
rect(81,81,120,120, lty=2)
rect(121,121,160,160, lty=2)
rect(161,161,200,200, lty=2)
rect(201,201,240,240, lty=2)
rect(241,241,280,280, lty=2)
rect(281,281,320,320, lty=2)
rect(361,1,400,400, lty=2)

}
}

\keyword{methods}

