\name{MSGLasso.cv}
\alias{MSGLasso.cv}

\title{Fit the MSGLasso for a series sets of tuning parameters and use the k-fold cross validation to select the optimal tunning parameter set.}

\description{Fit the MSGLasso for a series sets of tuning parameters and use the k-fold cross validation to select the optimal tunning parameter set.}

\usage{
MSGLasso.cv(X, Y, grpWTs, Pen.L, Pen.G, PQgrps, GRgrps, lam1.v, lamG.v, 
   fold = 10, seed = 1, Beta.ini = NULL, grp_Norm = NULL)
}

\arguments{
  \item{X}{numeric predictor matrix (n by p): columns correspond to predictor variables and rows correspond to samples. Missing values are not allowed.}
  \item{Y}{numeric predictor matrix (n by q): columns correspond to response variables and rows correspond to samples. Missing values are not allowed.}
  \item{grpWTs}{user specified adaptive group weighting matrix of g by r, for putting different penalization levels on different groups. Missing values are not allowed.}
  \item{Pen.L}{user specified single-entry level penalization indictor mateix of p by q. 1 for being penalized and 0 for not. Missing values are not allowed.}
  \item{Pen.G}{user specified group level penalization indictor mateix of g by r. 1 for being penalized and 0 for not. Missing values are not allowed.}
  \item{PQgrps}{the group attributing matrix of (p+q) by (gmax+1), where gmax is max number of different groups a single variable belongs to. Each row corresponds to a (predictor or response) varaible, and starts with group indexes the variable belongs to and followed by 999.}
  \item{GRgrps}{the variable attributing matrix of (g+r)*(cmax+1), where cmax is max number of variables a single group contains. Each row corresponds to a (predictor or response) group, and starts with variable indexes the group contains to and followed by 999.}
  \item{lam1.v}{lasso panelty parameter scaler.}
  \item{lamG.v}{group penalty parameter matrix (g by r).}
  \item{fold}{a positive integer for the corss validation fold. Default=5. }
  \item{seed}{a numeric scaler, specifying the seed of the random number generator in R for generating cross validation subset for each fold. Default=1.}
  \item{Beta.ini}{a numeric matrix of p by q, specifying the starting values of the input Beta matrix for each fold. Default using the zero matrix.}
  \item{grp_Norm}{a numeric matrix (g by r) containing starting L2 group norm values. Should be calculated from the Beta starting value matrix Beta.ini.}
}

\details{
    Performs a k-fold cross-validation for seaching the optimal tunning parameter associated with the minimal prediction error on a two-dimensional grid.
}

\value{
    A list with two components:
    \item{rss.cv }{a numeric matrix recording the cross validation scores based on the MSGLasso estimators for each pair of (lam1, lamG).}
    \item{lams.c }{a list of tuning parameter pairs corresponding the validation scores in the vectorized rss.cv.}
}

\references{
    Y. Li, B. Nan and J. Zhu (2015)
    Multivariate sparse group lasso for the multivariate multiple linear regression with an arbitrary group structure. Biometrics. DOI: 10.1111/biom.12292
}

\author{Yanming Li, Bin Nan, Ji Zhu}

\examples{

\dontrun{
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

################################################
## cross validation using the example data
################################################
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

lam1.v <- seq(1.0, 1.5, length=6) 
lamG.v <- seq(0.19, 0.25, length=7) 

try.cv<- MSGLasso.cv(X.m, Y.m, grpWTs, Pen_L, Pen_G, PQgrps, GRgrps, 
   lam1.v, lamG.v, fold=5, seed=1)
MSGLassolam1 <- try.cv$lams.c[which.min(as.vector(try.cv$rss.cv))][[1]]$lam1
MSGLassolamG  <- try.cv$lams.c[which.min(as.vector(try.cv$rss.cv))][[1]]$lam3
MSGLassolamG.m <- matrix(rep(MSGLassolamG, G*R),G,R,byrow=TRUE)
 
system.time(try <-MSGLasso(X.m, Y.m, grpWTs, Pen_L, Pen_G, PQgrps, GRgrps, 
    grp_Norm0, MSGLassolam1, MSGLassolamG.m, Beta0=NULL))

}
}

\keyword{methods}

