\name{Cal_grpWTs}
\alias{Cal_grpWTs}

\title{An auxilary function calculating the group weighting matrix \code{grpWTs} required when calling the MSGlasso function.}

\description{An auxilary function calculating the group weighting matrix \code{grpWTs} required when calling the MSGlasso function.}

\usage{
Cal_grpWTs(P, Q, G, R, gmax, PQ.grps)
}

\arguments{
  \item{P}{a positive interger indicating number of predictor variables}
  \item{Q}{a positive interger indicating number of response variables}
  \item{G}{a positive interger indicating number of predictor groups}
  \item{R}{a positive interger indicating number of response groups}
  \item{gmax}{a positive interger indicating the max number of different groups a single variable (either a predictor or response variable) belongs to.}
  \item{PQ.grps}{a matrix of (p+q) by (gmax+1), with each row starting with group indicators that row variable belongs to, and followed by 999's till the row is filled.}
}

\details{Generates the required input group weighting matrix \code{grpWTs} when calling the main MSGlasso function. The \code{grpWTs} is a g by r matrix containing the adaptive weighting scores for each group. \code{MSGLasso.grpWTs} use the square root of the group size (number of entries the group contains) as the weight for each group.}

\value{
    A list with one components:
    \item{grpWTs }{the \code{grpWTs} matrix generated}
}

\references{
    Y. Li, B. Nan and J. Zhu (2015)
    Multivariate sparse group lasso for the multivariate multiple linear regression with an arbitrary group structure. Biometrics. DOI: 10.1111/biom.12292
}

\author{Yanming Li, Bin Nan, Ji Zhu}

\examples{

###########################################################################
## generating the grp.WTs matrix for an overlapping group structure
###########################################################################

P <- 200
Q <- 200
G <- 10
R <- 10
gmax <- 1

GarrStarts <-c(0,20,40,60,80,100,120,140,160,180)
GarrEnds <-c(19,39,59,79,99,119,139,159,179,199)
RarrStarts <-c(0,20,40,60,80,100,120,140,160,180)
RarrEnds <-c(19,39,59,79,99,119,139,159,179,199)

tmp <- FindingPQGrps(P, Q, G, R, gmax, GarrStarts, GarrEnds, RarrStarts, RarrEnds)
PQ.grps <- tmp$PQgrps

tmp1 <- Cal_grpWTs(P, Q, G, R, gmax, PQ.grps)
grp.WTs <- tmp1$grpWTs
}

\keyword{ auxilaries }


