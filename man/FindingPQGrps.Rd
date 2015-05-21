\name{FindingPQGrps}
\alias{FindingPQGrps}

\title{An auxilary function calculating the group attribution index for each predictor (or response) variable. It generates the required input \code{PQ.grps} matrix when calling the MSGlasso function.}

\description{An auxilary function calculating the group attribution index for each predictor (or response) variable. It generates the required input \code{PQ.grps} matrix when calling the MSGlasso function.}

\usage{
FindingPQGrps(P, Q, G, R, gmax, G.Starts, G.Ends, R.Starts, R.Ends)
}

\arguments{
  \item{P}{a positive interger indicating number of predictor variables}
  \item{Q}{a positive interger indicating number of response variables}
  \item{G}{a positive interger indicating number of predictor groups}
  \item{R}{a positive interger indicating number of response groups}
  \item{gmax}{a positive interger indicating the max number of different groups a single variable (either a predictor or response variable) belongs to.}
  \item{G.Starts}{a vector of starting coordinates for the predictor groups.}
  \item{G.Ends}{a vector of ending coordinates for the predictor groups.}
  \item{R.Starts}{a vector of starting coordinates for the response groups.}
  \item{R.Ends}{a vector of ending coordinates for the response groups.}
}

\details{Generates the required input \code{PQgrps} matrix when calling the main MSGlasso function, when user provide the starting and ending coordinates for each of the predictor and response groups. The \code{PQgrps} is a matrix of (p+q) by (gmax+1), with each row starting with group indicators that row variable belongs to, and followed by 999's till the row is filled.}

\value{
    A list with one components:
    \item{PQgrps }{the \code{PQgrps} matrix generated}
}

\references{
    Y. Li, B. Nan and J. Zhu (2015)
    Multivariate sparse group lasso for the multivariate multiple linear regression with an arbitrary group structure. Biometrics. DOI: 10.1111/biom.12292
}

\author{Yanming Li, Bin Nan, Ji Zhu}

\examples{

###########################################################################
## generating the PQ.grps matrix for an overlapping group structure
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

}

\keyword{ auxilaries }

