\name{FindingGRGrps}
\alias{FindingGRGrps}

\title{An auxilary function calculating the containing variable index for each predictor (or response) group. It generates the required input \code{GR.grps} matrix when calling the MSGlasso function.}

\description{An auxilary function calculating the containing variable index for each predictor (or response) group. It generates the required input \code{GR.grps} matrix when calling the MSGlasso function.}

\usage{
FindingGRGrps(P, Q, G, R, cmax, G.Starts, G.Ends, R.Starts, R.Ends)
}

\arguments{
  \item{P}{a positive interger indicating number of predictor variables}
  \item{Q}{a positive interger indicating number of response variables}
  \item{G}{a positive interger indicating number of predictor groups}
  \item{R}{a positive interger indicating number of response groups}
  \item{cmax}{a positive interger indicating the max number of variables a single group (either a predictor or response group) contains.}
  \item{G.Starts}{a vector of starting coordinates for the predictor groups.}
  \item{G.Ends}{a vector of ending coordinates for the predictor groups.}
  \item{R.Starts}{a vector of starting coordinates for the response groups.}
  \item{R.Ends}{a vector of ending coordinates for the response groups.}
}

\details{Generates the required input \code{GRgrps} matrix when calling the main MSGlasso function, when user provide the starting and ending coordinates for each of the predictor and response groups. The \code{GRgrps} is a matrix of (g+r) by (cmax+1), with each row starting with variable indicators that row group contains, and followed by 999's till the row is filled.}

\value{
    A list with one components:
    \item{GRgrps }{the \code{GRgrps} matrix generated}
}

\references{
    Y. Li, B. Nan and J. Zhu (2015)
    Multivariate sparse group lasso for the multivariate multiple linear regression with an arbitrary group structure. Biometrics. DOI: 10.1111/biom.12292
}

\author{Yanming Li, Bin Nan, Ji Zhu}

\examples{

###########################################################################
## generating the GR.grps matrix for an overlapping group structure
###########################################################################

P <- 200
Q <- 200
G <- 10
R <- 10
cmax <- 400

GarrStarts <-c(0,20,40,60,80,100,120,140,160,180)
GarrEnds <-c(19,39,59,79,99,119,139,159,179,199)
RarrStarts <-c(0,20,40,60,80,100,120,140,160,180)
RarrEnds <-c(19,39,59,79,99,119,139,159,179,199)

tmp <- FindingGRGrps(P, Q, G, R, cmax, GarrStarts, GarrEnds, RarrStarts, RarrEnds)
GR.grps <- tmp$GRgrps

}

\keyword{ auxilaries }

