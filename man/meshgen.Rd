% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/meshgen.R
\name{meshgen}
\alias{meshgen}
\title{Mesh generation}
\usage{
meshgen(rad, zrange, mxbnd, mxnode, bndd, mxn, mzn, mxb, mze, mxe, zcod,
  nzin, nzot, zcalc = TRUE, layer = c(0.175, 0.295, 0.4, 0.5, 0.55))
}
\arguments{
\item{rad}{numeric vectors. Shape parameter vector (<1.0).}

\item{zrange}{numeric. Length of the model.}

\item{mxbnd}{numeric. Maximum bandwidth.}

\item{mxnode}{numeric. Number of nodes.}

\item{bndd}{numeric.}

\item{mxn}{numeric.}

\item{mzn}{numeric.}

\item{mxb}{numeric.}

\item{mze}{numeric.}

\item{mxe}{numeric.}

\item{zcod}{numeric.}

\item{nzin}{numeric.}

\item{nzot}{numeric.}

\item{zcalc}{numeric.}

\item{layer}{numeric vector.}
}
\description{
Mesh generation
}
