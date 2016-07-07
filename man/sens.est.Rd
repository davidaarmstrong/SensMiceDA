\name{sens.est}
\alias{sens.est}
\title{Multivariate Imputation by Chained Equations (Iteration Step for sensitivity analysis).}
\usage{sens.est(mids.obj, vars_vals, digits=2)} 
\arguments{
    \item{mids.obj}{An object of class \code{mids}, typically produces by a previous call to \code{mice()} or \code{mice.mids()}}

    \item{vars_vals}{A named list identifying the variables to be imputed under the non-ignorable (NI) mechanism and the values that represent the differences between responders and non-responders (see \code{\link{sens.mice}} for more details). The names must be as they appear in \code{mids.obj[["data"]]}, but not all variables have to be included, only the one(s) being imputed under the hypothetical NI mechanism. }

    \item{digits}{Number of digits to round the values when they get used to name the elements of the resulting data list.}
}
\description{This is a convenience wrapper function to \code{\link{sens.mice}} that makes specifying the sensitivity analysis slightly easier and provide output that can be used easily with other convenience functions included here.}

\value{Returns a list of objects of class \code{mids} (multiply imputed data set) with usual components}

\author{Noemie Resseguier, with contributions of Roch Giorgi, 
David Hajage, Yann De Rycke, Xavier Paoletti and Dave Armstrong}

\references{Resseguier, N., Giorgi, R. and Paoletti, X. (submitted) \emph{How to perform a senstivity analysis exploring the impact
of missing not at random data under different sceanrios of non response mechanism with the R software.} 

Rubin, D.B. \emph{Multiple Imputation for Nonresponse in Surveys.} New York: John Wiley & Sons, 1987.

van Buuren, S., Groothuis-Oudshoorn, K. \emph{MICE: Multivariate Imputation by Chained Equations in R.}}
\seealso{\code{\link{mice}}, \code{\link{mids}}, \code{\link{sens.mice}}}

