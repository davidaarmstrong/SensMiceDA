\name{sens.pool}
\alias{sens.pool}
\title{Pooling of Sensitivity Analysis Results}
\usage{sens.pool(obj, sensData, impData, ...)} 
\arguments{
    \item{obj}{A model estimated on the original, list-wise deleted data.  The model will be re-estimated on the complete data from the sensitivity analysis and the original mice procedure.  Must be a model for which \code{update} will work.}

    \item{sensData}{Output from \code{\link{sens.est}} that includes multiply imputed datasets under a hypothetical non-ignorable mechanism.}
    
    \item{impData}{An object of class \code{mids} resulting from a call to \code{mice} that gives the imputed data assuming ignorability}
    
    \item{...}{other arguments passed down to the model function, currently not implemented}
}
\description{Uses \code{\link{MIcombine}} to pool results from models estimated on the sensitivity imputed datasets.}

\value{Returns a data frame of model estimates for the original imputation under ignorability (using \code{impData}) and under the non-ignorable mechanism (using \code{sensData}).}

\author{Noemie Resseguier, with contributions of Roch Giorgi, 
David Hajage, Yann De Rycke, Xavier Paoletti and Dave Armstrong}

\references{Resseguier, N., Giorgi, R. and Paoletti, X. (submitted) \emph{How to perform a senstivity analysis exploring the impact
of missing not at random data under different sceanrios of non response mechanism with the R software.} 

Rubin, D.B. \emph{Multiple Imputation for Nonresponse in Surveys.} New York: John Wiley & Sons, 1987.

van Buuren, S., Groothuis-Oudshoorn, K. \emph{MICE: Multivariate Imputation by Chained Equations in R.}}
\seealso{\code{\link{mice}}, \code{\link{mids}}, \code{\link{sens.mice}}, \code{\link{sens.est}}}

