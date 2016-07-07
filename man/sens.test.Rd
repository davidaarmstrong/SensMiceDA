\name{sens.test}
\alias{sens.test}
\title{Analysis of Deviance for Pooled Sensitivity Analysis Results}
\usage{sens.test(obj, var, sensData, impData, digits=3, ...)} 
\arguments{
    \item{obj}{A model estimated on the original, list-wise deleted data.  The model will be re-estimated on the complete data from the sensitivity analysis and the original mice procedure.  Must be a model for which \code{update} will work.}

    \item{var}{Character string identifying the variable for which the test is to be done.}

    \item{sensData}{Output from \code{\link{sens.est}} that includes multiply imputed datasets under a hypothetical non-ignorable mechanism.}
    
    \item{impData}{An object of class \code{mids} resulting from a call to \code{mice} that gives the imputed data assuming ignorability}

    \item{digits}{Number of digits to be printed in output}
    
    \item{...}{other arguments passed down to the model function, currently not implemented}
}
\description{Uses \code{\link{MIcombine}} to pool results from models estimated on the sensitivity imputed datasets.}

\details{The function performs an incremental chi-squared test for the exclusion of a variable (\code{var}).  The deviance is calculated for both the full and restricted models estimated on each completed dataset.  The chi-squared statistic is calculated for each difference in deviance and then the p-values are averaged across all of the \code{m} iterations of the imputations according to Resseguier et a.}

\value{Returns a matrix of average chi-squared statistics and average p-values}

\author{Noemie Resseguier, with contributions of Roch Giorgi, 
David Hajage, Yann De Rycke, Xavier Paoletti and Dave Armstrong}

\references{Resseguier, N., Giorgi, R. and Paoletti, X. (submitted) \emph{How to perform a senstivity analysis exploring the impact
of missing not at random data under different sceanrios of non response mechanism with the R software.} 

Rubin, D.B. \emph{Multiple Imputation for Nonresponse in Surveys.} New York: John Wiley & Sons, 1987.

van Buuren, S., Groothuis-Oudshoorn, K. \emph{MICE: Multivariate Imputation by Chained Equations in R.}}
\seealso{\code{\link{mice}}, \code{\link{mids}}, \code{\link{sens.mice}}, \code{\link{sens.est}}}

