\name{sens.mice}
\alias{sens.mice}
\title{Multivariate Imputation by Chained Equations (Iteration Step for sensitivity analysis).}
\usage{sens.mice(IM, ListMethod=ListMethod, SupPar=SupPar)} 
\arguments{
    \item{IM}{An object of class \code{mids}, typically produces by a previous call to \code{mice()} or \code{mice.mids()}}

    \item{ListMethod}{ is a vector of strings with length \code{ncol(IM[["data"]])}, specifying the column(s) in data which
    has (have) to be imputed with an imputation model accounting for MNAR data. For variables which have to be imputed with
    a modified imputation model, the method has to be specified as \code{"MyFunc"}. For variables which do not have to be
    imputed with a different imputation model, the method has to be be specified as \code{""}.}

    \item{SupPar}{ is a vector of numbers, specifying the supplementary parameters to be added to the imputation model(s).
    For \code{logreg} and \code{polyreg} methods, the supplementary parameters are expressed as odds-ratios (correponding to
    the excess of risk to present the modality of interest for non responders as compared to responders). The value for the
    reference need not be specified. For \code{pmm} and \code{norm} methods, the supplementary parameters are the difference
    between the expected values in responders and non responders.}
}
\description{Takes a \code{mids} object, and produces a new object of class \code{mids} with hypothetical non-ignorability included.}
\details{This function is based on the \code{mice} function's principle, which uses the MICE algorithm. \code{mice} allows
to generate multiple imputations for incomplete multivariate data by Gibbs sampling. The algorithm imputes an incomplete
column (the target column) by generating \sQuote{plausible} synthetic values given other columns in the data. Each
incomplete column must act as a target column, and has its own set of predictors. The default set of predictors for a given
target consists of all others columns in the data. For predictors that are incomplete themselves, the most recently
generated imputations are used to complete the predictors prior to imputation of the target column.

Built-in elementary imputation methods are: 
\code{pmm} Bayesian linear regression (numeric)
\code{norm} Predictive mean matching (numeric)
\code{norm.nob} Linear regression ignoring model error (numeric)
\code{mean} Unconditional mean imputation (numeric)
\code{2l.norm} Two-level normal imputation (numeric)
\code{logreg} Logistic regression (factor, 2 categories)
\code{polyreg} Polytomous (unordered) logistic regression (factor, >=2 categories)
\code{lda} Linear discriminant analysis (factor, >=2 categories)
\code{sample} Random sample from the observed values (any)

This function enables to impute missing values under the hypothesis of MNAR data, for one or more variable(s). A
sensitivity analysis can be performed using the \code{mids} object returned by the function.

It goes through 3 steps:

1. Estimate the parameters of the imputation model under ignorable missing data hypothesis using the function \code{mice}.
The fitted imputation model depends on the type of the variable which contains missing values, i.e. bayesian linear
regression for numeric variables, logistic regression for binary variables, polytomous unordered logistic regression for
categorical variables.

2. Modify the imputation model according to the explored scenario by specifying supplementary parameters as arguments for
the \code{sens.mice} function. Without fitting a mixture model, we appeal to its principle as proposed by Rubin. Indeed, the
addition of this(these) supplementary parameter(s) allows to specify that the distribution of the variable of interest is
different among subjects with missing value and among subjects witouht missing value, conditionnally on all variables
included in the imputation model. The direction and the size of this difference is expressed by the supplementary parameter
in the imputation model.

3. Impute the missing data using the function \code{sens.mice}, resulting into a \code{mids} object that contains m newly
imputed data sets under MNAR hypothesis. For this step, missing values are imputed with the usual MICE algorithm, using the
previously modified imputation model.

The sensitivity analysis can then be performed on the returned \code{mids} object, which contains the newly imputed data,
according to the assumption that values for the variabl(s) imputed with the \code{sens.mice} function come from two
different distributions (for responders and incomplete responders).

Built-in elementary imputation methods are: 

- \code{pmm}: Predictive mean matching (numeric)

- \code{norm}: Bayesian linear regression (numeric) 

- \code{logreg}: Logistic regression (factor, 2 categories) 

- \code{polyreg}: Polytomous (unordered) regression (factor, >= 2 categories)   

For \code{logreg} and \code{polyreg} methods, the supplementary parameters are expressed as odds-ratios (correponding to the
excess of risk to present the modality of interest for non responders as compared to responders). The value for the
reference need not be specified. For \code{pmm} and \code{norm} methods, the supplementary parameters are the difference
between the expected values in responders and non responders. The name associated with the value identifies the target
column(s) of the imputed data matrix. }

\value{Returns an object of class \code{mids} (multiply imputed data set) with
usual components}

\author{Noemie Resseguier, with contributions of Roch Giorgi, 
David Hajage, Yann De Rycke, Xavier Paoletti and Dave Armstrong}

\references{Resseguier, N., Giorgi, R. and Paoletti, X. (submitted) \emph{How to perform a senstivity analysis exploring the impact
of missing not at random data under different sceanrios of non response mechanism with the R software.} 

Rubin, D.B. \emph{Multiple Imputation for Nonresponse in Surveys.} New York: John Wiley & Sons, 1987.

van Buuren, S., Groothuis-Oudshoorn, K. \emph{MICE: Multivariate Imputation by Chained Equations in R.}}
\seealso{\code{\link{mice}}, \code{\link{mids}}}

