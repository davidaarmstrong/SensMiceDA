##' Takes a \code{mids} object, and produces a new object of class \code{mids}
##'
##' This function is based on the \code{mice} function's principle, which uses the MICE algorithm. \code{mice} 
##' allows to generate multiple imputations for incomplete multivariate data by Gibbs sampling. The algorithm
##' imputes an incomplete colum (the target column) by generating 'plausible' synthetic values given other columns in the data.
##' Each incomplete column must act as a target column, and has its own set of predictors. The default set of predictors
##' for a given target consists of all others colums in the data. For predictors that are incomplete themselves, the most recently
##' generated imputations are used to complete the predictors prior to imputation of the target column. 
##'
##' Built-in elementary imputation methods are: 
##' \code{pmm} Bayesian linear regression (numeric)
##' \code{norm} Predictive mean matching (numeric)
##' \code{norm.nob} Linear regression ignoring model error (numeric)
##' \code{mean} Unconditional mean imputation (numeric)
##' \code{2l.norm} Two-level normal imputation (numeric)
##' \code{logreg} Logistic regression (factor, 2 categories)
##' \code{polyreg} Polytomous (unordered) logistic regression (factor, >=2 categories)
##' \code{lda} Linear discriminant analysis (factor, >=2 categories)
##' \code{sample} Random sample from the observed values (any)
##'
##' This function enables to impute missing values under the hypothesis of MNAR data, for one or more variable(s). A
##' sensitivity analysis can be performed using the \code{mids} object returned by the function.
##'
##' It goes through 3 steps:
##'
##' - 1. Estimate the parameters of the imputation model under ignorable missing data hypothesis using the function \code{mice}.
##' The fitted imputation model depends on the type of the variable which contains missing values, i.e. bayesian linear regression
##' for numeric variables, logistic regression for binary variables, polytomous unordered logistic regression for categorical
##' variables.
##'
##' - 2. Modify the imputation model according to the explored scenario by specifying supplementary parameters as arguments for the \code{sens.mice} function.
##' Without fitting a mixture model, we appeal to its principle as proposed by Rubin. Indeed, the addition of this(these) supplementary parameter(s) allows to specify
##' that the distribution of the variable of interest is different among subjects with missing value and among subjects witouht missing value, 
##' conditionnally on all variables included in the imputation model. The direction and the size of this difference is expressed by the
##' supplementary parameter in the imputation model.
##'
##' - 3. Impute the missing data using the function \code{sens.mice}, resulting into a \code{mids} object that contains m newly imputed data sets under MNAR hypothesis.
##' For this step, missing values are imputed with the usual MICE algorithm, using the previously modified imputation model.
##'
##' The sensitivity analysis can then be performed on the returned \code{mids} object, which contains the newly imputed data, according to
##' the assumption that values for the variabl(s) imputed with the \code{sens.mice} function come from two different distributions 
##' (for responders and incomplete responders).
##'
##' Built-in elementary imputation methods are: 
##'
##'     - \code{pmm}: Predictive mean matching (numeric)
##'    
##'     - \code{norm}: Bayesian linear regression (numeric) 
##'
##'     - \code{logreg}: Logistic regression (factor, 2 categories) 
##'
##'     - \code{polyreg}: Polytomous (unordered) regression (factor, >= 2 categories)   
##'
##' \code{ListMethod} is a vector of strings with length \code{ncol(IM$data)}, specifying 
##' the column(s) in data which has (have) to be imputed with an imputation model accounting for MNAR data. 
##' For variables which have to be imputed with a modified imputation model, 
##' the method has to be specified as \code{"MyFunc"}.
##' For variables which do not have to be imputed with a different imputation
##' model, the method has to be be specified as \code{""}.
##'
##' \code{SupPar} is a vector of numbers, specifying the supplementary parameters 
##' to be added to the imputation model(s). For \code{logreg} and \code{polyreg} methods,  
##' the supplementary parameters are expressed as odds-ratios (correponding to the excess of risk 
##' to present the modality of interest for non responders as compared to responders). 
##' The value for the reference need not be specified.
##' For \code{pmm} and \code{norm} methods, the supplementary parameters
##' are the difference between the expected values in responders and non responders.
##'
##' @title Multivariate Imputation by Chained Equations (Iteration Step for sensitivity analysis).
##' @export
##' @param IM An object of class \code{mids}, typically produces by a previous call
##' to \code{mice()} or \code{mice.mids()}
##' @param ListMethod A vector of strings with length \code{ncol(IM$data)}, specifying 
##' the colum(s) in data which has (have) to be imputed with a different imputation model
##' @param SupPar A vector of numbers, specifying the supplementary parameters 
##' to be added to the imputation model(s)
##' 
##' @return Returns an object of class \code{mids} (multiply imputed data set) with
##' usual components
##' @author Noemie Resseguier, with contributions of Roch Giorgi, 
##' David Hajage, Yann De Rycke and Xavier Paoletti
##' @references 
##' Resseguier, N., Giorgi, R. and Paoletti, X. (submitted) \emph{How to perform a senstivity analysis exploring the impact
##' of missing not at random data under different sceanrios of non response mechanism with the R software.} 
##' 
##' Rubin, D.B. \emph{Multiple Imputation for Nonresponse in Surveys.} New York: John Wiley & Sons, 1987.
##'
##' van Buuren, S., Groothuis-Oudshoorn, K. \emph{MICE: Multivariate Imputation by Chained Equations in R.}
##'
##' @seealso \code{\link{mice}}, \code{\link{mids}}
##' @examples 
##'
##' ###### In a descriptive context : data set : popmis
##'
##' # Do multiple imputation 
##' data(popmis)
##' MYpopmis <- popmis[, c("popular", "sex", "texp", "teachpop")]
##' IM <- mice(MYpopmis, method=c("pmm", "", "", ""), seed=13, maxit=5, m=5)
##'
##' # Describe popular
##'
##' # Construction of the completed imputed data sets
##' temp1 <- complete(IM, action=1)
##' temp2 <- complete(IM, action=2)
##' temp3 <- complete(IM, action=3)
##' temp4 <- complete(IM, action=4)
##' temp5 <- complete(IM, action=5)
##'
##' # Description of the popularity score
##' library(mitools)
##' miData <- imputationList(list(temp1, temp2, temp3, temp4, temp5))
##' MEAN <- with(miData, mean(popular))
##' VAR <- with(miData, var(popular))
##' MIcombine(MEAN, VAR)
##'
##' # Imputation with a supplementary parameter -0.5 on the popular variable
##' IMHyp0.5 <- sens.mice(IM, ListMethod = c("MyFunc", "", "", ""), SupPar = c(-0.5))
##'
##' # Describe popular
##'
##' # Construction of the completed imputed data sets
##' temp1Hyp0.5 <- complete(IMHyp0.5, action=1)
##' temp2Hyp0.5 <- complete(IMHyp0.5, action=2)
##' temp3Hyp0.5 <- complete(IMHyp0.5, action=3)
##' temp4Hyp0.5 <- complete(IMHyp0.5, action=4)
##' temp5Hyp0.5 <- complete(IMHyp0.5, action=5)
##'
##' # Description of the popularity score
##' library(mitools)
##' miDataHyp0.5 <- imputationList(list(temp1Hyp0.5, temp2Hyp0.5, temp3Hyp0.5,
##' temp4Hyp0.5, temp5Hyp0.5))
##' MEANHyp0.5 <- with(miDataHyp0.5, mean(popular))
##' VARHyp0.5 <- with(miDataHyp0.5, var(popular))
##' MIcombine(MEANHyp0.5, VARHyp0.5)
##'
##' ###### In a causal context : data set CHAIN
##'
##' library(mi)
##' data(CHAIN)
##'
##' MyData <- CHAIN[apply(CHAIN, 1, function(x) !all(is.na(x))), ]                                                                                         
##' 
##' ### VL2cat = binary variable (<400 c/mL vs >= 400 c/mL)
##' MyData$VL2cat <- NA
##' MyData$VL2cat <- ifelse((MyData$h39b.W1 < log(400) & !is.na(MyData$h39b.W1)),
##' 0, MyData$VL2cat)
##' MyData$VL2cat <- ifelse((MyData$h39b.W1 >= log(400) & !is.na(MyData$h39b.W1)),
##' 1, MyData$VL2cat)
##'
##' # Do multiple imputation
##'
##' MyDataVL2cat <- MyData[, c("age.W1", "c28.W1", "pcs.W1", "mcs37.W1", "b05.W1",
##' "haartadhere.W1", "VL2cat")]
##' MyDataVL2cat$mcs37.W1 <- as.factor(MyDataVL2cat$mcs37.W1)
##' MyDataVL2cat$b05.W1 <- as.factor(MyDataVL2cat$b05.W1)
##' MyDataVL2cat$haartadhere.W1 <- as.factor(MyDataVL2cat$haartadhere.W1)
##' MyDataVL2cat$VL2cat <- as.factor(MyDataVL2cat$VL2cat)
##' MyDataVL2catMI <- mice(MyDataVL2cat, method = c("pmm", "pmm", "pmm", "logreg",
##' "polyreg", "polyreg", "logreg"), seed = 13)
##' 
##' # Construction of the completed imputed data sets
##' mi1VL2cat <- complete(MyDataVL2catMI, action=1)
##' mi2VL2cat <- complete(MyDataVL2catMI, action=2)
##' mi3VL2cat <- complete(MyDataVL2catMI, action=3)
##' mi4VL2cat <- complete(MyDataVL2catMI, action=4)
##' mi5VL2cat <- complete(MyDataVL2catMI, action=5)
##' miDataVL2cat <- imputationList(list(mi1VL2cat, mi2VL2cat, mi3VL2cat,
##' mi4VL2cat, mi5VL2cat))
##'
##' # Description of VL2cat and fit of the model
##' library(mitools)
##' tablesVL2cat <- with(miDataVL2cat, prop.table(table(VL2cat)))
##'
##' fitVL2cat <- with(miDataVL2cat, glm(mcs37.W1 ~ as.numeric(age.W1) +  
##' as.numeric(pcs.W1) + as.factor(haartadhere.W1) + as.factor(VL2cat) +                                                
##' as.numeric(c28.W1) + as.factor(b05.W1), family=binomial))                                              
##'                                                
##' coefsVL2cat <- MIextract(fitVL2cat, fun=coef)
##' varsVL2cat <- MIextract(fitVL2cat, fun=vcov)
##' resVL2cat <- summary(MIcombine(coefsVL2cat, varsVL2cat))
##' devfitVL2cat <- MIextract(fitVL2cat, fun=deviance)
##'
##' fitVL2catPVAL <- with(miDataVL2cat, glm(mcs37.W1 ~ as.numeric(age.W1) +
##' as.numeric(pcs.W1) + as.factor(haartadhere.W1) +  as.numeric(c28.W1) +                                                 
##' as.factor(b05.W1), family=binomial))                                                                                                    
##' devfitVL2catPVAL <- MIextract(fitVL2catPVAL, fun=deviance)
##' 
##' pvalVL2cat <- mean(1-pchisq(unlist(devfitVL2catPVAL)-unlist(devfitVL2cat), 1))
##'
##' # Imputation with a supplementary parameter 1.2 on the VL2cat variable
##' MyDataVL2catMI.SCEN1 <- sens.mice(MyDataVL2catMI, ListMethod = c("", "", "", 
##' "", "", "", "MyFunc"), SupPar = c(1.2))
##' 
##' # Construction of the completed imputed data sets
##' mi1.SCEN1VL2cat <- complete(MyDataVL2catMI.SCEN1, action=1)
##' mi2.SCEN1VL2cat <- complete(MyDataVL2catMI.SCEN1, action=2)
##' mi3.SCEN1VL2cat <- complete(MyDataVL2catMI.SCEN1, action=3)
##' mi4.SCEN1VL2cat <- complete(MyDataVL2catMI.SCEN1, action=4)
##' mi5.SCEN1VL2cat <- complete(MyDataVL2catMI.SCEN1, action=5)
##' miData.SCEN1VL2cat <- imputationList(list(mi1.SCEN1VL2cat, mi2.SCEN1VL2cat,
##' mi3.SCEN1VL2cat, mi4.SCEN1VL2cat, mi5.SCEN1VL2cat))
##'
##' # Description of VL2cat and fit of the model 
##' tablesSCEN1VL2cat <- with(miData.SCEN1VL2cat, prop.table(table(VL2cat)))
##' 
##' fitVL2cat.SCEN1 <- with(miData.SCEN1VL2cat, glm(mcs37.W1 ~ as.numeric(age.W1) + 
##' as.numeric(pcs.W1) + as.factor(haartadhere.W1) + as.factor(VL2cat) +
##' as.numeric(c28.W1) + as.factor(b05.W1), family=binomial))                                                           
##'                                                           
##' coefs.SCEN1VL2cat <- MIextract(fitVL2cat.SCEN1, fun=coef)
##' vars.SCEN1VL2cat <- MIextract(fitVL2cat.SCEN1, fun=vcov)
##' res.SCEN1VL2cat <- summary(MIcombine(coefs.SCEN1VL2cat, vars.SCEN1VL2cat))
##' devfitVL2cat.SCEN1 <- MIextract(fitVL2cat.SCEN1, fun=deviance)
##'
##' fitVL2catPVAL.SCEN1 <- with(miData.SCEN1VL2cat, glm(mcs37.W1 ~ as.numeric(age.W1) + 
##' as.numeric(pcs.W1) + as.factor(haartadhere.W1) + as.numeric(c28.W1) +
##' as.factor(b05.W1), family=binomial))                                                               
##'                                                                
##' devfitVL2catPVAL.SCEN1 <- MIextract(fitVL2catPVAL.SCEN1, fun=deviance)
##' 
##' pvalVL2cat.SCEN1 <- mean(1-pchisq(unlist(devfitVL2catPVAL.SCEN1)-
##' unlist(devfitVL2cat.SCEN1), 1))






sens.mice <- function(IM, ListMethod = ListMethod, SupPar = SupPar){
if(length(ListMethod) > length(names(IM$data))){
  stop("You have specified too much new methods to be applied.")
} 
if(length(ListMethod) < length(names(IM$data))){
  stop("You have not specified enough new methods to be applied.")
}
cpt <- 0
for(i in 1:length(ListMethod)){
  if(ListMethod[i]=="MyFunc"){
    if(IM$method[i] == "norm" | IM$method[i] == "pmm" | IM$method[i] == "logreg"){
    cpt <- cpt + 1
    }
    if(IM$method[i] == "polyreg"){
    cpt <- cpt + dim(table(IM$data[i]))-1
    }
  }
  if(ListMethod[i]==""){
    cpt <- cpt
  }
}
if(length(SupPar) > cpt){
  stop("You have specified too much supplementary parameters to be applied.")
} 
if(length(SupPar) < cpt){
  stop("You have not specified enough supplementary parameters to be applied.")
} 
for(i in 1:length(ListMethod)){
  if(ListMethod[i] != "MyFunc" & ListMethod[i] != ""){
    stop("Values available for ListMethod are ''MyFunc'' and '' ''.")
  }
  if(ListMethod[i]=="Myfunc" & (IM$method[i] =="norm.nob")){
    stop("norm.nob is not an available method for the function sens.mice.")
  }
  if(ListMethod[i]=="Myfunc" & (IM$method[i] =="mean")){
    stop("mean is not an available method for the function sens.mice.")
  }
  if(ListMethod[i]=="Myfunc" & (IM$method[i] =="2l.norm")){
    stop("2l.norm is not an available method for the function sens.mice.")
  }    
  if(ListMethod[i]=="Myfunc" & (IM$method[i] =="lda")){
    stop("lda is not an available method for the function sens.mice.")
  }
  if(ListMethod[i]=="Myfunc" & (IM$method[i] =="sample")){
    stop("sample is not an available method for the function sens.mice.")
  }  
} 
j <- 0
cpt <- 0
IMinit <- IM
MyMethod <- IM$method
listvar <- names(IM$data)
SumPr <- matrix(NA, nrow=length(ListMethod[ListMethod=="MyFunc"]), ncol=3)
for(ii in 1:length(listvar)){
  if(ListMethod[ii]=="MyFunc"){
    j <- j + 1
    cpt <- cpt + 1  
    if(MyMethod[ii]=="pmm" | MyMethod[ii]=="norm" | MyMethod[ii]=="logreg"){
      SumPrtemp <- c(listvar[ii], MyMethod[ii], SupPar[j])
      cat(SumPrtemp)   
    }
    if(MyMethod[ii]=="logreg"){
      if(SupPar[j] < 0){
        stop("Value for odds ratio can't be negative.") 
      }
    }
    if(MyMethod[ii]=="polyreg"){
      tempSupPar <- c(SupPar[j : (j + dim(table(IM$data[ii]))-2)])
      SumPrtemp <- c(listvar[ii], MyMethod[ii], tempSupPar)
      cat(SumPrtemp)
      for(m in 1:length(tempSupPar)){
        if(tempSupPar[m] < 0){
          stop("Value for odds ratio can't be negative.")  
        }
      }
    } 
    IMtemp <- IM
    IMtemp$pad$method <- c(rep("", length(names(IM$data))), rep("dummy", length(IM$pad$method) - length(names(IM$data))))
    IMtemp$pad$method[ii] <- "MyFunc"  
    laps <- lapply(1:IM$m, function(x)complete(IM, x))
    temp <- sapply(laps, function(x)mice.impute.MyFunc(IM$data[,ii], !is.na(IM$data[,ii]), 
        model.matrix(~., data=x[,-which(names(x) == listvar[ii])]), SupPar, MyMethod, ii, j))
    IMinit$imp[[listvar[ii]]] <- temp
  if(MyMethod[ii]=="polyreg"){
      j <- (j + dim(table(IM$data[ii]))-2)
  }
  cat("\n")
  SumPr[cpt, 1] <- SumPrtemp[1]
  SumPr[cpt, 2] <- SumPrtemp[2]
  SumPr[cpt, 3] <- SumPrtemp[3]
  if(length(SumPrtemp) > 3){
    temp <- SumPrtemp[3]
    for(l in 4 : length(SumPrtemp)){
      temp <- paste(temp, SumPrtemp[l], sep=" ; ")  
    }
    SumPr[cpt, 3] <- temp
  }         
  }
  if(ListMethod[ii]==""){
    IMinit$imp[listvar[ii]] <- IMinit$imp[listvar[ii]]
  }
}
dimnames(SumPr)[[2]] <- c("Variable", "Method", "SupPar") 
cat("Summary :")
cat("\n")
cat("\n")
print(SumPr)
IMfinal <- IMinit
}




##' Internal function for sens.mice
##'
##' @title Internal function for sens.mice
##' @export
##' @keywords internal
mice.impute.MyFunc <- function(y, ry, x, suppar , Mymethod, i, j){
  if(Mymethod[i]=="pmm"){
    x <- cbind(1, as.matrix(x))
    parm <- .norm.draw(y, ry, x)
    parm$beta[1] <- parm$beta[1] + suppar[j]
    yhatobs <- x[ry, ] %*% parm$coef
    yhatmis <- x[!ry, ] %*% parm$beta
    return(apply(as.array(yhatmis), 1, .pmm.match, yhat = yhatobs,
        y = y[ry]))
  }
  if(Mymethod[i]=="norm"){
    x <- cbind(1, as.matrix(x))
    parm <- .norm.draw(y, ry, x)
    parm$beta[1] <- parm$beta[1] + suppar[j]
    return(x[!ry, ] %*% parm$beta + rnorm(sum(!ry)) * parm$sigma)
  }
  if(Mymethod[i]=="logreg"){
    aug <- augment(y, ry, x)
    x <- as.matrix(aug$x)
    y <- aug$y
    ry <- aug$ry
    w <- aug$w
    suppressWarnings(fit <- glm.fit(x[ry, ], y[ry], family = binomial(link = logit),
        weights = w[ry]))
    fit.sum <- summary.glm(fit)
    beta <- coef(fit)
    beta[1] <- beta[1] + log(suppar[j])
    rv <- t(chol(fit.sum$cov.unscaled))
    beta.star <- beta + rv %*% rnorm(ncol(rv))
    p <- 1/(1 + exp(-(x[!ry, ] %*% beta.star)))
    vec <- (runif(nrow(p)) <= p)
    vec[vec] <- 1
    if (is.factor(y)) {
        vec <- factor(vec, c(0, 1), levels(y))
    }
    return(vec)
  }
  if(Mymethod[i]=="polyreg"){
    x <- as.matrix(x)
    aug <- augment(y, ry, x)
    x <- aug$x
    y <- aug$y
    ry <- aug$ry
    w <- aug$w   
## check whether this works instead of the assign
    tmpData <- cbind.data.frame(y, x)
    fit <- multinom(formula(tmpData), data = tmpData[ry, ], weights = w[ry],
        maxit = 200, trace = FALSE) 
    temp <- matrix(fit$wts, nrow=nlevels(y), byrow=T)  
    for(k in 2:nlevels(y)){
      temp[k, 2] <- temp[k, 2] + log(suppar[j])
      j <- j + 1
    }       
    temp <- t(temp)
    fit$wts <- c(temp[, ])
    post <- predict(fit, tmpData[!ry, ], type = "probs")
    if (sum(!ry) == 1)
        post <- matrix(post, nrow = 1, ncol = length(post))
    fy <- as.factor(y)
    nc <- length(levels(fy))
    un <- rep(runif(sum(!ry)), each = nc)
    if (is.vector(post))
        post <- matrix(c(1 - post, post), ncol = 2)
    draws <- un > apply(post, 1, cumsum)
    idx <- 1 + apply(draws, 2, sum)
    return(levels(fy)[idx])     
  }
}


sens.est <- function(mids.obj, vars_vals, digits=2){
    if(!all(names(vars_vals) %in% names(mids.obj$data))){
        stop("Some variables not in imputed data set")
    }
    if(any(sapply(vars_vals, is.matrix))){
        nums <- lapply(1:length(vars_vals), function(x)1:length(vars_vals))
        eg.nums <- do.call("expand.grid", nums)
        eg <-  cn <- NULL
        for(i in 1:ncol(eg.nums)){
            if(is.matrix(vars_vals[[i]])){
                eg <- cbind(eg, vars_vals[[i]][eg.nums[,i], , drop=F])
                cn <- c(cn, paste(names(vars_vals)[i], 1:ncol(vars_vals[[i]]), sep=""))
            }
            else{
                eg <- cbind(eg, vars_vals[[i]][eg.nums[,i]])
                cn <- c(cn, names(vars_vals)[i])
            }
        }
        colnames(eg) <- cn
    }
    else{
        eg <- do.call(expand.grid, vars_vals)[,,drop=FALSE]
        if(!is.matrix(eg)){
            eg <- matrix(eg[[1]], ncol=1)[,,drop=F]
        }    
        colnames(eg) <- names(vars_vals)
    }
    out <- list()
    for(i in 1:nrow(eg)){                                             
        ListMethod <- ifelse(names(mids.obj$data) %in% names(vars_vals), "MyFunc", "")
        SupPar <- c(unlist(eg[i, ,drop=F]))
        out[[i]] <- sens.mice(mids.obj, ListMethod, SupPar)
    }
    nms <- NULL
    for(i in 1:ncol(eg)){
        nms <- cbind(nms, paste(colnames(eg)[i], round(eg[,i], digits), sep=": "))
    }
    names(out) <- apply(nms, 1, paste, collapse=", ")
    out
}

sens.pool <- function(obj, sensData, impData, ...){
    nconds <- length(sensData)
    condlist <- list()
    j <- 1
    for(l in 1:nconds){
        condlist[[j]] <- list()
        for(i in 1:sensData[[j]]$m){
            condlist[[j]][[i]] <- complete(sensData[[j]], i)
        }
        j <- j+1
    }
    condlist[[j]] <- list()
    for(i in 1:impData$m){
        condlist[[(j)]][[i]] <- complete(impData, i)
    }
   { if(length(names(sensData)) > 0){
        names(condlist) <- c(names(sensData), "mice")
    }                                                
    else{
        names(condlist) <- c(as.character(1:nconds), "mice")
    }}
    cond.mods <- list()
    for(i in 1:length(condlist)){
        cond.mods[[i]] <- list()
        for(j in 1:length(condlist[[i]])){
            tmp <- obj
            attr(tmp$terms, ".Environment") <- environment()
            cond.mods[[i]][[j]] <- update(tmp, . ~ ., data=condlist[[i]][[j]])
        }
    }
    comb.mods <- invisible(lapply(cond.mods, MIcombine))
    sum.mods <- invisible(lapply(comb.mods, summary))
    names(comb.mods) <- names(sum.mods) <- names(condlist) 
    sub <- as.data.frame(rbind(do.call(rbind, sum.mods)))
    varnames <- gsub("mice.", "", grep("^mice", rownames(sub), value=T) , fixed=T)
    sub$vars <- as.factor(rep(varnames, length(condlist)))
    sub$conds <- factor(c(rep(names(condlist), each = length(varnames))), levels=names(condlist))
    rownames(sub) <- NULL
    class(sub) <- c("sens.pool", "data.frame")
    sub
}         

plot.sens.pool <- function(x, ...){
    p <- xyplot(results ~ conds | vars , 
        data=x, scales=list(x=list(rot=45), y=list(relation="free", rot=90)), pch=16, col="black", 
        lower=x[["(lower"]], upper=x[["upper)"]], 
        xlab = "", ylab = "Coefficients with 95% Confidence Intervals",
        prepanel=function (x, y, subscripts, lower, upper,...){
            list(ylim = range(c(lower[subscripts], upper[subscripts]), finite = TRUE))},
         panel=function(x,y,lower,upper,subscripts,...){
            panel.xyplot(x, y, ...)
            panel.segments(x, lower[subscripts], x, upper[subscripts], ...)  
            panel.abline(h=0, lty=3)
        })
    p
}                                             

sens.test <- function(obj, var, sensData, impData, digits=3, ...){
    nconds <- length(sensData)
    condlist <- list()
    j <- 1
    for(l in 1:nconds){
        condlist[[j]] <- list()
        for(i in 1:sensData[[j]]$m){
            condlist[[j]][[i]] <- complete(sensData[[j]], i)
        }
        j <- j+1
    }
    condlist[[j]] <- list()
    for(i in 1:impData$m){
        condlist[[(j)]][[i]] <- complete(impData, i)
    }
    if(length(names(sensData)) > 0){
        names(condlist) <- c(names(sensData), "mice")
    }                                                
    else{
        names(condlist) <- c(as.character(1:nconds), "mice")
    }
    cond.mods <- list()
    for(i in 1:length(condlist)){
        cond.mods[[i]] <- list()
        for(j in 1:length(condlist[[i]])){
            tmp <- obj
            attr(tmp$terms, ".Environment") <- environment()
            cond.mods[[i]][[j]] <- update(tmp, . ~ ., data=condlist[[i]][[j]])
        }
    }
    restr.mods <- list()
    for(i in 1:length(condlist)){
        restr.mods[[i]] <- list()
        for(j in 1:length(condlist[[i]])){
            tmp <- obj
            attr(tmp$terms, ".Environment") <- environment()
            restr.mods[[i]][[j]] <- update(tmp, paste0(". ~ .-", var), data=condlist[[i]][[j]])
        }
    }
    full.devs <- lapply(cond.mods, function(x)sapply(x, deviance))
    restr.devs <- lapply(restr.mods, function(x)sapply(x, deviance))
    names(full.devs) <- names(restr.devs) <- names(condlist)
    df.diff <- df.residual(restr.mods[[1]][[1]]) - df.residual(cond.mods[[1]][[1]])
    chisqs <- lapply(1:length(full.devs), function(i)restr.devs[[i]] - full.devs[[i]])
    out <- sapply(1:length(full.devs), function(i)mean(1-pchisq(chisqs[[i]], df.diff)))
    fmt <- paste0("%.", digits, "f")
    out <- cbind(sprintf(fmt, sapply(chisqs, mean)), sprintf(fmt, out))
    colnames(out) <- c("Average X2", "p-value")
    rownames(out) <- names(full.devs)
    cat("Test for exclusion of ", var, "(", df.diff, " degrees of freedom)\n")
    print(noquote(out))
}


augment <- function (y, ry, x, maxcat = 50, ...) {
    # augment comes from mice v. 2.25.  It was not exported from
    # the namespace so could not be imported from that package here
    # I copied the function in its entirety in the interest of continued compatability
    icod <- sort(unique(unclass(y)))
    k <- length(icod)
    if (k > maxcat) 
        stop(paste("Maximum number of categories (", maxcat, 
            ") exceeded", sep = ""))
    p <- ncol(x)
    if (p == 0) 
        return(list(y = y, ry = ry, x = x, w = rep(1, length(y))))
    if (sum(!ry) == 1) 
        return(list(y = y, ry = ry, x = x, w = rep(1, length(y))))
    mean <- apply(x, 2, mean)
    sd <- sqrt(apply(x, 2, var))
    minx <- apply(x, 2, min)
    maxx <- apply(x, 2, max)
    nr <- 2 * p * k
    a <- matrix(mean, nrow = nr, ncol = p, byrow = TRUE)
    b <- matrix(rep(c(rep(c(0.5, -0.5), k), rep(0, nr)), length = nr * 
        p), nrow = nr, ncol = p, byrow = FALSE)
    c <- matrix(sd, nrow = nr, ncol = p, byrow = TRUE)
    d <- a + b * c
    d <- pmax(matrix(minx, nrow = nr, ncol = p, byrow = TRUE), 
        d)
    d <- pmin(matrix(maxx, nrow = nr, ncol = p, byrow = TRUE), 
        d)
    e <- rep(rep(icod, each = 2), p)
    dimnames(d) <- list(paste("AUG", 1:nrow(d), sep = ""), dimnames(x)[[2]])
    xa <- rbind.data.frame(x, d)
    if (is.factor(y)) 
        ya <- as.factor(levels(y)[c(y, e)])
    else ya <- c(y, e)
    rya <- c(ry, rep(TRUE, nr))
    wa <- c(rep(1, length(y)), rep((p + 1)/nr, nr))
    return(list(y = ya, ry = rya, x = xa, w = wa))
}
