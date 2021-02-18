# ----------------------------------------------------------------------------
# Author: Shamim Mollah  
# Created: 10-12-2016
# Last updated: 1-12-2017
#
# Histone prediction using PLSR
#-----------------------------------------------------------------------------

# A modified version of pls.cv
pls.cv<-function (X, y, k = 10, groups=NULL, m = ncol(X), use.kernel = FALSE, 
                  compute.covariance = FALSE,method.cor="pearson") 
{
    cross <- nrow(X) 
    if (cross < 10 & cross > 4)
    {
        k = cross
    }
    n <- nrow(X)
    p <- ncol(X)
    if (is.null(groups) == FALSE) {
        f = as.factor(groups)
        k = length(levels(f))
        my.names = levels(f)
        all.folds <- split(1:n, f)
    }
    if (is.null(groups) == TRUE) {
        f <- rep(1:k, length = n)
        my.names <- 1:k
        all.folds <- split(sample(1:n), f)
    }
    
    ntrain = vector(length = k)
    for (i in 1:k) {
        ntrain[i] = n - length(all.folds[[i]])
    }
    ntrain.min = min(ntrain)
    m = min(m, ntrain.min - 1, p)
    cv.error.matrix = matrix(0, k, m + 1)
    rownames(cv.error.matrix) = my.names
    colnames(cv.error.matrix) = 0:m
    cor.error.matrix<-cv.error.matrix
    for (i in seq(k)) {
        omit <- all.folds[[i]]
        Xtrain = X[-omit, , drop = FALSE]
        ytrain = y[-omit]
        Xtest = X[omit, , drop = FALSE]
        ytest = y[omit]
        pls.object <- plsdof::pls.model(Xtrain, ytrain, m = m, Xtest = Xtest, 
                                ytest = ytest, compute.DoF = FALSE, use.kernel = use.kernel,method.cor=method.cor)
        cv.error.matrix[i, ] <- pls.object$mse
        cor.error.matrix[i, ] <- pls.object$cor
    }
    cv.error = apply(cv.error.matrix, 2, mean)
    cor.error<-apply(cor.error.matrix,2,mean)
    # change by smollah 11-13-2016
    
    k=2
    lst=""
    iter=length(pls.object$RSS) -1
    for (i in 1:iter) {
        if ((pls.object$RSS[i] - pls.object$RSS[k]) <= 0.05) {
            lst=append(lst,i)
        }
        k=k+1  
    }
    m.opt = as.numeric(lst[2]) 
    # plot RSS to choose the optimal component
    plot(pls.object$RSS,main=m.opt)
    #m.opt <- which.min(cv.error) - 1
    m.opt.cor<-which.max(cor.error) - 1
    if (compute.covariance == TRUE) {
        use.kernel = FALSE
    }
    pls.object <- plsdof::pls.model(X, y, m = max(m.opt, m.opt.cor,1), use.kernel = use.kernel, 
                            compute.DoF = compute.covariance, compute.jacobian = compute.covariance)
    intercept <- pls.object$intercept[m.opt + 1]
    
    coefficients <- pls.object$coefficients[, m.opt + 1]
    covariance <- pls.object$covariance
    # edited by SMollah on 10-31-2016
    DoF <- pls.object$DoF
    intercept.cor <- pls.object$intercept[m.opt.cor + 1]
    coefficients.cor <- pls.object$coefficients[, m.opt.cor + 1]
    if (compute.covariance == TRUE) {
        #covariancve.cor<-covariance[m.opt.cor + 1, , ]
        covariance <- covariance[m.opt + 1, , ]
    }
    outlist = list(cv.error.matrix = cv.error.matrix, cor.error.matrix=cor.error.matrix,cv.error = cv.error, cor.error=cor.error,
                   m.opt = m.opt, m.opt.cor=m.opt.cor,covariance = covariance, DoF = DoF, intercept = intercept, intercept.cor=intercept.cor,
                   coefficients = coefficients,coefficients.cor=coefficients.cor)
    class(outlist) = "plsdof"
    return(outlist)
}

#' Identify significant PLSR coefficients for pairs of left and right objects.
#' 
#' This function computes PLSR coefficients and p-values between the rows of the
#' two matrices provided ("left" and "right"). This regression is performed separately 
#' for each left item; i.e. in each iteration, a separate PLS regression is performed
#' which explains 
#' 
#' The maximum number of components for PLSR is the number of columns in the right 
#' matrix or the number of rows minus one, whichever is less.
#' 
#' @param left Matrix of left objects
#' @param right Matrix of right objects
#' @param p_val_threshold Threshold for significant p-values in PLSR.
#' @param verbose Whether to print output.
#'
#' @return A list containing the following named values:
#'         p_list - A list containing left-right pairs with significant values (above p_val_threshold)
#'         non_p_list - A list containing left-right pairs with non-significant values
#'         comp - the maximum number of PLS components.

generate_plsr <- function(left, right, p_val_threshold=.0001, verbose=FALSE) {

    cli::cli_alert_success("Running PLSR")
    XNormZ<-scale(as.matrix(right[,]))

    # m = maximal number of Partial Least Squares components. Default is m=min(ncol(X),nrow(X)-1)
    comp=min(ncol(XNormZ),nrow(XNormZ)-1)

    left_num = dim(left)[2]    
    right_num = dim(right)[2]  
    ii=1
    cc=0
    #for capturing error
    #the following two lines should be uncommented for debugging
    #zz <- file("all.Rout", open="wt")
    #sink(zz, type="message")

    left_names <- list()
    right_names <- list()
    abs_coefs <- list()
    pvals <- list()
    directions <- list()
    significant <- list()
    p <- list()
    non_p <- list()

    pdf_path = file.path("measured_predicted.pdf")
    pdf(file=pdf_path) 
    par(mfrow=c(2,3))

    all_left_names = c()
    all_right_names = c()
    all_abs_coefs = c()
    all_pvals = c()
    all_directions = c()
    all_significant = c()


    pb <- progress::progress_bar$new(
      format = "   Regressing... [:bar] :percent eta: :eta",
      total = left_num * right_num, clear = FALSE, width= 60)

    for (i in 1:left_num) {
        set.seed(1234)
        yCent=scale(as.vector(left[,i]), scale = FALSE)
        mypls1=plsdof::pls.model(XNormZ,yCent, m=comp, compute.DoF=TRUE)
        #mypls3_k <- 10
        #if (nrow(XNormZ) < mypls3_k) {
        #    mypls3_k <- nrow(XNormZ)
        #}
        mypls3=pls.cv(XNormZ,yCent,compute.covariance=TRUE,m=comp)
        my.vcov=vcov(mypls3)
        my.sd=sqrt(diag(my.vcov)) # standard deviation of the regression coefficients

        index= mypls3$m.opt +1
        myvec = mypls3$coefficients
        mat=XNormZ%*%myvec

        #plot(yCent,mat,xlab="measured",ylab="predicted(mycalc)", ylim=c(-2,2), xlim=c(-2,2),main=names(left)[i])
        plot(yCent,mat,xlab="measured",ylab="predicted(mycalc)", xlim=c(min(yCent),max(yCent)), ylim=c(min(mat),max(mat)),main=names(left)[i])
        ## add naive estimate
        lines(-2:2,-2:2,lwd=3)
        # add naive estimate
        #plot(yCent,mypls1$Yhat[,index],xlab="measured",ylab="predicted(Yhat)", ylim=c(-2,2), xlim=c(-2,2),main=names(left)[i])
        plot(yCent,mypls1$Yhat[,index],xlab="measured",ylab="predicted(Yhat)", ylim=c(-2,2), xlim=c(-2,2),main=names(left)[i])
        lines(-2:2,-2:2,lwd=3)

        for (k in 1:right_num) {
            pval=dt(mypls3$coefficients[k]/my.sd[k], (mypls3$m.opt))

            #if (is.na(pval)) {
            #    significant <- FALSE
            #}
            #else {
            if (pval < p_val_threshold) {
                significant = TRUE
            }
            else {
                significant = FALSE
            }
            #}
            if (mypls3$coefficients[k] < 0){
                direction = "neg"
            } 
            else {
                direction = "pos" 
            } # end of 2nd if else

            left_name <- names(left)[i]
            right_name <- names(right)[k]
            abs_coefs <- (abs(mypls3$coefficients[k]))
            all_left_names = append(all_left_names, left_name)
            all_right_names = append(all_right_names, right_name)
            all_abs_coefs = append(all_abs_coefs, abs_coefs)
            all_pvals = append(all_pvals, pval)
            all_directions = append(all_directions, direction)
            all_significant = append(all_significant, significant)
            pb$tick()
        }
    }

    # Combine attributes of each pair in a data.frame
    out <- data.frame(
        left_names = all_left_names,
        right_names = all_right_names,
        abs_coefs = all_abs_coefs,
        pvals = all_pvals,
        directions = all_directions,
        significant = all_significant
    )

    # Separate significant and non-significant pairs
    p_list <- out[unlist(out[["significant"]]),]
    non_p_list <- out[!unlist(out[["significant"]]),]

    p_list$significant <- NULL
    non_p_list$significant <- NULL

    ## Display the log file

    #sink()

    ## Display the log file
    #error <- readLines("all.Rout")

    dev.off()

    list(
        p_list=p_list,
        non_p_list=non_p_list,
        comp=comp
    )
}
