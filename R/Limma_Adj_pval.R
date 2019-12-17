#' Returns adjusted pvalue from decideTests function, instead of -1, 0 or 1 
#' @inheritParams limma::decideTests
#' @importFrom stats p.adjust
#' @export

Limma_Adj_pval <- function(object,method="separate",adjust.method="BH",p.value=0.05,lfc=0)
    #   This was copied directly from decideTests function
    #	Accept or reject hypothesis tests across genes and contrasts
    #	Gordon Smyth
    #	17 Aug 2004. Last modified 24 June 2007.
{
    if(!methods::is(object,"MArrayLM")) stop("Need MArrayLM object")
    if(is.null(object$t)) object <- limma::eBayes(object)
    method <- match.arg(method,c("separate","global","hierarchical","nestedF"))
    adjust.method <- match.arg(adjust.method,c("none","bonferroni","holm","BH","fdr","BY"))
    if(adjust.method=="fdr") adjust.method <- "BH"
    switch(method,separate={
        p <- as.matrix(object$p.value)
        tstat <- as.matrix(object$t)
        for (j in 1:ncol(p)) {
            o <- !is.na(p[,j])
            p[o,j] <- p.adjust(p[o,j],method=adjust.method)
        }
        results <- methods::new("TestResults",p)
    },global={
        p <- as.matrix(object$p.value)
        tstat <- as.matrix(object$t)
        o <- !is.na(p)
        p[o] <- p.adjust(p[o],method=adjust.method)
        results <- methods::new("TestResults",p)
    },hierarchical={
        if(any(is.na(object$F.p.value))) stop("Can't handle NA p-values yet")
        sel <- p.adjust(object$F.p.value,method=adjust.method) < p.value
        i <- sum(sel,na.rm=TRUE)
        n <- sum(!is.na(sel))
        a <- switch(adjust.method,
                    none=1,
                    bonferroni=1/n,
                    holm=1/(n-i+1),
                    BH=i/n,
                    BY=i/n/sum(1/(1:n))
        )
        results <- methods::new("TestResults",array(0,dim(object$t)))
        dimnames(results) <- dimnames(object$coefficients)
        if(any(sel)) results[sel,] <- limma::classifyTestsP(object[sel,],p.value=p.value*a,method=adjust.method)
    },nestedF={
        if(any(is.na(object$F.p.value))) stop("nestedF method can't handle NA p-values",call.=FALSE)
        sel <- p.adjust(object$F.p.value,method=adjust.method) < p.value
        i <- sum(sel,na.rm=TRUE)
        n <- sum(!is.na(sel))
        a <- switch(adjust.method,
                    none=1,
                    bonferroni=1/n,
                    holm=1/(n-i+1),
                    BH=i/n,
                    BY=i/n/sum(1/(1:n))
        )
        results <- methods::new("TestResults",array(0,dim(object$t)))
        dimnames(results) <- dimnames(object$coefficients)
        if(any(sel)) results[sel,] <- limma::classifyTestsF(object[sel,],p.value=p.value*a)
    })
    if(lfc>0) {
        if(is.null(object$coefficients))
            warning("lfc ignored because coefficients not found")
        else
            results@.Data <- results@.Data * (abs(object$coefficients)>lfc)
    }
    results
}
