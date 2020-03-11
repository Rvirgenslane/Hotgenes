#' Returns adjusted pvalue from decideTests function, instead of -1, 0 or 1 
#' @inheritParams limma::decideTests
#' @importFrom stats p.adjust pt
#' @importFrom methods new
#' @export

Limma_Adj_pval <- function(object,method="separate",adjust.method="BH",p.value=0.05,lfc=0)

{
.classifyTestsP <- function(object,df=Inf,p.value=0.05,method="holm") {
#	TestResults by rows for a matrix t-statistics using adjusted p-values
#	Gordon Smyth
#	12 July 2003.  Last modified 23 March 2004.

#	Method intended for MArrayLM objects but accept unclassed lists as well
if(is.list(object)) {
if(is.null(object$t)) stop("tstat cannot be extracted from object")
tstat <- object$t
if(!is.null(object$df.residual)) df <- object$df.residual
if(!is.null(object$df.prior)) df <- df+object$df.prior
} else {
tstat <- object
}
if(is.null(dim(tstat))) dim(tstat) <- c(1,length(tstat))
ngenes <- nrow(tstat)
P <- 2*pt(-abs(tstat),df=df)
result <- tstat
for (i in 1:ngenes) {
P[i,] <- p.adjust(P[i,],method=method)
result[i,] <- sign(tstat[i,])*(P[i,]<p.value)
}
new("TestResults",result)
}


#   This was copied directly from decideTests function
#	Accept or reject hypothesis tests across genes and contrasts
#	Gordon Smyth
#	17 Aug 2004. Last modified 24 June 2007.

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
if(any(sel)) results[sel,] <- .classifyTestsP(object[sel,],p.value=p.value*a,method=adjust.method)
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
