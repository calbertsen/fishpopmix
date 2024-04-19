##' @method simulate baseline_fit
##' @export
simulate.baseline_fit <- function(object, nsim, seed = NULL, mixture, sample.name = "mixture", ploidy = 2){
    afL <- object$alleleFrequencies
    nuL <- object$alleleNu
    nLoci <- ncol(afL[[1]])
    simIndi <- function(s) do.call("cbind",lapply(seq_len(nLoci),function(j) object$rGeno(afL[[s]][,j],nuL[s,j], ploidy)))
    simStock <- sample(seq_along(afL), nsim, prob = mixture, replace = TRUE)
    aMat <- simplify2array(lapply(simStock, function(s){
        simIndi(s)
    }))
    dimnames(aMat) <- list(1:nrow(afL[[1]]), colnames(afL[[1]]), paste(sample.name,seq_len(dim(aMat)[3]),sep="_"))
    class(aMat) <- c("gen","array")
    attr(aMat,"population") <- rep(sample.name,dim(aMat)[3])
    attr(aMat,"true_population") <- names(afL)[simStock]
    attr(aMat,"alleleNames") <- if(is.null(rownames(afL[[1]]))){ seq_len(dim(afL[[1]])[1]) }else {rownames(afL[[1]])}
    ## attr(aMat,"genotypes") <- NULL
    aMat
}


## TODO: option to add dirichlet/Maxwell-Conway parameter
## TODO: option to add missingness
##' @method simulate alleleFrequencyList
##' @export
simulate.alleleFrequencyList <- function(object, nsim, seed = NULL, mixture, sample.name = "mixture", ploidy = 2){
    afL <- object$Populations
    simIndi <- function(af) do.call("cbind",apply(af,2,function(p) rmultinom(1,ploidy,p), simplify=FALSE))
    simStock <- sample(seq_along(afL), nsim, prob = mixture, replace = TRUE)
    aMat <- simplify2array(lapply(simStock, function(s){
        simIndi(afL[[s]])
    }))
    dimnames(aMat) <- list(1:nrow(afL[[1]]), colnames(afL[[1]]), paste(sample.name,seq_len(dim(aMat)[3]),sep="_"))
    class(aMat) <- c("gen","array")
    attr(aMat,"population") <- rep(sample.name,dim(aMat)[3])
    attr(aMat,"true_population") <- names(afL)[simStock]
    attr(aMat,"alleleNames") <- if(is.null(rownames(afL[[1]]))){ seq_len(dim(afL[[1]])[1]) }else {rownames(afL[[1]])}
    ## attr(aMat,"genotypes") <- NULL
    aMat
}
