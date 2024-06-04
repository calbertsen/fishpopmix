##' @method simulate baseline_fit
##' @export
simulate.baseline_fit <- function(object, nsim, seed = NULL, mixture, sample.name = "mixture", ploidy = 2){
    afL <- object$alleleFrequencies$Populations
    nuL <- object$alleleNu
    simIndi <- function(s) {
        r <- lapply(seq_along(afL[[s]]),function(l){
            p <- afL[[s]][[l]]
            if(any(is.na(p))){
                a <- numeric(ploidy)
            }else{
                a <- object$rGeno(p,nuL[s,l], ploidy)
            }
            names(a) <- names(p)
            a
        })
        names(r) <- names(afL[[s]])
        class(r) <- "Genotype"
        attr(r,"ploidy") <- ploidy
        r
    }
    simStock <- sample(seq_along(afL), nsim, prob = mixture, replace = TRUE)
    res <- lapply(simStock, function(s){
        simIndi(s)
    })
    class(res) <- "gen"
    attr(res,"population") <- rep(sample.name,length(res))
    attr(res,"true_population") <- names(afL)[simStock] 
    res
}


## TODO: option to add dirichlet/Maxwell-Conway parameter
## TODO: option to add missingness
##' @method simulate alleleFrequencyList
##' @export
simulate.alleleFrequencyList <- function(object, nsim, seed = NULL, mixture, sample.name = "mixture", ploidy = 2){
    afL <- object$Populations
    simIndi <- function(s) {
        r <- lapply(seq_along(afL[[s]]),function(l){
            p <- afL[[s]][[l]]
            if(any(is.na(p))){
                a <- numeric(ploidy)
            }else{
                a <- rmultinom(1,ploidy,p)
            }
            names(a) <- names(p)
            a
        })
        names(r) <- names(afL[[s]])
        class(r) <- "Genotype"
        attr(r,"ploidy") <- ploidy
        r
    }
    simStock <- sample(seq_along(afL), nsim, prob = mixture, replace = TRUE)
    res <- lapply(simStock, function(s){
        simIndi(s)
    })
    class(res) <- "gen"
    attr(res,"population") <- rep(sample.name,length(res))
    attr(res,"true_population") <- names(afL)[simStock] 
    res
}
