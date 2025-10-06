##' @export
classify <- function(baseline, samples, prior, ...){
    UseMethod("classify")
}

##' @method classify alleleFrequencyList
##' @export
classify.alleleFrequencyList <- function(baseline, samples, prior, ...){
    alleleFrequencies <- baseline$Populations
    isAllMissing <- lapply(alleleFrequencies, function(af) which(sapply(af,function(x) any(is.na(x)))))
    nPop <- length(alleleFrequencies)
    nAllele <- num_allele(samples)
    nLoci <- num_loci(samples)
    nMixIndi <- length(samples)
    if(missing(prior))
        prior <- rep(1/nPop,nPop)

    ## Assign
    r <- lapply(seq_along(samples), function(i){ ## Over individuals
        ## P(Genotype | Population) * P(Population)
        lapply(seq_len(nPop), function(j){ ## Over populations
            Reduce("+",lapply(setdiff(seq_len(nLoci),isAllMissing[[j]]), function(l){ ## Over loci
                if(any(is.na(alleleFrequencies[[j]][[l]])))
                    return(0)
                dmultinom(samples[[i]][[l]],sum(samples[[i]][[l]]),alleleFrequencies[[j]][[l]],log = TRUE)                    
            })) + log(prior)[j]
        })            
    })
    ## Likelihood
    lli <- lapply(r,function(x) Reduce(logspace_add,x,init = -Inf))
    ## "Posterior"
    post <- exp(do.call(cbind,lapply(seq_along(r), function(i) RTMBconvenience::simplify(r[[i]]) - lli[[i]])))
    class(post) <- "fit_posterior"
    attr(post,"sample") <- attr(samples, "population")
    attr(post,"groups") <- names(alleleFrequencies)
    rownames(post) <- names(alleleFrequencies)
    r <- list(class = factor(names(alleleFrequencies)[apply(post,2,which.max)]),
         posterior = post)
    class(r) <- "classification"
    r
}


##' @method classify baseline_fit
##' @export
classify.baseline_fit <- function(baseline, samples, prior, ...){
    alleleFrequencies <- baseline$alleleFrequencies$Populations
    isAllMissing <- lapply(alleleFrequencies, function(af) which(sapply(af,function(x) any(is.na(x)))))
    alleleNu <- baseline$alleleNu    
    nPop <- length(alleleFrequencies)
    nAllele <- num_allele(samples)
    nLoci <- num_loci(samples)
    nMixIndi <- length(samples)

    if(missing(prior))
        prior <- rep(1/nPop,nPop)
    dGeno <- baseline$dGeno

     r <- lapply(seq_along(samples), function(i){ ## Over individuals
        ## P(Genotype | Population) * P(Population)
        lapply(seq_len(nPop), function(j){ ## Over populations
            Reduce("+",lapply(setdiff(seq_len(nLoci),isAllMissing[[j]]), function(l){ ## Over loci
                dGeno(samples[[i]][[l]], alleleFrequencies[[j]][[l]], alleleNu[j,l])
            })) + log(prior)[j]
        })            
    })
    ## Likelihood
    lli <- lapply(r,function(x) Reduce(logspace_add,x))
    ## "Posterior"
    post <- exp(do.call(cbind,lapply(seq_along(r), function(i) RTMBconvenience::simplify(r[[i]]) - lli[[i]])))
    class(post) <- "fit_posterior"
    attr(post,"sample") <- attr(samples, "population")
    attr(post,"groups") <- names(alleleFrequencies)
    rownames(post) <- names(alleleFrequencies)
    r <- list(class = factor(names(alleleFrequencies)[apply(post,2,which.max)]),
              posterior = post)
    class(r) <- "classification"
    r    
}

##' @method classify genetic_dapc
##' @export
classify.genetic_dapc <- function(baseline, samples, prior, ...){
    r <- predict(baseline, samples, prior)[c("class","posterior","X")]
    class(r) <- "classification"
    r$posterior <- t(r$posterior)
    class(r$posterior) <- "fit_posterior"
    attr(r$posterior,"sample") <- attr(samples, "population")
    attr(r$posterior,"groups") <- levels(r$class)
    r
}

##' @method plot classification
##' @export
plot.classification <- function(x, ...){
    plot(x$posterior,...)
}
