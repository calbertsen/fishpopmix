##' @export
classify <- function(baseline, samples, prior, ...){
    UseMethod("classify")
}

##' @method classify alleleFrequencyList
##' @export
classify.alleleFrequencyList <- function(baseline, samples, prior, ...){
    alleleFrequencies <- baseline$Populations
    isAllMissing <- lapply(alleleFrequencies,function(af) which(apply(af,2,function(x)any(is.na(x)))))  
    genotypeList <- gen2list(samples) ##lapply(split(samples,slice.index(samples,3)),matrix,nrow=dim(samples)[1],ncol=dim(samples)[2])
    nPop <- length(alleleFrequencies)
    nAllele <- dim(samples)[1]
    nLoci <- dim(samples)[2]
    nMixIndi <- dim(samples)[3]
    if(missing(prior))
        prior <- rep(1/nPop,nPop)

    ## Assign
    r <- lapply(seq_along(genotypeList), function(i){ ## Over individuals
        ## P(Genotype | Population) * P(Population)
        lapply(seq_len(nPop), function(j){ ## Over populations
            Reduce("+",lapply(setdiff(seq_len(nLoci),isAllMissing[[j]]), function(l){ ## Over loci                
                dmultinom(genotypeList[[i]][,l],sum(genotypeList[[i]][,l]),alleleFrequencies[[j]][,l],log = TRUE)                    
            })) + log(prior)[j]
        })            
    })
    ## Likelihood
    lli <- lapply(r,function(x) Reduce(logspace_add,x))
    ## "Posterior"
    post <- exp(do.call(cbind,lapply(seq_along(r), function(i) simplify(r[[i]]) - lli[[i]])))
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
    alleleFrequencies <- baseline$alleleFrequencies
    isAllMissing <- lapply(alleleFrequencies,function(af) which(apply(af,2,function(x)any(is.na(x)))))  
    alleleNu <- baseline$alleleNu
    genotypeList <- gen2list(samples) #lapply(split(samples,slice.index(samples,3)),matrix,nrow=dim(samples)[1],ncol=dim(samples)[2])
    nPop <- length(alleleFrequencies)
    nAllele <- dim(samples)[1]
    nLoci <- dim(samples)[2]
    nMixIndi <- dim(samples)[3]
    if(missing(prior))
        prior <- rep(1/nPop,nPop)
    dGeno <- baseline$dGeno

     r <- lapply(seq_along(genotypeList), function(i){ ## Over individuals
        ## P(Genotype | Population) * P(Population)
        lapply(seq_len(nPop), function(j){ ## Over populations
            Reduce("+",lapply(setdiff(seq_len(nLoci),isAllMissing[[j]]), function(l){ ## Over loci
                dGeno(genotypeList[[i]][,l], alleleFrequencies[[j]][,l], alleleNu[j,l])
            })) + log(prior)[j]
        })            
    })
    ## Likelihood
    lli <- lapply(r,function(x) Reduce(logspace_add,x))
    ## "Posterior"
    post <- exp(do.call(cbind,lapply(seq_along(r), function(i) simplify(r[[i]]) - lli[[i]])))
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
