##' @export
empirical_allele_frequencies <- function(genotypes, populations){
    genotypeList <- gen2list(genotypes) ## lapply(split(genotypes,slice.index(genotypes,3)),matrix,nrow=dim(genotypes)[1],ncol=dim(genotypes)[2])
    A <- lapply(genotypeList, function(x) apply(x,2,function(y){
        if(sum(y)==0)
            return(rep(NA,length(y)))
        y/sum(y)
    }, simplify = TRUE))

    getAF <- function(x){
        r <- apply(x,1:2,mean, na.rm=TRUE)
        r[is.nan(r)] <- NA
        r
    }
    AFtotal <- getAF(simplify2array(A))
    if(!missing(populations)){
        AL <- lapply(split(A, populations), simplify2array)
        AFsub <- lapply(AL, getAF)
    }else{
        AFsub <- NULL
    }
    r <- list(Total = AFtotal,
              Populations = AFsub)
    class(r) <- "alleleFrequencyList"
    r
}

##' @export
read_alleleFrequency_csv <- function(file, pop.col=1, skip=0, byrow=TRUE, ...){
    d <- read.csv(file, ...)
    if(!byrow)
        d <- t(d)
    nms <- d[,pop.col]
    af <- d[,-setdiff(c(pop.col,skip),0)]
    total <- colMeans(af)
    totaf <- rbind(total,1-total)
    rownames(totaf) <- 1:2
    pops <- split(af,seq_len(nrow(af)))
    names(pops) <- nms
    r <- list(Total = totaf,
              Populations = lapply(pops,function(x) rbind(x,1-x)))
    class(r) <- "alleleFrequencyList"
    r
}

##' @export
combine_alleleFrequencies <- function(x, new.pop){
    ii <- split(seq_along(x$Populations), new.pop)
    r <- list(Total = x$Total,
              Populations = lapply(ii, function(i) Reduce("+",x$Populations[i]) / length(i)))
    class(r) <- "alleleFrequencyList"
    r   
}
