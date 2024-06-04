##' @export
empirical_allele_frequencies <- function(genotypes, populations){
    ##genotypeList <- gen2list(genotypes) ## lapply(split(genotypes,slice.index(genotypes,3)),matrix,nrow=dim(genotypes)[1],ncol=dim(genotypes)[2])
    nAllele <- num_allele(genotypes, error_fun = error)
    nLoci <- num_loci(genotypes, error_fun = error)
    getAF <- function(gt){
        gbyl <- lapply(seq_len(nLoci), extract_locus, x=gt)
        ploidy <- get_ploidy(gt)
        ploidy <- as.numeric(names(ploidy))[which.max(ploidy)]
        v <- lapply(gbyl, function(x){
            r <- colSums(x,na.rm=TRUE) / sum(x,na.rm=TRUE)
            r[!is.finite(r)] <- NA
            isFull <- rowSums(x)==ploidy
            if(sum(isFull)==0){
                attr(r,"heterozygote") <- rep(NA_real_,length(r))
                attr(r,"homozygote") <- rep(NA_real_,length(r))
            }else{
                het <- colMeans(x[isFull,,drop=FALSE] > 0 & x[isFull,,drop=FALSE] < ploidy, na.rm=TRUE)
                hom <- colMeans(x[isFull,,drop=FALSE]==ploidy)
                attr(r,"heterozygote") <- het
                attr(r,"homozygote") <- hom
            }
            r
        })
        names(v) <- names(genotypes[[1]])
        class(v) <- "alleleFrequency"
        attr(v,"samples") <- length(gt)
        attr(v,"ploidy") <- ploidy
        v
    }
    AFtotal <- getAF(genotypes)
    if(!missing(populations)){
        AFsub <- lapply(split(seq_along(genotypes), populations), function(i) getAF(genotypes[i]))
    }else{
        AFsub <- NULL
    }
    r <- list(Total = AFtotal,
              Populations = AFsub)
    class(r) <- "alleleFrequencyList"
    r
}

##' @export
read_alleleFrequency_csv <- function(file, pop.col=1, skip=0, byrow=TRUE, ploidy=2,samplesize, ...){
    d <- read.csv(file, ...)
    if(!byrow)
        d <- t(d)
    if(pop.col > 0){
        nms <- d[,pop.col]
    }else{
        nms <- rep("Population",ncol(d))
    }
    if(missing(samplesize))
        samplesize <- rep(Inf,length(nms))
    af <- d[,-setdiff(c(pop.col,skip),0)]
    ## Total
    total <- colSums(af * samplesize / sum(samplesize))
    totalaf <- lapply(split(total,seq_along(total)),function(x) c("REF"=unname(x),"ALT"=1-unname(x)))
    names(totalaf) <- colnames(af)
    class(totalaf) <- "alleleFrequency"
    attr(totalaf,"samples") <- sum(samplesize)
    attr(totalaf,"ploidy") <- ploidy
    ## By population
    pops <- split(seq_len(nrow(af)),nms)
    popaf <- lapply(pops, function(ii){
        afp <- colSums(af[ii,] * samplesize[ii]/ sum(samplesize[ii]))
        r <- lapply(split(afp,seq_along(afp)),function(x){
            r <- c("REF"=unname(x),"ALT"=1-unname(x))
            he <- r * (1-r)
            ho <- r^2
            attr(r,"heterozygote") <- he
            attr(r,"homozygote") <- ho
            r
        })
        names(r) <- colnames(af)
        class(r) <- "alleleFrequency"
        attr(r,"samples") <- sum(samplesize[ii])
        attr(r,"ploidy") <- ploidy
        r
    })
    ## Return
    r <- list(Total = totalaf,
              Populations = popaf)
    class(r) <- "alleleFrequencyList"
    r
}

##' @export
combine_alleleFrequencies <- function(x, new.pop, weighted=TRUE){
    ii <- split(seq_along(x$Populations), new.pop)
    newPopCombi <- function(jj){
        y <- x$Populations[jj]
        pN <- sapply(y,attr,which = "samples")
        if(all(is.na(p))){
            pN[] <- 1
        }else{
            pN[is.na(p)] <- mean(p,na.rm=TRUE)           
        }
        p <- pmin(1e6,pmax(2,pN))
        p <- p / sum(p)
        nLoci <- num_loci(genotypes, error_fun = error)
        r <- lapply(seq_len(nLoci), function(l){
            r <- colSums(do.call("rbind",lapply(y,function(v) v[[l]])) * p)
            ## Handle heterozygote/homozygote!
            r
        })
        names(r) <- names(y[[1]])
        class(r) <- "alleleFrequency"
        attr(r,"samples") <- sum(pN)
        tp <- lapply(y, attr, which = "ploidy")
        tpN <- matrix(NA,length(tp), max(sapply(tp,length)),dimnames=list(NULL,unique(sort(Reduce(c,lapply(tp,names))))))
        for(i in seq_along(tp)){
            tpN[i,match(names(tp[[i]]),colnames(tpN))] <- tp[[i]] * p[i]
        }
        attr(r,"ploidy") <- colSums(tpN)
        r
    }
    r <- list(Total = x$Total,
              Populations = lapply(ii, function(i) newPopCombi(i))
              )
    class(r) <- "alleleFrequencyList"
    r   
}
