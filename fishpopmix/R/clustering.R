##' @export
cluster_admixture <- function(genotypes, nPop, HWD = 0, categorical = TRUE){
    if(nPop <= 1) stop("Number of clusters must be at least two.")
    if(HWD != 0 && categorical) warning("HWD currently has no effect when categorical is TRUE")
    nAllele <- num_allele(genotypes, error_fun = stop)
    nLoci <- num_loci(genotypes, error_fun = stop)
    nIndi <- length(genotypes)
    ploidy <- get_ploidy(genotypes)
    ploidy <- as.numeric(names(ploidy)[which.max(ploidy)])
    
    eaf <- empirical_allele_frequencies(genotypes)
    Ohomo <- sapply(eaf$Total,function(x) sum(attr(x,"homozygote")))
    Phomo <- sapply(eaf$Total,function(x)sum(x^2))
    xx <- (Phomo - Ohomo) / Phomo
    jOffset <- which.max(ifelse(!is.finite(xx),1,xx))

    af_locIndx <- rep(seq_len(nLoci),nAllele-1)
    af_popIndx <- rep(seq_len(nPop), each = length(af_locIndx))
    ap2af <- function(x){
        AFpl <- AFplIn <- lapply(split(x, af_popIndx),function(y) split(y,af_locIndx))
        for(p in seq_len(nPop)){
            for(l in seq_len(nLoci)){
                a <- AFplIn[[p]][[l]]
                if(p > 1 && l == jOffset){
                    a[1] <- exp(a[1]) + AFplIn[[p-1]][[l]][1]
                    AFplIn[[p]][[l]] <- a
                }
                a <- exp(c(a,0.0)) / sum(exp(c(a,0.0)))
                AFpl[[p]][[l]] <- a
            }
        }
        AFpl
    }
    ## ap2af <- function(x){
    ##     A0 <- array(x,dim = c(nAllele-1,nLoci,nPop))
    ##     A1L <- lapply(split(A0, slice.index(A0,3)),matrix, nrow=nAllele-1,ncol=nLoci)
    ##     ## WHEN Clustering
    ##     for(i in tail(seq_along(A1L),-1)){
    ##             A1L[[i]][1,jOffset] <- exp(A1L[[i]][1,jOffset]) + A1L[[i-1]][1,jOffset]
    ##     }
    ##     lapply(A1L, function(A1) do.call(safe_cbind,safe_apply(safe_rbind(A1,0), 2, function(a) exp(a) / sum(exp(a)))))
    ## }

    map <- list()
    alleleNu <- matrix(1, nPop, nLoci)
    map$alleleNu <- factor(seq_along(alleleNu))
    if(categorical){
        map$alleleNu <- factor(rep(NA,length(alleleNu)))
        dGeno <- NULL
        rGeno <- NULL
    }else{
        if(HWD == 0){
            dGeno <- function(x, p, nu) dmultinom(x,sum(x),p,log = TRUE)
            rGeno <- function(p, nu, ploidy){
                if(any(is.na(p) | is.na(nu)))
                    return(numeric(ploidy))
                rmultinom(1,ploidy, p)
            }
            map$alleleNu <- factor(rep(NA,length(alleleNu)))
        }else{        
            dGeno <- function(x, p, nu) dConwayMaxwellMultinomial(x,sum(x),p,nu,log = TRUE)
            rGeno <- function(p, nu, ploidy){
                if(any(is.na(p) | is.na(nu)))
                    return(numeric(ploidy))
                rConwayMaxwellMultinomial(1,ploidy, p, nu)
            }
            if(HWD == 1){
                map$alleleNu <- factor(row(alleleNu))
            }else if(HWD == 2){

            }
        }
    }

    getU <- function(uIn) do.call(safe_cbind,safe_apply(safe_rbind(uIn,0), 2, function(a) exp(a) / sum(exp(a))))
    
    nll_admixture <- function(par){
        loadNamespace("RTMB")
        alleleFrequencies <- ap2af(par$alleleFrequenciesIn)
        u <- getU(par$uIn)

        if(categorical){
            ll <- Reduce("+",lapply(seq_along(genotypes), function(i){ ## Over individuals
                Reduce("+",lapply(seq_len(nLoci), function(l){ ## Over loci
                    ## Convert from number per allele e.g. (1,1)
                    gC <- genotypes[[i]][[l]]
                    ## To genotype ("001002")
                    xNum <- rep(1:length(gC), times = gC)
                    if(length(xNum) == 0) return(0)
                    ## Sum over allele copies
                    Reduce("+", lapply(xNum, function(a){
                        Reduce(logspace_add,lapply(seq_len(nPop), function(j){ ## Over populations
                            log(alleleFrequencies[[j]][[l]][a]) + log(u[j,i])
                        }))
                    }))                
                }))
            }))
        }else{
            ll <- Reduce("+",lapply(seq_along(genotypes), function(i){ ## Over individuals
                Reduce("+",lapply(seq_len(nLoci), function(l){ ## Over loci
                    Reduce(logspace_add,lapply(seq_len(nPop), function(j){ ## Over populations
                        ## NOTE: Need to sum over alleles ??
                        dGeno(genotypes[[i]][,l],alleleFrequencies[[j]][[l]], par$alleleNu[j,l]) + log(u[j,i])
                    }))
                }))
            }))
        }
        -ll
    }

    af0 <- unlist(lapply(seq_len(nPop), function(p){
        unlist(lapply(seq_along(eaf$Total), function(l){
            y <- eaf$Total[[l]]
            y <- squeeze(y,1e-5)
            v <- log(head(y,-1)) - log(tail(y,1))            
            v[v < -10] <- -10
            v[v > 10] <- 10
            v[is.nan(v)] <- NA
            #if(l == jOffset){                
            v[1] <- v[1] - 1 + (p-1) * (1-(-1))/(nPop-1)                
            #}
            names(v) <- paste0(p,"_",l,"_",names(v))
            v            
        }))
    }))    
    par0 <- list(alleleFrequenciesIn = af0,
                 uIn = matrix(0, nPop-1, nIndi),
                 alleleNu = alleleNu)
    obj <- RTMB::MakeADFun(nll_admixture, par0, map = map)#, profile="uIn")
    obj$fn()
    
    low <- obj$par; low[] <- -Inf
    low[names(low) == "alleleFrequenciesIn"] <- -10
    low[names(low) == "uIn"] <- -10
    up <- obj$par; up[] <- Inf
    up[names(up) == "alleleFrequenciesIn"] <- 10
    up[names(up) == "uIn"] <- 10
    opt <- stats::nlminb(obj$par, obj$fn, obj$gr,#obj$he,
                         control = list(eval.max = 100000, iter.max = 100000),
                         lower = low,
                         upper = up)

    ##rp <- obj$report(obj$env$last.par.best)
    pl <- obj$env$parList(obj$env$last.par.best)

    posterior <- getU(pl$u)
    attr(posterior,"sample") <- attr(genotypes, "population")
    attr(posterior,"groups") <- paste0("P",1:nPop)
    attr(posterior,"type") <- "admixture"
    class(posterior) <- "fit_posterior"

    rp <- list(alleleFrequencies = ap2af(pl$alleleFrequenciesIn),
               alleleNu = pl$alleleNu)
    afO <- lapply(seq_along(rp$alleleFrequencies), function(k){
        gt <- rp$alleleFrequencies[[k]]
        v <- lapply(seq_along(gt), function(l){
            A <- gt[[l]]
            names(A) <- names(genotypes[[1]][[l]])
            if(categorical){
                attr(A,"heterozygote") <- A * (1-A)
                attr(A,"homozygote") <- A^2
            }else{
                Gteo <- nexcom(ploidy,length(A))            
                pOut <- exp(sapply(Gteo,dGeno, p = A, nu = rp$alleleNu[k,l]))
                GteoM <- do.call(rbind,Gteo)
                isHomo <- GteoM == ploidy
                isHete <- GteoM > 0 & GteoM < ploidy
                pHomo <- sapply(split(isHomo,col(isHomo)),function(jj) sum(pOut[jj]))
                pHete <- sapply(split(isHete,col(isHete)),function(jj) sum(pOut[jj]))
                attr(A,"heterozygote") <- pHete
                attr(A,"homozygote") <- pHomo
            }
            A
        })
        names(v) <- names(genotypes[[1]])
        class(v) <- "alleleFrequency"
        attr(v,"samples") <- Inf
        attr(v,"ploidy") <- ploidy
        v
    })
    names(afO) <- paste0("P",1:nPop)
    afOut <- list(Total = NULL,
                  Populations = afO)
    class(afOut) <- "alleleFrequencyList"

    res <- list(alleleFrequencies = afOut,
                alleleNu = pl$alleleNu,
                IndividualClusterProportions = posterior,
                HWD = HWD,
                dGeno = dGeno,
                rGeno = rGeno,
                opt = opt,
                obj = obj,
                nobs = length(genotypes),
                sampleNames = attr(genotypes, "population"),
                AIC = 2 * opt$objective + 2 * length(opt$par),
                Type = "admixture")
    class(res) <- "clustering_fit"
    res
}


##' @export
cluster_populations <- function(genotypes, nPop, HWD = 0){
    if(nPop <= 1) stop("Number of clusters must be at least two.")

     nAllele <- num_allele(genotypes, error_fun = stop)
    nLoci <- num_loci(genotypes, error_fun = stop)
    nIndi <- length(genotypes)
    ploidy <- get_ploidy(genotypes)
    ploidy <- as.numeric(names(ploidy)[which.max(ploidy)])

    
    eaf <- empirical_allele_frequencies(genotypes)
    Ohomo <- sapply(eaf$Total,function(x) sum(attr(x,"homozygote")))
    Phomo <- sapply(eaf$Total,function(x)sum(x^2))
    xx <- (Phomo - Ohomo) / Phomo
    jOffset <- which.max(ifelse(!is.finite(xx),1,xx))

    af_locIndx <- rep(seq_len(nLoci),nAllele-1)
    af_popIndx <- rep(seq_len(nPop), each = length(af_locIndx))
    ap2af <- function(x){
        AFpl <- AFplIn <- lapply(split(x, af_popIndx),function(y) split(y,af_locIndx))
        for(p in seq_len(nPop)){
            for(l in seq_len(nLoci)){
                a <- AFplIn[[p]][[l]]
                if(p > 1 && l == jOffset){
                    a[1] <- exp(a[1]) + AFplIn[[p-1]][[l]][1]
                    AFplIn[[p]][[l]] <- a
                }
                a <- exp(c(a,0.0)) / sum(exp(c(a,0.0)))
                AFpl[[p]][[l]] <- a
            }
        }
        AFpl
    }

    map <- list()
    alleleNu <- matrix(1, nPop, nLoci)
    map$alleleNu <- factor(seq_along(alleleNu))
    if(HWD == 0){
        dGeno <- function(x, p, nu) dmultinom(x,sum(x),p,log = TRUE)
           rGeno <- function(p, nu, ploidy){
            if(any(is.na(p) | is.na(nu)))
                return(numeric(ploidy))
            rmultinom(1,ploidy, p)
        }
        map$alleleNu <- factor(rep(NA,length(alleleNu)))
    }else{        
        dGeno <- function(x, p, nu) dConwayMaxwellMultinomial(x,sum(x),p,nu,log = TRUE)
        rGeno <- function(p, nu, ploidy){
            if(any(is.na(p) | is.na(nu)))
                return(numeric(ploidy))
            rConwayMaxwellMultinomial(1,ploidy, p, nu)
        }
        if(HWD == 1){
            map$alleleNu <- factor(row(alleleNu))
        }else if(HWD == 2){

        }
    }
    
    nll_populations <- function(par){
        loadNamespace("RTMB")
        alleleFrequencies <- ap2af(par$alleleFrequenciesIn)
        RTMB::REPORT(alleleFrequencies)

        u <- exp(c(par$uIn,0)) / sum(exp(c(par$uIn,0)))
        RTMB::REPORT(u)
        
        r <- lapply(seq_along(genotypes), function(i){ ## Over individuals
                lapply(seq_len(nPop), function(j){ ## Over populations
                    Reduce("+",lapply(seq_len(nLoci), function(l){ ## Over loci
                        dGeno(genotypes[[i]][[l]],alleleFrequencies[[j]][[l]], par$alleleNu[j,l])
                })) + log(u[j])
            })
        })
        lli <- lapply(r,function(x) Reduce(logspace_add,x))
        ## "Posterior"
        post <- exp(do.call(safe_cbind,lapply(seq_along(r), function(i) simplify(r[[i]]) - lli[[i]])))
        REPORT(post)
        -Reduce("+",lli)
    }
    
    af0 <- unlist(lapply(seq_len(nPop), function(p){
        unlist(lapply(seq_along(eaf$Total), function(l){
            y <- eaf$Total[[l]]
            y <- squeeze(y,1e-5)
            v <- log(head(y,-1)) - log(tail(y,1))            
            v[v < -10] <- -10
            v[v > 10] <- 10
            v[is.nan(v)] <- NA
                                        #if(l == jOffset){                
            v[1] <- v[1] - 1 + (p-1) * (1-(-1))/(nPop-1)                
                                        #}
            names(v) <- paste0(p,"_",l,"_",names(v))
            v            
        }))
    }))    
    par0 <- list(alleleFrequenciesIn = af0,
                 uIn = matrix(0, nPop-1, nIndi),
                 alleleNu = alleleNu)
    obj <- RTMB::MakeADFun(nll_populations, par0, map = map)
    obj$fn()
    
      low <- obj$par; low[] <- -Inf
    low[names(low) == "alleleFrequenciesIn"] <- -10
    low[names(low) == "uIn"] <- -10
    up <- obj$par; up[] <- Inf
    up[names(up) == "alleleFrequenciesIn"] <- 10
    up[names(up) == "uIn"] <- 10
   opt <- stats::nlminb(obj$par, obj$fn, obj$gr,
                         control = list(eval.max = 100000, iter.max = 100000),
                         lower = low,
                         upper = up)

    rp <- obj$report(obj$env$last.par.best)

    posterior <- rp$post
    attr(posterior,"sample") <- attr(genotypes, "population")
    attr(posterior,"groups") <- 1:nPop
    attr(posterior,"type") <- "posterior"
    class(posterior) <- "fit_posterior"

    afO <- lapply(seq_along(rp$alleleFrequencies), function(k){
        gt <- rp$alleleFrequencies[[k]]
        v <- lapply(seq_along(gt), function(l){
            A <- gt[[l]]
            names(A) <- names(genotypes[[1]][[l]])
            Gteo <- nexcom(ploidy,length(A))            
            pOut <- exp(sapply(Gteo,dGeno, p = A, nu = rp$alleleNu[k,l]))
            GteoM <- do.call(rbind,Gteo)
            isHomo <- GteoM == ploidy
            isHete <- GteoM > 0 & GteoM < ploidy
            pHomo <- sapply(split(isHomo,col(isHomo)),function(jj) sum(pOut[jj]))
            pHete <- sapply(split(isHete,col(isHete)),function(jj) sum(pOut[jj]))
            attr(A,"heterozygote") <- pHete
            attr(A,"homozygote") <- pHomo            
            A
        })
        names(v) <- names(genotypes[[1]])
        class(v) <- "alleleFrequency"
        attr(v,"samples") <- Inf
        attr(v,"ploidy") <- ploidy
        v
    })
    names(afO) <- paste0("P",1:nPop)
    afOut <- list(Total = NULL,
                  Populations = afO)
    class(afOut) <- "alleleFrequencyList"

    
    res <- list(alleleFrequencies = afOut,
                alleleNu = rp$alleleNu,
                IndividualClusterProportions = posterior,
                HWD = HWD,
                dGeno = dGeno,
                rGeno = rGeno,
                opt = opt,
                obj = obj,
                nobs = length(genotypes),
                sampleNames = attr(genotypes, "population"),
                AIC = 2 * opt$objective + 2 * length(opt$par),
                Type = "populations")
    class(res) <- "clustering_fit"
    res
}

##' @method plot fit_posterior
##' @export
plot.fit_posterior <- function(x, col, sort_samples = TRUE, ...){
    if(is.null(attr(x,"sample"))){
        attr(x,"sample") <- rep(1,ncol(x))
    }
    nGrp <- nlevels(factor(attr(x,"sample")))
    
    if(missing(col))
        col <- colors(nrow(x))

    nIndi <- length(attr(x,"sample"))
    i2grp <- split(seq_along(attr(x,"sample")),attr(x,"sample"))
    grpMeanCP <- do.call("rbind",lapply(i2grp, function(ii) rowMeans(x[,ii,drop=FALSE])))
    grpNames <- levels(factor(attr(x,"sample")))
    if(sort_samples){
        grpOrder <- do.call(order,c(split(grpMeanCP,slice.index(grpMeanCP,2)),list(decreasing=TRUE)))
        i2grp <- i2grp[grpOrder]
        grpNames <- grpNames[grpOrder]
    }
    plot.new()
    #rect(0,0,1,1)
    dx <- 1/(nIndi + nGrp - 1)
    npg <- sapply(i2grp, length)
    for(i in seq_along(i2grp)){
        xGroup <- sum(c(0,npg[seq_len(i-1)]))+i-1
        uInGrp <- x[,i2grp[[i]], drop = FALSE]
        if(sort_samples){
            inGrpOrder <- do.call(order,c(split(uInGrp,slice.index(uInGrp,1))[order(grpMeanCP[[i]],decreasing=TRUE)],list(decreasing=TRUE)))
            uInGrp <- uInGrp[,inGrpOrder]
        }
        for(j in seq_along(i2grp[[i]])){
            nj <- xGroup+j
            uu <- cumsum(c(0,uInGrp[,j]))
            for(k in head(seq_along(uu),-1))
                rect((nj-1)*dx,uu[k],(nj)*dx,uu[k+1], col = col[k], border=NA)
        }
    }
    xGs <- sapply(1:nGrp, function(i) sum(c(0,npg[seq_len(i-1)]))+i-1) * dx
    xGe <- c(tail(xGs,-1),1)
    text((xGe + xGs)/2,0,grpNames,pos = 1, font = 2, ...)
    rect(0,0,1,1)    
}

##' @method plot clustering_fit
##' @export
plot.clustering_fit <- function(x, col, sort_samples = TRUE, ...){
    plot(x$IndividualClusterProportions,col=col,sort_samples=sort_samples,...)
}
