##' @importFrom RTMB dmultinom REPORT MakeADFun matrix 
##' @importFrom stats nlminb
##' @export
baseline_model <- function(genotypes, HWD = 0, HWD_density = c("CMM","DM")){
    HWD_density <- match.arg(HWD_density)
    genotypeList <- gen2list(genotypes)
    #lapply(split(genotypes,slice.index(genotypes,3)),matrix,nrow=dim(genotypes)[1],ncol=dim(genotypes)[2])
    Population <- factor(attr(genotypes,"population"))
    nPop <- nlevels(Population)
    nAllele <- dim(genotypes)[1]
    nLoci <- dim(genotypes)[2]
    nIndi <- dim(genotypes)[3]
    
    ap2af <- function(x){
        A0 <- array(x,dim = c(nAllele-1,nLoci,nPop))
        A1L <- lapply(split(A0, slice.index(A0,3)),RTMB::matrix, nrow=nAllele-1,ncol=nLoci)
        lapply(A1L, function(A1) do.call(safe_cbind,safe_apply(safe_rbind(A1,0), 2, function(a) exp(a) / sum(exp(a)))))
    }

    map <- list()
    alleleNu <- matrix(1, nPop, nLoci)
    map$alleleNu <- factor(seq_along(alleleNu))
    if(HWD == 0){        
        dGeno <- function(x, p, nu){
            dmultinomGen(x,sum(x),p,log = TRUE)
        }
        rGeno <- function(p, nu, ploidy){
            if(any(is.na(p) | is.na(nu)))
                return(numeric(ploidy))
            rmultinom(1,ploidy, p)
        }
        map$alleleNu <- factor(rep(NA,length(alleleNu)))
    }else{
        if(HWD_density == "CMM"){
            dGeno <- function(x, p, nu){
                dConwayMaxwellMultinomial(x,sum(x),p,nu,log = TRUE)
            }
            rGeno <- function(p, nu, ploidy){
                if(any(is.na(p) | is.na(nu)))
                    return(numeric(ploidy))
                rConwayMaxwellMultinomial(1,ploidy, p, nu)
            }
        }else if(HWD_density == "DM"){
            dGeno <- function(x, p, nu){
                dDirichletMultinomial(x,sum(x),p,exp(nu),log = TRUE)
            }
            rGeno <- function(p, nu, ploidy){
                if(any(is.na(p) | is.na(nu)))
                    return(numeric(ploidy))
                rDirichletMultinomial(1,ploidy, p, exp(nu))
            }
        }
        if(HWD == 1){
            map$alleleNu <- factor(row(alleleNu))
        }else if(HWD == 2){
          
        }
    }
    
    nll_baseline <- function(par){
        loadNamespace("RTMB")
        alleleFrequencies <- ap2af(par$alleleFrequenciesIn)
        RTMB::REPORT(alleleFrequencies)
        alleleNu <- par$alleleNu
        RTMB::REPORT(alleleNu)
        ll <- Reduce("+",lapply(seq_along(genotypeList), function(i){ ## Over individuals
            Reduce("+",lapply(seq_len(nLoci), function(l){ ## Over loci
                if(isAllMissing[[as.integer(Population)[i]]][l])
                    return(0)
                dGeno(genotypeList[[i]][,l],alleleFrequencies[[as.integer(Population)[i]]][,l], alleleNu[Population[i],l])                
            }))
        }))        
        -ll
    }

    ## Initialize using empirical allele frequencies
    eaf <- empirical_allele_frequencies(genotypes,Population)$Population
    af0 <- as.vector(simplify2array(lapply(eaf, function(x){
        v <- t(log(t(x[-nrow(x),,drop=FALSE])/x[nrow(x),]))
        v[v < -10] <- -10
        v[is.nan(v) | v > 10] <- 10
        v
    })))    
    map$alleleFrequenciesIn <- seq_along(af0)
    map$alleleFrequenciesIn[is.na(af0)] <- NA
    map$alleleFrequenciesIn <- factor(map$alleleFrequenciesIn)
    isAllMissing <- lapply(eaf,function(af) apply(af,2,function(x)any(is.na(x))))
    ## map$alleleNu[which(do.call("rbind",isAllMissing))] <- NA
    ## map$alleleNu <- factor(map$alleleNu)
    ##af0 <- numeric(nLoci * (nAllele-1) * nPop)
    par0 <- list(alleleFrequenciesIn = af0,
                 alleleNu = alleleNu)
    obj <- RTMB::MakeADFun(nll_baseline, par0, map = map)
    obj$fn()

    low <- obj$par; low[] <- -Inf
    low[names(low) == "alleleFrequenciesIn"] <- -10   
    up <- obj$par; up[] <- Inf
    up[names(up) == "alleleFrequenciesIn"] <- 10
    ## if(HWD_density == "DM"){
    ##     low[names(low) == "alleleNu"] <- -10
    ##     up[names(up) == "alleleNu"] <- 10
    ## }
    opt <- stats::nlminb(obj$par, obj$fn, obj$gr,
                  control = list(eval.max = 10000, iter.max = 10000),
                  lower = low,
                  upper = up)
    rp <- obj$report(obj$env$last.par.best)
    rp$alleleFrequencies <- lapply(rp$alleleFrequencies, function(x){
        colnames(x) <- dimnames(genotypes)[[2]]
        x
    })
    names(rp$alleleFrequencies) <- levels(Population)
    rownames(rp$alleleNu) <- levels(Population)
    colnames(rp$alleleNu) <- dimnames(genotypes)[[2]]
    res <- list(alleleFrequencies = rp$alleleFrequencies,
                alleleNu = rp$alleleNu,
                opt = opt,
                obj = obj,
                nobs = length(genotypeList),
                sampleNames = attr(genotypes, "population"),
                AIC = 2 * opt$objective + 2 * length(opt$par),
                HWD = HWD,
                dGeno = dGeno,
                rGeno = rGeno,
                Type = "baseline")
    class(res) <- "baseline_fit"
    res
}
