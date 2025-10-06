
##' @importFrom RTMB dmultinom REPORT MakeADFun matrix 
##' @importFrom stats nlminb
##' @export
baseline_model <- function(genotypes, HWD = 0, HWD_density = c("CMM","DM")){
    HWD_density <- match.arg(HWD_density)
    ##genotypeList <- gen2list(genotypes)
    #lapply(split(genotypes,slice.index(genotypes,3)),matrix,nrow=dim(genotypes)[1],ncol=dim(genotypes)[2])
    Population <- factor(attr(genotypes,"population"))
    nPop <- nlevels(Population)
    nAllele <- num_allele(genotypes, error_fun = stop)
    nLoci <- num_loci(genotypes, error_fun = stop)
    nIndi <- length(genotypes)
    ploidy <- get_ploidy(genotypes)
    ploidy <- as.numeric(names(ploidy)[which.max(ploidy)])

    af_locIndx <- rep(seq_len(nLoci),nAllele-1)
    af_popIndx <- rep(seq_len(nPop), each = length(af_locIndx))
    ap2af <- function(x){
        lapply(split(x, af_popIndx),function(y){
            lapply(split(y,af_locIndx),function(a){
                exp(c(a,0.0)) / sum(exp(c(a,0.0)))
            })
        })
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
        ##RTMB::REPORT(alleleFrequencies)
        alleleNu <- par$alleleNu
        ##RTMB::REPORT(alleleNu)
        ll <- Reduce("+",lapply(seq_along(genotypes), function(i){ ## Over individuals
            Reduce("+",lapply(seq_len(nLoci), function(l){ ## Over loci
                if(isAllMissing[[as.integer(Population)[i]]][l])
                    return(0)
                dGeno(genotypes[[i]][[l]],alleleFrequencies[[as.integer(Population)[i]]][[l]], alleleNu[Population[i],l])                
            }))
        }))        
        -ll
    }

    ## Initialize using empirical allele frequencies
    eaf <- empirical_allele_frequencies(genotypes,Population)$Populations
    af0 <- unlist(lapply(eaf, function(x){
        unlist(lapply(x, function(y){
            y <- squeeze(y,1e-5)
            v <- log(head(y,-1)) - log(tail(y,1))
            v[v < -10] <- -10
            v[v > 10] <- 10
            v[is.nan(v)] <- NA
            ##names(v) <- paste0(p,"_",l,"_",names(v))
            v
        }))
    }))    
    map$alleleFrequenciesIn <- seq_along(af0)
    map$alleleFrequenciesIn[is.na(af0)] <- NA
    map$alleleFrequenciesIn <- factor(map$alleleFrequenciesIn)
    isAllMissing <- lapply(eaf,function(af) sapply(af,function(x)any(is.na(x))))
    ## map$alleleNu[which(do.call("rbind",isAllMissing))] <- NA
    ## map$alleleNu <- factor(map$alleleNu)
    ##af0 <- numeric(nLoci * (nAllele-1) * nPop)
    par0 <- list(alleleFrequenciesIn = unname(af0),
                 alleleNu = alleleNu)
    obj <- RTMB::MakeADFun(nll_baseline, par0, map = map)
    obj$fn()

    low <- obj$par; low[] <- -Inf
    low[names(low) == "alleleFrequenciesIn"] <- -10   
    up <- obj$par; up[] <- Inf
    up[names(up) == "alleleFrequenciesIn"] <- 10
    ## if(HWD_density == "DM"){
    low[names(low) == "alleleNu"] <- -10
    up[names(up) == "alleleNu"] <- 10
    ## }
    opt <- stats::nlminb(obj$par, obj$fn, obj$gr,
                  control = list(eval.max = 10000, iter.max = 10000),
                  lower = low,
                  upper = up)
    ##rp <- obj$report(obj$env$last.par.best)
    pl <- obj$env$parList(obj$env$last.par.best)
    rp <- list(alleleFrequencies = ap2af(pl$alleleFrequenciesIn),
               alleleNu = pl$alleleNu)
    nSamp <- table(Population)
    afO <- lapply(seq_along(rp$alleleFrequencies), function(k){
        gt <- rp$alleleFrequencies[[k]]
        v <- lapply(seq_along(gt), function(l){
            A <- gt[[l]]
            names(A) <- names(eaf[[k]][[l]])
            if(any(is.na(A))){
                attr(A,"heterozygote") <- A*NA
                attr(A,"homozygote") <- A*NA
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
        names(v) <- names(eaf[[k]])
        class(v) <- "alleleFrequency"
        attr(v,"samples") <- nSamp[k]
        attr(v,"ploidy") <- ploidy
        v
    })
    names(afO) <- levels(Population)
    afOut <- list(Total = NULL,
                  Populations = afO)
    class(afOut) <- "alleleFrequencyList"
    rownames(rp$alleleNu) <- levels(Population)
    colnames(rp$alleleNu) <- names(genotypes[[1]])
    res <- list(genotype = genotypes,
                alleleFrequencies = afOut,
                alleleNu = rp$alleleNu,
                opt = opt,
                obj = obj,
                nobs = length(genotypes),
                sampleNames = attr(genotypes, "population"),
                AIC = 2 * opt$objective + 2 * length(opt$par),
                HWD = HWD,
                dGeno = dGeno,
                rGeno = rGeno,
                Type = "baseline")
    class(res) <- "baseline_fit"
    res
}
