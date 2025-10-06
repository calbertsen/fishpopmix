##' @export
mixture_model <- function(baseline, samples, formulaProportion, data, ...){
    UseMethod("mixture_model")
}

##' @method mixture_model alleleFrequencyList
##' @export
mixture_model.alleleFrequencyList <- function(baseline, samples, formulaProportion = ~1, data,...){
    ## Baseline is fixed, assuming HWE
    alleleFrequencies <- baseline$Populations
                                        #genotypeList <- lapply(split(samples,slice.index(samples,3)),matrix,nrow=dim(samples)[1],ncol=dim(samples)[2])    
    nPop <- length(alleleFrequencies)
    Population <- names(alleleFrequencies)
    ## nAllele <- dim(samples)[1]
    ## nLoci <- dim(samples)[2]
    ## nMixIndi <- dim(samples)[3]
    nAllele <- num_allele(samples, error_fun = stop)
    nLoci <- num_loci(samples, error_fun = stop)
    nMixIndi <- length(samples)
    ploidy <- get_ploidy(samples)
    ploidy <- as.numeric(names(ploidy)[which.max(ploidy)])

    alleleNu <- matrix(1, nPop, nLoci)
    rownames(alleleNu) <- levels(Population)
    colnames(alleleNu) <- names(baseline[[1]])

    dGeno <- function(x, p, nu){
        dmultinomGen(x,sum(x),p,log = TRUE)
    }
    rGeno <- function(p, nu, ploidy){
        if(any(is.na(p) | is.na(nu)))
            return(numeric(ploidy))
        rmultinom(1,ploidy, p)
    }
    
    if(missing(data))
        data <- data.frame(ID = seq_len(nMixIndi))

    ## Setup formula for proportions
    mfTheta <- model.frame(lme4::nobars(formulaProportion), data,
                           na.action = na.pass)
    XTheta <- Matrix::sparse.model.matrix(terms(mfTheta),
                                          data = mfTheta,
                                          transpose = FALSE,
                                          row.names = FALSE,
                                          drop.unused.levels = drop.unused.levels)
    ## if(!inherits(XTheta,"dgTMatrix"))
    XTheta <- as(XTheta,"TsparseMatrix")

    if(is.null(lme4::findbars(formulaProportion))){
        ZT <- list()
        UT <- array(0,dim = c(0))
        UT_rdim <- UT_cdim <- integer(0)
        UTcor <- array(0, dim = c(0))
        attr(UTcor,"rdim") <- integer(0)
        attr(UTcor,"cdim") <- integer(0)        
        UTlogSd <- array(0, dim = c(0))
        attr(UTlogSd,"rdim") <- integer(0)
        attr(UTlogSd,"cdim") <- integer(0)
        Uvec2list <- function(x,...) x
    }else{
        lc <- lme4::lmerControl(check.nobs.vs.rankZ = "ignore",
                                check.nobs.vs.nlev = "ignore",
                                check.nlev.gtreq.5 = "ignore",
                                check.nlev.gtr.1 = "ignore",
                                check.nobs.vs.nRE= "ignore",
                                check.rankX = "ignore",
                                check.scaleX = "ignore",
                                check.formula.LHS = "ignore")

        rtZT <- lme4::lFormula(formulaProportion,data, na.action = na.pass, control = lc)$reTrms
        ZT <- lapply(rtZT$Ztlist,function(xx){
            as(xx,"TsparseMatrix")
        })
        ZTnms <- rtZT$cnms
        ZTrdim <- sapply(rtZT$cnms,length)
        ZTcdim <- sapply(rtZT$flist,nlevels)
        UT <- numeric(sum(ZTrdim * ZTcdim * (nPop-1)))
        UT_rdim <- as.integer(ZTrdim * ZTcdim)
        UT_cdim <- as.integer(rep(nPop-1, length(ZT)))
        n <- ZTrdim * (ZTrdim - 1) / 2
        UTcor <- array(0, dim = sum(n*(nPop-1)))
        UTC_rdim <- as.integer(n)
        UTC_cdim <- as.integer(rep(nPop-1, length(ZT)))
        UTlogSd <- array(2, dim = sum(ZTrdim * (nPop-1)))
        UTS_rdim <- as.integer(ZTrdim)
        UTS_cdim <- as.integer(rep(nPop-1, length(ZT)))

        Uvec2list <- function(U, rdim, cdim){
            res <- unname(split(U,factor(rep(seq_along(cdim),times = rdim * cdim), levels = seq_along(cdim))))
            lapply(1:length(res), function(i) RTMB::matrix(res[[i]], rdim[i], cdim[i]))
        }
    }

    isAllMissing <- lapply(alleleFrequencies,function(af) sapply(af,function(x)any(is.na(x))))
    
    par0 <- list(UTheta = UT,
                 UThetaCor = UTcor,
                 UThetalogSd = UTlogSd,
                 betaTheta = matrix(0, ncol(XTheta), nPop-1))
    
    nll_mix <- function(par){
        
        ## Random effects (u)
        if(length(par$UTheta) > 0){
            U_list <- Uvec2list(par$UTheta,UT_rdim, UT_cdim)
            U_cor_list <- Uvec2list(par$UThetaCor,UTC_rdim, UTC_cdim)
            U_logsd_list <- Uvec2list(par$UThetalogSd,UTS_rdim, UTS_cdim)
            ## RE Likelihood
            llRE <- Reduce("+",lapply(seq_along(U_list), function(ii){
                nsd <- nrow(U_logsd_list[[ii]])
                if(nsd == 1){
                    corfn <- function(x) RTMB::matrix(1,1,1)
                }else{
                    corfn <- RTMB::unstructured(nsd)$corr
                }
                Reduce("+",lapply(seq_len(nPop-1), function(g){
                    sdv <- exp(U_logsd_list[[ii]][,g])
                    suppressWarnings(CorMat <- corfn(U_cor_list[[ii]][,g]))
                    DSV <- RTMB::diag(sdv,length(sdv),length(sdv))
                    Sigma <- DSV %*% CorMat %*% DSV                    
                    Uv <- RTMB::matrix(U_list[[ii]][,g],nrow=length(sdv))
                    sum(RTMB:::dmvnorm(t(Uv), mu = 0, Sigma = Sigma, log = TRUE))#, scale = sdv))
                }))
            }))
            ## Prediction
            RE_theta <- lapply(seq_along(ZT), function(i) Matrix::t(t(U_list[[i]]) %*% ZT[[i]]))
            RE <- Reduce("+",RE_theta)                
        }else{
            llRE <- 0
            RE <- matrix(0,nrow(XTheta),nPop-1)
        }
      
        ## Theta 
        logTheta <-  toRowLogPropMatrix(XTheta %*% par$betaTheta + RE) ## + RE
  
        ## Observations
        r <- lapply(seq_along(samples), function(i){ ## Over individuals
            ## P(Genotype | Population) * P(Population)
            lapply(seq_len(nPop), function(j){ ## Over populations
                Reduce("+",lapply(seq_len(nLoci), function(l){ ## Over loci
                    if(isAllMissing[[j]][l])
                        return(0)
                    dGeno(samples[[i]][[l]],alleleFrequencies[[j]][[l]], alleleNu[j,l])
                })) + logTheta[i,j]
            })            
        })
        ## Likelihood
        lli <- lapply(r,function(x) Reduce(logspace_add,x))
        ## "Posterior"
        post <- exp(do.call(cbind,lapply(seq_along(r), function(i) simplify(r[[i]]) - lli[[i]])))
        RTMB::REPORT(post)
        -Reduce("+",lli) - llRE        
    }

    obj <- RTMB::MakeADFun(nll_mix, par = par0, random = "UTheta")
    opt <- nlminb(obj$par, obj$fn, obj$gr)
    rep <- obj$report()
    res <- list(baseline_genotypes = NULL,
                samples_genotypes = samples,
                alleleFrequencies = baseline,
                alleleNu = alleleNu,
                opt = opt,
                obj = obj,
                AIC = 2 * opt$objective + 2 * length(opt$par),
                HWD = 0,
                dGeno = dGeno,
                rGeno = rGeno,
                call = match.call(),
                Type = "Mixture")
    class(res) <- "mixture_fit"
    res
}



##' @method mixture_model baseline_fit
##' @export
mixture_model.baseline_fit <- function(baseline, samples, formulaProportion = ~1, data, conditional = FALSE, ...){
    
    alleleFrequencies <- baseline$alleleFrequencies$Populations
    baselinePopulation <- factor(attr(baseline$genotype,"population"))
                                        #genotypeList <- lapply(split(samples,slice.index(samples,3)),matrix,nrow=dim(samples)[1],ncol=dim(samples)[2])    
    nPop <- length(alleleFrequencies)
    Population <- names(alleleFrequencies)
    ## nAllele <- dim(samples)[1]
    ## nLoci <- dim(samples)[2]
    ## nMixIndi <- dim(samples)[3]
    nAllele <- num_allele(samples, error_fun = stop)
    nLoci <- num_loci(samples, error_fun = stop)
    nMixIndi <- length(samples)
    ploidy <- get_ploidy(samples)
    ploidy <- as.numeric(names(ploidy)[which.max(ploidy)])

    alleleNu <- baseline$alleleNu

    af_locIndx <- rep(seq_len(nLoci),nAllele-1)
    af_popIndx <- rep(seq_len(nPop), each = length(af_locIndx))
    ap2af <- function(x){
        lapply(split(x, af_popIndx),function(y){
            lapply(split(y,af_locIndx),function(a){
                exp(c(a,0.0)) / sum(exp(c(a,0.0)))
            })
        })
    }
    
    
    dGeno <- baseline$dGeno
    rGeno <- baseline$rGeno

    map <- list()
    map$alleleNu <- factor(seq_along(alleleNu))
    if(baseline$HWD == 0){
        map$alleleNu <- factor(rep(NA,length(alleleNu)))
    }else if(baseline$HWD == 1){
        map$alleleNu <- factor(row(alleleNu))
    }else if(baseline$HWD == 2){
        ## Do nothing
    }else{
        stop("HWD type not implemented")
    }
  
    
    if(missing(data))
        data <- data.frame(ID = seq_len(nMixIndi))

    ## Setup formula for proportions
    mfTheta <- model.frame(lme4::nobars(formulaProportion), data,
                           na.action = na.pass)
    XTheta <- Matrix::sparse.model.matrix(terms(mfTheta),
                                          data = mfTheta,
                                          transpose = FALSE,
                                          row.names = FALSE,
                                          drop.unused.levels = drop.unused.levels)
    ## if(!inherits(XTheta,"dgTMatrix"))
    XTheta <- as(XTheta,"TsparseMatrix")

    if(is.null(lme4::findbars(formulaProportion))){
        ZT <- list()
        UT <- array(0,dim = c(0))
        UT_rdim <- UT_cdim <- integer(0)
        UTcor <- array(0, dim = c(0))
        attr(UTcor,"rdim") <- integer(0)
        attr(UTcor,"cdim") <- integer(0)        
        UTlogSd <- array(0, dim = c(0))
        attr(UTlogSd,"rdim") <- integer(0)
        attr(UTlogSd,"cdim") <- integer(0)
        Uvec2list <- function(x,...) x
    }else{
        lc <- lme4::lmerControl(check.nobs.vs.rankZ = "ignore",
                                check.nobs.vs.nlev = "ignore",
                                check.nlev.gtreq.5 = "ignore",
                                check.nlev.gtr.1 = "ignore",
                                check.nobs.vs.nRE= "ignore",
                                check.rankX = "ignore",
                                check.scaleX = "ignore",
                                check.formula.LHS = "ignore")

        rtZT <- lme4::lFormula(formulaProportion,data, na.action = na.pass, control = lc)$reTrms
        ZT <- lapply(rtZT$Ztlist,function(xx){
            as(xx,"TsparseMatrix")
        })
        ZTnms <- rtZT$cnms
        ZTrdim <- sapply(rtZT$cnms,length)
        ZTcdim <- sapply(rtZT$flist,nlevels)
        UT <- numeric(sum(ZTrdim * ZTcdim * (nPop-1)))
        UT_rdim <- as.integer(ZTrdim * ZTcdim)
        UT_cdim <- as.integer(rep(nPop-1, length(ZT)))
        n <- ZTrdim * (ZTrdim - 1) / 2
        UTcor <- array(0, dim = sum(n*(nPop-1)))
        UTC_rdim <- as.integer(n)
        UTC_cdim <- as.integer(rep(nPop-1, length(ZT)))
        UTlogSd <- array(2, dim = sum(ZTrdim * (nPop-1)))
        UTS_rdim <- as.integer(ZTrdim)
        UTS_cdim <- as.integer(rep(nPop-1, length(ZT)))

        Uvec2list <- function(U, rdim, cdim){
            res <- unname(split(U,factor(rep(seq_along(cdim),times = rdim * cdim), levels = seq_along(cdim))))
            lapply(1:length(res), function(i) RTMB::matrix(res[[i]], rdim[i], cdim[i]))
        }
    }

    af0 <- unlist(lapply(alleleFrequencies, function(x){
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
    if(conditional){
        map$alleleFrequenciesIn[] <- NA
        map$alleleNu[] <- NA
    }
    map$alleleFrequenciesIn <- factor(map$alleleFrequenciesIn)
    isAllMissing <- lapply(alleleFrequencies,function(af) sapply(af,function(x)any(is.na(x))))
    
    
    par0 <- list(alleleFrequenciesIn = unname(af0),
                 alleleNu = alleleNu,
                 UTheta = UT,
                 UThetaCor = UTcor,
                 UThetalogSd = UTlogSd,
                 betaTheta = matrix(0, ncol(XTheta), nPop-1))
    
    nll_mix <- function(par){
        alleleFrequencies <- ap2af(par$alleleFrequenciesIn)
        alleleNu <- par$alleleNu
        RTMB::REPORT(alleleFrequencies)
        RTMB::REPORT(alleleNu)
        ## Random effects (u)
        if(length(par$UTheta) > 0){
            U_list <- Uvec2list(par$UTheta,UT_rdim, UT_cdim)
            U_cor_list <- Uvec2list(par$UThetaCor,UTC_rdim, UTC_cdim)
            U_logsd_list <- Uvec2list(par$UThetalogSd,UTS_rdim, UTS_cdim)
            ## RE Likelihood
            llRE <- Reduce("+",lapply(seq_along(U_list), function(ii){
                nsd <- nrow(U_logsd_list[[ii]])
                if(nsd == 1){
                    corfn <- function(x) RTMB::matrix(1,1,1)
                }else{
                    corfn <- RTMB::unstructured(nsd)$corr
                }
                Reduce("+",lapply(seq_len(nPop-1), function(g){
                    sdv <- exp(U_logsd_list[[ii]][,g])
                    suppressWarnings(CorMat <- corfn(U_cor_list[[ii]][,g]))
                    DSV <- RTMB::diag(sdv,length(sdv),length(sdv))
                    Sigma <- DSV %*% CorMat %*% DSV                    
                    Uv <- RTMB::matrix(U_list[[ii]][,g],nrow=length(sdv))
                    sum(RTMB:::dmvnorm(t(Uv), mu = 0, Sigma = Sigma, log = TRUE))#, scale = sdv))
                }))
            }))
            ## Prediction
            RE_theta <- lapply(seq_along(ZT), function(i) Matrix::t(t(U_list[[i]]) %*% ZT[[i]]))
            RE <- Reduce("+",RE_theta)                
        }else{
            llRE <- 0
            RE <- matrix(0,nrow(XTheta),nPop-1)
        }
      
        ## Theta 
        logTheta <-  toRowLogPropMatrix(XTheta %*% par$betaTheta + RE) ## + RE

        ## Baseline
        if(!conditional){
            llBaseline <- Reduce("+",lapply(seq_along(baseline$genotype), function(i){ ## Over individuals
                Reduce("+",lapply(seq_len(nLoci), function(l){ ## Over loci
                    if(isAllMissing[[as.integer(baselinePopulation)[i]]][l])
                        return(0)
                    dGeno(baseline$genotype[[i]][[l]],alleleFrequencies[[as.integer(baselinePopulation)[i]]][[l]], alleleNu[baselinePopulation[i],l])                
                }))
            })) 
        }else{
            llBaseline <- 0
        }
        ## Observations
        r <- lapply(seq_along(samples), function(i){ ## Over individuals
            ## P(Genotype | Population) * P(Population)
            lapply(seq_len(nPop), function(j){ ## Over populations
                Reduce("+",lapply(seq_len(nLoci), function(l){ ## Over loci
                    if(isAllMissing[[j]][l])
                        return(0)
                    dGeno(samples[[i]][[l]],alleleFrequencies[[j]][[l]], alleleNu[j,l])
                })) + logTheta[i,j]
            })            
        })
        ## Likelihood
        lli <- lapply(r,function(x) Reduce(logspace_add,x))
        ## "Posterior"
        post <- exp(do.call(cbind,lapply(seq_along(r), function(i) simplify(r[[i]]) - lli[[i]])))
        RTMB::REPORT(post)
        -Reduce("+",lli) - llRE - llBaseline        
    }
    
    obj <- RTMB::MakeADFun(nll_mix, par = par0, map = map, random = "UTheta")

    low <- obj$par; low[] <- -Inf
    low[names(low) == "alleleFrequenciesIn"] <- -10   
    up <- obj$par; up[] <- Inf
    up[names(up) == "alleleFrequenciesIn"] <- 10
    ## if(HWD_density == "DM"){
    low[names(low) == "alleleNu"] <- -10
    up[names(up) == "alleleNu"] <- 10
 
    opt <- nlminb(obj$par, obj$fn, obj$gr,
                  control = list(eval.max = 10000, iter.max = 10000),
                  lower = low,
                  upper = up)
    pl <- obj$env$parList(par = obj$env$last.par.best)
    rp <- obj$report(par = obj$env$last.par.best)

    sdr <- sdreport(obj, par.fixed = opt$par)
    
    nSamp <- table(baselinePopulation)
    
    afO <- lapply(seq_along(rp$alleleFrequencies), function(k){
        gt <- rp$alleleFrequencies[[k]]
        v <- lapply(seq_along(gt), function(l){
            A <- gt[[l]]
            names(A) <- names(alleleFrequencies[[k]][[l]])
            if(any(is.na(A))){
                attr(A,"heterozygote") <- A*NA
                attr(A,"homozygote") <- A*NA
            }else{                
                Gteo <- nexcom(ploidy,length(A))            
                pOut <- exp(sapply(Gteo,dGeno, p = A, nu = pl$alleleNu[k,l]))
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
        names(v) <- names(alleleFrequencies[[1]])
        class(v) <- "alleleFrequency"
        attr(v,"samples") <- nSamp[k]
        attr(v,"ploidy") <- ploidy
        v
    })
    names(afO) <- levels(baselinePopulation)
    afOut <- list(Total = NULL,
                  Populations = afO)
    class(afOut) <- "alleleFrequencyList"
    rownames(pl$alleleNu) <- levels(baselinePopulation)
    colnames(pl$alleleNu) <- names(baseline$genotypes[[1]])

  
    res <- list(baseline_genotypes = baseline$genotypes,
                samples_genotypes = samples,
                conditional = conditional,
                alleleFrequencies = afOut,
                alleleNu = pl$alleleNu,
                opt = opt,
                obj = obj,
                sdr = sdr,
                AIC = 2 * opt$objective + 2 * length(opt$par),
                HWD = baseline$HWD,
                dGeno = dGeno,
                rGeno = rGeno,
                classifications = rp$post,
                call = match.call(),
                data = data,
                mf = mfTheta,
                terms = terms(mfTheta),
                xlevels = .getXlevels(terms(mfTheta), mfTheta),
                all_vars = get_all_vars(terms(mfTheta), data),
                formula = formulaProportion,
                Type = "Mixture")
    class(res) <- "mixture_fit"
    res
}


##' @export
getGroupProportion <- function(fit, newdata, ...){
    UseMethod("getGroupProportion")
}

##' @method getGroupProportion mixture_fit
##' @export
getGroupProportion.mixture_fit <- function(fit, newdata, randEff=FALSE, ...){

    if(missing(newdata))
        newdata <- fit$mf

    XTheta <- Matrix::sparse.model.matrix(fit$terms,
                                          data = newdata,
                                          xlev = fit$xlevels,
                                          transpose = FALSE,
                                          row.names = FALSE,
                                          drop.unused.levels = drop.unused.levels)
    XTheta <- as(XTheta,"TsparseMatrix")

    if(!randEff || is.null(lme4::findbars(fit$formulaProportion))){
        ZT <- list()
        UT <- array(0,dim = c(0))
        UT_rdim <- UT_cdim <- integer(0)
        UTcor <- array(0, dim = c(0))
        attr(UTcor,"rdim") <- integer(0)
        attr(UTcor,"cdim") <- integer(0)        
        UTlogSd <- array(0, dim = c(0))
        attr(UTlogSd,"rdim") <- integer(0)
        attr(UTlogSd,"cdim") <- integer(0)
        Uvec2list <- function(x,...) x
    }else{
        stop("Not implemented with random effects yet")
        lc <- lme4::lmerControl(check.nobs.vs.rankZ = "ignore",
                                check.nobs.vs.nlev = "ignore",
                                check.nlev.gtreq.5 = "ignore",
                                check.nlev.gtr.1 = "ignore",
                                check.nobs.vs.nRE= "ignore",
                                check.rankX = "ignore",
                                check.scaleX = "ignore",
                                check.formula.LHS = "ignore")

        rtZT <- lme4::lFormula(fit$formulaProportion,newdata, na.action = na.pass, control = lc)$reTrms
        ZT <- lapply(rtZT$Ztlist,function(xx){
            as(xx,"TsparseMatrix")
        })
        ZTnms <- rtZT$cnms
        ZTrdim <- sapply(rtZT$cnms,length)
        ZTcdim <- sapply(rtZT$flist,nlevels)
        UT <- numeric(sum(ZTrdim * ZTcdim * (nPop-1)))
        UT_rdim <- as.integer(ZTrdim * ZTcdim)
        UT_cdim <- as.integer(rep(nPop-1, length(ZT)))
        n <- ZTrdim * (ZTrdim - 1) / 2
        UTcor <- array(0, dim = sum(n*(nPop-1)))
        UTC_rdim <- as.integer(n)
        UTC_cdim <- as.integer(rep(nPop-1, length(ZT)))
        UTlogSd <- array(2, dim = sum(ZTrdim * (nPop-1)))
        UTS_rdim <- as.integer(ZTrdim)
        UTS_cdim <- as.integer(rep(nPop-1, length(ZT)))

        Uvec2list <- function(U, rdim, cdim){
            res <- unname(split(U,factor(rep(seq_along(cdim),times = rdim * cdim), levels = seq_along(cdim))))
            lapply(1:length(res), function(i) RTMB::matrix(res[[i]], rdim[i], cdim[i]))
        }
    }

    predfn <- function(p){
        betaTheta <- RTMB::matrix(p,nrow = ncol(XTheta))
        logTheta <-  toRowLogPropMatrix(XTheta %*% betaTheta) ## + RE
        for(i in 1:nrow(logTheta))
            for(j in 1:ncol(logTheta))
                logTheta[i,j] <- logTheta[i,j] - logspace_sub(0.0,logTheta[i,j])
        logTheta
    }
    p0 <- fit$opt$par[names(fit$opt$par) %in% "betaTheta"]
    tp <- RTMB::MakeTape(predfn, p0)

    lfitted <- as.matrix(predfn(p0))
    fitted <- plogis(lfitted)
    G <- tp$jacobian(p0)
    Sigma <- fit$sdr$cov.fixed[names(fit$opt$par) %in% "betaTheta",names(fit$opt$par) %in% "betaTheta",drop=FALSE]
    Sdval <- sqrt(diag(G %*% Sigma %*% t(G)))
    CIlow <- plogis(lfitted - 2 * Sdval)
    CIhigh <- plogis(lfitted + 2 * Sdval)
    rownames(fitted) <- rownames(CIlow) <- rownames(CIhigh) <- rownames(newdata)
    colnames(fitted) <- colnames(CIlow) <- colnames(CIhigh) <- names(fit$alleleFrequencies$Populations)
    list(fit = fitted,
         CI_low = CIlow,
         CI_high = CIhigh)
}
