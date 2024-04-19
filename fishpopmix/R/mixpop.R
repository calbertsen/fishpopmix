##' @export
mixture_model <- function(baseline, samples, formula, ...){
    UseMethod("mixture_fit")
}

##' @method mixture_model alleleFrequencyList
##' @export
mixture_model.alleleFrequencyList <- function(baseline, samples, formulaProportion = ~1, data){
    ## Baseline is fixed, assuming HWE
    alleleFrequencies <- baseline$Populations
    genotypeList <- lapply(split(samples,slice.index(samples,3)),matrix,nrow=dim(samples)[1],ncol=dim(samples)[2])
    nPop <- length(alleleFrequencies)
    nAllele <- dim(samples)[1]
    nLoci <- dim(samples)[2]
    nMixIndi <- dim(samples)[3]

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
        attr(UT,"rdim") <- integer(0)
        attr(UT,"cdim") <- integer(0)
        attr(UT,"adim") <- integer(0)
        UTcor <- array(0, dim = c(0))
        attr(UTcor,"rdim") <- integer(0)
        attr(UTcor,"cdim") <- integer(0)        
        UTlogSd <- array(0, dim = c(0))
        attr(UTlogSd,"rdim") <- integer(0)
        attr(UTlogSd,"cdim") <- integer(0)        
    }else{
        rtZT <- lme4::lFormula(formulaProportion,data, na.action = na.pass, control = lc)$reTrms
        ZT <- lapply(rtZT$Ztlist,function(xx){
            as(xx,"TsparseMatrix")
        })
        ZTnms <- rtZT$cnms
        ZTrdim <- sapply(rtZT$cnms,length)
        ZTcdim <- sapply(rtZT$flist,nlevels)
        UT <- array(numeric(sum(ZTrdim * ZTcdim * (nPop-1))))
        attr(UT,"rdim") <- as.integer(ZTrdim * ZTcdim)
        attr(UT,"cdim") <- as.integer(rep(nPop-1, length(ZT)))
        n <- ZTrdim * (ZTrdim - 1) / 2
        UTcor <- array(0, dim = sum(n*(nPop-1)))
        attr(UTcor,"rdim") <- as.integer(n)
        attr(UTcor,"cdim") <- as.integer(rep(nPop-1, length(ZT)))
        UTlogSd <- array(2, dim = sum(ZTrdim * (nPop-1)))
        attr(UTlogSd,"rdim") <- as.integer(ZTrdim)
        attr(UTlogSd,"cdim") <- as.integer(rep(nPop-1, length(ZT)))
    }

    par0 <- list(UTheta = UT,
                 UThetaCor = UTcor,
                 UThetalogSd = UTlogSd,
                 betaTheta = matrix(0, ncol(XTheta), nPop-1))
    
    nll_mix <- function(par){

        ## Random effects (u)

        ## Theta 
        logTheta <-  toRowLogPropMatrix(XTheta %*% par$betaTheta) ## + RE
  
        ## Observations
        r <- lapply(seq_along(genotypeList), function(i){ ## Over individuals
            ## P(Genotype | Population) * P(Population)
            lapply(seq_len(nPop), function(j){ ## Over populations
                Reduce("+",lapply(seq_len(nLoci), function(l){ ## Over loci                
                    dmultinom(genotypeList[[i]][,l],sum(genotypeList[[i]][,l]),alleleFrequencies[[j]][,l],log = TRUE)                    
                })) + logTheta[j]
            })            
        })
        ## Likelihood
        lli <- lapply(r,function(x) Reduce(logspace_add,x))
        ## "Posterior"
        post <- exp(do.call(cbind,lapply(seq_along(r), function(i) simplify(r[[i]]) - lli[[i]])))
        REPORT(post)
        -Reduce("+",lli)        
    }

    obj <- RTMB::MakeADFun(nll_mix, par = par0)
    opt <- nlminb(obj$par, obj$fn, obj$gr)

}
