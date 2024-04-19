colors <- function(n){
    ## Simplified version of scales::hue_pal
    hcl(h = seq(15, 375 - 360/n, length = n), l = 65, c = 100)
}


simplify <- function(x){
    if(!is.list(x)) return(x)
    if(length(x) == 1) return(x[[1]])
    r <- RTMB:::magic(numeric(length(x)))
    for(i in seq_along(x))
        r[i] <- RTMB:::magic(x[[i]])
    r
}

Value <- function(x){
    if(is(x,"advector"))
        return(Im(unclass(x)))
    x
}

##' @importFrom RTMB cbind.advector
safe_cbind <- function(...){
    args <- list(...)
    isAD <- sapply(args,function(x) is(x,"advector"))
    if(any(isAD))
        return(do.call(RTMB::cbind.advector,args))
    return(do.call(cbind,args))
}

##' @importFrom RTMB rbind.advector
safe_rbind <- function(...){
    args <- list(...)
    isAD <- sapply(args,function(x) is(x,"advector"))
    if(any(isAD))
        return(do.call(RTMB::rbind.advector,args))
    return(do.call(rbind,args))
}


safe_apply <- function(x, MARGIN, FUN, ...){
    indx <- slice.index(x, MARGIN)
    xl <- split(x, indx)
    lapply(xl, FUN, ...)
}

safe_aggregate <- function(x, by, FUN){
    r0 <- lapply(split(x,do.call("paste",c(by,list(sep=":")))),FUN)
    cl <- lapply(by,class)
    ii <- as.list(as.data.frame(do.call("rbind",strsplit(names(r0),":"))))
    ii <- lapply(seq_along(ii), function(qq) as(ii[[qq]],cl[[qq]]))
    if(!is.null(names(by)))
        names(ii) <- names(by)
    oo <- do.call(order,rev(ii))
    c(list(x = simplify(r0[oo])),lapply(ii,function(y)y[oo]))
}


squeeze <- function(u, eps = .Machine$double.eps){
    (1.0 - eps) * (u - 0.5) + 0.5
}


svd_solve <- function(x){
    ss <- svd(x)
    ss$v %*% diag(1/ss$d, length(ss$d), length(ss$d)) %*% t(ss$u)
}


lgamma <- function(x) UseMethod("lgamma")    
lgamma.default <- function(x){
    (.Primitive("lgamma"))(x)    
}
lgamma.advector <- function(x){
    RTMB:::Math1(x,"lgamma")
}


toRowLogPropMatrix <- function(x){
    y <- cbind(x,0)
    ys <- rowLogSumExp(y)
    y - ys[row(y)]
}


### Logspace functions

logspace_sub_num <- function(x,y){
    x + ifelse(x > -log(2), log(-expm1(y-x)), log1p(-exp(y-x)))
}
logspace_sub_ad <- function(x,y){
    pars <- cbind(RTMB:::magic(x),RTMB:::magic(y))
    pl <- split(pars,row(pars))
    simplify(lapply(pl, .tape_logspace_sub))
    ##.tape_logspace_sub(c(x,y))
}

##' @importFrom RTMB advector
setGeneric("logspace_sub", function(x, y){
    standardGeneric("logspace_sub")
})
setMethod("logspace_sub", signature(x="numeric",y="numeric"), function(x,y){
    logspace_sub_num(x,y)
})
setMethod("logspace_sub", signature(x="advector",y="numeric"), function(x,y){
    logspace_sub_ad(x,y)
})
setMethod("logspace_sub", signature(x="numeric",y="advector"), function(x,y){
    logspace_sub_ad(x,y)
})
setMethod("logspace_sub", signature(x="advector",y="advector"), function(x,y){
    logspace_sub_ad(x,y)
})

logspace_add_num <- function(x,y){
    ifelse(x < y, y + log1p(exp(x-y)), x + log1p(exp(y-x)))
}

logspace_add_ad <- function(x,y){
    pars <- cbind(RTMB:::magic(x),RTMB:::magic(y))
    pl <- split(pars,row(pars))
    simplify(lapply(pl, .tape_logspace_add))
    ##.tape_logspace_add(c(x,y))
}


setGeneric("logspace_add", function(x, y){
    standardGeneric("logspace_add")
})
setMethod("logspace_add", signature(x="numeric",y="numeric"), function(x,y){
    logspace_add_num(x,y)
})
setMethod("logspace_add", signature(x="advector",y="numeric"), function(x,y){
    logspace_add_ad(x,y)
})
setMethod("logspace_add", signature(x="numeric",y="advector"), function(x,y){
    logspace_add_ad(x,y)
})
setMethod("logspace_add", signature(x="advector",y="advector"), function(x,y){
    logspace_add_ad(x,y)
})


logSumExp <- function(x){
    Reduce(logspace_add, x[-1], init = x[1])
}
colLogSumExp <- function(x){
    xL <- split(x,col(x))
    simplify(lapply(xL, logSumExp))
}
rowLogSumExp <- function(x){
    xL <- split(x,row(x))
    simplify(lapply(xL, logSumExp))
}

### pnorm


log_ipnorm <- function(x,y,mu,sd){
    pars <- cbind(RTMB:::magic(x),RTMB:::magic(y),RTMB:::magic(mu),RTMB:::magic(sd))
    pl <- split(pars,row(pars))
    simplify(lapply(pl, .tape_log_ipnorm))
                                        #.tape_log_ipnorm(c(x,y,mu,sd))
}
pnorm2 <- function(x, mu, sd, lower.tail = TRUE, log.p = FALSE){
    pars <- cbind(RTMB:::magic(x),RTMB:::magic(mu),RTMB:::magic(sd))
    pl <- split(pars,row(pars))
    if(lower.tail){
        log_v <- simplify(lapply(pl, .tape_log_low_pnorm))
    }else{
        log_v <- simplify(lapply(pl, .tape_log_up_pnorm))
    }
    if(log.p)
        return(log_v)
    exp(log_v)
}

### PDE schemes

## pde_scheme_ps <- function(Pf){
##     simplify(lapply(Pf,.tape_schemes_ps))
## }



lospec_colors <- c("#201923", "#fcff5d", "#7dfc00", "#0ec434", "#228c68", "#8ad8e8",
                  "#235b54", "#29bdab", "#3998f5", "#37294f", "#277da7", "#3750db", "#f22020",
                  "#991919", "#ffcba5", "#e68f66", "#c56133", "#96341c", "#632819", "#ffc413",
                  "#f47a22", "#2f2aa0")
