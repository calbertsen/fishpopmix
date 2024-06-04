dmultinomGen <- function(x, size, prob, log) {
    ## Modified from from RTMB::dmultinom, which modified from stats::dmultinom
    K <- length(prob)
    if (length(x) != K)
        stop("x[] and prob[] must be equal length vectors.")
    s <- sum(prob)
    prob <- prob / s
    if (missing(size))
        size <- sum(x)
    r <- lgamma(size + 1) + sum(x * log(prob) - lgamma(x + 1))
    if (log)
        r
    else exp(r)
}

rConwayMaxwellMultinomial <- function(n, size, p, nu){
    p <- p / sum(p)
    possible <- do.call(rbind,nexcom(size,length(p)))
    ## possible <- do.call("expand.grid",lapply(seq_along(p), function(i) seq_len(size+1)-1))
    ## possible <- possible[rowSums(possible)==size,]
    getOne <- function(k, log = FALSE){
        if(log){
            return(nu * (lgamma(size+1) - sum(lgamma(k+1))) + sum(k * log(p)))
        }
        (gamma(size+1) / prod(gamma(k+1)))^nu * prod(p^k)

    }
    v2 <- sapply(safe_apply(possible,1, getOne, log = TRUE),exp)
    i <- sample(seq_len(nrow(possible)),n,prob = v2, replace = TRUE)
    res <- possible[i,]
    unname(res)
    t(res)
}

dConwayMaxwellMultinomial <- function(x, size, p, nu, log = FALSE){
    p <- p / sum(p)
    possible <- do.call(rbind,nexcom(size,length(p)))
    ## possible <- as.matrix(do.call("expand.grid",lapply(seq_along(p), function(i) seq_len(size+1)-1)))
    ## possible <- possible[rowSums(possible)==size,]
    getOne <- function(k, log = FALSE){
        if(log){
            return(nu * (lgamma(size+1) - sum(lgamma(k+1))) + sum(k * log(p)))
        }
        (gamma(size+1) / prod(gamma(k+1)))^nu * prod(p^k)

    }
    v2 <- Reduce(logspace_add,safe_apply(possible,1, getOne, log = TRUE))
    v1 <- getOne(x, TRUE)
    if(log)
        return(v1 - v2)
    return(exp(v1 - v2))
}


rDirichletMultinomial <- function(n, size, p, alpha, log = FALSE){    
    p0 <- apply(matrix(rgamma(n*length(p),shape = p*alpha, rate = 1),length(p),n),2,function(x)x/sum(x))
    apply(p0,2, function(pp) rmultinom(1,size,pp))
}
    
dDirichletMultinomial <- function(x, size, p, alpha, log = FALSE){
    r <- lgamma(alpha) + lgamma(size + 1) - lgamma(size + alpha) + Reduce("+",lgamma(x+p*alpha)-lgamma(p*alpha)-lgamma(x+1))
    if(!log)
        return(exp(r))
    r
}
