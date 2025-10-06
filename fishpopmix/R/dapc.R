
##' Prepare allele matrix for PCA analysis
##'
##' @param x output from \\link{read.gen}
##' @param alleleMeans Mean allele frequencies to impute. Is calculated if missing.
##' @return allele data for PCA analysis
##' @author Christoffer Moesgaard Albertsen
gen2PCA <- function(x, alleleMeans, skipFirst = TRUE){
    x2 <- do.call("rbind",lapply(x, function(gt){
        unlist(lapply(gt, function(y){
            if(sum(y)==0) return(rep(NA_real_,length(y)-ifelse(skipFirst,1,0)))
            if(skipFirst)
                return(y[-1] / sum(y))
            y / sum(y)
        }))
    }))
    if(missing(alleleMeans))
        alleleMeans <- colMeans(x2, na.rm=TRUE)
    x2[is.na(x2)] <- alleleMeans[col(x2)[is.na(x2)]]
    r <- x2 - alleleMeans[col(x2)]  
    attr(r,"alleleMeans") <- alleleMeans
    r
}




##' @export
pca <- function(genotypes, nPC = Inf, skipFirst = TRUE){
    A2 <- gen2PCA(genotypes, skipFirst = skipFirst)
    allMissing <- apply(A2,2,function(x) all(x==0))
    A2 <- A2[,!allMissing]
    CC <- scale(A2,TRUE,TRUE)
    ss <- svd(CC)
    loadings <- ss$v[,seq_len(pmin(nPC,ncol(ss$v)))]
    XX <- CC %*% loadings
    colnames(ss$v) <- colnames(XX) <- paste("PC",seq_len(ncol(XX)))
    rownames(ss$v) <- dimnames(genotypes)[[2]]
    rownames(XX) <- dimnames(genotypes)[[3]]
    
    r <- list(X = XX,
              loadings = loadings,
              center = attr(CC,"scaled:center"),
              scale = attr(CC,"scaled:scale"),
              alleleMeans = attr(A2,"alleleMeans"),
              skipFirst = skipFirst,
              data = genotypes
              )
    class(r) <- "genetic_pca"
    r
}

##' @importFrom stats predict
##' @method predict genetic_pca
##' @export
predict.genetic_pca <- function(object, newdata){
    if(missing(newdata))
        return(object$X)
    A2 <- gen2PCA(newdata, object$alleleMeans, skipFirst = object$skipFirst)
    CC <- scale(A2, object$center, object$scale)
    XX <- CC %*% object$loadings    
    XX
}

##' @export
dapc <- function(genotypes, prior, nPC = Inf, skipFirst = TRUE){
    v <- pca(genotypes, nPC = nPC, skipFirst = skipFirst)
    vG <- split(as.data.frame(v$X),attr(genotypes,"population"))
    nPG <- sapply(vG,nrow)
    muG <- lapply(vG,function(x) colMeans(as.matrix(x)))
    xc <- as.matrix(v$X) - do.call("rbind",muG)[attr(genotypes,"population"),]
    if(missing(prior))
        prior <- rep(1/length(muG),length(muG))
    prior <- prior / sum(prior)
    f1 <- sqrt(diag(var(xc)))
    scaling <- diag(1/f1, ncol(v$X), ncol(v$X))
    fac <- 1/(nrow(v$X) - length(muG))
    X <- sqrt(fac) * xc %*% scaling
    X.s <- svd(X, nu = 0L)
    rank <- sum(X.s$d > 1e-4)
    scaling <- scaling %*% X.s$v[,1:rank] %*% diag(1/X.s$d[1:rank],rank,rank)
    xbar <- colSums(prior %*% do.call("rbind",muG))
    fac <- 1 / (length(muG)-1)
    X <- sqrt((nrow(v$X) * prior) * fac) * scale(do.call("rbind",muG), center = xbar, scale = FALSE) %*% scaling
    X.s <- svd(X, nu = 0)
    rank <- sum(X.s$d > 1e-4 * X.s$d[1L])
    scaling <- scaling %*% X.s$v[, 1L:rank]
    r <- list(X = scale(as.matrix(v$X), xbar, FALSE) %*% scaling,
              population = attr(genotypes,"population"),
              scaling = scaling,
              group_means = do.call("rbind",muG),
              prior = prior,              
              pca = v)
    colnames(r$X) <- paste("Coordinate",1:ncol(r$X))
    class(r) <- "genetic_dapc"
    r
}

##' @importFrom stats predict
##' @method predict genetic_dapc
##' @export
predict.genetic_dapc <- function(object, newdata, prior){
    if(missing(newdata))
        newdata <- object$pca$data
    if(missing(prior))
        prior <- object$prior
    prior <- prior / sum(prior)
    X <- predict(object$pca, newdata)
    xbar <- colMeans(object$group_means)
    x <- scale(X,xbar,FALSE) %*% object$scaling
    colnames(x) <- paste("Coordinate",1:ncol(x))
    dm <- scale(object$group_means,xbar,FALSE) %*% object$scaling
    dist <- matrix(0.5 * rowSums(dm^2) - log(prior), nrow(x), length(prior), byrow = TRUE) - x %*% t(dm)
    dist <- exp(-(dist - apply(dist, 1L, min, na.rm = TRUE)))
    clas <- factor(rownames(object$group_means)[apply(dist, 1, which.max)], levels = rownames(object$group_means))
    list(class = clas,
         posterior = dist / rowSums(dist),
         X = x)
}

print.genetic_dapc <- function(x,...){

}

##' @method plot genetic_dapc
##' @export
plot.genetic_dapc <- function(x, pch = 16, col, legend_pos = "topleft", bty="n", ...){
      if(missing(col))
          col <- colors(length(x$prior))
      plot(x$X, pch=pch, col = col[factor(x$population)], ...)
      legend(legend_pos, legend = levels(factor(x$population)), col = col, pch=pch,bty=bty)
}
