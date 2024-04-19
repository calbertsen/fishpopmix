
kde2d <- function (X, bw = covafillr::suggestBandwith(X, -1),
                   npred = 100,
                   from = min(X), 
                   to = max(X)) 
{
    d <- ifelse(is.matrix(X), ncol(X), 1)
    n <- ifelse(is.matrix(X), nrow(X), length(X))
    cf <- covafillr::covafill(coord = X, obs = rep(1, n), p = -1L, h = bw)
    I <- d * (d + 2) * gamma(d/2)/(4 * pi^(d/2))
    if (length(npred) < d) 
        npred <- rep(npred, length = d)
    if (length(from) < d) 
        from <- rep(from, length = d)
    if (length(to) < d) 
        to <- rep(to, length = d)
    coords <- expand.grid(lapply(as.list(1:d), function(i) seq(from[i], 
        to[i], length = npred[i])))
    dens <- apply(coords, 1, function(x) I * cf$predict(matrix(x, 
        1)))
    if (!is.null(colnames(X))) 
        colnames(coords) <- colnames(X)
    return(list(coord = coords, density = dens, predict = function(coords){apply(coords, 1, function(x) I * cf$predict(matrix(x, 1)))}))
}


hdr <- function(X, q = c(0.05,0.5),
                npred = 100,
                from = min(X), 
                to = max(X)){
    usr <- par("usr")
    Z <- kde2d(X, npred=npred, from = c(usr[1],usr[3]), to = c(usr[2],usr[4]))
    ##contour(unique(Z$coord[,1]),unique(Z$coord[,2]),matrix(Z$density,100,100, byrow=FALSE),add=TRUE, col = caMisc::addTint("#000000",0.5), lwd = 1, nlevels = 10,labels=NA)
    ## Highest density region
    ## DOI: 10.1080/00031305.1996.10474359
    falph <- quantile(Z$predict(X),q)
    cnt <- contourLines(unique(Z$coord[,1]),unique(Z$coord[,2]),matrix(Z$density,100,100, byrow=FALSE), levels = falph)
    cnt
}
