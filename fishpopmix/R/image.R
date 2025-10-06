

toGreyScale <- function(im){
    if (!is.matrix(im) & dim(im)[3] >= 3)
        im <- round((0.2989 * im[, , 1] + 0.587 * im[, , 2] + 0.114 * im[, , 3]) * 255)
    im    
}

#' @importFrom jpeg readJPEG
#' @importFrom png readPNG
#' @importFrom tiff readTIFF
getPixelMatrix <- function (file, grey = FALSE) {
    if (grepl("\\.jp(e)*g$", file, ignore.case = TRUE)) {
        im <- jpeg::readJPEG(file)
    }
    else if (grepl("\\.png$", file, ignore.case = TRUE)) {
        im <- png::readPNG(file)
    }
    else if (grepl("\\.tif(f)*$", file, ignore.case = TRUE)) {
        im <- tiff::readTIFF(file)
    }
    else {
        stop(sprintf("Unsupported file format: %s. Only JPEG, PNG, and TIFF are supported.", 
            tail(strsplit(file, ".")[[1]], 1)))
    }
    if (grey & length(dim(im)) == 3) {
        if (dim(im)[3] >= 3) 
            imOut <- toGreyScale(im)
    }
    else {
        imOut <- round(im * 255)
    }
    return(imOut)
}

watershed <- function(seeds, priority){
    storage.mode(seeds) <- "integer"
    .Call(C_watershed,seeds,priority,PACKAGE="fishpopmix")
}

read_otolith_image <- function(im, blur_sigma = 2, blur_k = 3){
    pic <- getPixelMatrix(im, TRUE)
    edge <- sqrt(3 * otoclass:::sobel(otoclass:::extendPic(otoclass:::gaussian_blur(otoclass:::extendPic(pic,blur_k),blur_k,blur_sigma),3)))
    priority <- 1 / (1 + edge)
    ## Check for number of otoliths
    ## Seeds (assume 2 otoliths for now)
    cm <- colMeans(pic)
    rm <- rowMeans(pic)
    midrange <- which.max(cm[1:floor(length(cm)/2)]):floor(length(cm)/2)+which.max(cm[ceiling(length(cm)/2):length(cm)])
    xv <- quantile(cm,c(0.25,0.75))
    yv <- quantile(rm,c(0.25,0.75))
    seed <- matrix(0,nrow(pic),ncol(pic))
    ## Background
    isBG <- col(pic) %in% which(cm < xv[1]) | row(pic) %in% which(rm < yv[1]) | col(pic) %in% (min(midrange) + which.min(cm[midrange]) - 1)
    seed[isBG] <- 1
    ## Left otolith
    isLeft <- col(pic) %in% which(cm > xv[2]) & row(pic) %in% which(rm > yv[2]) & col(pic) < (min(midrange) + which.min(cm[midrange]) - 1)
    seed[isLeft] <- 2
    ## Right otolith
    isRight <- col(pic) %in% which(cm > xv[2]) & row(pic) %in% which(rm > yv[2]) & col(pic) > (min(midrange) + which.min(cm[midrange]) - 1)
    seed[isRight] <- 3
    ## watershed
    wt <- watershed(seed, priority)
    list(pic=pic,watershed=wt)
}


downsample <- function(pic,n=2) matrix(pic[row(pic)%%n==0 & col(pic)%%n==0],floor(nrow(pic)/n),floor(ncol(pic)/n))


otsu_threshold <- function(pic){
    pic <- toGreyScale(pic)
    h <- table(factor(floor(pic),0:255)) / length(pic)
    ii <- 1:256
    getSig <- function(t){
        w0 <- sum(h[head(ii,t)])
        m0 <- sum((head(ii,t) - 1) * h[head(ii,t)]) / w0
        w1 <- sum(h[tail(ii,-t)])
        m1 <- sum((tail(ii,-t)-1) * h[tail(ii,-t)]) / w1
        w0*w1*(m0-m1)^2
    }
    ii <- which.max(sapply(ii,getSig))-1
    picx <- 0 * pic
    picx[pic >= ii] <- 1
    picx
}


otsu_2D_threshold <- function(pic, n = 3){
    pic <- toGreyScale(pic)
    mPic <- otoclass:::meanfilter(otoclass:::extendPic(pic,n),n)
    h <- table(factor(pic,0:255), factor(mPic,0:255)) / length(pic)
    ii <- 1:256
    getSig <- Vectorize(function(t,s){
        h0 <- h[head(ii,t),head(ii,s),drop=FALSE]
        w0 <- sum(h0)
        m0 <- c(sum((row(h0)-1)*h0)/w0,sum((col(h0)-1)*h0)/w0)
        h1 <- h[tail(ii,-t),tail(ii,-s),drop=FALSE]
        w1 <- sum(h1)
        m1 <- c(sum((t+row(h1)-1)*h1)/w1,sum((s+col(h1)-1)*h1)/w1)
        mT <- w0*m0 + w1*m1
        w0*sum((m0-mT)^2)+w1*sum((m1-mT)^2)
    })
    v <- outer(ii,ii,getSig)
    ii <- which(v==max(v,na.rm=TRUE),TRUE)-1
    picx <- 0*pic+1
    picx[pic < min(ii[1]) & mPic < min(ii[2])] <- 0
    picx
}


otsu_2DD_threshold <- function(pic, n = 3){
    pic <- toGreyScale(pic)
    dx<- otoclass:::sobel_x(otoclass:::extendPic(pic,n))
    dy <- otoclass:::sobel_y(otoclass:::extendPic(pic,n))
    Theta <- round(atan2(dy,dx) * 180/pi)
    m <- sqrt(dx^2 + dy^2)
    mPic <- Theta #otoclass:::meanfilter(otoclass:::extendPic(pic,n),n)
    h <- table(factor(pic,0:255), factor(mPic,-180:180)) / length(pic)
    ii <- 0:255
    jj <- -180:180
    getSig <- Vectorize(function(ti,si){
        t <- match(ti,ii)
        s <- match(ti,jj)
        h0 <- h[head(ii,t),head(ii,s),drop=FALSE]
        w0 <- sum(h0)
        m0 <- c(sum((row(h0)-1)*h0)/w0,sum((col(h0)-1)*h0)/w0)
        h1 <- h[tail(ii,-t),tail(ii,-s),drop=FALSE]
        w1 <- sum(h1)
        m1 <- c(sum((t+row(h1)-1)*h1)/w1,sum((s+col(h1)-1)*h1)/w1)
        mT <- w0*m0 + w1*m1
        w0*sum((m0-mT)^2)+w1*sum((m1-mT)^2)
    })
    v <- outer(seq,ii,getSig)
    ii <- which(v==max(v,na.rm=TRUE),TRUE)-1
    picx <- 0*pic
    picx[pic >= ii[1] & mPic >= ii[2]] <- 1
    picx
}

## MOVE TO C++
connectedComponentLabeling <- function(l){
    ## https://enpc.hal.science/hal-03651336/documents
    cc <- matrix(-1,nrow(l),ncol(l))
    n <- 1
    H <- list()
    for(p in seq_along(l)){
        if(cc[p] == -1){
            cc[p] <- n
            C <- p
            H[n] <- 0
            while(length(C) > 0){
                q <- C[1]
                C <- C[-1]
                H[[n]] <- H[[n]] + 1
                ## find neighbours
                qi <- (q-1) %% nrow(l) + 1
                qj <- (q-1) %/% nrow(l) + 1
                windx <- outer(-1:1,-1:1,function(r,c) qi+r + (qj+c-1) * nrow(l))
                windx <- windx[windx >= 1 & windx <= length(l)]
                for(r in windx){
                    if(l[r] == l[q] & cc[r] == -1){
                        cc[r] <- n
                        C <- c(C,r)
                    }
                }
            }
            n <- n + 1
        }
    }
    list(cc=cc, H=H)    
}

## MOVE TO C++
findSuperPixels <- function(pic, k, m = 20, tol = 0.001, max.pic = 255, onlyImage = FALSE,splitFinalPixels=FALSE){
    if(m < 1 || m>40)
        stop("m should be in the range [1,40]")
    if(k < 1)
        stop("there should be at least one super pixel")
    if(is.matrix(pic))
        pic <- simplify2array(list(pic,pic,pic))
    
    nRow <- dim(pic)[1]
    nCol <- dim(pic)[2]
    nPix <- prod(dim(pic)[1:2])
    S <- round(sqrt(nPix/k))
    Cxy0 <- round(as.matrix(expand.grid(seq(1+S/3,nRow,S),seq(1+S/3,nCol,S))))
    ## Image gradient
    G <- simplify2array(apply(pic,3,function(p) otoclass:::sobel(otoclass:::extendPic(p,3)), simplify=FALSE))
    ## Move centers
    Cxy <- t(apply(Cxy0, 1, function(xy){
        mv <- as.matrix(expand.grid(-1:1,-1:1))
        g0 <- apply(mv,1,function(ij) sqrt(sum((G[xy[[1]]+ij[[1]],xy[[2]]+ij[[2]],])^2)))
        xy + mv[which.min(g0),]
    }))
    ## Colors at centers
    Crgb <- t(apply(Cxy,1,function(xx) pic[xx[1],xx[2],]))
    Clab <- grDevices::convertColor(Crgb,"sRGB","Lab", scale.in=max.pic)
    C <- cbind(Clab,Cxy)
    ## Initialize labels and distances
    labels <- matrix(-1,dim(pic)[1],dim(pic)[2])
    dist <- matrix(Inf,dim(pic)[1],dim(pic)[2])
    Index <- matrix(seq_len(nPix),dim(pic)[1],dim(pic)[2])
    Row <- row(labels)
    Col <- col(labels)
    colLab <- grDevices::convertColor(cbind(as.vector(pic[,,1]),as.vector(pic[,,2]),as.vector(pic[,,3])),"sRGB","Lab", scale.in=max.pic)
    E <- Inf
    tol <- 0.1
    poti <- expand.grid(-S:S,-S:S)
    iter <- 0
    while(E > tol){
        iter <- iter + 1
        for(j in seq_len(nrow(C))){
           # cat(j,"\n")
            Ck <- C[j,]
            ## poti_k <- t(poti[poti[,1]+Ck[4] >= 1 & poti[,1]+Ck[4] <= nRow & poti[,2]+Ck[5] >= 1 & poti[,2]+Ck[5] <= nCol,]) + Ck[4:5]
            ## iAll <-  apply(poti_k,2,function(xx)Index[xx[1],xx[2]])
            i0 <- poti[,1]+Ck[4]
            j0 <- poti[,2]+Ck[5]
            toUse <- i0 >= 1 & i0 <= nRow & j0 >= 1 & j0 <= nCol
            i0 <- i0[toUse]
            j0 <- j0[toUse]
            iAll <- i0 + (j0-1) * nRow
            col_lab <- t(colLab[iAll,])
            dc <- sqrt(colSums((col_lab - Ck[1:3])^2))
            ds <- sqrt(colSums((rbind(i0,j0) - Ck[4:5])^2))
            D <- sqrt(dc^2 + (ds/S)^2*m^2)
            indx0 <- which(D < dist[iAll])            
            dist[iAll[indx0]] <- D[indx0]
            labels[iAll[indx0]] <- j
        }
        ## New centers
        iList <- split(Index,labels)
        Cnew <- do.call("rbind",lapply(iList,function(i) colMeans(cbind(colLab[i,],Row[i],Col[i]))))
        E <- mean(sqrt(rowSums((Cnew - C)^2)))
        C <- Cnew
        ## Label connectivity
        
        cat("Iteration:",iter,"E:",E,"\n")
    }
    finalRGB <- grDevices::convertColor(C[,1:3],"Lab","sRGB")
    superPic <- simplify2array(list(matrix(finalRGB[labels,1],nrow(labels),ncol(labels)),
                                    matrix(finalRGB[labels,2],nrow(labels),ncol(labels)),
                                    matrix(finalRGB[labels,3],nrow(labels),ncol(labels))))
    
    if(onlyImage)
        return(superPic)
    if(splitFinalPixels){
        C <- do.call("rbind",lapply(seq_len(nrow(C)), function(ki){
            clk <- contourLines(1:nrow(labels),1:ncol(labels),t(labels==ki),levels=1)
            do.call("rbind",lapply(clk,function(xx){
                pa <- otoclass:::polygon_area(cbind(xx$x,xx$y))
                if(pa == 0)
                    return(NULL)
                c(C[ki,1:3],X=attr(pa,"Cx"),Y=attr(pa,"Cy"),label=ki)
            }))
        }))
    }
  
    list(superpixels = C,
         RGB = finalRGB,
         labels = labels,
         distance = dist,
         pic = superPic)
}



##########################################################################################
## Local polynomial smoothing
##########################################################################################

localPolynomialSmoothing <- function(pic, n = 3, degree = 2){
    picE <- otoclass:::extendPic(toGreyScale(pic),(n)*2+1) / 255    
    H <- n+0.5
    if(degree == 0){
        ff <- ~1
    }else{
        ff <- as.formula(sprintf("~(%s)*(%s)",paste(sapply(seq_len(degree),function(d)sprintf("I(R^%d)",d)),collapse="+"),
                                 paste(sapply(seq_len(degree),function(d)sprintf("I(C^%d)",d)),collapse="+")))
    }
    d0 <- expand.grid(R=seq(-n,n,1),C=seq(-n,n,1))
    X0 <- model.matrix(ff,d0)
    K <- function(u) pmax(0,3 * (1 - u^2) / 4)
    W0 <- K(apply(d0,1,function(x) sqrt(sum((x/H)^2))))
    d <- d0[W0>0,]
    X <- X0[W0>0,]
    W <- diag(W0[W0>0],sum(W0>0))
    toBeta <- solve(t(X)%*%W%*%X)%*%t(X)%*%W

    toIndex<-function(i,j){
        ii <- d[,1] + i + n
        jj <- d[,2] + j + n
        Indx <- ii + (jj-1) * nrow(picE)    
        Indx
    }

    picEL <- convertColor(cbind(as.vector(picE),as.vector(picE),as.vector(picE)),"sRGB","Lab")[,1]

    picO <- outer(seq_len(nrow(pic)),
                  seq_len(ncol(pic)),
                  Vectorize(function(i,j){
                      ## g <- picE[toIndex(i,j)]
                      ## Lab <- convertColor(cbind(g,g,g),"sRGB","Lab")
                      ## Lstar <- (toBeta %*% Lab[,"L"])[1,1]
                      ## convertColor(cbind(pmax(0,pmin(100,Lstar)),0,0),"Lab","sRGB")[1,1]
                      (toBeta %*% picEL[toIndex(i,j)])[1,1]
                  }))
    picORGB <- convertColor(cbind(as.vector(picO),0,0),"Lab","sRGB")[,1]
    matrix(round(picORGB*255),nrow(pic),ncol(pic))
}

## pic <- getPixelMatrix("~/dtu/otolith_images/experiments/bg_remove_test/examples/8505362_ALA_RLX_XX.jpg",TRUE)
## mpic <- downsample(otoclass:::whiteBalance(pic),5)
## pp <- localPolynomialSmoothing(mpic,5,3)
## vx <- otoclass:::to01(otoclass:::floodfill(otoclass:::sobel(otoclass:::to01(otoclass:::extendPic(pp,3)))*255,0.1))
## sv <- which(vx > 200/255,TRUE)
## for(i in seq_len(nrow(sv)))
##     vx <- 1-otoclass:::floodfill((1-vx)*255,0.5,sv[i,1]-1,sv[i,2]-1)/255
## caMisc::imagePlot(otoclass:::to01(1-vx),"contain")

## image(vx < 0.05)

## caMisc::imagePlot(pp/255,"contain")
## vy <- mpic
## caMisc::imagePlot(mpic/255,"contain")

## caMisc::imagePlot(vx,"contain")

## hist(pp,0:255)


##########################################################################################

## rm <- rowMeans(mpic)
## rq70 <- quantile(rm,0.5)
## rq30 <- quantile(rm,0.3)
## cm <- colMeans(mpic)
## cq70 <- quantile(cm,0.7)
## cq30 <- quantile(cm,0.3)

## isBG <- row(mpic) %in% which(rm<rq30) | col(mpic) %in% which(cm < cq30) | col(mpic) %in% (round(ncol(mpic)/2) + seq(-round(ncol(mpic)/100),round(ncol(mpic)/100),1))

## mpic2 <- mpic
## mpic2[isBG] <- 0

## sp <- findSuperPixels(pic,300,20,max.pic=255,onlyImage=FALSE)

## caMisc::imagePlot(toGreyScale(sp$pic)/255,"contain")


## Data Lstar, u index, isKnownBG
##
fnSPClust <- function(par){
    #phi <- 1/(1+exp(-par$phiIn))
    mu1 <- cumsum(c(par$mu[1,1],exp(tail(par$mu[,1],-1))))
    mu2 <- cumsum(c(par$mu[1,2],exp(tail(par$mu[,2],-1))))
    ##cumsum(exp(par$mu[,2])) + par$mu[1,1]    
    p1 <- exp(c(par$p[,1],0)) / sum(exp(c(par$p[,1],0)))
    p2 <- exp(c(par$p[,2],0)) / sum(exp(c(par$p[,2],0)))
    REPORT(mu1)
    REPORT(mu2)
    REPORT(p1)
    REPORT(p2)
    ## f1 <- function(x) dautoreg(x, phi=phi,log=TRUE)
    ## nll1 <- -dseparable(f1,f1)(par$u - par$umu)
    nu <- exp(par$logNu)
    rho <- exp(par$logRho) #* nrow(pic0)
    C <- diag(1,nrow(D))
    ## ## AR
    C <- exp(par$sigma)^2 / (2 * rho) * exp(-rho * D)
    ## Matern
    ## C[lower.tri(D)] <- 2^(1-nu) / gamma(nu) * (D[lower.tri(D)]/rho*sqrt(2*nu))^nu *  besselK(x = D[lower.tri(D)]/rho*sqrt(2*nu), nu = rho)
    ## C[upper.tri(D)] <- 2^(1-nu) / gamma(nu) * (D[upper.tri(D)]/rho*sqrt(2*nu))^nu *  besselK(x = D[upper.tri(D)]/rho*sqrt(2*nu), nu = rho)
    ## C <- C * exp(par$sigma)^2
    nll1 <- -dmvnorm(par$u,par$umu,C,log=TRUE)
    ## Prior
    Lmod <- qlogis(squeeze(Lstar/100,1e-4))
    vB <- dmixnorm(Lmod,mu1,exp(par$logSd[,1]),p1)
    logpF <- -logspace_add(0,-par$u)
    logpB <- -par$u + logpF
    ## 1/(1+exp(-par$u)
    vF <- dmixnorm(Lmod,mu2,exp(par$logSd[,2]),p2)
    ## Known background
    nll2 <- vB + logpB #[isBG])+1e-8)
    ## Known foreground
    nll3 <- vF + logpF
    ## Unknown
    ##nll3 <- log(vB[!isBG] * (1-pF[!isBG]) + vF[!isBG] * pF[!isBG] +1e-8)
    nll4 <- logspace_add(vB+logpB,vF+logpF)
    ##sum(nll1) +
    posterior <- nll3 - nll4
    REPORT(posterior)   
    sum(nll1) + sum(-nll2[isBG]) + sum(-nll3[isFG]) + sum(-nll4[!isBG & !isFG])
}

##
## pic0 <- toGreyScale(sp$pic)
## ii <- sample(seq_len(nrow(sp$superpixels)
## Lstar <- sp$superpixels[,1]
## uIndx <- floor(sp$superpixels[,4]) + (floor(sp$superpixels[,5])-1) * nrow(pic0)
## D <- as.matrix(dist(sp$superpixels[,4:5])) / nrow(pic0)


## rm <- rowMeans(pic)
## rq70 <- quantile(rm,0.5)
## rq30 <- quantile(rm,0.3)
## cm <- colMeans(pic)
## cq70 <- quantile(cm,0.7)
## cq30 <- quantile(cm,0.3)
## isMid <- which(cm < cm[round(ncol(pic0)/2)] & seq_along(cm) > which.max(cm[1:round(ncol(pic0)/2)]) & seq_along(cm) < which.max(cm[round(ncol(pic0)/2):ncol(pic0)])+round(ncol(pic0)/2-1))
## isMid <- min(isMid):max(isMid)


## isBG <- uIndx %in% which(row(pic0) %in% which(rm<rq30) | col(pic0) %in% which(cm < cq30) | col(pic0) %in% (isMid))
## isFG <- uIndx %in% which(row(pic0) %in% which(rm>rq70) & col(pic0) %in% which(cm > cq70)) & !isBG

## plot(sp$superpixels[,4:5],col=ifelse(isBG,"red",ifelse(isFG,"green","black")),ylim=c(1,ncol(pic)),xlim=c(1,nrow(pic)))

## nM <- 5
## par <- list(mu = matrix(-2,nM,2),
##             logSd=matrix(0,nM,2),
##             p = matrix(0,nM-1,2),
##             umu = 0,
##             u=numeric(nrow(D)),#*pic0,
##             sigma = -2,
##             logNu = 5,
##             logRho = -2
##                                         #phiIn=-5
##             )
## obj <- MakeADFun(fnSPClust, par,random="u"
##                  )
##                  #)

## obj$fn()
## opt <- nlminb(obj$par,obj$fn,obj$gr, control = list(trace=1, iter.max=1000,eval.max=1000))


## spFP <- as.numeric(exp(obj$report(obj$env$last.par.best)$posterior) > 0.6)
## table(spFP)
## image(matrix(spFP[sp$labels],nrow(pic0),ncol(pic0)))

## table(spFP,isBG)
## table(spFP,isFG)

## table(isFG,isBG,spFP)


## ##########################################################################################

## library(OpenImageR)


## im = OpenImageR::readImage("~/dtu/otolith_images/bg_remove_test/examples/8505362_ALA_RLX_XX.jpg")


## #--------------
## # "slic" method
## #--------------
## pic2 <- otoclass:::whiteBalance(mpic)
## res_slic = superpixels(input_image = pic2,
##                        method = "slic",
##                        superpixel = 1000, 
##                        compactness = 20,
##                        return_slic_data = TRUE,
##                        return_labels = TRUE, 
##                        write_slic = "", 
##                        verbose = TRUE)


## plot_slic = OpenImageR::NormalizeObject(res_slic$slic_data)
## plot_slic = grDevices::as.raster(plot_slic)
## graphics::plot(plot_slic)

## ##########################################################################################


## pic <- downsample(toGreyScale(getPixelMatrix("~/dtu/otolith_images/bg_remove_test/examples/8505362_ALA_RLX_XX.jpg")/255),10)
## caMisc::imagePlot(pic/255)


## sp <- findSuperPixels(otoclass:::whiteBalance(pic),1000,20,onlyImage=FALSE,splitFinalPixels=TRUE)
## graphics::plot(as.raster(sp$pic))

## sp$superpixels


## ## Fit

## fnSPClust <- function(par){
##     #phi <- 1/(1+exp(-par$phiIn))
##     mu1 <- cumsum(c(par$mu[1,1],exp(tail(par$mu[,1],-1))))
##     mu2 <- cumsum(c(par$mu[1,2],exp(tail(par$mu[,2],-1))))
##     ##cumsum(exp(par$mu[,2])) + par$mu[1,1]    
##     p1 <- exp(c(par$p[,1],0)) / sum(exp(c(par$p[,1],0)))
##     p2 <- exp(c(par$p[,2],0)) / sum(exp(c(par$p[,2],0)))
##     REPORT(mu1)
##     REPORT(mu2)
##     REPORT(p1)
##     REPORT(p2)
##     kap <- advector(exp(as.vector(Xkappa %*% par$kappa)))
##     Qk <- new("adsparse", x=kap, i=Qd@i, p=Qd@p, Dim=Qd@Dim)
##     Q2 <- (Q + Qk) * exp(par$sigma) ##  %*% Qd
##     nll1 <- dgmrf(par$u,Q=Q2,log=TRUE)
##     umat <- (matrix(par$u,nrow(pic0),ncol(pic0)) + par$umu) #* exp(par$sigma)
##     REPORT(umat)
##     ## f1 <- function(x) dautoreg(x, phi=phi,log=TRUE)
##     ## nll1 <- dseparable(f1,f1)(par$u)
##     ## nu <- exp(par$logNu)
##     ## rho <- exp(par$logRho) #* nrow(pic0)
##     ## C <- diag(1,nrow(D))
##     ## ## AR
##     ## C <- exp(par$sigma)^2 / (2 * rho) * exp(-rho * D)
##     ## Matern
##     ## C[lower.tri(D)] <- 2^(1-nu) / gamma(nu) * (D[lower.tri(D)]/rho*sqrt(2*nu))^nu *  besselK(x = D[lower.tri(D)]/rho*sqrt(2*nu), nu = rho)
##     ## C[upper.tri(D)] <- 2^(1-nu) / gamma(nu) * (D[upper.tri(D)]/rho*sqrt(2*nu))^nu *  besselK(x = D[upper.tri(D)]/rho*sqrt(2*nu), nu = rho)
##     ## C <- C * exp(par$sigma)^2
##     ## nll1 <- -dmvnorm(par$u,par$umu,C,log=TRUE)
##     ## Prior
##     Lmod <- qlogis(squeeze(Lstar/100,1e-4))
##     vB <- dmixnorm(Lmod,mu1,exp(par$logSd[,1]),p1)
##     logpF <- -logspace_add(0,-umat[uIndx])
##     logpB <- -umat[uIndx] + logpF
##     ## 1/(1+exp(-par$u)
##     vF <- dmixnorm(Lmod,mu2,exp(par$logSd[,2]),p2)
##     ## Known background
##     nll2 <- vB + logpB #[isBG])+1e-8)
##     ## Known foreground
##     nll3 <- vF + logpF
##     ## Unknown
##     ##nll3 <- log(vB[!isBG] * (1-pF[!isBG]) + vF[!isBG] * pF[!isBG] +1e-8)
##     nll4 <- logspace_add(vB+logpB,vF+logpF)
##     ##sum(nll1) +
##     posterior <- nll3 - nll4
##     REPORT(posterior)   
##     sum(-nll1) + sum(-nll2[isBG]) + sum(-nll3[isFG]) + sum(-nll4[!isBG & !isFG])
## }

## ##
## pic0 <- toGreyScale(sp$pic)
## ii <- sample(seq_len(nrow(sp$superpixels)),3000)
## Lstar <- sp$superpixels[ii,1]
## uIndx <- floor(sp$superpixels[ii,4]) + (floor(sp$superpixels[ii,5])-1) * nrow(pic0)
## D <- as.matrix(dist(sp$superpixels[ii,4:5])) / nrow(pic0)
## Tri <- do.call("rbind",lapply(1:length(pic0), function(k){
##     ## Row
##     i <- (k-1) %% nrow(pic0) + 1
##     ## Column
##     j <- (k-1) %/% nrow(pic0) + 1    
##     nei <- c(
##     #(i-1) + (j-1-1)*nrow(pic0),
##     (i-1) + (j-0-1)*nrow(pic0),
##     #(i-1) + (j+1-1)*nrow(pic0),
##     (i) + (j-1-1)*nrow(pic0),
##     #(i) + (j-0-1)*nrow(pic0),
##     (i) + (j+1-1)*nrow(pic0),
##     #(i+1) + (j-1-1)*nrow(pic0),
##     (i+1) + (j-0-1)*nrow(pic0)#,
##     #(i+1) + (j+1-1)*nrow(pic0)
##     )
##     nei <- nei[nei > 0 & nei <= length(pic0)]
##     rbind(cbind(k,nei,-1),
##           c(k,k,length(nei)))
## }))
## Q <- Matrix::sparseMatrix(Tri[,1],Tri[,2],x=Tri[,3],repr="C")
## ##Qd <- Matrix::sparseMatrix(1:nrow(Q),1:nrow(Q),x=1,repr="C")
## smoothPic <- otoclass:::whiteBalance(otoclass:::to01(otoclass:::unsharpmask(otoclass:::extendPic(pic0,11),11,11))*255)
## Lval <- 1 - convertColor(cbind(as.vector(smoothPic),as.vector(smoothPic),as.vector(smoothPic)),"sRGB","Lab",scale.in=255)[,1] / 100
## Qd <- Matrix::sparseMatrix(1:nrow(Q),1:nrow(Q),x=1,repr="C")
## Xkappa <- model.matrix(~poly(R,2)+poly(C,2),data.frame(R=as.vector(row(pic0)), C=as.vector(col(pic0))))

## rm <- rowMeans(pic)
## rq70 <- quantile(rm,0.5)
## rq30 <- quantile(rm,0.2)
## cm <- colMeans(pic)
## cq70 <- quantile(cm,0.7)
## cq30 <- quantile(cm,0.2)
## isMid <- which(cm < cm[round(ncol(pic0)/2)] & seq_along(cm) > which.max(cm[1:round(ncol(pic0)/2)]) & seq_along(cm) < which.max(cm[round(ncol(pic0)/2):ncol(pic0)])+round(ncol(pic0)/2-1))
## isMid <- round(quantile(isMid,0.1)):round(quantile(isMid,0.9))


## isBG <- uIndx %in% which(row(pic0) %in% which(rm<rq30) | col(pic0) %in% which(cm < cq30) | col(pic0) %in% (isMid))
## isFG <- uIndx %in% which(row(pic0) %in% which(rm>rq70) & col(pic0) %in% which(cm > cq70)) & !isBG

## plot(sp$superpixels[ii,4:5],col=ifelse(isBG,"red",ifelse(isFG,"green","black")),ylim=c(1,ncol(pic)),xlim=c(1,nrow(pic)))

## nM <- 3
## par <- list(mu = matrix(-2,nM,2),
##             logSd=matrix(0,nM,2),
##             p = matrix(0,nM-1,2),
##             umu = 0,
##             u=as.vector(0*pic0),
##             sigma = 0,
##             ## logNu = 5,
##             kappa = numeric(ncol(Xkappa))
##             )
## obj <- MakeADFun(fnSPClust, par,random="u"
##                  )
##                                         #)


## obj$fn()
## obj$gr()
## opt <- nlminb(obj$par,obj$fn,obj$gr, control = list(trace=1, iter.max=1000,eval.max=1000))

## image(matrix(as.vector(exp(Xkappa %*% obj$env$parList()$kappa)),nrow(pic0),ncol(pic0)))


## spFP <- as.numeric(exp(obj$report(obj$env$last.par.best)$posterior) > 0.5)
## table(spFP)
## image(obj$report(par=obj$env$last.par.best)$umat>=0)

## pic2 <- pic * (obj$report(par=obj$env$last.par.best)$umat>=0)
## caMisc::imagePlot(pic2/255)

## table(spFP,isBG)
## table(spFP,isFG)

## table(isFG,isBG,spFP)



kmeans_SLIC <- function(pic, k = 1000, m = 20, g = 3L, connect = TRUE){
    if(length(dim(pic)) != 3 || (length(dim(pic)) == 3 && dim(pic)[3] != 3))
        pic <- array(pic[seq_len(prod(dim(pic)[1:2]))],dim = c(dim(pic)[1:2],3))
    .Call(C_kmeans_SLIC,pic=pic,k=as.integer(k),m=as.integer(m),g=as.integer(g),connect=as.logical(connect))
}
