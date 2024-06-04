##' @export
genetic_distance <- function(x, ...){
    UseMethod("genetic_distance")
}

##' @method genetic_distance list
##' @export
genetic_distance.list <- function(x, type = c("Theta","D","DR","Euclid","Gst","Fst","Chord","ChordTime"), ...){
    type <- match.arg(type)
    if(type == "Theta")
        return(genotype_distance_Theta(x))
    if(type == "D")
        return(genotype_distance_D(x))
    if(type == "DR")
        return(genotype_distance_DR(x))
    if(type == "Euclid")
        return(genotype_distance_Euclid(x))
    if(type == "Gst")
        return(genotype_distance_Gst(x))
    if(type == "Chord")
        return(genotype_distance_Chord(x))
    if(type == "ChordTime")
        return(genotype_distance_ChordTime(x))
    stop("Wrong distance type")
}

##' @method genetic_distance alleleFrequencyList
##' @export
genetic_distance.alleleFrequencyList <- function(x, ...){   
   genetic_distance(x$Populations, ...)
}


##' @method genetic_distance gen
##' @export
genetic_distance.gen <- function(x, ..., individual = FALSE){
    if(individual){
        pop <- dimnames(x)[[3]]
    }else{
        pop <- attr(x,"population")
    }
   genetic_distance(empirical_allele_frequencies(x, pop), ...)
}


## Nei (1972) DOI:10.1086/282771  D
genotype_distance_D <- function(alleleFrequencies){
    ##AL <- empirical_allele_frequencies(genotypes, attr(genotypes,"population"))$Population
    AL <- alleleFrequencies
    D <- outer(seq_along(AL),seq_along(AL), Vectorize(function(i,j) {
        X <- AL[[i]] #[-1,,drop=FALSE]
        Y <- AL[[j]] #[-1,,drop=FALSE]
        ##jx <- colSums(X^2)
        jx <- sapply(X,function(x) sum(x[-1]^2))
        ##jy <- colSums(Y^2)
        jy <- sapply(Y,function(x) sum(x[-1]^2))
        ##jxy <- colSums(X*Y)
        jxy <- sapply(seq_along(X), function(l) sum(X[[l]][-1] * Y[[l]][-1]))
        allNum <- !is.na(jx) & !is.na(jy) & !is.na(jxy)
        Jx <- mean(jx[allNum])
        Jy <- mean(jy[allNum])
        Jxy <- mean(jxy[allNum])
        -log(Jxy / sqrt(Jx * Jy))
        ##mean(-log(jxy/sqrt(jx*jy)), na.rm=TRUE)
    }))
    attr(D, "population") <- names(AL)    
    attr(D,"type") <- expression(D[Nei])
    class(D) <- "genotype_distance"
    D
}

## Roger's distance 1972
genotype_distance_DR <- function(alleleFrequencies){
    AL <- alleleFrequencies
    D <- outer(seq_along(AL),seq_along(AL), Vectorize(function(i,j) {
        X <- AL[[i]]#[-1,,drop=FALSE]
        Y <- AL[[j]]#[-1,,drop=FALSE]
        v <- sapply(seq_along(X), function(l){
            p1 <- X[[l]]
            p2 <- Y[[l]]
            if(any(is.na(p1)) || any(is.na(p1)))
                return(NA)
            sqrt(sum((p1-p2)^2) / attr(X,"ploidy"))
        })
        mean(v, na.rm = TRUE)
    }))
    attr(D, "population") <- names(AL)    
    attr(D,"type") <- expression(D[Roger])
    class(D) <- "genotype_distance"
    D
}

## B. S. Weir, C. Clark Cockerham (1984) DOI: 10.2307/2408641. Implementaiton of theta[w]
genotype_distance_Theta <- function(alleleFrequencies){
    ## Truncate negative to zero Roesti et al (2012)
    AL <- alleleFrequencies
    D <- outer(seq_along(AL),seq_along(AL), Vectorize(function(i,j) {
        X <- AL[[i]]#[-1,,drop=FALSE]
        Y <- AL[[j]]#[-1,,drop=FALSE]
        n1 <- attr(X,"samples")
        if(is.na(n1)) n1 <- 1
        n2 <- attr(Y,"samples")
        if(is.na(n2)) n2 <- 1
        n1 <- pmin(pmax(2,n1),1e6)
        n2 <- pmin(pmax(2,n2),1e6)
        r <- 2
        nbar <- (n1 + n2) / r
        nc <- (r*nbar - (n1^2 + n2^2)/(r*nbar)) /(r-1)      
        numer <- sapply(seq_along(X), function(l) {
            h1 <- attr(X[[l]],"heterozygote")[-1]
            h2 <- attr(Y[[l]],"heterozygote")[-1]
            p1 <- as.numeric(X[[l]])[-1]
            p2 <- as.numeric(Y[[l]])[-1]
            if(any(is.na(p1)) || any(is.na(p1)))
                return(NA)
            hbar <- (n1*h1 + n2*h2) / (r * nbar)
            pbar <- (n1 * p1 + n2 * p2) / (r * nbar)
            ssq <-  (n1 / nbar * (p1-pbar)^2 + n2 / nbar * (p2-pbar)^2) / (r -1)
            aa <- nbar/nc * (ssq - (1/(nbar-1)) * (pbar*(1-pbar) - (r-1)/r * ssq - hbar/4))
            sum(pmax(0,aa))
        })
        denom <- (sapply(seq_along(X), function(l) {
            h1 <- attr(X[[l]],"heterozygote")[-1]
            h2 <- attr(Y[[l]],"heterozygote")[-1]
            p1 <- as.numeric(X[[l]])[-1]
            p2 <- as.numeric(Y[[l]])[-1]
            if(any(is.na(p1)) || any(is.na(p1)))
                return(NA)
            hbar <- (n1*h1 + n2*h2) / (r * nbar)
            pbar <- (n1 * p1 + n2 * p2) / (r * nbar)
            ssq <-  (n1 * (p1-pbar)^2 + n2 * (p2-pbar)^2) / ((r -1) * nbar)
            aa <- nbar/nc * (ssq - (1/(nbar-1)) * (pbar*(1-pbar) - (r-1)/r * ssq - hbar/4))
            bb <- nbar / (nbar - 1) * ( pbar * (1 - pbar) - (r-1)/r * ssq  - (2*nbar-1)/(4*nbar) * hbar )
            cc <- hbar/2
            sum(pmax(0,aa), na.rm=TRUE) + sum(pmax(0,bb), na.rm=TRUE) + sum(pmax(0,cc), na.rm=TRUE)
        }))
        sum(numer) / sum(denom)
    }))
    attr(D, "population") <- names(AL)    
    attr(D,"type") <- expression(theta)
    class(D) <- "genotype_distance"
    D
}

## Nei '73
genotype_distance_Gst <- function(alleleFrequencies){
    ## https://gsajournals.figshare.com/articles/software/Supplemental_Material_for_Kitada_Nakamichi_and_Kishino_2021/14813490?file=28505820
    ## We only have allele frequencies, so we assume they are known without error (e.g, no bias correction)
    ##AL <- empirical_allele_frequencies(genotypes, attr(genotypes,"population"))$Population
    AL <- alleleFrequencies
    D <- outer(seq_along(AL),seq_along(AL), Vectorize(function(i,j) {
        ## G[ST] = (H[T] - H[S]) / H[T]
        X <- AL[[i]]#[-1,,drop=FALSE]
        Y <- AL[[j]]#[-1,,drop=FALSE]
        Hs <- sapply(seq_along(X), function(l){
            p1 <- as.numeric(X[[l]])
            p2 <- as.numeric(Y[[l]])
            if(any(is.na(p1)) || any(is.na(p1)))
                return(NA)
            1 - sum(p1^2 + p2^2) / 2
        })
        Ht <- sapply(seq_along(X), function(l){
            p1 <- as.numeric(X[[l]])
            p2 <- as.numeric(Y[[l]])
            pbar <- (p1 + p2) / 2
            if(any(is.na(p1)) || any(is.na(p1)))
                return(NA)
            ps <- (p1+p2)/2
            1 - sum(ps^2)
        })
        sum(Ht-Hs,na.rm=TRUE) / sum(Ht,na.rm=TRUE)
    }))
    attr(D, "population") <- names(AL)    
    attr(D,"type") <- expression(G[st])
    class(D) <- "genotype_distance"
    D
}



## Cavalii-Sforza and Edwards (1967) chord distance
## Heuch 1975 DOI: 10.2307/2529552
genotype_distance_Chord <- function(alleleFrequencies){
    ##AL <- empirical_allele_frequencies(genotypes, attr(genotypes,"population"))$Population
    AL <- alleleFrequencies
    D <- outer(seq_along(AL),seq_along(AL), Vectorize(function(i,j) {
        X <- AL[[i]]#[-1,,drop=FALSE]
        Y <- AL[[j]]#[-1,,drop=FALSE]
        cosTheta <- sapply(seq_along(X), function(l){
            p1 <- X[[l]]
            p2 <- Y[[l]]
            if(any(is.na(p1)) || any(is.na(p1)))
                return(NA)
            sum(sqrt(p1*p2))
        })
        sqrt(sum(2 * (1 - cosTheta), na.rm = TRUE))
    }))
    attr(D, "population") <- names(AL)    
    attr(D,"type") <- "Chord distance"
    class(D) <- "genotype_distance"
    D
}

## Cavalii-Sforza and Edwards (1967) chord distance
## Heuch 1975 DOI: 10.2307/2529552
genotype_distance_ChordTime <- function(alleleFrequencies){
    ##AL <- empirical_allele_frequencies(genotypes, attr(genotypes,"population"))$Population
    AL <- alleleFrequencies
    D <- outer(seq_along(AL),seq_along(AL), Vectorize(function(i,j) {
        X <- AL[[i]]#[-1,,drop=FALSE]
        Y <- AL[[j]]#[-1,,drop=FALSE]
         cosTheta <- sapply(seq_along(X), function(l){
            p1 <- X[[l]]
            p2 <- Y[[l]]
            if(any(is.na(p1)) || any(is.na(p1)))
                return(NA)
            sum(sqrt(p1*p2))
        })
        f <- 4 * sum(1-cosTheta) / (sum(lociUse) * (nrow(X)-1))
        -2 * log(1 - f)
    }))
    attr(D, "population") <- names(AL)    
    attr(D,"type") <- "Time"
    class(D) <- "genotype_distance"
    D
}


genotype_distance_Euclid <- function(alleleFrequencies){
    ##AL <- empirical_allele_frequencies(genotypes, attr(genotypes,"population"))$Population
    AL <- alleleFrequencies
      D <- outer(seq_along(AL),seq_along(AL), Vectorize(function(i,j) {
        X <- AL[[i]]#[-1,,drop=FALSE]
        Y <- AL[[j]]#[-1,,drop=FALSE]
        v <- sapply(seq_along(X), function(l){
            p1 <- X[[l]]
            p2 <- Y[[l]]
            if(any(is.na(p1)) || any(is.na(p1)))
                return(NA)
            sum((p1-p2)^2)
        })
        sum(sqrt(v), na.rm = TRUE)
    }))   
    attr(D, "population") <- names(AL)    
    attr(D,"type") <- "Euclidian distance"
    class(D) <- "genotype_distance"
    D
}

buildTree_upgma <- function(D){
    ## UPGMA
### Number of leaves
    nL <- nrow(D)
### Number of internal nodes
    nI <- nL - 1
### Total number of nodes
    nT <- 2*nI+1
    ## ALuse <- vector("list",nT)
    ## ALuse[1:nL] <- AL[1:nL]
    Joins <- as.list(1:nL)
    nInJoin <- rep(NA,nT)
    nInJoin[1:nL] <- 1
    isActive <- c(rep(TRUE,nL),rep(FALSE,nT-nL))
    branchLength <- numeric(nT)    
    Di <- D
    ## LOOP
    while(nrow(Di) > 1){
        step <- max(which(isActive))
        ## Clustering
        delta <- min(Di[lower.tri(Di)])
        nn <- as.vector(which(Di == delta & row(Di) < col(Di), TRUE)[1,])
        Joins[[step+1]] <- nnNodes <- which(isActive)[nn]
        nInJoin[step+1] <- sum(nInJoin[nnNodes])
        ## Branch length
        branchLength[step+1] <- delta/2
        ## New distance
        Dx <- Di[-nn,-nn,drop=FALSE]
        if(nrow(Di) > 2){
        Dnew <- colSums(Di[nn,-nn,drop=FALSE] * matrix(nInJoin[nnNodes] / nInJoin[step+1],length(nn),ncol(Di)-length(nn)))
        Di <- cbind(rbind(Dx,c(Dnew)),c(Dnew,0))
        }else{
            Di <- matrix(0,1,1)
        }
        ## Prep
        isActive[nnNodes] <- FALSE
        isActive[step+1] <- TRUE
    }    
    isLeaf <- which(sapply(Joins,length)==1)
   
    nodeOrder <- c(length(Joins))
    for(p in rev(seq_along(Joins))){
        if(is.na(match(p,isLeaf))){
            pii <- which(nodeOrder==p)
            ch <- Joins[[p]]
            ch <- ch[order(nInJoin[ch])]
            nodeOrder <- c(nodeOrder[seq_len(pii-1)],ch[1],p,ch[2],nodeOrder[setdiff(seq_along(nodeOrder),seq_len(pii))])
        }
    }
    nodeY <- branchLength / max(branchLength)    
    nodeX <- rep(NA,length(Joins))  
    leafPos <- seq(0,1,len=length(isLeaf))
    leafOrder <- order(intersect(nodeOrder,isLeaf))
    nodeX[isLeaf] <- leafPos[leafOrder]
    for(i in setdiff(seq_along(Joins),isLeaf)){            
        nodeX[i] <- mean(nodeX[Joins[[i]]])
    }
    list(Joins = Joins,branchLength = branchLength, nodeY=nodeY, nodeX=nodeX,         
         isLeaf = isLeaf, names = attr(D,"population"), distanceType = attr(D,"type"))        
}


buildTree_nj <- function(D, rooted = TRUE){
    ## Neighbour joining
### Number of leaves
    nL <- nrow(D)
### Number of internal nodes
    nI <- nL - 1
### Total number of nodes
    nT <- 2*nI+1
    ## ALuse <- vector("list",nT)
    ## ALuse[1:nL] <- AL[1:nL]
    Joins <- as.list(1:nL)
    nInJoin <- rep(NA,nT)
    nInJoin[1:nL] <- 1
    isActive <- c(rep(TRUE,nL),rep(FALSE,nT-nL))
    branchLength <- numeric(nT)
    Di <- D
    D2Q <- function(Dt) (nrow(Dt)-2) * Dt - rowSums(Dt)[row(Dt)] - colSums(Dt)[col(Dt)]
    Qi <- D2Q(Di)
    ## LOOP
    while(nrow(Di) > 2){
        step <- max(which(isActive))
        ## Clustering
        nn <- as.vector(which(Qi == min(Qi[lower.tri(Qi)]) & row(Qi) > col(Qi), TRUE)[1,])
        Joins[[step+1]] <- nnNodes <- which(isActive)[nn]
        nInJoin[step+1] <- sum(nInJoin[nnNodes])
        ## Branch lengths to parent
        branchLength[nnNodes[1]] <- 0.5 * Di[nn[1],nn[2]] + 1/(2 * (nrow(Di)-2)) * (sum(Di[nn[1],]) - sum(Di[nn[2],]))
        branchLength[nnNodes[2]] <- Di[nn[1],nn[2]] - branchLength[nnNodes[1]]
        if(!rooted & nrow(Di) == 3){
            nnx <- setdiff(1:3,nn)
            branchLength[nnNodes[2]] <- Di[nn[1],nnx] - branchLength[nnNodes[1]]
            Joins[[step+1]] <- c(Joins[[step+1]],which(isActive)[nnx])
        }
        ## Avoid negative branch lengths
        branchLength[nnNodes] <- branchLength[nnNodes] - pmin(0,min(branchLength[nnNodes]))
        ## New distance
        Dx <- Di[-nn,-nn,drop=FALSE]
        if(nrow(Di) > 2){
            Dnew <- sapply(setdiff(seq_len(nrow(Di)), nn), function(i) 0.5 * (sum(Di[nn,i]) - Di[nn[1],nn[2]]))           
            Di <- cbind(rbind(Dx,c(Dnew)),c(Dnew,0))
        }else{
            Di <- matrix(0,1,1)
        }
        Qi <- D2Q(Di)
        ## Prep
        isActive[nnNodes] <- FALSE
        isActive[step+1] <- TRUE
    }
    if(rooted){
        ## Join the last two
        step <- max(which(isActive))
        Joins[[step+1]] <- nnNodes <- which(isActive)
        nInJoin[step+1] <- sum(nInJoin[nnNodes])
        branchLength[nnNodes[1]] <- 0.5 * Di[1,2] + 0.5 * (sum(Di[1,]) - sum(Di[2,]))
        branchLength[nnNodes[2]] <- Di[1,2] - branchLength[nnNodes[1]]
        ## Avoid negative branch lengths
        branchLength[nnNodes] <- branchLength[nnNodes] - pmin(0,min(branchLength[nnNodes]))
    }
    isLeaf <- which(sapply(Joins,length)==1)

    if(rooted){
        treeDepthAbove <- numeric(0)
        nodeOrder <- c(length(Joins))
        pcc <- do.call("rbind",lapply(setdiff(seq_along(Joins),isLeaf),function(p) c(p,Joins[[p]])))
        for(p in rev(seq_along(Joins))){
            fp <- which(pcc[,2]==p | pcc[,3]==p)
            if(length(fp)==0){
                treeDepthAbove[p] <- 0
            }else{
                treeDepthAbove[p] <- sum(treeDepthAbove[pcc[fp,1]]) + branchLength[p]
            }
            if(is.na(match(p,isLeaf))){
                pii <- which(nodeOrder==p)
                ch <- Joins[[p]]
                ch <- ch[order(nInJoin[ch])]
                nodeOrder <- c(nodeOrder[seq_len(pii-1)],ch[1],p,ch[2],nodeOrder[setdiff(seq_along(nodeOrder),seq_len(pii))])
            }
        }
        nodeY <- 1 - treeDepthAbove / max(treeDepthAbove)    
        nodeX <- rep(NA,length(Joins))  
        leafPos <- seq(0,1,len=length(isLeaf))
        leafOrder <- order(intersect(nodeOrder,isLeaf))
        nodeX[isLeaf] <- leafPos[leafOrder]
        for(i in setdiff(seq_along(Joins),isLeaf)){            
            nodeX[i] <- mean(nodeX[Joins[[i]]])
        }
    }else{
        stop("Not ready")
        nodeX <- rep(NA,length(Joins))
        nodeY <- rep(NA,length(Joins))
        anglesIn <- rep(NA,length(Joins))
        ## First one
        p <- length(Joins)
        nodeX[p] <- 0
        nodeY[p] <- 0
        anglesIn[p] <- 0
        ## Set children
        for(p in rev(seq_along(Joins))){
            ## Find self
            x0 <- nodeX[p]
            y0 <- nodeY[p]
            angle0 <- anglesIn[p] + pi
            ## angles for children
            anglesIn[Joins[[p]]] <- (angle0 + 2*pi/(length(Joins[[p]])+1) * 1:(length(Joins[[p]]))) %% (2*pi)
            ## Set childen
            nodeX[Joins[[p]]] <- x0 + branchLength[Joins[[p]]] * cos(anglesIn[Joins[[p]]])
            nodeY[Joins[[p]]] <- y0 + branchLength[Joins[[p]]] * sin(anglesIn[Joins[[p]]])      
        }
    }
    list(Joins = Joins,branchLength = branchLength, nodeY=nodeY, nodeX=nodeX,         
         isLeaf = isLeaf, names = attr(D,"population"), distanceType = attr(D,"type"))
}

pcoa <- function(D){
    D2 <- D^2
    C <- diag(1,nrow(D2),ncol(D2)) - matrix(1/nrow(D2),nrow(D2),ncol(D2))
    B <- -0.5 * C %*% D2 %*% C
    ei <- eigen(B)
    r <- sum(ei$values > 0)
    X <- ei$vectors[,1:r] %*% diag(sqrt(ei$values[1:r]),r)
    w <- ei$values[1:r] / sum(ei$values[1:r])
    colnames(X) <- paste("Coord.",1:r)
    list(X=X,w=w, names = attr(D,"population"))
}

##' @method plot genotype_distance
##' @export
plot.genotype_distance <- function(D, type = c("cladogram","dendrogram","MDS","PCoA"), treeType = c("nj_root","upgma"), horizontal=TRUE, plot = TRUE, point_cex = 1.5, point_pch = 16, lwd = 5, legend.pos="topright", text_col = "black",text_cex=1,text_font=1, text_srt = ifelse(horizontal,0,90), text_adj = if(horizontal){NULL}else{1.25},text_pos = if(horizontal){4}else{NULL}, axes = FALSE, xlab,...){
    type <- match.arg(type)
    treeType <- match.arg(treeType)
    if(type == "cladogram" || type == "dendrogram"){        
        ## Plot the tree
        if(plot)
            plot.new()
                                        #leafOrder <- intersect(unlist(Joins[-isLeaf]),isLeaf)
        usr <- par("usr")
        x1 <- ifelse(horizontal,-1,0)#usr[1] + 0.05 * diff(usr[1:2])
        x2 <- ifelse(horizontal,0,1)#usr[2] - 0.05 * diff(usr[1:2])
        y1 <- 0#usr[3] + 0.05 * diff(usr[3:4])
        y2 <- 1#usr[4] - 0.05 * diff(usr[3:4])

        if(treeType == "upgma"){
            tree <- buildTree_upgma(D)
        }else if(treeType == "nj_root"){
            tree <- buildTree_nj(D,TRUE)
        }

        if(type %in% c("cladogram","dendrogram")){
            if(horizontal){
                wx <- max(-tree$nodeY) + max(strwidth(tree$names,cex=text_cex,font=text_font)*1.1)
                hx <- 0
            }else{
                wx <- 0
                hx <- min(tree$nodeY) - max(strwidth(tree$names,cex=text_cex,font=text_font)*1.1)
            }
        }else{
            wx <- hx <- 0
        }
        
        par(usr = c(x1-0.1,pmax(x2,wx)+0.1,pmin(y1,hx)-0.1,y2+0.1))
        if(horizontal){
            tmpY <- tree$nodeY
            tmpX <- tree$nodeX
            tree$nodeX <- -tmpY
            tree$nodeY <- tmpX
        }
        if(plot){      
            if(type == "cladogram"){
                if(sum(!is.na(match(tree$Joins[[length(tree$Joins)]],tree$isLeaf))) == 1){
                    vv <- !is.na(match(tree$Joins[[length(tree$Joins)]],tree$isLeaf))
                    if(horizontal){
                        tree$nodeY[length(tree$Joins)] <- 2/3 * tree$nodeY[tree$Joins[[length(tree$Joins)]][vv]] + 1/3 * tree$nodeY[tree$Joins[[length(tree$Joins)]][!vv]]
                    }else{
                        tree$nodeX[length(tree$Joins)] <- 2/3 * tree$nodeX[tree$Joins[[length(tree$Joins)]][vv]] + 1/3 * tree$nodeX[tree$Joins[[length(tree$Joins)]][!vv]]
                    }
                }
                ## Plot leafs
                points(tree$nodeX,tree$nodeY, col = "black", cex = point_cex, pch = point_pch)
                for(i in seq_along(tree$Joins)){
                    if(length(tree$Joins) > 1){
                        j <- tree$Joins[[i]][1]
                        segments(tree$nodeX[j],tree$nodeY[j],tree$nodeX[i],tree$nodeY[i], lwd = lwd)
                        k <- tree$Joins[[i]][2]
                        segments(tree$nodeX[k],tree$nodeY[k],tree$nodeX[i],tree$nodeY[i], lwd = lwd)
                    }
                }
                lapply(seq_along(tree$names), function(i) text(tree$nodeX[tree$isLeaf[i]],tree$nodeY[tree$isLeaf[i]],tree$names[[i]],pos=text_pos,font=2,cex=1, col = text_col[(i-1)%%length(text_col)+1], srt = text_srt,adj=text_adj))
                ## prettyNum(unique(branchLength))                                
            }else if(type == "dendrogram"){
                points(tree$nodeX,tree$nodeY, col = "black", cex = point_cex, pch = point_pch)
                 for(i in seq_along(tree$Joins)){
                     if(length(tree$Joins) > 1){                         
                         j <- tree$Joins[[i]][1]
                         k <- tree$Joins[[i]][2]
                         if(horizontal){
                             segments(tree$nodeX[j],tree$nodeY[j],tree$nodeX[i],tree$nodeY[j], lwd = lwd)
                             segments(tree$nodeX[i],tree$nodeY[j],tree$nodeX[i],tree$nodeY[i], lwd = lwd)
                             segments(tree$nodeX[k],tree$nodeY[k],tree$nodeX[i],tree$nodeY[k], lwd = lwd)
                             segments(tree$nodeX[i],tree$nodeY[k],tree$nodeX[i],tree$nodeY[i], lwd = lwd)
                         }else{
                             segments(tree$nodeX[j],tree$nodeY[j],tree$nodeX[j],tree$nodeY[i], lwd = lwd)
                             segments(tree$nodeX[j],tree$nodeY[i],tree$nodeX[i],tree$nodeY[i], lwd = lwd)
                             segments(tree$nodeX[k],tree$nodeY[k],tree$nodeX[k],tree$nodeY[i], lwd = lwd)
                             segments(tree$nodeX[k],tree$nodeY[i],tree$nodeX[i],tree$nodeY[i], lwd = lwd)
                         }
                     }
                 }
                 lapply(seq_along(tree$names), function(i) text(tree$nodeX[tree$isLeaf[i]],tree$nodeY[tree$isLeaf[i]],tree$names[[i]],pos=text_pos,font=text_font,cex=text_cex, col = text_col[(i-1)%%length(text_col)+1], srt = text_srt,adj=text_adj))
            }
            tt <- axisTicks((par("usr")[1:2]-c(0,wx)) * max(tree$branchLength),FALSE)
            if(axes){
                axis(ifelse(horizontal,1,2), tt / max(tree$branchLength),ifelse(horizontal,-1,1) * tt, lwd = 5, font = 2)
            }
            if(missing(xlab))
               xlab = tree$distanceType
            mtext(xlab,ifelse(horizontal,1,2),line=3,font=2)
        }
        return(invisible(tree))
    }else if(type == "MDS" || type == "PCoA"){
        v <- pcoa(D)
        if(plot){
            labs <- paste(colnames(v$X)[1:2],sprintf("(%.1f%%)",v$w[1:2]*100))
            plot(v$X[,1:2], xlab=labs[1], ylab = labs[2], col = factor(v$names), cex = point_cex, pch = point_pch, axes = axes, ...)
            box()
            legend(legend.pos, legend = v$names, col = 1:nlevels(factor(v$names)), pch = point_pch)
        }
        return(v)
    }else{
        stop("Wrong plot type")
    }
}

##' @method print genotype_distance
##' @export
print.genotype_distance <- function(x,...){
    x2 <- x
    diag(x2) <- 0
    tab <- as.matrix(x2)
    rownames(tab) <- colnames(tab) <- attr(x,"population")
    tab[upper.tri(tab)] <- NA
    cat(paste("Genetic distance metric:",attr(x,"type")),"\n\n")
    xx <- format(unclass(x2), digits = 2, justify = "right")
    nn <- max(nchar(xx))
    print.table(tab,digits = 2, zero.print = paste(rep("-",nn),collapse=""))
}
