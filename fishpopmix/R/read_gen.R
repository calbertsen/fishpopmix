
na2 <- function(x, to){
    ifelse(is.na(x), to, x)
}

toAlleleMatrix <- function(G, giveNames = FALSE){
    if(is.numeric(G))
        G <- sprintf("%06d",na2(G,0))
    xx <- do.call("rbind",lapply(strsplit(G,""),function(x){
        if(length(x) %in% c(4,6)){
            n <- length(x)
            v <- c(paste0(x[seq_len(n/2)],collapse=""),paste0(x[n/2 + seq_len(n/2)],collapse=""))
            return(sprintf("%03d",as.numeric(v)))
        ## }else if(length(x) <= 3){
        ##     ##stop("Haploid data are not currently supported")
        ##     v <- paste0(x,collapse="")
        ##     return(c("000",sprintf("%03d",as.numeric(v))))
        }else{
            stop("Wrong number of allele digits")
        }
    }
    ))
    lvls <- sort(unique(c("000",as.vector(xx))))
    if(giveNames)
        return(lvls[-1])
    apply(matrix(as.numeric(factor(xx,levels=lvls))-1,ncol=2),1,function(x) as.numeric(table(factor(x,1:(length(lvls)-1)))))
}


##' Prepare allele matrix for PCA analysis
##'
##' @param x output from \\link{read.gen}
##' @param alleleMeans Mean allele frequencies to impute. Is calculated if missing.
##' @return allele data for PCA analysis
##' @author Christoffer Moesgaard Albertsen
gen2PCA <- function(x, alleleMeans){
    x2 <- apply(x,c(2,3),function(x) x/sum(x))
    if(missing(alleleMeans))
        alleleMeans <- apply(x2,1:2,mean, na.rm=TRUE)
    md <- unique(which(is.nan(x2), arr.ind = TRUE)[,2:3])
    for(i in seq_len(nrow(md))){
        x2[,md[i,1],md[i,2]] <- alleleMeans[,md[i,1]]
    }
    for(i in seq_len(dim(x2)[3]))
        x2[,,i] <- x2[,,i] - alleleMeans
    r <- x2[seq_len(dim(x2)[1]-1),,]
    attr(r,"alleleMeans") <- alleleMeans
    r
}

##' Read genepop data files
##'
##' @param f file name
##' @param pop.names population names. If missing, the ID of the last individual is used
##' @param sort.loci Should loci be sorted by names?
##' @param sort.individuals Should individuals be sorted by id?
##' @return an allele array
##' @author Christoffer Moesgaard Albertsen
##' @export
read_gen <- function(f, pop.names, sort.loci = FALSE, sort.individuals = FALSE, NAlleleKeep = NA){
    l <- readLines(f, warn = FALSE)
    l <- l[!grepl("^[[:space:]]*$",l)]
    l <- gsub("[[:space:]]+"," ",l)
    popPlace <- grep("^[[:space:]]*pop[[:space:]]*$",tolower(l))
    if(length(popPlace) == 0)
        stop("The file must have a POP")
    if(popPlace[1] < 2)
        stop("Loci names must come before the first population.")
    popStart <- popPlace+1
    popEnd <- unique(c(popPlace[-1] - 1, length(l)))
    if(missing(pop.names)){
        pop.names <- gsub("[[:space:]]+$","",gsub("(^.+)(,[[:space:]].+$)","\\1",l[popEnd]))
    }else if(is.function(pop.names)){
        pop.names <- pop.names(gsub("[[:space:]]+$","",gsub("(^.+)(,[[:space:]].+$)","\\1",l[popEnd])))
    }
    if(length(pop.names) != length(popStart))
        stop("Number of pop.names does not match the number of POP")
    ## Get data
    NinPop <- popEnd - popStart + 1
    inPop <- rep(pop.names, times = NinPop)
    dataIndx <- unlist(sapply(seq_along(popStart), function(i) seq(popStart[i], popEnd[i], by = 1)))
    indiId <- gsub("(^.+)(,[[:space:]].+$)","\\1",l[dataIndx])
    dataList <- strsplit(gsub("^.+,[[:space:]]*","",l[dataIndx]),"[[:space:]]+")
    genoMat <- do.call("rbind",dataList)
    ## Get colnames before converting to frequencies
    ## Case 1, one line with many colnames
    if(popPlace[1] == 2 || popPlace[1] == 3){
        cn <- unlist(strsplit(l[popPlace[1]-1],","))
    }else{## Case 2, many lines with one colname
        ## Find number of columns
        cn <- tail(l[seq_len(popPlace[1]-1)],ncol(genoMat))
    }
    alleleNames <- lapply(as.list(as.data.frame(genoMat, stringsAsFactors = FALSE)),toAlleleMatrix, giveNames = TRUE)
    Nallele <- sapply(alleleNames, length)
    if(!all(Nallele == Nallele[1])){
        if(is.na(NAlleleKeep))
            NAlleleKeep <- as.numeric(names(which.max(table(Nallele))))
        aIndx <- unname(which(Nallele == NAlleleKeep))
        warning(sprintf("Only loci with the same number of alleles are currently supported. %s %s removed.",paste(cn[Nallele != NAlleleKeep],collapse=", "),ifelse(sum(Nallele != NAlleleKeep)==1,"was","were")))
        ## Remove
        cn <- cn[aIndx]
        genoMat <- genoMat[,aIndx,drop=FALSE]
        alleleNames <- alleleNames[aIndx]
        Nallele <- Nallele[aIndx]
        dataList <- dataList[aIndx]
    }
    names(alleleNames) <- cn
    dimnames(genoMat) <- list(indiId, cn)
    aMat <- aperm(simplify2array(lapply(as.list(as.data.frame(genoMat, stringsAsFactors = FALSE)),toAlleleMatrix)),c(1,3,2))
    dimnames(aMat) <- list(NULL, cn, indiId)
    if(sort.loci){
        aMat <- aMat[,order(cn),,drop=FALSE]
        alleleNames <- alleleNames[order(cn)]
        genoMat <- genoMat[,order(cn),drop=FALSE]
    }
    if(sort.individuals){
        aMat <- aMat[,,order(indiId),drop=FALSE]
        genoMat <- genoMat[order(indiId),,drop=FALSE]
        inPop <- inPop[order(indiId)]
    }
    class(aMat) <- c("gen","array")
    attr(aMat,"population") <- inPop
    attr(aMat,"alleleNames") <- alleleNames
    ##attr(aMat,"genotypes") <- genoMat
    aMat
}

gen2list <- function(x){
    d <- dim(x)
    dnm <- dimnames(x)
    class(x) <- "array"
    xL <- lapply(split(x, slice.index(x,3)), matrix, nrow = d[[1]], ncol = d[[2]], dimnames=dnm[1:2])
    names(xL) <- dnm[[3]]
    xL
}

##' @method c gen
##' @export
c.gen <- function(..., dropLoci=FALSE){
    x <- list(...)
    pop <- do.call("c", lapply(x, function(y) attr(y,"population")))
    res <- simplify2array(do.call("c",lapply(x, gen2list)))
    class(res) <- c("gen","array")
    attr(res,"population") <- pop
    attr(res,"alleleNames") <- attr(x[[1]],"alleleNames")
    otherA <- setdiff(unique(sapply(x,function(y)names(attributes(y)))),c("dim","dimnames","class","population","alleleNames"))
    for(aa in otherA){
        attr(res,aa) <- lapply(x,function(y) attr(y,aa))
    }
    res
}

##' @method `[[` gen
##' @export
`[[.gen` <- function(x, i, l, ...){
    if(missing(l))
        l <- seq_len(dim(x)[2])
    class(x) <- "array"
    y <- as.array(x)[,l,i,drop=FALSE]
    class(y) <- c("gen","array")
    attr(y,"population") <- attr(x,"population")[i]
    attr(y,"alleleNames") <- attr(x,"alleleNames")
    y
}

##' @method `[` gen
##' @export
`[.gen` <- function(x,i,j,k,...){    
    y <- x
    class(y) <- "array"
    y <- y[i,j,k,drop=FALSE]
    class(y) <- c("gen","array")
    if(!missing(k)){
        attr(y,"population") <- attr(x,"population")[k]
    }else{
        attr(y,"population") <- attr(x,"population")
    }
    if(!missing(i)){
        attr(y,"alleleNames") <- attr(x,"alleleNames")[i]
    }else{
        attr(y,"alleleNames") <- attr(x,"alleleNames")
    }
    y    
}

write.gen <- function(x, file){


## toGen <- function(am, file, P, info = ""){
##     if(!is.list(am))
##         am <- list(am)
##     v <- sprintf("P_%s_=%f",names(P),P)
##     cat(paste(c("Simulated data",v,info), collapse = "; "),"\n", file = file)
##     cat(colnames(am[[1]]), file = file, sep = "\n", append = TRUE)
##     for(pp in seq_along(am)){
##         cat("POP", file = file, sep = "\n", append = TRUE)
##         rn <- rownames(am[[pp]])
##         a <- sapply(seq_len(nrow(am[[pp]])), function(i){
##             cat(paste0(rn[i],", ", paste(am[[pp]][i,], collapse = " "),"\n"), file = file, append = TRUE)
##         })
##     }
## }

    
}
