
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
    if(length(lvls) == 1){ ## Only missing
        lvls <- c(lvls,"REF","ALT")
    }else if(length(lvls) == 2){ ## All the same
        lvls <- c(lvls,"ALT")
    }
    r <- apply(matrix(as.numeric(factor(xx,levels=lvls))-1,ncol=2),1,function(x) as.numeric(table(factor(x,1:(length(lvls)-1)))))
    rownames(r) <- lvls[-1]   
    r
}


## S3 class: Genotype
## names list of Locus
## attr ploidy: ploidy of the individual
##' @export
##' @method `[[` Genotype
`[[.Genotype` <- function(x,j,...){
    if(missing(j)) return(x)
    if(length(j) == 1){
        y <- unclass(x)[[j]]
        attr(y,"ploidy") <- attr(x,"plodiy")        
        ## class(y) <- "Genotype"
        return(y)
    }        
    y <- unclass(x)[j]
    attr(y,"ploidy") <- attr(x,"ploidy")
    class(y) <- "Genotype"
    y
}

##' @export
##' @method `[` Genotype
`[.Genotype` <- function(x,i,j,...){
    ## i: allele number
    if(missing(i)){
        y <- unclass(x)[j]
        attr(y,"ploidy") <- attr(x,"ploidy")
        class(y) <- "Genotype"
        return(y)
    }
    ## j: locus number
    lapply(unclass(x)[j],function(y) y[i])   
}
##' @export
##' @method print Genotype
print.Genotype <- function(x,...){
    isHomo <- sapply(x,function(y) any(y==attr(x,"ploidy")))
    isMis <- sapply(x,function(y) sum(y) < attr(x,"ploidy"))
    cat(sprintf("Genotype. Ploidy: %d. Loci: %d. Homozygous: %.1f%%. Missingness: %.1f%%.",attr(x,"ploidy"),length(x),sum(isHomo)/length(x)*100,sum(isMis)/length(x)*100),"\n")
}

num_allele <- function(x,...){
    UseMethod("num_allele")
}

##' @export
##' @method num_allele Genotype
num_allele.Genotype <- function(x,...){
    sapply(x,length)
}

##' @export
##' @method num_allele gen
num_allele.gen <- function(x,error_fun=warning,...){
    r <- do.call("rbind",lapply(x,num_allele))
    isOK <- apply(r,2,function(v) all(as.integer(v)==as.integer(mean(v))))
    if(!all(isOK)) error_fun("Individuals have different number of alleles for the same locus.")
    as.integer(apply(r,2,mean))
}


num_loci <- function(x,...){
    UseMethod("num_loci")
}

##' @export
##' @method num_loci Genotype
num_loci.Genotype <- function(x,...){
    length(x)
}

##' @export
##' @method num_loci gen
num_loci.gen <- function(x,error_fun=warning,...){
    r <- do.call("rbind",lapply(x,num_loci))
    isOK <- apply(r,2,function(v) all(as.integer(v)==as.integer(mean(v))))
    if(!all(isOK)) error_fun("Individuals have different number of loci.")
    as.integer(mean(r))
}

extract_locus <- function(x, l){
    UseMethod("extract_locus")
}

##' @export
##' @method extract_locus gen
extract_locus.gen <- function(x, l){
    r <- do.call("rbind",lapply(x,function(y) y[[l]]))
    rownames(r) <- names(x)
    r
}

add_locus <- function(x, lname, aname, warn=TRUE){
    UseMethod("add_locus")
}

##' @export
##' @method add_locus gen
add_locus.Genotype <- function(x, lname, aname, warn=TRUE){
    if(lname %in% names(x)){
        if(warn) warning("Locus already part of genotype.")
        return(x)
    }
    v <- numeric(length(aname))
    names(v) <- aname
    x[[lname]] <- v
    x
}


##' @export
##' @method add_locus gen
add_locus.gen <- function(x, lname, aname, warn=TRUE){
   lapply(x, add_locus, lname=lname, aname=aname,warn=warn)
}


get_ploidy <- function(x){
    UseMethod("get_ploidy")
}

##' @export
##' @method get_ploidy Genotype
get_ploidy.Genotype <- function(x){
    attr(x,"ploidy")
}

##' @export
##' @method get_ploidy gen
get_ploidy.gen <- function(x){
    v <- sapply(x,get_ploidy)
    table(v,dnn=NULL) / length(v)
}



## S3 class: gen
## list of Genotype

##' @method c gen
##' @export
c.gen <- function(..., lociAction=c("extend","keep","drop")){
    lociAction <- match.arg(lociAction)
    x <- list(...)
    if(length(x) == 1)
        return(x[[1]])   
    pop <- do.call("c", lapply(x, function(y) attr(y,"population")))
    nms <- do.call("c", lapply(x, names))
    if(any(sapply(x,function(y) !is.null(attr(y,"true_population"))))){
        tpop <- do.call("c", lapply(x, function(y){
            if(is.null(attr(y,"true_population"))){
                return(rep(NA,length(y)))
            }
            attr(y,"true_population")
        }))
    }else{
        tpop <- NULL
    }
    res <- do.call("c",lapply(x,unclass))
    class(res) <- "gen"
    attr(res,"population") <- pop
    attr(res,"true_population") <- tpop
    ## attr(res,"alleleNames") <- attr(x[[1]],"alleleNames")
    otherA <- setdiff(unique(sapply(x,function(y)names(attributes(y)))),c("class","population","true_population","names"))
    if(length(otherA) > 0)
        for(aa in otherA){
            attr(res,aa) <- lapply(x,function(y) attr(y,aa))
        }
    if(lociAction == "drop"){
        loci2use <- Reduce(intersect,lapply(res,names))
        res <- res[,,loci2use]
        attr(res,"true_population") <- tpop       
    }else if(lociAction == "extend"){
        loci2use <- Reduce(union,lapply(res,names))
        alleleNames <- lapply(loci2use,function(ll) unique(unlist(lapply(res,function(x) names(x[[ll]])))))
        res <- lapply(res, function(xx){
            newloci <- match(setdiff(loci2use,names(xx)),loci2use)
            for(i in newloci){
                xx <- add_locus(xx,loci2use[i],alleleNames[[i]])
            }
            xx
        })
        class(res) <- "gen"
        attr(res,"population") <- pop
        attr(res,"true_population") <- tpop
        ## attr(res,"alleleNames") <- attr(x[[1]],"alleleNames")
        res <- res[,,loci2use]
        attr(res,"true_population") <- tpop
        otherA <- setdiff(unique(sapply(x,function(y)names(attributes(y)))),c("class","population","true_population","names"))
        if(length(otherA) > 0)
            for(aa in otherA){
                attr(res,aa) <- lapply(x,function(y) attr(y,aa))
            }
    }
    names(res) <- nms
    res
}

##' @method `[[` gen
##' @export
`[[.gen` <- function(x, i, ...){
    unclass(x)[[i]]
}


##' @method `[` gen
##' @export
`[.gen` <- function(x,i,j,k,...){
    if(missing(i))
        i <- seq_along(x)
    if(missing(j) && missing(k)){
        y <- unclass(x)[i]
        class(y) <- c("gen")
        attr(y,"population") <- attr(x,"population")[i]
        if(!is.null(attr(x,"true_population")))
            attr(y,"true_population") <- attr(x,"true_population")[i]
        return(y)
    }
    if(!missing(j) && missing(k))
        r <- (lapply(x[i],function(y) y[j,]))
    if(missing(j)&& !missing(k))
        r <- (lapply(x[i],function(y) y[,k]))
    if(!missing(j)&& !missing(k))
        r <- (lapply(x[i],function(y) y[j,k]))
    class(r) <- c("gen")
    attr(r,"population") <- attr(x,"population")[i]
    if(!is.null(attr(x,"true_population")))
        attr(r,"true_population") <- attr(x,"true_population")[i]
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
    ## TODO: Move to C file (?)
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
    ##alleleNames <- lapply(as.list(as.data.frame(genoMat, stringsAsFactors = FALSE)),toAlleleMatrix, giveNames = TRUE)
    alleleMatList <- lapply(as.list(as.data.frame(genoMat, stringsAsFactors = FALSE)),toAlleleMatrix)
    names(alleleMatList) <- cn    
    Nallele <- sapply(alleleMatList, nrow)
    Gchar <- apply(genoMat,1,function(x) max(nchar(x)))
    ploidy <- unname(sapply(as.character(Gchar),function(x)switch(x,`2`=1,`3`=1,`4`=2,`6`=2,NA)))
    res <- lapply(seq_along(indiId),function(indi){
        r <- lapply(alleleMatList,function(y) y[,indi])
        names(r) <- cn
        class(r) <- "Genotype"
        attr(r,"ploidy") <- ploidy[indi]
        r
    })
    names(res) <- indiId
   
    ## aMat <- aperm(simplify2array(lapply(as.list(as.data.frame(genoMat, stringsAsFactors = FALSE)),toAlleleMatrix)),c(1,3,2))
    ## dimnames(aMat) <- list(NULL, cn, indiId)
    ## if(sort.loci){
    ##     aMat <- aMat[,order(cn),,drop=FALSE]
    ##     alleleNames <- alleleNames[order(cn)]
    ##     genoMat <- genoMat[,order(cn),drop=FALSE]
    ## }
    ## if(sort.individuals){
    ##     aMat <- aMat[,,order(indiId),drop=FALSE]
    ##     genoMat <- genoMat[order(indiId),,drop=FALSE]
    ##     inPop <- inPop[order(indiId)]
    ## }
    class(res) <- "gen"
    attr(res,"population") <- inPop
    ## attr(aMat,"alleleNames") <- alleleNames
    ##attr(aMat,"genotypes") <- genoMat
    res
}

to_txt <- function(x, ...){
    UseMethod("to_txt")
}

#' @export
to_txt.Genotype <- function(x, name="", ...){
    paste0(name,", ", paste(sapply(x,function(y){ v <- tail(c(0,0,rep(c(1,2),times=y)),2); sprintf("%03d%03d",v[1],v[2]) }), collapse = " "))
}

#' @export
to_txt.gen <- function(x, ...){
    nm <- names(x)
    if(is.null(nm))
        nm <- sprintf(sprintf("ID%%0%dd",ceiling(log10(length(x)))),seq_along(x))
    do.call(c,lapply(seq_along(x),function(i) to_txt(x[[i]],nm[i])))
}

#' @export
to_txt.list <- function(x, ...){
    ## Make sure they have the same size
    x2 <- do.call(c,x)
    xL <- split(x2, rep(seq_along(x), times = sapply(x,length)))
    r <- unname(do.call(c,lapply(xL, function(y) c("POP",to_txt(y)))))
    attr(r, "loci") <- names(x2[[1]])
    attr(r,"population") <- names(x)
    r
}

##' @export
write_gen <- function(x, file, info = "", ...){
    if(is(x,"gen"))
        x <- list(x)
    v <- to_txt(x)
    if(!is.null(attr(v,"population")))
        info <- paste0(info,ifelse(info=="","","; "),"Populations: ",paste(attr(v,"population"),collapse=", "))
    cat(info,"\n", file = file)
    cat(attr(v,"loci"), file = file, sep = "\n", append = TRUE)
    cat(v, file = file, sep = "\n", append = TRUE)
}
    
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

    

