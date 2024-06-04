##' @export
ld_matrix <- function(x, measure, structured, ...){
    UseMethod("ld_matrix")
}

##' If structured = TRUE: Nei & Li (1973) DOI: 10.1093/genetics/75.1.213
##' @method ld_matrix gen
##' @export
ld_matrix.gen <- function(x,measure = c("r2","r2hap","D'","D'hap","Chisq","Chisqdf"), structured = TRUE, ...){
    ploidy <- get_ploidy(x)
    if(structured){
        popi <- attr(x,"population")
    }else{
        popi <- rep("Total",length(x))
    }
    if(as.numeric(names(ploidy[which.max(ploidy)]))!=2 || any(ploidy != 1))
        stop("ld_matrix is urrently only implemented for diploid individuals.")
    measure <- match.arg(measure)
    if(measure == "r2"){
        FUN <- ld_metric_r2
    }else if(measure == "r2hap"){
        FUN <- ld_metric_r2hap
    }else if(measure == "D'"){
        FUN <- ld_metric_D
    }else if(measure == "D'hap"){
        FUN <- ld_metric_Dhap
    }else if(measure == "Chisq"){
        FUN <- ld_metric_Chisq
    }else if(measure == "Chisqdf"){
        FUN <- ld_metric_Chisqdf
    }else{
        stop("Unknown LD metric.")
    }
    l1 <- seq_len(num_loci(x))
    XX <- lapply(l1, function(i) extract_locus(x,i))
    LDM <- outer(l1,l1, Vectorize(function(i,j){
        if(j >= i) return(NA)
        Xi <- XX[[i]]
        Xj <- XX[[j]]
        FUN(Xi,Xj, popi)
    }))
    dimnames(LDM) <- list(names(x[[1]]),names(x[[1]]))
    class(LDM) <- c("ld_matrix")
    attr(LDM,"structured") <- structured
    LDM
}

ld_test_chisq <- function(Xi, Xj, pop){
    optI <- sort(sapply(fishpopmix:::nexcom(2,ncol(Xi)), function(x) paste(rep(LETTERS[1:length(x)],times=x),collapse="")))
    optJ <- sort(sapply(fishpopmix:::nexcom(2,ncol(Xj)), function(x) paste(rep(LETTERS[1:length(x)],times=x),collapse="")))
    MakeOne <- function(popi){
        iuse <- rowSums(Xi) == 2 & rowSums(Xj) == 2 & seq_len(nrow(Xi)) %in% popi
        Gi <- apply(Xi[iuse,,drop=FALSE],1,function(x) paste(rep(LETTERS[1:length(x)],times=x),collapse=""))
        Gj <- apply(Xj[iuse,,drop=FALSE],1,function(x) paste(rep(LETTERS[1:length(x)],times=x),collapse=""))
        tryCatch({chisq.test(table(factor(Gi,optI),factor(Gj,optJ)))$p.value},error=function(e)NA)
    }
    popl <- split(seq_len(nrow(Xi)), pop)
    v <- log(sapply(popl, MakeOne))
    X <- -2 * sum(v,na.rm=TRUE)
    1-pchisq(X,2*sum(!is.na(v)))
}

ld_metric_r2 <- function(Xi, Xj, pop){
    MakeOne <- function(popi){
        iuse <- rowSums(Xi) == 2 & rowSums(Xj) == 2 & seq_len(nrow(Xi)) %in% popi
        r <- sum(outer(seq_len(ncol(Xi)),seq_len(ncol(Xj)), Vectorize(function(a,b){
            ## Only for diploid at the moment
            tab <- table(factor(Xi[iuse,a],0:2),factor(Xj[iuse,b],0:2))
            tab <- tab / sum(tab)
            p_A <- 0.5 * sum(tab[2,]) + sum(tab[3,])
            p_B <- 0.5 * sum(tab[,2]) + sum(tab[,3])
            x_AB <- tab[3,3] + tab[2,3]/2 + tab[3,2] / 2 + tab[2,2] / 4
            x_aB <- p_B - x_AB
            x_Ab <- p_A - x_AB
            x_ab <- 1 - x_AB - x_Ab - x_aB
            D <- (x_AB * x_ab - x_Ab * x_aB)
            p_AB <- x_AB
            ((D)^2 / (squeeze(p_A) * squeeze(1 - p_A) * squeeze(p_B) * squeeze(1-p_B))) * p_A * p_B             
        })))
        attr(r,"n") <- sum(iuse)
        r
    }
    popl <- split(seq_len(nrow(Xi)), pop)
    v0 <- lapply(popl, MakeOne)
    n <- sapply(v0, attr, which="n")
    v <- unlist(v0)
    sum(v[is.finite(v)] * n[is.finite(v)] / sum(n[is.finite(v)]))
}


ld_metric_r2hap <- function(Xi, Xj, pop){  
    MakeOne <- function(popi){
        iuse <- rowSums(Xi) == 2 & rowSums(Xj) == 2 & seq_len(nrow(Xi)) %in% popi
        r <- sum(outer(seq_len(ncol(Xi)),seq_len(ncol(Xj)), Vectorize(function(a,b){
            ## Only for diploid at the moment
            tab <- table(factor(Xi[iuse,a],0:2),factor(Xj[iuse,b],0:2))
            tab <- tab / sum(tab)
            p_A <- 0.5 * sum(tab[2,]) + sum(tab[3,])
            p_B <- 0.5 * sum(tab[,2]) + sum(tab[,3])
            x_AB <- tab[3,3] + tab[2,3]/2 + tab[3,2] / 2 + tab[2,2] / 4
            x_aB <- p_B - x_AB
            x_Ab <- p_A - x_AB
            x_ab <- 1 - x_AB - x_Ab - x_aB
            D <- (x_AB * x_ab - x_Ab * x_aB)
            p_AB <- x_AB
            (D^2 / (squeeze(p_A) * squeeze(1 - p_A) * squeeze(p_B) * squeeze(1-p_B))) * p_AB
        })))
        attr(r,"n") <- sum(iuse)
        r
    }
    popl <- split(seq_len(nrow(Xi)), pop)
    v0 <- lapply(popl, MakeOne)
    n <- sapply(v0, attr, which="n")
    v <- unlist(v0)
    sum(v[is.finite(v)] * n[is.finite(v)] / sum(n[is.finite(v)]))
}


ld_metric_D <- function(Xi, Xj, pop){  
    MakeOne <- function(popi){
        iuse <- rowSums(Xi) == 2 & rowSums(Xj) == 2 & seq_len(nrow(Xi)) %in% popi
        r <- sum(outer(seq_len(ncol(Xi)),seq_len(ncol(Xj)), Vectorize(function(a,b){
            ## Only for diploid at the moment
            tab <- table(factor(Xi[iuse,a],0:2),factor(Xj[iuse,b],0:2))
            tab <- tab / sum(tab)
            p_A <- 0.5 * sum(tab[2,]) + sum(tab[3,])
            p_B <- 0.5 * sum(tab[,2]) + sum(tab[,3])
            x_AB <- tab[3,3] + tab[2,3]/2 + tab[3,2] / 2 + tab[2,2] / 4
            x_aB <- p_B - x_AB
            x_Ab <- p_A - x_AB
            x_ab <- 1 - x_AB - x_Ab - x_aB
            D <- (x_AB * x_ab - x_Ab * x_aB)
            p_AB <- x_AB
            if(D >= 0){
                Dmax <- pmin(squeeze(p_A)*squeeze(1-p_B),squeeze(1-p_A)*squeeze(p_B))
            }else{
                Dmax <- pmin(squeeze(p_A)*squeeze(p_B),squeeze(1-p_A)*squeeze(1-p_B))
            }
            abs(D / Dmax) * p_A * p_B
        })))
        attr(r,"n") <- sum(iuse)
        r
    }
    popl <- split(seq_len(nrow(Xi)), pop)
    v0 <- lapply(popl, MakeOne)
    n <- sapply(v0, attr, which="n")
    v <- unlist(v0)
    sum(v[is.finite(v)] * n[is.finite(v)] / sum(n[is.finite(v)]))
}


ld_metric_Dhap <- function(Xi, Xj, pop){  
    MakeOne <- function(popi){
        iuse <- rowSums(Xi) == 2 & rowSums(Xj) == 2 & seq_len(nrow(Xi)) %in% popi
        r <- sum(outer(seq_len(ncol(Xi)),seq_len(ncol(Xj)), Vectorize(function(a,b){
            ## Only for diploid at the moment
            tab <- table(factor(Xi[iuse,a],0:2),factor(Xj[iuse,b],0:2))
            tab <- tab / sum(tab)
            p_A <- 0.5 * sum(tab[2,]) + sum(tab[3,])
            p_B <- 0.5 * sum(tab[,2]) + sum(tab[,3])
            x_AB <- tab[3,3] + tab[2,3]/2 + tab[3,2] / 2 + tab[2,2] / 4
            x_aB <- p_B - x_AB
            x_Ab <- p_A - x_AB
            x_ab <- 1 - x_AB - x_Ab - x_aB
            D <- (x_AB * x_ab - x_Ab * x_aB)
            p_AB <- x_AB
            if(D >= 0){
                Dmax <- pmin(squeeze(p_A)*squeeze(1-p_B),squeeze(1-p_A)*squeeze(p_B))
            }else{
                Dmax <- pmin(squeeze(p_A)*squeeze(p_B),squeeze(1-p_A)*squeeze(1-p_B))
            }
            abs(D / Dmax) * p_AB
        })))
        attr(r,"n") <- sum(iuse)
        r
    }
    popl <- split(seq_len(nrow(Xi)), pop)
    v0 <- lapply(popl, MakeOne)
    n <- sapply(v0, attr, which="n")
    v <- unlist(v0)
    sum(v[is.finite(v)] * n[is.finite(v)] / sum(n[is.finite(v)]))
}


ld_metric_Chisq <- function(Xi, Xj, pop){  
    MakeOne <- function(popi){
        iuse <- rowSums(Xi) == 2 & rowSums(Xj) == 2 & seq_len(nrow(Xi)) %in% popi
        r <- sum(outer(seq_len(ncol(Xi)),seq_len(ncol(Xj)), Vectorize(function(a,b){
            ## Only for diploid at the moment
            tab <- table(factor(Xi[iuse,a],0:2),factor(Xj[iuse,b],0:2))
            tab <- tab / sum(tab)
            p_A <- 0.5 * sum(tab[2,]) + sum(tab[3,])
            p_B <- 0.5 * sum(tab[,2]) + sum(tab[,3])
            x_AB <- tab[3,3] + tab[2,3]/2 + tab[3,2] / 2 + tab[2,2] / 4
            x_aB <- p_B - x_AB
            x_Ab <- p_A - x_AB
            x_ab <- 1 - x_AB - x_Ab - x_aB
            D <- (x_AB * x_ab - x_Ab * x_aB)
            p_AB <- x_AB
            2 * sum(iuse) * D^2 / (squeeze(p_A) * squeeze(p_B)) 
        })))
        attr(r,"n") <- sum(iuse)
        r
    }
    popl <- split(seq_len(nrow(Xi)), pop)
    v0 <- lapply(popl, MakeOne)
    n <- sapply(v0, attr, which="n")
    v <- unlist(v0)
    sum(v[is.finite(v)] * n[is.finite(v)] / sum(n[is.finite(v)]))
}

ld_metric_Chisqdf <- function(Xi, Xj, pop){  
    MakeOne <- function(popi){
        iuse <- rowSums(Xi) == 2 & rowSums(Xj) == 2 & seq_len(nrow(Xi)) %in% popi
        r <- sum(outer(seq_len(ncol(Xi)),seq_len(ncol(Xj)), Vectorize(function(a,b){
            ## Only for diploid at the moment
            tab <- table(factor(Xi[iuse,a],0:2),factor(Xj[iuse,b],0:2))
            tab <- tab / sum(tab)
            p_A <- 0.5 * sum(tab[2,]) + sum(tab[3,])
            p_B <- 0.5 * sum(tab[,2]) + sum(tab[,3])
            x_AB <- tab[3,3] + tab[2,3]/2 + tab[3,2] / 2 + tab[2,2] / 4
            x_aB <- p_B - x_AB
            x_Ab <- p_A - x_AB
            x_ab <- 1 - x_AB - x_Ab - x_aB
            D <- (x_AB * x_ab - x_Ab * x_aB)
            p_AB <- x_AB
            2 * sum(iuse) * D^2 / (squeeze(p_A) * squeeze(p_B))
        }))) / (2 * sum(iuse) * (col(Xi)-1) * (col(Xj)-1))
        attr(r,"n") <- sum(iuse)
        r
    }
    popl <- split(seq_len(nrow(Xi)), pop)
    v0 <- lapply(popl, MakeOne)
    n <- sapply(v0, attr, which="n")
    v <- unlist(v0)
    sum(v[is.finite(v)] * n[is.finite(v)] / sum(n[is.finite(v)]))
}


##' @method print ld_matrix
##' @export
print.ld_matrix <- function(x, ...){
    cat("Linkage disequlibrium measured by","\n")
    x2 <- x
    x2[upper.tri(x2)] <- NA
    diag(x2) <- NA
    tab <- round(as.matrix(x2),2)
    xx <- format(unclass(x2), digits = 2, justify = "right",scientific=FALSE)
    nn <- max(nchar(xx))
    print.table(tab, digits = 2) ##  zero.print = paste(rep("-", nn)
}
