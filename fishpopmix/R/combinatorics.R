## get length k subsets of length n vector
nexksb <- function(n,k){
    ## Nijenhuis & Wild (1978) DOI: 10.1016/C2013-0-11243-3
    i <- 1
    r <- list()
    m <- 0; h <- k
    a <- numeric(k)
    repeat{
        ##cat(i,"/",ncombi,"\n")
        for(j in 1:h){
            a[k+j-h] <- m+j            
        }
        r[[length(r)+1]] <- a
        if(a[1] == n-k+1){
            break;
        }
        if(m < n-h){
           h <- 0
        }
        h <- h + 1
        m <- a[k+1-h]
        i <- i + 1
        if(i > ncombi+10)
            break;
    }
    r
}

## Get length k components that sum to n
nexcom <- function(n, k){
    ## Nijenhuis & Wild (1978) DOI: 10.1016/C2013-0-11243-3
    rlist <- vector("list",choose(n+k-k,k))
    r <- numeric(k)
    r[1] <- n
    t <- n
    h <- 0
    i <- 1
    repeat{
        rlist[[i]] <- r
        if(r[k] == n){
            break;
        }
        if(t > 1){
           h <- 0
        }
        h <- h + 1
        t <- r[h]
        r[h] <- 0
        r[1] <- t-1
        r[h+1] <- r[h+1]+1
        i <- i+1
    }
    rlist
}

