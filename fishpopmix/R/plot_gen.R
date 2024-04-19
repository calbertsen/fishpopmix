##' @method plot gen
##' @export
plot.gen <- function(x, ...){
    indx <- split(seq_len(dim(x)[3]),attr(x,"population"))
    x2 <- apply(apply(x, c(1,2), sum),2,which.max)
    x3 <- ifelse(x2==1,2,1)
    xClass <- matrix("white",dim(x)[3], dim(x)[2])
    for(i in seq_len(dim(x)[2])){
        v <- apply(x[,i,],2:3,function(y){
            if(any(is.na(y))) return(palette()[4])
            if(y[1]==2 && y[2]==0) return(palette()[1])
            if(y[1]==0 && y[2]==2) return(palette()[3])
            if(y[1]==1 && y[2]==1) return(palette()[2])
            return(palette()[4])
        })
        xClass[,i] <- v
    }

    d <- split(as.data.frame(xClass),attr(x,"population"))
    d2 <- vector("list",length(d)+length(d)-1)
    d2[seq_along(d2) %% 2 == 1] <- d
    d2[seq_along(d2) %% 2 == 0] <- replicate(length(d)-1,matrix("white",round(dim(x)[3]/100),dim(x)[2]), simplify = FALSE)

    xClass2 <- do.call("rbind",d2)

    vx <- sapply(d2,nrow)
    px <- c(0,cumsum(vx/sum(vx)))

    plot.new()
    rasterImage(as.raster(as.matrix(xClass2)), par("usr")[1] + strheight("C",cex=2), par("usr")[3]+ strheight("C",cex=4), par("usr")[2], par("usr")[4],
                interpolate = FALSE)
    for(s in seq_along(d))
        text(par("usr")[1] + strheight("C",cex=2),
             par("usr")[4] - mean(px[(s-1)*2+1:2]) * (diff(par("usr")[3:4]) - strheight("C",cex=4)),
             names(d)[s],pos=2,srt=0,cex=1,xpd=TRUE)
    legend("bottom",legend=c("Homozygot Reference", "Heterozygot", "Homozygot Alternative","Missing"),
           fill = 1:4,ncol=4,border=NA,box.col=NA)

}
