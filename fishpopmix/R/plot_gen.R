##' @method plot gen
##' @export
plot.gen <- function(x, col = colors(4), ...){
    ## indx <- split(seq_along(x),attr(x,"population"))
    nAllele <- num_allele(x)
    if(any(nAllele != 2))
        stop("The plot is currenly only implemented for biallelic genotypes")
    ploidy <- get_ploidy(x)
    if(length(ploidy) != 1 && names(ploidy)[1] != 2)
        stop("The plot is currently only implemented for diploid genotypes")
    if(length(col) != 4)
        stop("col should have length 4")
    nLoci <- num_loci(x)
    nIndi <- length(x)
    xClass <- matrix("white",nIndi, nLoci)
    for(l in seq_len(nLoci)){
        gt <- extract_locus(x,l)
        gs <- rowSums(gt)
        g <- gt[,1]
        xClass[gs < 2,l] <- col[4]
        xClass[gs == 2 & g == 2,l] <- col[1]
        xClass[gs == 2 & g == 1,l] <- col[2]
        xClass[gs == 2 & g == 0,l] <- col[3]        
    }
    d <- split(as.data.frame(xClass),attr(x,"population"))
    d2 <- vector("list",length(d)+length(d)-1)
    d2[seq_along(d2) %% 2 == 1] <- d
    d2[seq_along(d2) %% 2 == 0] <- replicate(length(d)-1,matrix("white",round(nIndi/100),nLoci), simplify = FALSE)

    xClass2 <- do.call("rbind",d2)

    vx <- sapply(d2,nrow)
    px <- c(0,cumsum(vx/sum(vx)))

    plot.new()
    strw <- max(strwidth(names(d),cex=1))
    rasterImage(as.raster(as.matrix(xClass2)), par("usr")[1] + strw , par("usr")[3]+ strheight("C",cex=4), par("usr")[2], par("usr")[4],
                interpolate = FALSE)
    for(s in seq_along(d))
        text(par("usr")[1] + strw ,
             par("usr")[4] - mean(px[(s-1)*2+1:2]) * (diff(par("usr")[3:4]) - strheight("C",cex=4)),
             names(d)[s],pos=2,srt=0,cex=1,xpd=TRUE)
    legend("bottom",legend=c("Homozygot Reference", "Heterozygot", "Homozygot Alternative","Missing"),
           fill = col,ncol=4,border=NA,box.col=NA)

}
