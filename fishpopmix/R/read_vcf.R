
##' @export
read_vcf <- function(file, n = 10000, simplified = TRUE, pooled = FALSE){
    ## TODO:
    ## Pooled VCF -> AD gives the allele frequencies
    ## Non-pooled -> GT gives genotype
    ## Should extract allele names
    ## Convert to gen (pooled = FALSE) class or alleleFrequencyList (pooled=TRUE)
    f <- base::file(file)
    open(f)

    lFirst <- readLines(f,1)
    if(!grepl("^##fileformat=VCFv4\\.",lFirst))
        stop("The file does not appear to be a valid VCFv4 file")
    header <- NA
    tabDat <- list()
    nobs <- 0
    iter <- 0
    repeat {
        iter <- iter + 1
        cat("Iteration",iter,"Lines",1+n*(iter-1),"-",n*iter,"Processed data:",nobs,"\n")
        lBlock <- readLines(f, n)
        if(length(lBlock) == 0)
            break;
        ## Remove meta data
        ik <- !grepl("^##",lBlock)
        if(sum(ik) == 0)
            next;
        if(any(grepl("^#CHROM",lBlock))){
            header <- lBlock[grepl("^#CHROM",lBlock)]
            ik <- !grepl("^#",lBlock)
            if(sum(ik) == 0)
                next;        
        }    
        parseOne <- function(l){
            x <- strsplit(l,"[[:space:]]+")[[1]]
            name <- gsub("_\\.$","",sprintf("chr%s_%s_%s",x[1],x[2],x[3]))
            i <- match("GT",strsplit(x[9],":")[[1]])
            ni <- length(x) - 9
            if(is.na(i)){
                return(c(name,x[7],rep("000000",ni)))
            }
            c(name, x[7], sapply(strsplit(x[-(1:9)],":"),function(y) paste(sprintf("%03d",as.numeric(strsplit(gsub("\\.","-1",y[i]),"[|/]")[[1]])+1),collapse="")))
        }
        newDat <- lapply(lBlock[ik], parseOne)
        nobs <- nobs + length(newDat)
        tabDat <- list(tabDat, newDat)
        if(length(lBlock) < n)
            break;
    }
    close(f)

    datMat <- matrix(unlist(tabDat),2+length(strsplit(header,"[[:space:]]+")[[1]])-9,byrow=TRUE)
    datMat
}
