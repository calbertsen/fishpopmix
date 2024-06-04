##' @export
read_ped <- function(file, map = gsub("\\.ped$",".map",file),
                     nIndividuals = 1000,nLoci = 10000,
                     missing_genotype = "0",
                     no_fid = FALSE,
                     no_parents = FALSE,
                     no_sex = FALSE,
                     no_pheno = FALSE,
                     progress = TRUE){
    
    if(length(file) != 1 || !is.character(file))
        stop("file should be a length one character vector")
    if(length(map) != 1 || !is.character(map))
        stop("map should be a length one character vector")
    if(length(missing_genotype) != 1 || !is.character(missing_genotype))
        stop("missing_genotype should be a length one character vector")
    if(length(no_fid) != 1 || !is.logical(no_fid))
        stop("no_fid should be a length one logical vector")

    .Call("read_ped",file,map,nIndividuals,nLoci, !no_fid, !no_parents, !no_sex, !no_pheno,progress,missing_genotype, PACKAGE="fishpopmix")
    
}
