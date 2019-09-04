#' @import GenomicAlignments
#' @import GenomicRanges



#' @title transfer bed file to Granges
#' @description This function is used to transfer bed to Granges
#' @param bed_file bed file
#' @return 
#' \item{gr}{Granges object of bed file}
#' @keywords data structure
#' @examples 
#' \dontrun{
#' gr<-bed2granges('path/to/bed/file/peak.bed')
#' }
#' @export

bed2granges <- function(bed_file){
    df <- utils::read.table(bed_file,
                     header=F,
                     stringsAsFactors=F)
    if(length(df) > 6){
        df <- df[,-c(7:length(df))]
    }
    header <- c('chr','start','end','id','score','strand')
    names(df) <- header[1:length(names(df))]
    if('strand' %in% colnames(df)){
        df$strand <- gsub(pattern="[^+-]+", replacement = '*', x = df$strand)
    }
    if(length(df)==3){
        gr <- with(df, GRanges(chr, IRanges(start, end)))
    } else if (length(df)==4){
        gr <- with(df, GRanges(chr, IRanges(start, end), id=id))
    } else if (length(df)==5){
        gr <- with(df, GRanges(chr, IRanges(start, end), id=id, score=score))
    } else if (length(df)==6){
        gr <- with(df, GRanges(chr, IRanges(start, end), id=id, score=score, strand=strand))
    }
    return(gr)
}



#' @title get count matrix from peak file and bam files
#' @description This function is used to get count matrix from bam files and peak file.
#' @param bed_file peak file
#' @param bam_file_path path to bam files
#' @return 
#' \item{countmatrix}{count matrix, peak by sample}
#' @keywords count matrix
#' @examples 
#' \dontrun{
#' countmat<-getcount('path/to/bed/file/peak.bed','path/to/bam/files')
#' }
#' @export

getcount<-function(bed_file,bam_file_path){
    gr<-bed2granges(bed_file)
    start(gr)<-start(gr)-99
    end(gr)<-end(gr)+99
    
    f <- list.files(bam_file_path,pattern="\\.bam$")
    satac <- sapply(sapply(f,readGAlignmentPairs),GRanges)
    n <- names(satac)
    satac <- lapply(satac,function(i) {
        start(i) <- end(i) <- round((start(i) + end(i))/2)
        i
    })
    names(satac) <- n
    
    countmatrix <- sapply(names(satac),function(sid) {
        tp<-countOverlaps(gr,satac[[sid]],ignore.strand=T)
        tp
    })
    return(countmatrix)
}




