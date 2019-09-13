#' @import GenomicAlignments
#' @import GenomicRanges
#' @import Rsamtools


#' @title do peak calling
#' @description This fuction is used to do peak calling.
#' @param bam_file_path path to single-cell bam files (The number of sc bam files must >= 2.)
#' @param genome genome assembly to be used
#' @return
#' \item{peak}{peak regions from merged single cells}
#' @keywords peak calling
#' @examples 
#' \dontrun{
#' peak<-peakcalling('path/to/bam/files',genome='hg19')
#' }
#' @export

peakcalling<-function(bam_file_path,genome=c('hg19','hg38','mm9','mm10')){
    f <- list.files(bam_file_path,pattern="\\.bam$",full.names = T)
    mergeBam(files=f,destination = paste0(bam_file_path,'/','merged_bam.bam'))

    utils::data(genomes)
    if(genome=='hg19'){
        assembly<-genomes$hg19
    }
    else if(genome=='hg38'){
        assembly<-genomes$hg38
    }
    else if(genome=='mm9'){
        assembly<-genomes$mm9
    }
    else if(genome=='mm10'){
        assembly<-genomes$mm10
    }
    else{
        stop(sprintf('Genome assembly \"%s\" ',genome),'is not found.') 
    }
    
    bk<-assembly
    GenomicRanges::start(bk)<-GenomicRanges::end(bk)<-round((GenomicRanges::start(bk) + GenomicRanges::end(bk))/2)
    GenomicRanges::start(bk)<-GenomicRanges::start(bk)-4999
    GenomicRanges::end(bk)<-GenomicRanges::end(bk)+5000
    
    satac<-GenomicAlignments::readGAlignmentPairs(paste0(bam_file_path,'/','merged_bam.bam'))
    satac<-GenomicRanges::GRanges(satac)
    GenomicRanges::start(satac)<-GenomicRanges::end(satac)<-round((GenomicRanges::start(satac) + GenomicRanges::end(satac))/2)
    count.templt<-GenomicRanges::countOverlaps(assembly,satac,ignore.strand=T)
    count.bk<-GenomicRanges::countOverlaps(bk,satac,ignore.strand=T)
    
    flag<-count.templt/count.bk
    flag.0<-flag
    flag.0[is.nan(flag.0)]<-0
    
    flag.n<-flag[!is.nan(flag)]
    set.seed(2019)
    km<-stats::kmeans(flag.n,centers=2,iter.max=100)
    rsd<-stats::sd(flag.0[which(km$cluster==which.max(km$centers))])
    
    cutoff<-max(km$centers)+rsd
    idx<-which(flag.0>cutoff)
    peak<-assembly[idx]
    return(peak)
}



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
    
    start(gr)<-start(gr)-99
    end(gr)<-end(gr)+99
    
    return(gr)
}



#' @title get count matrix from peak file and bam files
#' @description This function is used to get count matrix from bam files and peak file.
#' @param peak Granges format of peak regions 
#' @param bam_file_path path to bam files
#' @return 
#' \item{countmatrix}{count matrix, peak by sample}
#' @keywords count matrix
#' @examples 
#' \dontrun{
#' countmat<-getcount(bed2granges('path/to/bed/file/peak.bed'),'path/to/bam/files')
#' }
#' @export

getcount<-function(peak,bam_file_path){

    f <- list.files(bam_file_path,pattern="\\.bam$",full.names = T)
    satac <- sapply(sapply(f,readGAlignmentPairs),GRanges)
    n <- names(satac)
    satac <- lapply(satac,function(i) {
        GenomicRanges::start(i) <- GenomicRanges::end(i) <- round((GenomicRanges::start(i) + GenomicRanges::end(i))/2)
        i
    })
    names(satac) <- n
    
    countmatrix <- sapply(names(satac),function(sid) {
        tp<-GenomicRanges::countOverlaps(peak,satac[[sid]],ignore.strand=T)
        tp
    })
    return(countmatrix)
}




