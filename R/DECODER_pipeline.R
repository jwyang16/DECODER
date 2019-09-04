#' @import mgcv
#' @import lsei
#' @import preprocessCore
#' @import GenomicAlignments
#' @import GenomicRanges
#' @import umap
#' @importFrom mclust Mclust priorControl

#' @title DECODER: Deconvolution of DNA accessibility based on regression
#' @description This function is used to deconvolve bulk ATAC-seq samples.
#' @param bulk.dataset bulk ATAC-seq raw count matrix
#' @param ref the type of reference: user/PBMC/mouse pfc/hematopoietic
#' @param singlecell.dataset scATAC-seq raw count matrix. A required input when ref==user; equal to NULL when using built-in reference.
#' @param celltype.idx cell type labels of scATAC-seq. DECODER will cluster the scATAC-seq data when equal to NULL in 'user' mode; equal to NULL when using built-in reference.
#' @return 
#' \item{proportion.n}{estimated proportion matrix for bulk samples}
#' @keywords Deconvolution
#' @examples 
#' \dontrun{
#' prop.n<-DECODER(bulk.dataset,ref='user',singlecell.dataset,celltype.idx)
#' }
#' @export

DECODER<-function(bulk.dataset, ref=c('user','hematopoietic','pbmc','mouse_pfc'), singlecell.dataset=NULL,celltype.idx=NULL){
    if(ref=='user'){
        celltype.ref<-DECODER.user(singlecell.dataset,celltype.idx)
    }
    else if(ref=='hematopoietic'){
        #celltype.ref<-readRDS('hematopoietic.rds')
        utils::data(reference)
        celltype.ref<-reference$hematopoietic
    }
    else if(ref=='pbmc'){
        #celltype.ref<-readRDS('pbmc.rds')
        utils::data(reference)
        celltype.ref<-reference$pbmc
    }
    else if(ref=='mouse_pfc'){
        #celltype.ref<-readRDS('mouse_pfc.rds')
        utils::data(reference)
        celltype.ref<-reference$mouse_pfc
    }
    else{
        stop(sprintf('Reference for \"%s\" ',ref),'is not found in the built-in reference profiles.')
    }
    
    anchor<-select.anchor(celltype.ref,ths=0.95)
    anchor.idx<-anchor$anchor.idx
    anchor.locus<-anchor$anchor.locus
    variable.idx<-select.variable(celltype.ref,ths=0.05)
    
    bulk.dataset<-as.matrix(bulk.dataset)
    num.bulk<-ncol(bulk.dataset)
    num.celltype<-ncol(celltype.ref)
    proportion.n<-matrix(rep(0,num.celltype*num.bulk), num.celltype)
    rownames(proportion.n)<-colnames(celltype.ref)
    colnames(proportion.n)<-colnames(bulk.dataset)

    for(i in 1:num.bulk){
        bulk<-bulk.dataset[,i]
        bulk.n<-bulk.normalize(bulk,variable.idx,anchor.idx,anchor.locus)
        
        celltype.ref.n<-celltype.ref[variable.idx,] 
        celltype.ref.n<-normalize.quantiles.use.target(celltype.ref.n,target=as.vector(bulk.n))
        
        proportion.n[,i]<-pnnls(log2(1+celltype.ref.n),as.matrix(log2(1+bulk.n)))$x
        proportion.n[which(proportion.n[,i]<=0),i]<-0
        proportion.n[,i]<-proportion.n[,i]/sum(proportion.n[,i])
    }
    
    return(proportion.n)
}



