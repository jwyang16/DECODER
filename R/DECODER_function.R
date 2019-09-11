#' @import mgcv
#' @import lsei
#' @import preprocessCore
#' @import GenomicAlignments
#' @import GenomicRanges
#' @import umap
#' @import Rsamtools



#' @title Transfer raw count matrix to count-per-million (CPM)
#' @description This function is used to transfer the input of raw count matrix to count-per-million (CPM)
#' @param dataset raw count matrix
#' @param c the transfered sequencing depth constant. Default=1e6.
#' @return 
#' \item{dataset.t}{CPM matrix}
#' @keywords preprocessing
#' @examples 
#' \dontrun{
#' dataset.cpm<-count2cpm(dataset,c=1e6)
#' }
#' @export

count2cpm<-function(dataset,c=1e6){
    depth<-colSums(dataset)
    dataset.t<-sweep(dataset,2L,depth,'/')*c
    return(dataset.t)
}



#' @title Select constant anchor loci from single-cell ATAC-seq reference
#' @description This function is used to select constant anchor loci from scATAC-seq reference
#' @param celltype.ref scATAC-seq reference
#' @param ths threshold for choosing loci whose differences are lower than ths\% of total differences as constant anchors. Default=0.95
#' @return 
#' \item{anchor.idx}{the index of constant anchor loci}
#' \item{anchor.locus}{the profile of constant anchor loci}
#' @keywords gam-based anchor selection
#' @examples 
#' \dontrun{
#' anchor<-select.anchor(celltype.ref,ths=0.95)
#' }
#' @export

select.anchor<-function(celltype.ref,ths=0.95){
    celltype.mean<-apply(celltype.ref,1,mean)
    celltype.var<-apply(celltype.ref,1,stats::var)
    data_fit <- data.frame(X=log2(1+celltype.mean),Y=log2(1+celltype.var))
    fit_model <- gam(Y~s(X),data=data_fit)
    
    diff_var <- log2(1+celltype.var) - fit_model$fitted.values
    diff_var_sort <- sort(diff_var,decreasing=TRUE)
    diff_var_th <- diff_var_sort[round(length(diff_var_sort)*ths)]
    
    anchor.idx <- which(diff_var < diff_var_th)
    anchor.locus<-celltype.mean[anchor.idx]
    
    return(list(anchor.idx=anchor.idx, anchor.locus=anchor.locus))
}



#' @title Select highly variable loci from single-cell ATAC-seq reference
#' @description This function is used to select highly variable anchors from scATAC-seq reference
#' @param celltype.ref scATAC-seq reference
#' @param ths threshold for choosing loci whose differences are higher than ths\% of total differences as variable anchors. Default=0.05.
#' @return 
#' \item{variable.idx}{the index of highly variable loci}
#' @keywords gam-based anchor selection
#' @examples 
#' \dontrun{
#' variable.idx<-select.variable(celltype.ref,ths=0.05)
#' }
#' @export

select.variable<-function(celltype.ref,ths=0.05){
    celltype.mean<-apply(celltype.ref,1,mean)
    celltype.var<-apply(celltype.ref,1,stats::var)
    data_fit <- data.frame(X=log2(1+celltype.mean),Y=log2(1+celltype.var))
    fit_model <- gam(Y~s(X),data=data_fit)
    
    diff_var <- log2(1+celltype.var) - fit_model$fitted.values
    diff_var_sort <- sort(diff_var,decreasing=TRUE)
    diff_var_th <- diff_var_sort[round(length(diff_var_sort)*ths)]
    
    variable.idx <- which(diff_var > diff_var_th)
    return(variable.idx)
}



#' @title Normalize the bulk sample using anchors
#' @description This function is used to normalize the bulk sample using anchors.
#' @param bulk bulk sample to be deconvolved
#' @param variable.idx the index of highly variable loci
#' @param anchor.idx the index of constant anchor loci
#' @param anchor.locus the profile of constant anchor loci
#' @return 
#' \item{bulk.normalized}{the anchor-based normalized bulk sample}
#' @keywords normalization
#' @examples
#' \dontrun{
#' bulk.normalized<-bulk.normalize(bulk,variable.idx,anchor.idx,anchor.locus)
#' }
#' @export

bulk.normalize<-function(bulk,variable.idx,anchor.idx,anchor.locus){
    
    anchor.bulk<-bulk[anchor.idx]
    
    fv<-mspline(anchor.bulk,anchor.locus)
    bulk.normalized<-mpredict(fv,bulk[variable.idx])
    
    bulk.normalized[which(bulk.normalized<0)]<-0
    
    return(bulk.normalized)
}



#' @title Fit monotonous smooth spline
#' @description This function is used to fit monotonous smooth spline.
#' @param x the independent variables
#' @param y the dependent variables
#' @param k number of knots. Default=8
#' @param lower the lower bound of spline. Default=NA
#' @param upper the upper bound of spline. Default=NA
#' @return 
#' \item{sm}{smooth condition}
#' \item{p}{spline}
#' @keywords Smooth spline with monotonous restriction
#' @examples 
#' \dontrun{
#' fv<-mspline(anchor.bulk,anchor.locus,k=8)
#' }
#' @export

mspline<-function(x,y,k=8,lower=NA,upper=NA){
    dat<-data.frame(x=x,y=y)
    init_gam <- gam(y~s(x,k=k,bs="cr"))
    sm <- smoothCon(s(x,k=k,bs="cr"),dat,knots=NULL)[[1]]
    mc <- mono.con(sm$xp,lower=lower,upper=upper) 
    M <- list(X=sm$X,y=y,C=matrix(0,0,0),Ain=mc$A, bin=mc$b,sp=init_gam$sp, p=sm$xp,S=sm$S,w=y*0+1,off=0)
    p<-pcls(M)
    return(list(sm=sm,p=p))
}


#' @title Predict using the monotonous smooth spline
#' @description This function is used to predict values using the spline.
#' @param msp the monotonous smooth spline
#' @param x the new independent variables
#' @return 
#' \item{res}{the predicted value}
#' @keywords Prediction
#' @examples 
#' \dontrun{
#' bulk.normalized<-mpredict(fv,bulk)
#' }
#' @export

mpredict<-function(msp,x){
    res<-Predict.matrix(msp$sm,data.frame(x=x))%*%msp$p
    return(res)
}



#' @title Evaluate the result 
#' @description  This function is used to evaluate the results.
#' @param prop1 the first proportion 
#' @param prop2 the second proportion
#' @return 
#' \item{RMSD}{Root-mean-square-deviation between two proportions}
#' \item{pearson}{Pearson correlation between two proportions}
#' @keywords Evaluation
#' @examples 
#' \dontrun{
#' res<-evaluate(prop.real,prop.n)
#' }
#' @export

evaluate<-function(prop1,prop2){
    RMSD<-rep(0,ncol(prop1))
    pearson<-rep(0,ncol(prop1))

    for(i in 1:ncol(prop1)){
        a<-prop1[,i]-prop2[,i]
        RMSD[i]<-sqrt(mean(a^2))
        pearson[i]<-stats::cor(prop1[,i],prop2[,i],method="pearson")
    }
    
    return(list(RMSD=RMSD, pearson=pearson))
}



#' @title Construct cell type-specific reference when ref==user
#' @description This function is used to construct scATAC-seq reference when ref==user.
#' @param singlecell.dataset scATAC-seq raw count matrix.
#' @param celltype.idx Cell type labels of scATAC-seq. DECODER will cluster the scATAC-seq data when equal to NULL.
#' @return 
#' \item{res}{if celltype.idx==NULL, res is a list containing reference of optimal cluster number k and reference of each k from 2 to 50; otherwise res is a reference matrix.}
#' @keywords Reference construction
#' @examples
#' \dontrun{
#' res<-DECODER.user(singlecell.dataset,celltype.idx)
#' }
#' @export

DECODER.user<-function(singlecell.dataset,celltype.idx){
    if(is.null(singlecell.dataset)){
        stop('Single-cell raw count matrix is missing.')
    }
    else if(is.null(celltype.idx)){
        base::message("Cell type labels are missing.")
        base::message("DECODER will cluster the cells and return a list of proportions")
        allref<-list()
        
        y<-singlecell.dataset[rowMeans(singlecell.dataset>0)>0.1,]
        y<-normalize.quantiles(y)
        prres <- stats::prcomp(t(y),scale. = T,rank.=50)$x
        set.seed(1234)
        tsne <- umap(prres)$layout
        
        ds<-stats::dist(tsne,method='euclidean')
        hc<-stats::hclust(ds,method='complete')
        for(k in 2:min(50,ncol(singlecell.dataset))){
            celltype.idx<-stats::cutree(hc,k=k)
            celltype.unique<-unique(celltype.idx)
            len<-length(unique(celltype.idx))
            y<-matrix(data=NA,nrow=nrow(singlecell.dataset),ncol=len)
            for(i in 1:len){
                idx<-which(celltype.idx==celltype.unique[i])
                y[,i]<-rowMeans(as.matrix(singlecell.dataset[,idx]))
            }
            celltype.ref<-count2cpm(y)
            tmp<-paste0('k=',k)
            allref[[tmp]]<-celltype.ref
        }
        
        xaxis <- 2:min(50,ncol(singlecell.dataset))
        yaxis <- sapply(xaxis, function(k) {
            clu <- stats::cutree(hc,k)
            cluSS <- sum(sapply(unique(clu),function(i) {
                sum(rowSums(sweep(tsne[clu==i,,drop=F],2,colMeans(tsne[clu==i,,drop=F]),"-")^2))
            }))
            1-cluSS/sum((sweep(tsne,2,colMeans(tsne),"-"))^2)
        })      
        yaxis <- c(0,yaxis)
        xaxis <- c(1,xaxis)
        clunum <- which.min(sapply(xaxis, function(i) {
            x2 <- pmax(0,xaxis-i)
            sum(stats::lm(yaxis~xaxis+x2)$residuals^2)
        }))
        
        opt<-paste0('Optimal k=',clunum)
        res<-list()
        res[[opt]]<-allref[[paste0('k=',clunum)]]
        res[['All']]<-allref
        
    }
    else{
        celltype.unique<-unique(celltype.idx)
        len<-length(unique(celltype.idx))
        y<-matrix(data=NA,nrow=nrow(singlecell.dataset),ncol=len)
        for(i in 1:len){
            idx<-which(celltype.idx==celltype.unique[i])
            y[,i]<-rowMeans(as.matrix(singlecell.dataset[,idx]))
        }
        res<-count2cpm(y)
    }
    return(res)
}


