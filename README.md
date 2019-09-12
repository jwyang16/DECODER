# DECODER: DECOnvolution of DNA accEssibility using anchor-based Regression

DECODER is an anchor-based deconvolution method designed for ATAC-seq and other chromatin accessibility data. DECODER deconvolves bulk ATAC-seq samples using cell type-specific reference generated from scATAC-seq data. 


## DECODER Installation

DECODER can be installed via Github.
To install the latest version of DECODER package via Github, run following commands in R:

	if (!require("devtools"))
	   install.packages("devtools")
	devtools::install_github("jwyang16/DECODER")


## Getting Started

Users can load DECODER as follows:
	
	library(DECODER)

Users can choose either to use the built-in cell type reference or to construct reference by their own scATAC-seq data.

	data(demo)
	bulk.sample <- demo$bulk.sample
	singlecell.dataset <- demo$singlecell.dataset
	celltype.idx <- demo$celltype.idx
	
	prop <- DECODER(bulk.sample, ref='hematopoietic', singlecell.dataset=NULL, celltype.idx=NULL)
	prop <- DECODER(bulk.sample, ref='user', singlecell.dataset=singlecell.dataset, celltype.idx=celltype.idx)
	prop <- DECODER(bulk.sample, ref='user', singlecell.dataset=singlecell.dataset, celltype.idx=NULL)

Each mode returns the estimated cell proportions of cell types in the chosen reference.

Users can run the DECODER pipeline of reading in bam, peak calling and counting, and performing DECODER.
	
	prop <- DECODER.pipeline(bulk_bam_file_path, ref='user', sc_bam_file_path, celltype.idx, genome='hg19')
	prop <- DECODER.pipeline(bulk_bam_file_path, ref='hematopoietic', sc_bam_file_path=NULL, celltype.idx=NULL, genome=NULL)
