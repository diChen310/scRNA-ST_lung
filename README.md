# scRNA-ST_lung
scRNA_workflow.R:  The main analysis of scRNA-seq including the QC, integration, clustering, cell type definition as well as identification of the differences among multiple lesions of the same patients. 

stRNA_workflow.R: The analyses based on spatial transcriptomics (ST) or integration of both scRNA-seq and ST including the QC, integration, clustering, cell type prediction based on the scRNA-seq data, and pathological sub-type analysis. 

epiSubClustering.R: subtype analysis for the epithelial cells.

tCellSubClustering.R:subtype analysis for the T&NK cells.

BCellSubClustering.R:subtype analysis for the B cells.

endoCellSubClustering.R: subtype analysis for the endothelial cells.

myeloidCellSubClustering.R: subtype analysis for the myeloid cells.

inferCNV.R & inferCNVDownStreamAnalysis.R: cnv based analysis.

scDC_cellTypeCompositionDiffAnalysis: identification of cell composition differences between different samples.
