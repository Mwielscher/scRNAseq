## About this Repository  
These are Jupyternotebooks for scRNAseq analysis. Comments and more detailed explainations are in the notebooks.  
>* this [notebook](singleSAMPLE_scRNAseq.ipynb) takes output from [10xcell ranger](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/output/overview) and performs an initial exploratory analysis of the reads incl. doublet score calculation. Then it runs the Seurat SCT pipeline. The notebook als contains scirpt cells for differential gene expresssion analysis. It compares each cell population (cluster) against each remaining cluster (cell population) in the analysis. Resulting lists can be exported or can be further analysed with enrichR to gather information about possible cell type identities. Finally, cluster specfic markers are TIERed so you can see the top markers associtated to each cluster.
>* this [notebook](integreate_multiple10x_samples.ipynb) merges and integrates several 10x scRNAseq datasests together and then performs the same differential gene expression analysis as outlined above.  
