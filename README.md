# cluster_seurat3_10x_ms
An updated version of the cluster_pooled_10x_ms pipeline

1. InitializeVars - Import data from 10X, Make combined object
2. Cluster - Scale, run PCA, run UMAP, Find Neighbors, Find Clusters, Visualize Batch Effects, Identify Cell Type Markers
3. ComparePopCounts - compare population sizes between factors of interest
4. CompareDE - compare differential expression between factors of interest

S. Script to demultiplex on the basis of sex, and identify putative doublets using co-expression of Xist and Ymarker genes
C. Script for comparing cluster stability when clusters made using different methods?
