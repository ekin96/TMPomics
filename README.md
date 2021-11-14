# TMPomics
#R function that performs an automatic statistical analysis of metabolomic, proteomic and trascrittomic data matrixes:

- data = dataset name
- names_col = how columns are named in the dataset
- groups = n° of groups considered in the analysis (e.g. young and old)
- p_threshold = p-value threshold
- n4group = n° of individuals in each group
- PCA_dimensions = n° of dimensions considered in PCA analysis
- fold_threshold = Log2 fold change threshold value
- p_col = color assigned to transcripts/protein/metabolites below the threshold p-value in volcano plot 
- fold_col = color assigned to transcripts/protein/metabolites exceeding the threshold |log2 fold change| in volcano plot
- both_col = color assigned to transcripts/protein/metabolites exceeding the threshold |log2 fold change| and below the threshold p-value in volcano plot
- cex_volcanolab = size of labels in volcano plot
- x_lim_volcano = x dimension in volcano plot
- y_lim_volcano = y dimension in volcano plot
- pch_volcano = size of points in volcano plot 
