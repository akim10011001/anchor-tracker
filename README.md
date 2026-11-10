# anchor-tracker

To check whether the projection workflow is actually credible: not biased by a certain batch variable, and the projection on ref.umap and the prediction results match each other. Also for better visualizations.

Built on Seurat transfer anchors (FindTransferAnchors + MapQuery, reference PCA projection).

## Requirements

R with Seurat (>= 5), SeuratObject, dplyr, ggplot2, tidyr, tibble, Matrix.
Optional: MASS, DESeq2, ggalluvial, SummarizedExperiment, glmGamPoi, future.

## Example (Usage)

```r
source("R/anchor_projection_functions.R")

ref <- prep_reference(reference, assay = "RNA")
qry <- prep_query(query, reference = ref)
res <- run_anchors(qry, ref)

anchors <- update_anchors(res$anchorset, res$query, res$reference,
                          variables           = "CellClass",
                          query_meta_cols     = c("seurat_clusters", "protocol"),
                          reference_meta_cols = "CellClass")
anchors <- transfer_reduction(anchors, query = res$query, reference = res$reference)
anchors <- max_anchors(anchors)
```

`prep_reference()`, `prep_query()` and `run_anchors()` reuse what is already on the objects, never overwrite existing reductions or clusters, and return a `$report` of what was reused or built.

## Checks

Batch bias
```r
max_anchor_bar_plot(anchors, group.by = "protocol")
scatter_plot(anchors, group.by = "protocol", plot.by = "CellClass")
stacked_count_plot(anchors, group.by = "protocol", plot.by = "CellClass")
```

ref.umap position vs. prediction
```r
prediction_plot(anchors, group.by = "predicted.CellClass", reduction = "ref.umap")
density_plot(anchors)
max_anchor_plot(anchors, ref.reduction = ref$report$names$umap)
matrix_plot(anchors, ref_col = "CellClass")
prop_plot(anchors, ref_col = "CellClass")
```
