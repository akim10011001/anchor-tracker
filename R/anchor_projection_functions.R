# =============================================================================
# Anchor projection functions (Seurat 5)
# -----------------------------------------------------------------------------
# Project a query single-cell dataset onto a reference with Seurat transfer
# anchors (reference PCA projection + reference UMAP model), keep a light
# "anchors" object (metadata, anchor pairs, embeddings; no expression data),
# and plot anchor / prediction results.
#
# Typical workflow
#   ref     <- prep_reference(reference, assay = "RNA")
#   qry     <- prep_query(query, reference = ref)
#   res     <- run_anchors(qry, ref)
#   anchors <- update_anchors(res$anchorset, res$query, res$reference,
#                             variables           = "CellClass",
#                             query_meta_cols     = "seurat_clusters",
#                             reference_meta_cols = "CellClass")
#   anchors <- transfer_reduction(anchors, query = res$query, reference = res$reference)
#   anchors <- max_anchors(anchors)
#   max_anchor_plot(anchors, ref.reduction = ref$report$names$umap)
#
# prep_reference(), prep_query() and run_anchors() return their object(s)
# plus $report (input state, what was reused / built / skipped, names
# written, genes, parameters used). Nothing is assigned into the caller's
# environment, and nothing already on an object is overwritten unless
# overwrite = TRUE (new results go to new names such as "pca_2").
#
# Parallel: prep_reference(), prep_query() and run_anchors() take workers = <n>
# (future multisession) for Seurat's parallel steps (ScaleData regression,
# SCTransform, UMAP threads, ...); the caller's plan is restored afterwards.
# workers = NULL keeps the caller's own future plan.
#
# Required packages: Seurat (>= 5.0), SeuratObject, dplyr, ggplot2, tidyr,
#                    tibble, Matrix
# Optional packages: MASS (density_plot), DESeq2 (compute_pseudobulk),
#                    ggalluvial (sandkey_plot_seurat), SummarizedExperiment
#                    (transfer_metadata on SingleCellExperiment objects),
#                    glmGamPoi (faster SCTransform)
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(tibble)
  library(Matrix)
})

`%||%` <- function(a, b) if (is.null(a)) b else a


# ── internal helpers ──────────────────────────────────────────────────────────

.need <- function(pkg, what) {
  if (!requireNamespace(pkg, quietly = TRUE))
    stop(sprintf("package '%s' is required for %s -- install it first", pkg, what), call. = FALSE)
}

# progress line: function, stage i/N, % done, elapsed minutes
.say <- function(verbose, fn, i, n, t0, ...) {
  if (!isTRUE(verbose)) return(invisible(NULL))
  message(sprintf("[%s | %5.1f min] %s %d/%d (%3.0f%%) %s",
                  format(Sys.time(), "%H:%M:%S"),
                  as.numeric(difftime(Sys.time(), t0, units = "mins")),
                  fn, i, n, 100 * i / n, paste0(...)))
}

# do.call with user overrides merged onto defaults; `force` always wins
# Arguments are passed as symbols from a private environment, so Seurat's
# command log (LogSeuratCommand) records a short call instead of deparsing an
# inlined Seurat object (that deparse takes minutes on real data).
.call <- function(f, defaults, extra = list(), force = list()) {
  args <- utils::modifyList(defaults, extra %||% list())
  args <- utils::modifyList(args, force)
  env  <- new.env(parent = parent.frame())
  syms <- vector("list", length(args)); names(syms) <- names(args)
  for (i in seq_along(args)) {
    v <- paste0(".arg_", names(args)[i])
    assign(v, args[[i]], envir = env)
    syms[[i]] <- as.name(v)
  }
  do.call(f, syms, envir = env)
}

.log_step <- function(log, step, action, detail = "") {
  rbind(log, data.frame(step = step, action = action, detail = detail, stringsAsFactors = FALSE))
}

.set_workers <- function(workers, globals_max_gb) {
  if (is.null(workers)) return(NULL)
  .need("future", "workers = <n>")
  old_plan <- future::plan(future::multisession, workers = workers)
  old_opt  <- options(future.globals.maxSize = globals_max_gb * 1024^3)
  list(plan = old_plan, opt = old_opt)
}
.restore_workers <- function(old) {
  if (is.null(old)) return(invisible(NULL))
  future::plan(old$plan); options(old$opt); invisible(NULL)
}

.empty_log <- function() {
  data.frame(step = character(0), action = character(0), detail = character(0), stringsAsFactors = FALSE)
}

.sort_genes <- function(x) x[order(tolower(x), x)]

.free_name <- function(existing, name, overwrite) {
  if (isTRUE(overwrite) || !name %in% existing) return(name)
  i <- 2L
  while (paste0(name, "_", i) %in% existing) i <- i + 1L
  paste0(name, "_", i)
}

.key <- function(name) {
  k <- gsub("[^A-Za-z0-9]", "", name)
  if (!grepl("^[A-Za-z]", k)) k <- paste0("k", k)
  paste0(k, "_")
}

.layers <- function(obj, assay, search = NA) {
  suppressWarnings(tryCatch(SeuratObject::Layers(obj[[assay]], search = search), error = function(e) character(0)))
}

.n_counts_layers <- function(obj, assay) length(.layers(obj, assay, "counts"))

.join_layers <- function(obj, assay) {
  obj[[assay]] <- SeuratObject::JoinLayers(obj[[assay]])
  obj
}

.layer_values <- function(obj, assay, layer, n_cells = 200) {
  lyr <- .layers(obj, assay, layer)
  if (!length(lyr)) return(numeric(0))
  m <- tryCatch(SeuratObject::LayerData(obj, assay = assay, layer = lyr[1]), error = function(e) NULL)
  if (is.null(m) || !length(m) || ncol(m) == 0) return(numeric(0))
  m <- m[, seq_len(min(ncol(m), n_cells)), drop = FALSE]
  if (inherits(m, "dgCMatrix")) m@x else as.numeric(as.matrix(m))
}

# log-normalized data are not whole numbers; raw counts (or data == counts) are
.is_lognormalized <- function(obj, assay) {
  v <- .layer_values(obj, assay, "data")
  length(v) > 0 && any(abs(v - round(v)) > 1e-6)
}

.is_sct <- function(obj, assay) {
  !is.null(assay) && assay %in% Assays(obj) && inherits(obj[[assay]], "SCTAssay")
}

.sct_models <- function(obj, assay) {
  if (!.is_sct(obj, assay)) return(character(0))
  ml <- slot(obj[[assay]], "SCTModel.list")
  nm <- names(ml)
  if (is.null(nm)) nm <- rep("", length(ml))
  nm[!nzchar(nm)] <- paste0("model", which(!nzchar(nm)))
  nm
}

.scaled_features <- function(obj, assay) {
  m <- suppressWarnings(tryCatch(SeuratObject::LayerData(obj, assay = assay, layer = "scale.data"), error = function(e) NULL))
  if (is.null(m)) character(0) else (rownames(m) %||% character(0))
}

.red_info <- function(obj, name) {
  if (is.null(name) || !name %in% Reductions(obj)) return(list(present = FALSE))
  r <- obj[[name]]
  list(present    = TRUE,
       assay      = tryCatch(DefaultAssay(r), error = function(e) NA_character_),
       ndims      = ncol(Embeddings(r)),
       n_loadings = nrow(Loadings(r)),
       has_model  = length(r@misc$model) > 0)
}

# reduction a UMAP was built from, read from the object's command log (NA if not logged)
.umap_source <- function(obj, umap_name) {
  src <- NA_character_
  for (cmd in obj@commands) {
    if (!startsWith(slot(cmd, "name"), "RunUMAP")) next
    p <- slot(cmd, "params")
    if (identical(p$reduction.name %||% "umap", umap_name)) src <- p$reduction %||% NA_character_
  }
  src
}

# accepts a Seurat object or the list returned by prep_query() / prep_reference()
.unwrap <- function(x, what) {
  if (inherits(x, "Seurat")) return(list(obj = x, report = NULL))
  if (is.list(x) && inherits(x[[what]], "Seurat")) return(list(obj = x[[what]], report = x$report))
  stop(sprintf("expected a Seurat object or the list returned by prep_%s()", what), call. = FALSE)
}


# ── object_status ─────────────────────────────────────────────────────────────
# What is on a Seurat object: per assay (class, layers, normalized?, scaled
# genes, variable genes, SCT models) and per reduction (assay, dims, loadings,
# UMAP model, reduction a UMAP was built from).

object_status <- function(obj, assays = NULL, reductions = NULL) {
  if (!inherits(obj, "Seurat")) stop("obj must be a Seurat object")
  assays <- assays %||% Assays(obj)
  assay_df <- do.call("rbind", lapply(assays, function(a) {
    data.frame(
      assay          = a,
      class          = class(obj[[a]])[1],
      default        = identical(DefaultAssay(obj), a),
      n_features     = nrow(obj[[a]]),
      counts_layers  = length(.layers(obj, a, "counts")),
      data_layers    = length(.layers(obj, a, "data")),
      log_normalized = .is_lognormalized(obj, a),
      n_scaled       = length(.scaled_features(obj, a)),
      n_variable     = length(tryCatch(VariableFeatures(obj[[a]]), error = function(e) character(0))),
      sct_models     = length(.sct_models(obj, a)),
      stringsAsFactors = FALSE)
  }))
  reductions <- reductions %||% Reductions(obj)
  red_df <- if (length(reductions)) do.call("rbind", lapply(reductions, function(r) {
    i <- .red_info(obj, r)
    data.frame(reduction = r, assay = i$assay, ndims = i$ndims, n_loadings = i$n_loadings,
               has_umap_model = i$has_model, umap_built_from = .umap_source(obj, r),
               stringsAsFactors = FALSE)
  })) else data.frame()
  list(n_cells = ncol(obj), assays = assay_df, reductions = red_df,
       meta_cols = colnames(obj@meta.data))
}


# ── small helpers  ─────────────────────────────────────────────────

add_orig_id <- function(obj, col = "orig_id") {
  if (!col %in% colnames(obj@meta.data))
    obj <- AddMetaData(obj, metadata = colnames(obj), col.name = col)
  obj
}

kl_divergence <- function(p, q) {
  p[is.na(p)] <- 0; q[is.na(q)] <- 0
  p <- p / sum(p); if (sum(p) == 0) return(Inf)
  q <- q / sum(q); if (sum(q) == 0) return(Inf)
  mask <- p > 0 & q > 0
  sum(p[mask] * log(p[mask] / q[mask])) + sum(q[mask] * log(q[mask] / p[mask]))
}

kl_arrange <- function(mat, anchor) {
  if (!anchor %in% colnames(mat)) stop(sprintf("'%s' not found in columns", anchor))
  ref_vec <- mat[, anchor]
  other   <- setdiff(colnames(mat), anchor)
  if (length(other) == 0) return(anchor)
  divs    <- sapply(other, function(col) kl_divergence(ref_vec, mat[, col]))
  c(anchor, names(sort(divs)))
}

# Drop rows whose id is duplicated within one metadata data frame.
ignore_duplication <- function(df, id_col = "obj_id", label = "") {
  dup_ids <- names(table(df[[id_col]]))[table(df[[id_col]]) > 1]
  n_dropped_ids   <- length(dup_ids)
  dropped_rows    <- df[[id_col]] %in% dup_ids
  n_dropped_cells <- sum(dropped_rows)
  if (n_dropped_ids > 0) {
    message(sprintf("%s: %d %s value(s) duplicated, dropping %d cell(s)",
                    label, n_dropped_ids, id_col, n_dropped_cells))
  } else {
    message(sprintf("%s: no duplicated %s values", label, id_col))
  }
  list(df = df[!dropped_rows, , drop = FALSE], dup_ids = dup_ids, n_dropped = n_dropped_cells)
}


# ── prep_reference ────────────────────────────────────────────────────────────
# Gets a reference ready for anchor projection, reusing whatever is already
# valid on it. Six stages: layers -> normalization -> features/scaling ->
# PCA -> clusters -> UMAP model (return.model = TRUE, built on the PCA used
# for mapping).
#
# assay           assay to anchor on. normalization = NULL picks "SCT" for an
#                 SCTAssay (or assay = "SCT"), otherwise "LogNormalize".
# counts_assay    assay holding raw counts, used when an SCT assay is built.
# build_*         "auto" = build only what is missing / not valid,
#                 "always" = rebuild, "never" = use as is (stop if missing).
# sct_model       which SCT model to keep when the SCT assay holds several
#                 (FindTransferAnchors needs exactly one); name or index.
# features        gene list for the PCA (NULL = variable features).
# vars_to_regress e.g. "donor_id"; used by SCTransform / ScaleData.
# overwrite       FALSE = new results go to free names ("pca_2", "umap_2",
#                 "seurat_clusters_2"); TRUE = replace.
# *_args          extra arguments passed straight to the Seurat function.
#
# Returns list(reference = <Seurat>, report = <list>);
# report$names holds the assay / pca / umap / cluster names to use downstream.

prep_reference <- function(
    reference,
    assay             = "RNA",
    normalization     = NULL,
    counts_assay      = "RNA",
    build_assay       = c("auto", "always", "never"),
    sct_model         = NULL,
    join_layers       = TRUE,
    features          = NULL,
    nfeatures         = 3000,
    vars_to_regress   = NULL,
    vst.flavor        = "v2",
    scale.factor      = 10000,
    npcs              = 30,
    dims              = 1:30,
    build_pca         = c("auto", "always"),
    build_clusters    = c("auto", "always", "never"),
    build_umap        = c("auto", "always"),
    pca_name          = "pca",
    umap_name         = "umap",
    cluster_name      = "seurat_clusters",
    resolution        = 0.8,
    cluster_algorithm = 1,
    k.param           = 20,
    overwrite         = FALSE,
    seed              = 42,
    verbose           = TRUE,
    workers           = NULL,
    globals_max_gb    = 8,
    normalize_args    = list(),
    hvg_args          = list(),
    scale_args        = list(),
    sct_args          = list(),
    pca_args          = list(),
    neighbor_args     = list(),
    cluster_args      = list(),
    umap_args         = list()
) {
  fn <- "prep_reference"; N <- 6; t0 <- Sys.time()
  if (!inherits(reference, "Seurat")) stop("reference must be a Seurat object")
  build_assay    <- match.arg(build_assay)
  build_pca      <- match.arg(build_pca)
  build_clusters <- match.arg(build_clusters)
  build_umap     <- match.arg(build_umap)
  normalization  <- normalization %||%
    (if (.is_sct(reference, assay) || identical(assay, "SCT")) "SCT" else "LogNormalize")
  normalization  <- match.arg(normalization, c("LogNormalize", "SCT"))
  if (max(dims) > npcs) stop(sprintf("max(dims) = %d is larger than npcs = %d", max(dims), npcs))
  old_workers <- .set_workers(workers, globals_max_gb)
  on.exit(.restore_workers(old_workers), add = TRUE)

  input <- object_status(reference)
  log   <- .empty_log()

  # 1. layers ------------------------------------------------------------------
  .say(verbose, fn, 1, N, t0, "layers")
  layer_assay <- if (normalization == "LogNormalize") assay else counts_assay
  if (layer_assay %in% Assays(reference) && .n_counts_layers(reference, layer_assay) > 1) {
    n_split <- .n_counts_layers(reference, layer_assay)
    if (join_layers) {
      reference <- .join_layers(reference, layer_assay)
      log <- .log_step(log, "layers", "joined", sprintf("%s: %d counts layers -> 1", layer_assay, n_split))
    } else {
      log <- .log_step(log, "layers", "kept split", sprintf("%s: %d counts layers (join_layers = FALSE)", layer_assay, n_split))
    }
  } else {
    log <- .log_step(log, "layers", "ok", sprintf("%s: no split counts layers", layer_assay))
  }

  # 2. normalization -------------------------------------------------------------
  .say(verbose, fn, 2, N, t0, "normalization (", normalization, ")")
  if (normalization == "LogNormalize") {
    if (!assay %in% Assays(reference))
      stop(sprintf("assay '%s' not found; available: %s", assay, paste(Assays(reference), collapse = ", ")))
    normalized <- .is_lognormalized(reference, assay)
    if (build_assay == "always" || (build_assay == "auto" && !normalized)) {
      reference <- .call("NormalizeData",
                         list(object = reference, assay = assay, normalization.method = "LogNormalize",
                              scale.factor = scale.factor, verbose = verbose),
                         normalize_args)
      log <- .log_step(log, "normalization", "built",
                       sprintf("NormalizeData on %s (scale.factor %g)", assay, scale.factor))
    } else if (!normalized) {
      stop(sprintf("assay '%s' has no log-normalized data layer and build_assay = 'never'", assay))
    } else {
      log <- .log_step(log, "normalization", "reused", sprintf("%s data layer already log-normalized", assay))
    }
  } else {
    has_sct <- .is_sct(reference, assay)
    if (build_assay == "always" || (build_assay == "auto" && !has_sct)) {
      if (!counts_assay %in% Assays(reference))
        stop(sprintf("counts_assay '%s' not found; needed to build the SCT assay", counts_assay))
      new_assay <- .free_name(Assays(reference), assay, overwrite)
      reference <- .call("SCTransform",
                         list(object = reference, assay = counts_assay, new.assay.name = new_assay,
                              vst.flavor = vst.flavor, vars.to.regress = vars_to_regress,
                              variable.features.n = nfeatures, seed.use = seed, verbose = verbose),
                         sct_args)
      log <- .log_step(log, "normalization", "built",
                       sprintf("SCTransform %s -> %s (vst.flavor %s; vars.to.regress: %s)", counts_assay, new_assay,
                               vst.flavor, paste(vars_to_regress %||% "none", collapse = ", ")))
      assay <- new_assay
    } else if (!has_sct) {
      stop(sprintf("assay '%s' is not an SCT assay and build_assay = 'never'", assay))
    } else {
      log <- .log_step(log, "normalization", "reused", sprintf("SCT assay %s", assay))
    }
  }
  DefaultAssay(reference) <- assay

  # SCT model choice: checked now, applied at the end (GetResidual needs all models)
  models    <- .sct_models(reference, assay)
  model_idx <- NULL
  if (length(models) > 1) {
    if (is.null(sct_model))
      stop(sprintf("assay '%s' has %d SCT models (%s); FindTransferAnchors needs one -- set sct_model to a model name or index",
                   assay, length(models), paste(models, collapse = ", ")))
    model_idx <- if (is.numeric(sct_model)) as.integer(sct_model) else match(sct_model, models)
    if (length(model_idx) != 1 || anyNA(model_idx) || model_idx < 1 || model_idx > length(models))
      stop("sct_model must name or index exactly one model of: ", paste(models, collapse = ", "))
  }

  # 3. features + scaling ----------------------------------------------------------
  .say(verbose, fn, 3, N, t0, "features + scaling")
  pca_info <- .red_info(reference, pca_name)
  pca_reusable <- build_pca == "auto" && is.null(features) && isTRUE(pca_info$present) &&
    isTRUE(pca_info$n_loadings > 0) && isTRUE(pca_info$ndims >= max(dims)) &&
    identical(pca_info$assay, assay) &&
    all(rownames(Loadings(reference[[pca_name]])) %in% rownames(reference[[assay]]))

  dropped_features <- character(0)
  if (pca_reusable) {
    feats <- rownames(Loadings(reference[[pca_name]]))
    log <- .log_step(log, "features", "reused", sprintf("%d genes from %s loadings", length(feats), pca_name))
    log <- .log_step(log, "scaling", "skipped", "PCA reused")
  } else {
    if (!is.null(features)) {
      feats <- intersect(features, rownames(reference[[assay]]))
      dropped_features <- .sort_genes(setdiff(features, feats))
      log <- .log_step(log, "features", "user list",
                       sprintf("%d of %d genes present in %s", length(feats), length(features), assay))
    } else {
      feats <- tryCatch(VariableFeatures(reference[[assay]]), error = function(e) character(0))
      rebuild_hvg <- normalization == "LogNormalize" && (build_assay == "always" || !length(feats))
      if (rebuild_hvg) {
        reference <- .call("FindVariableFeatures",
                           list(object = reference, assay = assay, nfeatures = nfeatures, verbose = verbose),
                           hvg_args)
        feats <- VariableFeatures(reference[[assay]])
        log <- .log_step(log, "features", "built", sprintf("FindVariableFeatures: %d genes", length(feats)))
      } else if (!length(feats) && normalization == "SCT") {
        feats <- .scaled_features(reference, assay)
        log <- .log_step(log, "features", "from scale.data", sprintf("%d genes in SCT scale.data", length(feats)))
      } else {
        log <- .log_step(log, "features", "reused", sprintf("%d variable features on %s", length(feats), assay))
      }
    }
    if (!length(feats)) stop("no features left for the PCA")

    if (normalization == "LogNormalize") {
      reference <- .call("ScaleData",
                         list(object = reference, assay = assay, features = feats,
                              vars.to.regress = vars_to_regress, verbose = verbose),
                         scale_args)
      log <- .log_step(log, "scaling", "built",
                       sprintf("ScaleData %d genes (vars.to.regress: %s)", length(feats),
                               paste(vars_to_regress %||% "none", collapse = ", ")))
    } else {
      missing <- setdiff(feats, .scaled_features(reference, assay))
      if (length(missing)) {
        reference <- GetResidual(reference, features = missing, assay = assay, verbose = verbose)
        log <- .log_step(log, "scaling", "built", sprintf("GetResidual for %d genes", length(missing)))
      } else {
        log <- .log_step(log, "scaling", "reused", "all features in SCT scale.data")
      }
      feats <- intersect(feats, .scaled_features(reference, assay))
    }
  }

  # 4. PCA -----------------------------------------------------------------------------
  .say(verbose, fn, 4, N, t0, "PCA")
  if (pca_reusable) {
    log <- .log_step(log, "pca", "reused", sprintf("%s (%d dims, %d genes)", pca_name, pca_info$ndims, length(feats)))
  } else {
    pca_out <- .free_name(Reductions(reference), pca_name, overwrite)
    reference <- .call("RunPCA",
                       list(object = reference, assay = assay, features = feats, npcs = npcs,
                            reduction.name = pca_out, reduction.key = .key(pca_out), seed.use = seed,
                            verbose = verbose, ndims.print = 1, nfeatures.print = 5),
                       pca_args)
    log <- .log_step(log, "pca", "built", sprintf("%s: %d PCs on %d genes", pca_out, npcs, length(feats)))
    pca_name <- pca_out
  }

  # 5. clusters ------------------------------------------------------------------------
  .say(verbose, fn, 5, N, t0, "clusters")
  has_cl <- cluster_name %in% colnames(reference@meta.data)
  if (build_clusters == "never" || (build_clusters == "auto" && has_cl)) {
    log <- .log_step(log, "clusters", if (has_cl) "reused" else "skipped",
                     if (has_cl) sprintf("%s already in metadata", cluster_name) else "build_clusters = 'never'")
    if (!has_cl) cluster_name <- NA_character_
  } else {
    cl_out <- .free_name(colnames(reference@meta.data), cluster_name, overwrite)
    graphs <- c(.free_name(Graphs(reference), paste0(pca_name, "_nn"), overwrite),
                .free_name(Graphs(reference), paste0(pca_name, "_snn"), overwrite))
    keep_sc <- reference@meta.data[["seurat_clusters"]]
    keep_id <- Idents(reference)
    reference <- .call("FindNeighbors",
                       list(object = reference, reduction = pca_name, dims = dims, k.param = k.param,
                            graph.name = graphs, verbose = verbose),
                       neighbor_args)
    reference <- .call("FindClusters",
                       list(object = reference, graph.name = graphs[2], resolution = resolution,
                            algorithm = cluster_algorithm, random.seed = seed, verbose = verbose),
                       cluster_args)
    new_cl <- Idents(reference)
    if (cl_out != "seurat_clusters") reference@meta.data[["seurat_clusters"]] <- keep_sc
    reference@meta.data[[cl_out]] <- new_cl
    Idents(reference) <- keep_id
    log <- .log_step(log, "clusters", "built",
                     sprintf("%s: %d clusters (resolution %g, graph %s)", cl_out, nlevels(new_cl), resolution, graphs[2]))
    cluster_name <- cl_out
  }

  # 6. UMAP model ----------------------------------------------------------------------
  .say(verbose, fn, 6, N, t0, "UMAP model")
  um     <- .red_info(reference, umap_name)
  um_src <- if (isTRUE(um$present)) .umap_source(reference, umap_name) else NA_character_
  um_reusable <- build_umap == "auto" && isTRUE(um$present) && isTRUE(um$has_model) &&
    (identical(um_src, pca_name) || (is.na(um_src) && pca_reusable))
  if (um_reusable) {
    log <- .log_step(log, "umap", "reused",
                     sprintf("%s (model present; built from %s)", umap_name,
                             if (is.na(um_src)) "a reduction not in the command log" else um_src))
  } else {
    reason <- if (!isTRUE(um$present)) "absent" else if (!isTRUE(um$has_model)) "no UMAP model" else
      sprintf("built from %s, not %s", if (is.na(um_src)) "an unlogged reduction" else um_src, pca_name)
    um_out <- .free_name(Reductions(reference), umap_name, overwrite)
    reference <- .call("RunUMAP",
                       list(object = reference, reduction = pca_name, dims = dims, reduction.name = um_out,
                            reduction.key = .key(um_out), seed.use = seed, verbose = verbose),
                       umap_args,
                       force = list(return.model = TRUE))
    log <- .log_step(log, "umap", "built",
                     sprintf("%s on %s dims %d-%d, return.model = TRUE (previous %s: %s)",
                             um_out, pca_name, min(dims), max(dims), umap_name, reason))
    umap_name <- um_out
  }

  # keep one SCT model -------------------------------------------------------------------
  sct_kept <- NA_character_
  if (!is.null(model_idx)) {
    sct <- reference[[assay]]
    slot(sct, "SCTModel.list") <- slot(sct, "SCTModel.list")[model_idx]
    reference[[assay]] <- sct
    sct_kept <- models[model_idx]
    log <- .log_step(log, "sct model", "selected", sprintf("kept %s of %d models", sct_kept, length(models)))
  }

  report <- list(
    fn       = fn,
    names    = list(assay = assay, normalization = normalization, pca = pca_name, umap = umap_name,
                    clusters = cluster_name, sct_model = sct_kept),
    steps    = log,
    features = list(n = length(feats), genes = .sort_genes(feats), dropped = dropped_features),
    input    = input,
    output   = object_status(reference),
    params   = list(nfeatures = nfeatures, npcs = npcs, dims = dims, vars_to_regress = vars_to_regress,
                    vst.flavor = vst.flavor, scale.factor = scale.factor, resolution = resolution,
                    cluster_algorithm = cluster_algorithm, k.param = k.param, seed = seed,
                    build_assay = build_assay, build_pca = build_pca, build_clusters = build_clusters,
                    build_umap = build_umap, overwrite = overwrite, join_layers = join_layers,
                    workers = workers %||% "caller's plan"),
    seconds  = as.numeric(difftime(Sys.time(), t0, units = "secs"))
  )
  if (verbose) print(log, row.names = FALSE)
  list(reference = reference, report = report)
}


# ── prep_query ────────────────────────────────────────────────────────────────
# Gets a query ready for run_anchors(). Nothing on the query is overwritten:
# its own PCA / UMAP / clusters stay as they are.
#   LogNormalize : NormalizeData only if the data layer is not normalized yet
#                  (the query's own variable genes / scaling are not used by
#                  FindTransferAnchors, so they are not computed).
#   SCT          : only raw counts are needed; FindTransferAnchors computes the
#                  query residuals with the reference SCT model
#                  (run_anchors(recompute.residuals = TRUE)).
# reference      a Seurat object or the prep_reference() output; used to read
#                the normalization and to report which anchor genes the query
#                lacks.
#
# Returns list(query = <Seurat>, report = <list>).

prep_query <- function(
    query,
    reference       = NULL,
    assay           = "RNA",
    normalization   = NULL,
    reference_assay = NULL,
    pca_name        = NULL,
    join_layers     = TRUE,
    normalize       = c("auto", "always", "never"),
    scale.factor    = 10000,
    verbose         = TRUE,
    workers         = NULL,
    globals_max_gb  = 8,
    normalize_args  = list()
) {
  fn <- "prep_query"; N <- 3; t0 <- Sys.time()
  query     <- .unwrap(query, "query")$obj
  ref       <- if (!is.null(reference)) .unwrap(reference, "reference") else NULL
  normalize <- match.arg(normalize)
  old_workers <- .set_workers(workers, globals_max_gb)
  on.exit(.restore_workers(old_workers), add = TRUE)
  if (!assay %in% Assays(query))
    stop(sprintf("assay '%s' not found on the query; available: %s", assay, paste(Assays(query), collapse = ", ")))
  ref_assay <- reference_assay %||% ref$report$names$assay %||% (if (!is.null(ref)) DefaultAssay(ref$obj) else NULL)
  normalization <- normalization %||%
    (if (!is.null(ref) && .is_sct(ref$obj, ref_assay)) "SCT" else "LogNormalize")
  normalization <- match.arg(normalization, c("LogNormalize", "SCT"))

  input <- object_status(query)
  log   <- .empty_log()

  # 1. layers ----------------------------------------------------------------------
  .say(verbose, fn, 1, N, t0, "layers")
  n_split <- .n_counts_layers(query, assay)
  if (n_split > 1 && join_layers) {
    query <- .join_layers(query, assay)
    log <- .log_step(log, "layers", "joined", sprintf("%s: %d counts layers -> 1", assay, n_split))
  } else if (n_split > 1) {
    log <- .log_step(log, "layers", "kept split", sprintf("%s: %d counts layers (join_layers = FALSE)", assay, n_split))
  } else if (n_split == 0) {
    stop(sprintf("assay '%s' has no counts layer", assay))
  } else {
    log <- .log_step(log, "layers", "ok", sprintf("%s: one counts layer", assay))
  }

  # 2. normalization ---------------------------------------------------------------
  .say(verbose, fn, 2, N, t0, "normalization (", normalization, ")")
  if (normalization == "LogNormalize") {
    normalized <- .is_lognormalized(query, assay)
    if (normalize == "always" || (normalize == "auto" && !normalized)) {
      query <- .call("NormalizeData",
                     list(object = query, assay = assay, normalization.method = "LogNormalize",
                          scale.factor = scale.factor, verbose = verbose),
                     normalize_args)
      log <- .log_step(log, "normalization", "built", sprintf("NormalizeData on %s (scale.factor %g)", assay, scale.factor))
    } else if (!normalized) {
      stop(sprintf("assay '%s' is not log-normalized and normalize = 'never'", assay))
    } else {
      log <- .log_step(log, "normalization", "reused", sprintf("%s data layer already log-normalized", assay))
    }
  } else {
    log <- .log_step(log, "normalization", "deferred",
                     sprintf("SCT residuals from %s counts are computed in run_anchors() with the reference model", assay))
  }
  DefaultAssay(query) <- assay

  # 3. genes vs reference ------------------------------------------------------------
  .say(verbose, fn, 3, N, t0, "genes vs reference")
  genes <- list(n_anchor_genes = NA_integer_, n_present = NA_integer_, missing = character(0),
                n_query_only = NA_integer_)
  pca_name <- pca_name %||% ref$report$names$pca %||% "pca"
  if (!is.null(ref) && pca_name %in% Reductions(ref$obj)) {
    anchor_genes <- rownames(Loadings(ref$obj[[pca_name]]))
    present      <- intersect(anchor_genes, rownames(query[[assay]]))
    genes <- list(n_anchor_genes = length(anchor_genes), n_present = length(present),
                  missing = .sort_genes(setdiff(anchor_genes, present)),
                  n_query_only = length(setdiff(rownames(query[[assay]]), rownames(ref$obj))))
    log <- .log_step(log, "genes", "checked",
                     sprintf("%d of %d reference %s genes in query %s", length(present), length(anchor_genes), pca_name, assay))
  } else {
    log <- .log_step(log, "genes", "skipped", "no reference PCA given")
  }

  report <- list(
    fn      = fn,
    names   = list(assay = assay, normalization = normalization, reference_assay = ref_assay, pca = pca_name),
    steps   = log,
    genes   = genes,
    input   = input,
    output  = object_status(query),
    params  = list(scale.factor = scale.factor, join_layers = join_layers, normalize = normalize),
    seconds = as.numeric(difftime(Sys.time(), t0, units = "secs"))
  )
  if (verbose) print(log, row.names = FALSE)
  list(query = query, report = report)
}


# ── tag_ids / verify_ids ──────────────────────────────────────────────────────
# IDs are tagged once, inside run_anchors(), right after MapQuery. Each object
# gets <prefix>_id_backup (never touched again) and <prefix>_id_surface (used
# for matching everywhere downstream). verify_ids() forces id_surface back to
# id_backup if they ever drift apart.

tag_ids <- function(obj, prefix) {
  ids <- colnames(obj)
  obj <- AddMetaData(obj, metadata = ids, col.name = paste0(prefix, "_id_backup"))
  obj <- AddMetaData(obj, metadata = ids, col.name = paste0(prefix, "_id_surface"))
  obj
}

verify_ids <- function(obj, prefix) {
  backup_col  <- paste0(prefix, "_id_backup")
  surface_col <- paste0(prefix, "_id_surface")
  if (!backup_col %in% colnames(obj@meta.data) || !surface_col %in% colnames(obj@meta.data))
    stop(sprintf("'%s'/'%s' not found -- run_anchors() must tag ids before this step", backup_col, surface_col))
  if (!identical(obj@meta.data[[surface_col]], obj@meta.data[[backup_col]])) {
    warning(sprintf("'%s' out of sync with '%s' -- forcing '%s' to match backup", surface_col, backup_col, surface_col))
    obj@meta.data[[surface_col]] <- obj@meta.data[[backup_col]]
  }
  obj
}


# ── run_anchors ───────────────────────────────────────────────────────────────
# FindTransferAnchors + MapQuery (projection only; label transfer happens in
# update_anchors()). query / reference can be Seurat objects or the outputs of
# prep_query() / prep_reference(); names (assay, PCA, UMAP model) are then read
# from their reports.
#
# recompute.residuals  SCT only: TRUE = query residuals from the reference SCT
#                      model (query.assay = raw-count assay); FALSE = use the
#                      query's own SCT assay (set query.assay to it).
# features             NULL = reference PCA genes present in the query.
# new.reduction.name   name of the projected UMAP on the query ("ref.umap").
# integrated_name      name of MapQuery's integrated embedding (NULL = "ref.<pca>").
# anchor_args / mapquery_args  extra arguments for FindTransferAnchors / MapQuery.
#
# Returns list(anchorset, query, reference, report).

run_anchors <- function(
    query,
    reference,
    reference.assay      = NULL,
    query.assay          = NULL,
    normalization.method = NULL,
    recompute.residuals  = TRUE,
    reference.reduction  = NULL,
    reduction.model      = NULL,
    dims                 = 1:30,
    features             = NULL,
    k.anchor             = 5,
    k.filter             = NA,
    k.score              = 30,
    max.features         = 200,
    new.reduction.name   = "ref.umap",
    integrated_name      = NULL,
    verbose              = TRUE,
    workers              = NULL,
    globals_max_gb       = 8,
    anchor_args          = list(),
    mapquery_args        = list()
) {
  fn <- "run_anchors"; N <- 3; t0 <- Sys.time()
  q  <- .unwrap(query, "query");         qo <- q$obj
  r  <- .unwrap(reference, "reference"); ro <- r$obj

  reference.assay      <- reference.assay %||% r$report$names$assay %||% DefaultAssay(ro)
  query.assay          <- query.assay %||% q$report$names$assay %||% DefaultAssay(qo)
  normalization.method <- normalization.method %||% (if (.is_sct(ro, reference.assay)) "SCT" else "LogNormalize")
  reference.reduction  <- reference.reduction %||% r$report$names$pca %||% "pca"
  reduction.model      <- reduction.model %||% r$report$names$umap %||% "umap"
  integrated_name      <- integrated_name %||% paste0("ref.", reference.reduction)
  old_workers <- .set_workers(workers, globals_max_gb)
  on.exit(.restore_workers(old_workers), add = TRUE)

  if (!reference.assay %in% Assays(ro)) stop(sprintf("reference.assay '%s' not on the reference", reference.assay))
  if (!query.assay %in% Assays(qo))     stop(sprintf("query.assay '%s' not on the query", query.assay))
  if (!reference.reduction %in% Reductions(ro))
    stop(sprintf("reference.reduction '%s' not on the reference -- run prep_reference() first", reference.reduction))
  if (nrow(Loadings(ro[[reference.reduction]])) == 0)
    stop(sprintf("'%s' has no feature loadings; PCA projection needs a linear (orthogonal) reduction", reference.reduction))
  if (max(dims) > ncol(Embeddings(ro[[reference.reduction]])))
    stop(sprintf("max(dims) = %d but '%s' has %d dims", max(dims), reference.reduction, ncol(Embeddings(ro[[reference.reduction]]))))
  if (!is.null(reduction.model) &&
      (!reduction.model %in% Reductions(ro) || length(ro[[reduction.model]]@misc$model) == 0)) {
    warning(sprintf("reference has no UMAP model '%s' (RunUMAP return.model = TRUE) -- UMAP projection skipped", reduction.model))
    reduction.model <- NULL
  }
  if (normalization.method == "SCT") {
    n_models <- length(.sct_models(ro, reference.assay))
    if (n_models > 1)
      stop(sprintf("reference SCT assay has %d models; keep one with prep_reference(sct_model = ...)", n_models))
    if (!recompute.residuals && !.is_sct(qo, query.assay))
      stop("recompute.residuals = FALSE needs query.assay to be an SCT assay")
  }

  load_genes <- rownames(Loadings(ro[[reference.reduction]]))
  features   <- intersect(features %||% load_genes, load_genes)
  features   <- intersect(features, rownames(qo[[query.assay]]))
  if (!length(features)) stop("no reference PCA genes found in the query")
  dropped    <- .sort_genes(setdiff(load_genes, features))

  DefaultAssay(ro) <- reference.assay
  DefaultAssay(qo) <- query.assay

  .say(verbose, fn, 1, N, t0, sprintf("FindTransferAnchors (%s, %d genes, dims %d-%d)",
                                      normalization.method, length(features), min(dims), max(dims)))
  anchorset <- .call("FindTransferAnchors",
                     list(reference = ro, query = qo, normalization.method = normalization.method,
                          recompute.residuals = recompute.residuals, reference.assay = reference.assay,
                          query.assay = query.assay, reference.reduction = reference.reduction,
                          dims = dims, features = features, k.anchor = k.anchor, k.filter = k.filter,
                          k.score = k.score, max.features = max.features, verbose = verbose),
                     anchor_args)

  .say(verbose, fn, 2, N, t0, "MapQuery")
  qo <- .call("MapQuery",
              list(anchorset = anchorset, reference = ro, query = qo,
                   reference.reduction = reference.reduction, reduction.model = reduction.model,
                   new.reduction.name = integrated_name,
                   projectumap.args = list(reduction.name = new.reduction.name,
                                           reduction.key  = .key(new.reduction.name)),
                   verbose = verbose),
              mapquery_args)

  .say(verbose, fn, 3, N, t0, "tagging ids")
  qo <- tag_ids(qo, "query")
  ro <- tag_ids(ro, "ref")

  anc <- slot(anchorset, "anchors")
  report <- list(
    fn       = fn,
    names    = list(reference.assay = reference.assay, query.assay = query.assay,
                    normalization.method = normalization.method, reference.reduction = reference.reduction,
                    reduction.model = reduction.model %||% NA_character_, integrated = integrated_name,
                    ref_umap = if (is.null(reduction.model)) NA_character_ else new.reduction.name),
    features = list(n = length(features), dropped = dropped),
    anchors  = list(n_pairs = nrow(anc), n_query_cells = length(unique(anc[, "cell2"])),
                    n_ref_cells = length(unique(anc[, "cell1"])), n_query_total = ncol(qo)),
    params   = list(dims = dims, recompute.residuals = recompute.residuals, k.anchor = k.anchor,
                    k.filter = k.filter, k.score = k.score, max.features = max.features),
    seconds  = as.numeric(difftime(Sys.time(), t0, units = "secs"))
  )
  if (verbose) message(sprintf("run_anchors: %d anchor pairs, %d of %d query cells anchored",
                               report$anchors$n_pairs, report$anchors$n_query_cells, ncol(qo)))
  list(anchorset = anchorset, query = qo, reference = ro, report = report)
}


# ── update_anchors ────────────────────────────────────────────────────────────
# Keeps only light pieces of the anchorset / query / reference: metadata
# data.frames + an anchor pair table, plus TransferData predictions for
# `variables`. The result never holds the AnchorSet, a Seurat object or
# expression data:
#   $query$meta  query metadata + query_id_surface/backup + predicted.<var>(.score)
#   $ref$meta    reference metadata + ref_id_surface/backup
#   $anchor_df   one row per anchor pair: query_id_surface, ref_id_surface, score, ...

update_anchors <- function(
    anchorset,
    query,
    reference,
    variables           = NULL,   # ref meta cols to run TransferData on
    query_meta_cols     = NULL,   # query meta cols to attach to anchor_df
    reference_meta_cols = NULL,   # ref meta cols to attach to anchor_df
    dims                = 1:30,
    k.weight            = 50,
    verbose             = TRUE,
    transfer_args       = list()
) {
  query     <- verify_ids(.unwrap(query, "query")$obj, "query")
  reference <- verify_ids(.unwrap(reference, "reference")$obj, "ref")

  for (col in c(variables, reference_meta_cols))
    if (!col %in% colnames(reference@meta.data)) stop(sprintf("'%s' not found in reference metadata", col))
  for (col in query_meta_cols)
    if (!col %in% colnames(query@meta.data)) stop(sprintf("'%s' not found in query metadata", col))

  adf <- as.data.frame(slot(anchorset, "anchors"))
  adf$query_id_surface <- query@meta.data$query_id_surface[adf$cell2]
  adf$ref_id_surface   <- reference@meta.data$ref_id_surface[adf$cell1]
  adf$cell1 <- NULL
  adf$cell2 <- NULL

  for (col in query_meta_cols)
    adf[[paste0("query.", col)]] <- query@meta.data[[col]][match(adf$query_id_surface, query@meta.data$query_id_surface)]
  for (col in reference_meta_cols)
    adf[[paste0("ref.", col)]] <- reference@meta.data[[col]][match(adf$ref_id_surface, reference@meta.data$ref_id_surface)]

  q_meta <- query@meta.data
  r_meta <- reference@meta.data

  if (!is.null(variables)) {
    refdata <- setNames(lapply(variables, function(v) reference@meta.data[[v]]), variables)
    preds <- .call("TransferData",
                   list(anchorset = anchorset, refdata = refdata, dims = dims, k.weight = k.weight, verbose = verbose),
                   transfer_args)
    if (is.data.frame(preds)) preds <- setNames(list(preds), variables)   # one variable -> data.frame
    for (v in names(preds)) {
      q_meta[[paste0("predicted.", v)]]          <- preds[[v]]$predicted.id
      q_meta[[paste0("predicted.", v, ".score")]] <- preds[[v]]$prediction.score.max
    }
  }

  list(query = list(meta = q_meta), ref = list(meta = r_meta), anchor_df = adf)
}


# ── transfer_reduction ────────────────────────────────────────────────────────
# Attaches embeddings as plain data.frames (no expression data).
#   Path A - 'embeddings': named list of matrices/data.frames (groups x dims)
#            with rownames = ids, matched against anchors[[side]]$meta[[id_col]]
#            (e.g. pseudobulk coordinates per sample_id / donor_id).
#   Path B - 'query' / 'reference' Seurat objects: per-cell embeddings keyed by
#            query_id_surface / ref_id_surface (or id_col).

transfer_reduction <- function(
    anchors,
    query             = NULL,
    reference         = NULL,
    query_reductions  = NULL,   # NULL = all reductions on the query
    ref_reductions    = NULL,   # NULL = all reductions on the reference
    embeddings        = NULL,
    side              = c("auto", "query", "ref"),
    id_col            = NULL
) {
  side <- match.arg(side)

  if (!is.null(embeddings)) {
    if (is.null(id_col))
      stop("id_col must be specified when passing 'embeddings' directly -- which meta column do these ids match?")
    embeddings <- embeddings[!sapply(embeddings, is.null)]
    if (length(embeddings) == 0) stop("embeddings is empty after dropping NULL entries -- nothing to attach")
    ids_emb <- rownames(embeddings[[1]])
    if (is.null(ids_emb)) stop("embeddings must have rownames to match against id_col")
    if (side == "auto") {
      q_overlap <- if (!is.null(anchors$query$meta[[id_col]])) sum(anchors$query$meta[[id_col]] %in% ids_emb) else 0
      r_overlap <- if (!is.null(anchors$ref$meta[[id_col]]))   sum(anchors$ref$meta[[id_col]]   %in% ids_emb) else 0
      if (q_overlap == 0 && r_overlap == 0)
        stop(sprintf("none of embeddings' rownames matched anchors$query$meta$%s or anchors$ref$meta$%s", id_col, id_col))
      side <- if (q_overlap >= r_overlap) "query" else "ref"
    }
    meta_ids <- anchors[[side]]$meta[[id_col]]
    if (is.null(meta_ids)) stop(sprintf("anchors$%s$meta$%s not found", side, id_col))
    matched <- lapply(embeddings, function(emb) {
      emb <- as.data.frame(emb)
      emb[match(meta_ids, rownames(emb)), , drop = FALSE]
    })
    anchors[[side]]$emb <- utils::modifyList(anchors[[side]]$emb %||% list(), matched)
    return(anchors)
  }

  if (is.null(query) || is.null(reference))
    stop("either 'embeddings' or both 'query' and 'reference' must be supplied")
  query     <- verify_ids(.unwrap(query, "query")$obj, "query")
  reference <- verify_ids(.unwrap(reference, "reference")$obj, "ref")

  if (is.null(query_reductions)) query_reductions <- Reductions(query)
  if (is.null(ref_reductions))   ref_reductions   <- Reductions(reference)

  q_id_col <- id_col %||% "query_id_surface"
  r_id_col <- id_col %||% "ref_id_surface"

  .emb_list <- function(obj, reds, id_column) {
    ids <- obj@meta.data[[id_column]]
    setNames(lapply(reds, function(r) {
      emb <- as.data.frame(Embeddings(obj, reduction = r))
      colnames(emb) <- paste0(toupper(r), "_", seq_len(ncol(emb)))
      rownames(emb) <- ids
      emb
    }), reds)
  }

  anchors$query$emb <- .emb_list(query,     query_reductions, q_id_col)
  anchors$ref$emb   <- .emb_list(reference, ref_reductions,   r_id_col)
  anchors
}


# ── max_anchors ───────────────────────────────────────────────────────────────
# Best-scoring reference cell per query cell -> anchors$max_anchor_df.

max_anchors <- function(anchors) {
  adf    <- anchors$anchor_df
  sorted <- adf[order(adf$score, decreasing = TRUE), ]
  max_df <- sorted[!duplicated(sorted$query_id_surface), ]
  rownames(max_df) <- NULL
  anchors$max_anchor_df <- max_df
  anchors
}


# ── transfer_metadata ─────────────────────────────────────────────────────────
# Joins one numeric column from any object (Seurat, SingleCellExperiment /
# cell_data_set, data.frame / list, or a named numeric vector) into anchors by
# id. Side (query / ref) is picked from the id overlap unless `side` is set.

.extract_id_value <- function(obj, col = NULL, id_col = NULL) {
  if (is.atomic(obj) && is.null(dim(obj)) && !is.null(names(obj)))
    return(list(id = names(obj), value = unname(obj)))
  if (inherits(obj, "Seurat")) {
    if (is.null(col)) stop("col must be specified for a Seurat object")
    if (is.null(id_col)) {
      obj    <- add_orig_id(obj)
      id_col <- "orig_id"
    }
    md <- obj@meta.data
    return(list(id = md[[id_col]], value = md[[col]]))
  }
  if (methods::is(obj, "SummarizedExperiment") || inherits(obj, "cell_data_set")) {
    .need("SummarizedExperiment", "transfer_metadata() on SingleCellExperiment objects")
    if (is.null(col)) stop("col must be specified for a SingleCellExperiment/cell_data_set object")
    cd  <- SummarizedExperiment::colData(obj)
    ids <- if (!is.null(id_col)) cd[[id_col]] else rownames(cd)
    return(list(id = ids, value = cd[[col]]))
  }
  if (is.data.frame(obj) || is.list(obj)) {
    if (is.null(col)) stop("col must be specified for a data.frame/list object")
    ids <- if (!is.null(id_col)) obj[[id_col]] else rownames(obj)
    if (is.null(ids)) stop("could not determine ids: supply id_col")
    return(list(id = ids, value = obj[[col]]))
  }
  stop(sprintf("unsupported object class: %s", paste(class(obj), collapse = "/")))
}

transfer_metadata <- function(
    anchors,
    obj,
    col      = NULL,
    id_col   = NULL,
    side     = c("auto", "query", "ref"),
    out_name = NULL     # column name to write (default: col, or the vector's name)
) {
  side <- match.arg(side)
  extracted <- .extract_id_value(obj, col, id_col)
  ids  <- extracted$id
  vals <- extracted$value
  if (!is.numeric(vals)) stop(sprintf("'%s' is not numeric in obj", col %||% "<vector>"))

  if (side == "auto") {
    q_overlap <- sum(ids %in% anchors$query$meta$query_id_surface)
    r_overlap <- sum(ids %in% anchors$ref$meta$ref_id_surface)
    if (q_overlap == 0 && r_overlap == 0)
      stop("none of obj's ids matched query_id_surface or ref_id_surface")
    side <- if (q_overlap >= r_overlap) "query" else "ref"
  }

  id_field  <- if (side == "query") "query_id_surface" else "ref_id_surface"
  df_prefix <- if (side == "query") "query." else "ref."
  out_name  <- out_name %||% (if (is.null(col)) deparse(substitute(obj)) else col)
  out_col   <- paste0(df_prefix, out_name)

  m <- anchors[[side]]$meta
  m[[out_name]] <- vals[match(m[[id_field]], ids)]
  anchors[[side]]$meta <- m

  if (!is.null(anchors$anchor_df))
    anchors$anchor_df[[out_col]] <- vals[match(anchors$anchor_df[[id_field]], ids)]
  if (!is.null(anchors$max_anchor_df))
    anchors$max_anchor_df[[out_col]] <- vals[match(anchors$max_anchor_df[[id_field]], ids)]
  anchors
}


# ── filter_anchor_df ──────────────────────────────────────────────────────────
# Score window filter used by the plotting functions.

filter_anchor_df <- function(df, min_score = NULL, max_score = NULL) {
  if (!is.null(min_score)) df <- df[df$score >= min_score, ]
  if (!is.null(max_score)) df <- df[df$score <= max_score, ]
  df
}

.add_title <- function(p, title, default) {
  title <- title %||% default
  if (is.na(title) || !nzchar(title)) p else p + ggtitle(title)
}


# ── density_plot ──────────────────────────────────────────────────────────────
# Query cell density on the projected reference UMAP (needs transfer_reduction()).
# title: NULL = default text, NA or "" = no title.

density_plot <- function(
    anchors,
    anchored_only   = TRUE,
    query.reduction = "ref.umap",
    query_subset    = NULL,
    ref_subset      = NULL,
    min_score       = NULL,
    max_score       = NULL,
    n               = 200,
    pt.size         = 0.3,
    bg.pt.size      = pt.size,
    bg.color        = "lightgrey",
    viridis_option  = "D",
    title           = NULL
) {
  .need("MASS", "density_plot()")
  if (is.null(anchors$max_anchor_df)) stop("anchors$max_anchor_df not found -- run max_anchors() first")
  q_emb <- anchors$query$emb[[query.reduction]]
  if (is.null(q_emb)) stop(sprintf("anchors$query$emb$%s not found -- run transfer_reduction() first", query.reduction))

  q_meta <- anchors$query$meta
  adf    <- filter_anchor_df(anchors$max_anchor_df, min_score, max_score)

  all_coords <- q_emb
  colnames(all_coords)[1:2] <- c("UMAP_1", "UMAP_2")
  all_coords$query_id_surface <- rownames(all_coords)

  if (is.logical(query_subset)) query_subset <- q_meta$query_id_surface[query_subset]
  if (is.logical(ref_subset))   ref_subset   <- anchors$ref$meta$ref_id_surface[ref_subset]

  use_ids <- if (anchored_only) unique(adf$query_id_surface) else q_meta$query_id_surface
  if (!is.null(query_subset)) use_ids <- intersect(use_ids, query_subset)
  if (!is.null(ref_subset))   use_ids <- intersect(use_ids, adf$query_id_surface[adf$ref_id_surface %in% ref_subset])

  plot_coords <- all_coords[all_coords$query_id_surface %in% use_ids, ]
  dens <- MASS::kde2d(plot_coords$UMAP_1, plot_coords$UMAP_2, n = n,
                      lims = c(range(all_coords$UMAP_1), range(all_coords$UMAP_2)))
  ix <- pmin(pmax(findInterval(plot_coords$UMAP_1, dens$x), 1), nrow(dens$z))
  iy <- pmin(pmax(findInterval(plot_coords$UMAP_2, dens$y), 1), ncol(dens$z))
  plot_coords$density <- dens$z[cbind(ix, iy)]

  p <- ggplot() +
    geom_point(data = all_coords,  aes(UMAP_1, UMAP_2), color = bg.color, size = bg.pt.size) +
    geom_point(data = plot_coords, aes(UMAP_1, UMAP_2, color = density), size = pt.size) +
    scale_color_viridis_c(option = viridis_option) + theme_classic()
  .add_title(p, title, ifelse(anchored_only, "Anchored query cell density on reference UMAP",
                              "All query cell density on reference UMAP"))
}


# ── max_anchor_plot ───────────────────────────────────────────────────────────
# Highest anchor score per cell on the reference or query UMAP.

max_anchor_plot <- function(
    anchors,
    umap_source     = c("reference", "query"),
    query_subset    = NULL,
    ref_subset      = NULL,
    query.reduction = "umap",
    ref.reduction   = "umap",
    min_score       = NULL,
    max_score       = NULL,
    pt.size         = 0.3,
    bg.pt.size      = pt.size,
    bg.color        = "lightgrey",
    viridis_option  = "D",
    limits          = NULL,     # colour scale limits, e.g. c(0, 1)
    title           = NULL
) {
  umap_source <- match.arg(umap_source)
  if (is.null(anchors$max_anchor_df)) stop("anchors$max_anchor_df not found -- run max_anchors() first")

  q_meta <- anchors$query$meta
  adf    <- filter_anchor_df(anchors$max_anchor_df, min_score, max_score)
  if (is.logical(query_subset)) query_subset <- q_meta$query_id_surface[query_subset]
  if (is.logical(ref_subset))   ref_subset   <- anchors$ref$meta$ref_id_surface[ref_subset]
  if (!is.null(query_subset)) adf <- adf[adf$query_id_surface %in% query_subset, ]
  if (!is.null(ref_subset))   adf <- adf[adf$ref_id_surface   %in% ref_subset, ]

  if (umap_source == "reference") {
    bg <- anchors$ref$emb[[ref.reduction]]
    if (is.null(bg)) stop(sprintf("anchors$ref$emb$%s not found -- run transfer_reduction() first", ref.reduction))
    colnames(bg)[1:2] <- c("UMAP_1", "UMAP_2")
    bg$ref_id_surface   <- rownames(bg)
    max_scores          <- tapply(adf$score, adf$ref_id_surface, max)
    bg$max_anchor_score <- max_scores[bg$ref_id_surface]
    default_title       <- "Max anchor score on reference UMAP"
  } else {
    bg <- anchors$query$emb[[query.reduction]]
    if (is.null(bg)) stop(sprintf("anchors$query$emb$%s not found -- run transfer_reduction() first", query.reduction))
    colnames(bg)[1:2] <- c("UMAP_1", "UMAP_2")
    bg$query_id_surface <- rownames(bg)
    max_scores          <- tapply(adf$score, adf$query_id_surface, max)
    bg$max_anchor_score <- max_scores[bg$query_id_surface]
    default_title       <- "Max anchor score on query UMAP"
  }
  scored <- bg[!is.na(bg$max_anchor_score), ]

  p <- ggplot() +
    geom_point(data = bg,     aes(UMAP_1, UMAP_2), color = bg.color, size = bg.pt.size) +
    geom_point(data = scored, aes(UMAP_1, UMAP_2, color = max_anchor_score), size = pt.size) +
    scale_color_viridis_c(option = viridis_option, limits = limits) + theme_classic()
  .add_title(p, title, default_title)
}


# ── prediction_plot ───────────────────────────────────────────────────────────
# A query metadata column (e.g. predicted.<var>) on a query embedding.

prediction_plot <- function(
    anchors,
    group.by,
    query_subset  = NULL,
    ref_subset    = NULL,
    reduction     = "umap",
    anchored_only = FALSE,
    pt.size       = 0.3,
    bg.pt.size    = pt.size,
    bg.color      = "lightgrey",
    title         = NULL,
    ...
) {
  if ((anchored_only || !is.null(ref_subset)) && is.null(anchors$max_anchor_df))
    stop("anchors$max_anchor_df not found -- run max_anchors() first")
  q_meta <- anchors$query$meta
  q_emb  <- anchors$query$emb[[reduction]]
  if (is.null(q_emb)) stop(sprintf("anchors$query$emb$%s not found -- run transfer_reduction() first", reduction))
  adf    <- anchors$max_anchor_df

  if (is.logical(query_subset)) query_subset <- q_meta$query_id_surface[query_subset]
  if (is.logical(ref_subset))   ref_subset   <- anchors$ref$meta$ref_id_surface[ref_subset]

  use_ids <- if (anchored_only) unique(adf$query_id_surface) else q_meta$query_id_surface
  if (!is.null(query_subset)) use_ids <- intersect(use_ids, query_subset)
  if (!is.null(ref_subset))   use_ids <- intersect(use_ids, adf$query_id_surface[adf$ref_id_surface %in% ref_subset])

  if (!group.by %in% colnames(q_meta)) stop("'", group.by, "' not found in query metadata")

  bg <- q_emb
  colnames(bg)[1:2] <- c("UMAP_1", "UMAP_2")
  bg$query_id_surface <- rownames(bg)
  bg[[group.by]]      <- q_meta[[group.by]][match(bg$query_id_surface, q_meta$query_id_surface)]
  plot_coords <- bg[bg$query_id_surface %in% use_ids, ]

  p <- ggplot() +
    geom_point(data = bg,          aes(UMAP_1, UMAP_2), color = bg.color, size = bg.pt.size) +
    geom_point(data = plot_coords, aes(UMAP_1, UMAP_2, color = .data[[group.by]]), size = pt.size) +
    theme_classic()
  .add_title(p, title, paste0(ifelse(anchored_only, "Anchored", "All"), " cells -- ", group.by))
}


# ── .transfer_anchor_value ────────────────────────────────────────────────────
# Makes sure value_col is a column of anchors$query$meta: query-side values are
# used as is; reference-side values are pulled onto each query cell through its
# best anchor (max_anchor_df).

.transfer_anchor_value <- function(anchors, value_col) {
  if (value_col %in% colnames(anchors$query$meta)) return(anchors)
  if (value_col %in% colnames(anchors$ref$meta)) {
    ma <- anchors$max_anchor_df
    if (is.null(ma)) stop("anchors$max_anchor_df not found -- run max_anchors() first")
    ref_vals <- anchors$ref$meta[[value_col]][match(ma$ref_id_surface, anchors$ref$meta$ref_id_surface)]
    names(ref_vals) <- ma$query_id_surface
    anchors$query$meta[[value_col]] <- ref_vals[anchors$query$meta$query_id_surface]
    return(anchors)
  }
  stop(sprintf("'%s' not found in query or reference metadata", value_col))
}


# ── anchor_value_plot ─────────────────────────────────────────────────────────
# Any numeric value (query side, or reference side through the best anchor) on
# a query embedding.

anchor_value_plot <- function(
    anchors,
    value_col,
    reduction      = "ref.umap",
    query_subset   = NULL,
    ref_subset     = NULL,
    min_score      = NULL,
    max_score      = NULL,
    pt.size        = 0.3,
    bg.pt.size     = pt.size,
    bg.color       = "lightgrey",
    viridis_option = "D",
    limits         = NULL,
    title          = NULL
) {
  anchors <- .transfer_anchor_value(anchors, value_col)
  if (is.null(anchors$max_anchor_df)) stop("anchors$max_anchor_df not found -- run max_anchors() first")

  q_meta <- anchors$query$meta
  q_emb  <- anchors$query$emb[[reduction]]
  if (is.null(q_emb)) stop(sprintf("anchors$query$emb$%s not found -- run transfer_reduction() first", reduction))
  adf    <- filter_anchor_df(anchors$max_anchor_df, min_score, max_score)

  if (is.logical(query_subset)) query_subset <- q_meta$query_id_surface[query_subset]
  if (is.logical(ref_subset))   ref_subset   <- anchors$ref$meta$ref_id_surface[ref_subset]

  use_ids <- q_meta$query_id_surface[!is.na(q_meta[[value_col]])]
  if (!is.null(query_subset)) use_ids <- intersect(use_ids, query_subset)
  if (!is.null(ref_subset))   use_ids <- intersect(use_ids, adf$query_id_surface[adf$ref_id_surface %in% ref_subset])
  if (!is.null(min_score) || !is.null(max_score)) use_ids <- intersect(use_ids, unique(adf$query_id_surface))

  bg <- q_emb
  colnames(bg)[1:2] <- c("UMAP_1", "UMAP_2")
  bg$query_id_surface <- rownames(bg)
  bg[[value_col]]     <- q_meta[[value_col]][match(bg$query_id_surface, q_meta$query_id_surface)]
  plot_df <- bg[bg$query_id_surface %in% use_ids, ]

  p <- ggplot() +
    geom_point(data = bg,      aes(UMAP_1, UMAP_2), color = bg.color, size = bg.pt.size) +
    geom_point(data = plot_df, aes(UMAP_1, UMAP_2, color = .data[[value_col]]), size = pt.size) +
    scale_color_viridis_c(option = viridis_option, limits = limits, na.value = bg.color) +
    theme_classic() + labs(color = value_col)
  .add_title(p, title, value_col)
}


# ── matrix_plot ───────────────────────────────────────────────────────────────
# Heatmap of average max anchor score: query groups x reference groups.

matrix_plot <- function(
    anchors,
    ref_col      = NULL,
    query_col    = "seurat_clusters",
    query_subset = NULL,
    ref_subset   = NULL,
    transpose    = FALSE,
    arrange_by   = NULL,
    digits       = 2,
    title_label  = NULL,
    text.size    = 2.5
) {
  if (is.null(ref_col)) stop("ref_col must be specified")
  if (is.null(anchors$max_anchor_df)) stop("anchors$max_anchor_df not found -- run max_anchors() first")

  q_meta <- anchors$query$meta
  if (is.logical(query_subset)) query_subset <- q_meta$query_id_surface[query_subset]
  if (is.logical(ref_subset))   ref_subset   <- anchors$ref$meta$ref_id_surface[ref_subset]

  df <- anchors$max_anchor_df
  if (!is.null(query_subset)) df <- df[df$query_id_surface %in% query_subset, ]
  if (!is.null(ref_subset))   df <- df[df$ref_id_surface   %in% ref_subset, ]

  q_col <- paste0("query.", query_col)
  r_col <- paste0("ref.",   ref_col)
  if (!q_col %in% colnames(df)) stop(sprintf("'%s' not found in max_anchor_df -- add it via update_anchors(query_meta_cols)", q_col))
  if (!r_col %in% colnames(df)) stop(sprintf("'%s' not found in max_anchor_df -- add it via update_anchors(reference_meta_cols)", r_col))

  group_props <- df |>
    group_by(.data[[r_col]]) |>
    summarise(n_cells = n_distinct(query_id_surface), .groups = "drop") |>
    mutate(prop = n_cells / sum(n_cells),
           col_label = sprintf("%s\n(%.1f%%)", .data[[r_col]], prop * 100))
  label_map      <- setNames(group_props$col_label, group_props[[r_col]])
  clean_to_label <- label_map
  df[[r_col]]    <- label_map[as.character(df[[r_col]])]

  mat_long <- df |>
    group_by(.data[[q_col]], query_id_surface, .data[[r_col]]) |>
    summarise(max_score = max(score), .groups = "drop") |>
    group_by(.data[[q_col]], .data[[r_col]]) |>
    summarise(avg_max_score = mean(max_score), .groups = "drop")

  if (!is.null(arrange_by)) {
    anchor_label <- clean_to_label[[arrange_by]]
    if (is.na(anchor_label)) stop(sprintf("'%s' not found in %s", arrange_by, ref_col))
    mat_wide <- mat_long |>
      tidyr::pivot_wider(names_from = all_of(r_col), values_from = avg_max_score, values_fill = 0) |>
      tibble::column_to_rownames(q_col)
    col_order <- kl_arrange(as.matrix(mat_wide), anchor_label)
    mat_long[[r_col]] <- factor(mat_long[[r_col]], levels = col_order)
  }

  if (is.null(title_label))
    title_label <- if (!transpose) sprintf("Matrix: %s x %s", query_col, ref_col) else
      sprintf("Matrix: %s x %s", ref_col, query_col)

  p <- if (!transpose)
    ggplot(mat_long, aes(x = .data[[r_col]], y = factor(.data[[q_col]]), fill = avg_max_score)) +
    labs(title = title_label, x = ref_col, y = query_col)
  else
    ggplot(mat_long, aes(x = factor(.data[[q_col]]), y = .data[[r_col]], fill = avg_max_score)) +
    labs(title = title_label, x = query_col, y = ref_col)

  p +
    geom_tile(color = "white") +
    geom_text(aes(label = round(avg_max_score, digits)), size = text.size) +
    scale_fill_viridis_c(name = "Avg Max Anchor Score") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}


# ── prop_plot ─────────────────────────────────────────────────────────────────
# Heatmap: proportion of each query group's anchored cells per reference group.

prop_plot <- function(
    anchors,
    ref_col      = NULL,
    query_col    = "seurat_clusters",
    query_subset = NULL,
    ref_subset   = NULL,
    transpose    = FALSE,
    arrange_by   = NULL,
    digits       = 2,
    title_label  = NULL,
    text.size    = 2.5
) {
  if (is.null(ref_col)) stop("ref_col must be specified")
  if (is.null(anchors$max_anchor_df)) stop("anchors$max_anchor_df not found -- run max_anchors() first")

  q_meta <- anchors$query$meta
  if (is.logical(query_subset)) query_subset <- q_meta$query_id_surface[query_subset]
  if (is.logical(ref_subset))   ref_subset   <- anchors$ref$meta$ref_id_surface[ref_subset]

  df <- anchors$max_anchor_df
  if (!is.null(query_subset)) df <- df[df$query_id_surface %in% query_subset, ]
  if (!is.null(ref_subset))   df <- df[df$ref_id_surface   %in% ref_subset, ]

  q_col <- paste0("query.", query_col)
  r_col <- paste0("ref.",   ref_col)
  if (!q_col %in% colnames(df)) stop(sprintf("'%s' not found in max_anchor_df -- add it via update_anchors(query_meta_cols)", q_col))
  if (!r_col %in% colnames(df)) stop(sprintf("'%s' not found in max_anchor_df -- add it via update_anchors(reference_meta_cols)", r_col))

  group_props <- df |>
    group_by(.data[[r_col]]) |>
    summarise(n_cells = n_distinct(query_id_surface), .groups = "drop") |>
    mutate(prop = n_cells / sum(n_cells),
           col_label = sprintf("%s\n(%.1f%%)", .data[[r_col]], prop * 100))
  label_map      <- setNames(group_props$col_label, group_props[[r_col]])
  clean_to_label <- label_map
  df[[r_col]]    <- label_map[as.character(df[[r_col]])]

  prop_long <- df |>
    group_by(.data[[q_col]], .data[[r_col]]) |>
    summarise(n_cells = n_distinct(query_id_surface), .groups = "drop") |>
    group_by(.data[[q_col]]) |>
    mutate(prop = n_cells / sum(n_cells)) |>
    ungroup() |>
    dplyr::filter(prop > 0)

  if (!is.null(arrange_by)) {
    anchor_label <- clean_to_label[[arrange_by]]
    if (is.na(anchor_label)) stop(sprintf("'%s' not found in %s", arrange_by, ref_col))
    mat_wide <- prop_long |>
      tidyr::pivot_wider(id_cols = all_of(q_col), names_from = all_of(r_col),
                         values_from = prop, values_fill = 0) |>
      tibble::column_to_rownames(q_col)
    col_order <- kl_arrange(as.matrix(mat_wide), anchor_label)
    prop_long[[r_col]] <- factor(prop_long[[r_col]], levels = col_order)
  }

  if (is.null(title_label))
    title_label <- if (!transpose) sprintf("Proportion: %s x %s", query_col, ref_col) else
      sprintf("Proportion: %s x %s", ref_col, query_col)

  p <- if (!transpose)
    ggplot(prop_long, aes(x = .data[[r_col]], y = factor(.data[[q_col]]), fill = prop)) +
    labs(title = title_label, x = ref_col, y = query_col)
  else
    ggplot(prop_long, aes(x = factor(.data[[q_col]]), y = .data[[r_col]], fill = prop)) +
    labs(title = title_label, x = query_col, y = ref_col)

  p +
    geom_tile(color = "white") +
    geom_text(aes(label = round(prop, digits)), size = text.size) +
    scale_fill_viridis_c(name = "Proportion", limits = c(0, 1), na.value = "white") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}


# ── scatter_plot ──────────────────────────────────────────────────────────────
# For each (group.by, predicted plot.by) pair: x = proportion of the group's
# cells with that predicted label, y = average max anchor score
# (y_metric = "anchor_score") or average prediction score for plot.by
# (y_metric = "prediction_score").
# group.by: query metadata column (colour). plot.by: variable name given to
# update_anchors() (reads predicted.<plot.by>).

scatter_plot <- function(
    anchors,
    group.by,
    plot.by,
    query_subset  = NULL,
    ref_subset    = NULL,
    min_score     = NULL,
    max_score     = NULL,
    y_metric      = c("anchor_score", "prediction_score"),
    anchored_only = TRUE,
    label         = TRUE,
    pt.size       = 3,
    label.size    = 3,
    title         = NULL
) {
  y_metric <- match.arg(y_metric)
  q_meta <- anchors$query$meta
  if (is.logical(query_subset)) query_subset <- q_meta$query_id_surface[query_subset]
  if (is.logical(ref_subset))   ref_subset   <- anchors$ref$meta$ref_id_surface[ref_subset]
  if (!group.by %in% colnames(q_meta)) stop(sprintf("'%s' not found in query metadata", group.by))

  p_col <- paste0("predicted.", plot.by)
  if (!p_col %in% colnames(q_meta))
    stop(sprintf("'%s' not found in query metadata -- plot.by is the variable name given to update_anchors() (e.g. 'CellClass')", plot.by))
  g_col <- group.by

  if (y_metric == "anchor_score") {
    if (is.null(anchors$max_anchor_df)) stop("anchors$max_anchor_df not found -- run max_anchors() first")
    mdf <- filter_anchor_df(anchors$max_anchor_df, min_score, max_score)
    if (!is.null(query_subset)) mdf <- mdf[mdf$query_id_surface %in% query_subset, ]
    if (!is.null(ref_subset))   mdf <- mdf[mdf$ref_id_surface   %in% ref_subset, ]
    mdf[[g_col]] <- q_meta[[g_col]][match(mdf$query_id_surface, q_meta$query_id_surface)]
    mdf[[p_col]] <- q_meta[[p_col]][match(mdf$query_id_surface, q_meta$query_id_surface)]
    mdf$y_value <- mdf$score
    y_label <- "Average max anchor score"
    x_label <- sprintf("Proportion of %s's anchored cells", group.by)
  } else {
    score_col <- paste0(p_col, ".score")
    if (!score_col %in% colnames(q_meta))
      stop(sprintf("'%s' not found in query metadata -- run update_anchors() with '%s' in variables", score_col, plot.by))
    if (!is.null(ref_subset) && !anchored_only)
      warning("ref_subset is ignored when y_metric = 'prediction_score' and anchored_only = FALSE")
    mdf <- q_meta
    if (!is.null(query_subset)) mdf <- mdf[mdf$query_id_surface %in% query_subset, ]
    if (anchored_only) {
      if (is.null(anchors$max_anchor_df)) stop("anchors$max_anchor_df not found -- run max_anchors() first")
      anchored_ids <- unique(anchors$max_anchor_df$query_id_surface)
      if (!is.null(ref_subset))
        anchored_ids <- intersect(anchored_ids,
                                  anchors$max_anchor_df$query_id_surface[anchors$max_anchor_df$ref_id_surface %in% ref_subset])
      mdf <- mdf[mdf$query_id_surface %in% anchored_ids, ]
    }
    if (!is.null(min_score)) mdf <- mdf[mdf[[score_col]] >= min_score, ]
    if (!is.null(max_score)) mdf <- mdf[mdf[[score_col]] <= max_score, ]
    mdf$y_value <- mdf[[score_col]]
    y_label <- sprintf("Prediction score (%s)", plot.by)
    x_label <- sprintf("Proportion of %s's %scells", group.by, ifelse(anchored_only, "anchored ", ""))
  }

  group_totals <- mdf |>
    group_by(.data[[g_col]]) |>
    summarise(group_total = n_distinct(query_id_surface), .groups = "drop")
  plot_df <- mdf |>
    group_by(.data[[g_col]], .data[[p_col]]) |>
    summarise(avg_y_value = mean(y_value), n_cells = n_distinct(query_id_surface), .groups = "drop") |>
    left_join(group_totals, by = g_col) |>
    mutate(prop = n_cells / group_total)
  legend_map <- setNames(sprintf("%s (%d)", group_totals[[g_col]], group_totals$group_total), group_totals[[g_col]])
  plot_df[[g_col]] <- legend_map[as.character(plot_df[[g_col]])]

  p <- ggplot(plot_df, aes(x = prop, y = avg_y_value, color = .data[[g_col]])) +
    geom_point(size = pt.size) +
    labs(y = y_label, x = x_label, color = group.by) +
    xlim(0, 1) + ylim(0, 1) + theme_classic()
  if (label) p <- p + geom_text(aes(label = .data[[p_col]]), vjust = -0.8, size = label.size, show.legend = FALSE)
  .add_title(p, title, sprintf("%s x %s -- proportion vs %s", group.by, plot.by, tolower(y_label)))
}


# ── stacked_count_plot ────────────────────────────────────────────────────────
# Stacked bar of predicted.<plot.by> composition per query group.by level.

stacked_count_plot <- function(
    anchors,
    group.by,
    plot.by,
    query_subset = NULL,
    proportion   = TRUE,
    title_label  = NULL
) {
  q_meta <- anchors$query$meta
  if (is.logical(query_subset)) query_subset <- q_meta$query_id_surface[query_subset]
  if (!is.null(query_subset)) q_meta <- q_meta[q_meta$query_id_surface %in% query_subset, ]
  if (!group.by %in% colnames(q_meta)) stop(sprintf("'%s' not found in query metadata", group.by))
  p_col <- paste0("predicted.", plot.by)
  if (!p_col %in% colnames(q_meta))
    stop(sprintf("'%s' not found in query metadata -- plot.by is the variable name given to update_anchors()", plot.by))
  g_col <- group.by

  plot_df <- q_meta |>
    dplyr::filter(!is.na(.data[[p_col]])) |>
    dplyr::count(.data[[g_col]], .data[[p_col]], name = "n") |>
    dplyr::group_by(.data[[g_col]]) |>
    dplyr::mutate(prop = n / sum(n)) |>
    dplyr::ungroup()

  if (is.null(title_label)) title_label <- sprintf("%s composition by %s", plot.by, group.by)
  y_var   <- if (proportion) "prop" else "n"
  y_label <- if (proportion) "Proportion" else "Cell count"

  ggplot(plot_df, aes(x = as.character(.data[[g_col]]), y = .data[[y_var]], fill = .data[[p_col]])) +
    geom_bar(stat = "identity", position = "stack") +
    labs(title = title_label, x = group.by, y = y_label, fill = plot.by) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}


# ── stacked_count_plot_seurat ─────────────────────────────────────────────────
# Same as stacked_count_plot, read straight from a Seurat object's metadata.

stacked_count_plot_seurat <- function(
    object,
    group.by,
    plot.by,
    cell_subset  = NULL,
    proportion   = TRUE,
    title_label  = NULL
) {
  if (!inherits(object, "Seurat")) stop("'object' must be a Seurat object")
  meta <- object@meta.data
  if (!is.null(cell_subset)) {
    meta <- if (is.logical(cell_subset)) meta[cell_subset, , drop = FALSE] else
      meta[rownames(meta) %in% cell_subset, , drop = FALSE]
  }
  if (!group.by %in% colnames(meta)) stop(sprintf("'%s' not found in Seurat object metadata", group.by))
  if (!plot.by %in% colnames(meta))  stop(sprintf("'%s' not found in Seurat object metadata", plot.by))
  g_col <- group.by
  p_col <- plot.by

  plot_df <- meta |>
    dplyr::filter(!is.na(.data[[p_col]])) |>
    dplyr::count(.data[[g_col]], .data[[p_col]], name = "n") |>
    dplyr::group_by(.data[[g_col]]) |>
    dplyr::mutate(prop = n / sum(n)) |>
    dplyr::ungroup()

  if (is.null(title_label)) title_label <- sprintf("%s composition by %s", plot.by, group.by)
  y_var   <- if (proportion) "prop" else "n"
  y_label <- if (proportion) "Proportion" else "Cell count"

  ggplot(plot_df, aes(x = as.character(.data[[g_col]]), y = .data[[y_var]], fill = .data[[p_col]])) +
    geom_bar(stat = "identity", position = "stack") +
    labs(title = title_label, x = group.by, y = y_label, fill = plot.by) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}


# ── sandkey_plot_seurat ───────────────────────────────────────────────────────
# Sankey of how ref_group.by categories map onto query_group.by categories for
# cells matched across two Seurat objects by a shared id column (id_col).
# Ids duplicated within one object are dropped first (ignore_duplication).
# equal_size = TRUE gives every category equal visual weight.

sandkey_plot_seurat <- function(
    query,
    ref,
    query_subset_cell = NULL,
    ref_subset_cell   = NULL,
    query_group.by,
    ref_group.by,
    id_col            = "obj_id",
    title_label       = NULL,
    equal_size        = FALSE,
    text.size         = 3
) {
  .need("ggalluvial", "sandkey_plot_seurat()")
  if (!inherits(query, "Seurat") || !inherits(ref, "Seurat"))
    stop("'query' and 'ref' must both be Seurat objects")
  q_meta <- query@meta.data
  r_meta <- ref@meta.data
  if (!id_col %in% colnames(q_meta) || !id_col %in% colnames(r_meta))
    stop(sprintf("both objects need an '%s' metadata column holding the shared cell id", id_col))
  if (!query_group.by %in% colnames(q_meta)) stop(sprintf("'%s' not found in query metadata", query_group.by))
  if (!ref_group.by %in% colnames(r_meta))   stop(sprintf("'%s' not found in ref metadata", ref_group.by))

  if (!is.null(query_subset_cell))
    q_meta <- if (is.logical(query_subset_cell)) q_meta[query_subset_cell, , drop = FALSE] else
      q_meta[rownames(q_meta) %in% query_subset_cell, , drop = FALSE]
  if (!is.null(ref_subset_cell))
    r_meta <- if (is.logical(ref_subset_cell)) r_meta[ref_subset_cell, , drop = FALSE] else
      r_meta[rownames(r_meta) %in% ref_subset_cell, , drop = FALSE]

  q_df <- data.frame(obj_id = q_meta[[id_col]], query_label = q_meta[[query_group.by]], stringsAsFactors = FALSE)
  r_df <- data.frame(obj_id = r_meta[[id_col]], ref_label   = r_meta[[ref_group.by]],   stringsAsFactors = FALSE)
  q_df <- q_df[!is.na(q_df$obj_id), ]
  r_df <- r_df[!is.na(r_df$obj_id), ]
  q_df <- ignore_duplication(q_df, id_col = "obj_id", label = "query")$df
  r_df <- ignore_duplication(r_df, id_col = "obj_id", label = "ref")$df

  merged <- dplyr::inner_join(r_df, q_df, by = "obj_id")
  if (nrow(merged) == 0) stop(sprintf("no shared '%s' values between query and ref", id_col))
  merged <- merged[!is.na(merged$ref_label) & !is.na(merged$query_label), ]
  flow_df <- dplyr::count(merged, ref_label, query_label, name = "n")

  y_label <- "Cell count"
  if (equal_size) {
    flow_df <- flow_df |>
      dplyr::group_by(ref_label) |> dplyr::mutate(w_ref = n / sum(n)) |> dplyr::ungroup() |>
      dplyr::group_by(query_label) |> dplyr::mutate(w_query = n / sum(n)) |> dplyr::ungroup() |>
      dplyr::mutate(n = (w_ref + w_query) / 2) |>
      dplyr::select(-w_ref, -w_query)
    y_label <- "Relative flow (equal group size)"
  }
  if (is.null(title_label))
    title_label <- sprintf("%s -> %s (n=%d matched cells)", ref_group.by, query_group.by, nrow(merged))

  ggplot(flow_df, aes(axis1 = ref_label, axis2 = query_label, y = n)) +
    ggalluvial::geom_alluvium(aes(fill = ref_label), width = 1/8) +
    ggalluvial::geom_stratum(width = 1/8, fill = "grey90", color = "grey30") +
    geom_text(stat = ggalluvial::StatStratum, aes(label = after_stat(stratum)), size = text.size) +
    scale_x_discrete(limits = c(ref_group.by, query_group.by), expand = c(.05, .05)) +
    labs(title = title_label, y = y_label, fill = ref_group.by) +
    theme_classic()
}


# ── compute_pseudobulk ────────────────────────────────────────────────────────
# Pseudobulk counts per group_col level (summed in batches straight from the
# sparse counts), DESeq2 size-factor normalized, log1p. Returns genes x groups.

compute_pseudobulk <- function(B, group_col = "sample_id", batch_size = 30, assay = NULL, layer = "counts") {
  .need("DESeq2", "compute_pseudobulk()")
  if (!group_col %in% colnames(B@meta.data)) stop(sprintf("group_col '%s' not found in metadata", group_col))
  assay <- assay %||% DefaultAssay(B)
  if (length(.layers(B, assay, layer)) > 1) B <- .join_layers(B, assay)

  counts_full <- SeuratObject::LayerData(B, assay = assay, layer = layer)
  sample_vec  <- as.character(B@meta.data[[group_col]])
  names(sample_vec) <- colnames(B)
  sample_ids_all <- unique(sample_vec[!is.na(sample_vec)])
  batches <- split(sample_ids_all, ceiling(seq_along(sample_ids_all) / batch_size))

  pseudobulk_list <- vector("list", length(batches))
  for (b in seq_along(batches)) {
    batch_ids <- batches[[b]]
    col_idx   <- which(sample_vec %in% batch_ids)
    sub_mat   <- counts_full[, col_idx, drop = FALSE]
    idx_by_sample <- split(seq_along(col_idx), sample_vec[col_idx])
    pseudobulk_list[[b]] <- vapply(batch_ids, function(sid) Matrix::rowSums(sub_mat[, idx_by_sample[[sid]], drop = FALSE]),
                                   numeric(nrow(sub_mat)))
  }
  raw_counts <- do.call("cbind", pseudobulk_list)
  rownames(raw_counts) <- rownames(counts_full)
  raw_counts <- raw_counts[, sample_ids_all, drop = FALSE]
  rm(counts_full, pseudobulk_list); gc()

  size_factors <- DESeq2::estimateSizeFactorsForMatrix(as.matrix(raw_counts))
  log1p(sweep(as.matrix(raw_counts), 2, size_factors, "/"))
}


# ── project_pseudobulk_reductions ─────────────────────────────────────────────
# Projects a pseudobulk profile (genes x groups) onto the reference's EXISTING
# PCA (loadings reused, data scaled with the reference's per-gene mean / sd);
# optional harmony coordinates from a stored harmony model (harmony_dims).
# Returns list(pca = <groups x pca_dims>, harmony = <groups x harmony_dims> or NULL).

project_pseudobulk_reductions <- function(
    expr_profile,
    B,
    pca_dims     = c(1, 2),
    harmony_dims = NULL,
    pca_name     = "pca",
    harmony_name = "harmony",
    assay        = NULL
) {
  if (!pca_name %in% Reductions(B)) stop(sprintf("reduction '%s' not found", pca_name))
  assay <- assay %||% DefaultAssay(B[[pca_name]])
  pca_loadings_full <- Loadings(B, reduction = pca_name)
  common_genes      <- intersect(rownames(pca_loadings_full), rownames(expr_profile))
  if (!length(common_genes)) stop("no genes shared between expr_profile and the PCA loadings")

  ref_data <- SeuratObject::LayerData(B, assay = assay, layer = "data")[common_genes, , drop = FALSE]
  mean_vec <- Matrix::rowMeans(ref_data)
  mean_sq  <- Matrix::rowMeans(ref_data^2)
  n_cells  <- ncol(ref_data)
  sd_vec   <- sqrt((mean_sq - mean_vec^2) * n_cells / (n_cells - 1))
  sd_vec[sd_vec == 0 | is.na(sd_vec)] <- 1e-8

  .project_onto_pca <- function(dims) {
    loadings_d <- pca_loadings_full[common_genes, dims, drop = FALSE]
    x_scaled   <- (expr_profile[common_genes, , drop = FALSE] - mean_vec) / sd_vec
    t(t(loadings_d) %*% x_scaled)
  }
  pc_coords <- .project_onto_pca(pca_dims)
  colnames(pc_coords) <- paste0("PC_", pca_dims)
  if (is.null(harmony_dims)) return(list(pca = pc_coords, harmony = NULL))

  harmony_model <- B@reductions[[harmony_name]]@misc$harmony_model
  if (is.null(harmony_model)) stop(sprintf("harmony_dims was set but %s@misc$harmony_model not found", harmony_name))
  Y     <- harmony_model$Y
  sigma <- as.vector(harmony_model$sigma)
  pca_for_harmony <- .project_onto_pca(seq_len(nrow(Y)))

  dist_mat <- matrix(0, nrow = ncol(Y), ncol = nrow(pca_for_harmony))
  for (k in seq_len(ncol(Y))) dist_mat[k, ] <- rowSums(sweep(pca_for_harmony, 2, Y[, k])^2)
  R_sample <- exp(-dist_mat / sigma)
  R_sample <- sweep(R_sample, 2, colSums(R_sample), "/")

  Z_corrected <- pca_for_harmony
  for (n in seq_len(nrow(pca_for_harmony)))
    Z_corrected[n, ] <- pca_for_harmony[n, ] - as.vector(Y %*% R_sample[, n])
  harmony_coords <- Z_corrected[, harmony_dims, drop = FALSE]
  colnames(harmony_coords) <- paste0("harmony_", harmony_dims)
  list(pca = pc_coords, harmony = harmony_coords)
}


# ── transfer_matrix ───────────────────────────────────────────────────────────
# Stores any matrix (e.g. pseudobulk expression, genes x groups) on anchors[[side]].

transfer_matrix <- function(
    anchors,
    mat,
    side      = c("ref", "query"),
    slot_name = "pseudobulk_expr",
    group_col = "sample_id"
) {
  side <- match.arg(side)
  if (!is.list(anchors[[side]])) stop(sprintf("anchors$%s is not a list", side))
  anchors[[side]][[slot_name]] <- mat
  attr(anchors[[side]][[slot_name]], "group_col") <- group_col
  anchors
}


# ── summary_write ─────────────────────────────────────────────────────────────
# Text summary of anchored cells per query group and per reference group.

summary_write <- function(
    anchors,
    query_col    = "seurat_clusters",
    ref_col      = NULL,
    query_subset = NULL,
    ref_subset   = NULL,
    filename     = NULL,
    out_dir      = "."
) {
  if (is.null(query_col)) stop("query_col must be specified")
  if (is.null(anchors$max_anchor_df)) stop("anchors$max_anchor_df not found -- run max_anchors() first")

  q_meta <- anchors$query$meta
  if (is.logical(query_subset)) query_subset <- q_meta$query_id_surface[query_subset]
  if (is.logical(ref_subset))   ref_subset   <- anchors$ref$meta$ref_id_surface[ref_subset]
  df <- anchors$max_anchor_df
  if (!is.null(query_subset)) df <- df[df$query_id_surface %in% query_subset, ]
  if (!is.null(ref_subset))   df <- df[df$ref_id_surface   %in% ref_subset, ]

  filename <- filename %||% sprintf("summary_%s_%s.txt", query_col, ref_col %||% "ref")
  q_col <- paste0("query.", query_col)
  r_col <- if (!is.null(ref_col)) paste0("ref.", ref_col) else NULL
  if (!q_col %in% colnames(df)) stop(sprintf("'%s' not found in max_anchor_df", q_col))
  if (!is.null(r_col) && !r_col %in% colnames(df)) stop(sprintf("'%s' not found in max_anchor_df", r_col))

  cluster_counts <- table(q_meta[[query_col]])
  anchor_summary <- df |>
    group_by(.data[[q_col]]) |>
    summarise(n_anchored = n_distinct(query_id_surface), .groups = "drop") |>
    mutate(total_cells   = as.integer(cluster_counts[as.character(.data[[q_col]])]),
           prop_in_query = round(total_cells / sum(as.integer(cluster_counts)) * 100, 2),
           prop_anchored = round(n_anchored / total_cells * 100, 2))
  ord_num <- suppressWarnings(as.numeric(as.character(anchor_summary[[q_col]])))
  anchor_summary <- anchor_summary[if (!anyNA(ord_num)) order(ord_num) else order(as.character(anchor_summary[[q_col]])), ]

  group_col   <- r_col %||% q_col
  group_label <- ref_col %||% query_col
  summary_df <- df |>
    group_by(.data[[q_col]], query_id_surface, .data[[group_col]]) |>
    summarise(max_score = max(score), .groups = "drop") |>
    group_by(.data[[q_col]], .data[[group_col]]) |>
    summarise(n_uniq_cells = n(), avg_max_score = mean(max_score), .groups = "drop") |>
    group_by(.data[[q_col]]) |>
    mutate(prop = n_uniq_cells / sum(n_uniq_cells)) |>
    ungroup() |>
    arrange(.data[[q_col]], desc(avg_max_score))

  out <- c(sprintf("=== Cluster Overview: %s in query ===", query_col),
           sprintf("Generated: %s", format(Sys.time(), "%Y%m%d_%H%M%S")), "",
           sprintf("%-20s %12s %12s %16s %16s", q_col, "total_cells", "n_anchored", "prop_in_query(%)", "prop_anchored(%)"),
           sprintf("%-20s %12d %12d %16.2f %16.2f", as.character(anchor_summary[[q_col]]), anchor_summary$total_cells,
                   anchor_summary$n_anchored, anchor_summary$prop_in_query, anchor_summary$prop_anchored),
           "", sprintf("=== Detailed Summary: %s x %s ===", query_col, group_label), "",
           sprintf("%-28s %-28s %14s %14s %10s", q_col, group_label, "n_uniq_cells", "avg_max_score", "prop"),
           sprintf("%-28s %-28s %14d %14.4f %10.4f",
                   sprintf("%s (n=%d)", as.character(summary_df[[q_col]]),
                           as.integer(cluster_counts[as.character(summary_df[[q_col]])])),
                   as.character(summary_df[[group_col]]), summary_df$n_uniq_cells,
                   summary_df$avg_max_score, summary_df$prop))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  path <- file.path(out_dir, filename)
  writeLines(out, path)
  invisible(path)
}


# ── save_upload_ec2 ───────────────────────────────────────────────────────────
# Saves a ggplot (ggsave, white background) and, if rclone_remote and
# drive_dir are both given, copies the file there with rclone.

save_upload_ec2 <- function(
    plot_obj      = NULL,
    filename      = NULL,
    width         = 7,
    height        = 6,
    dpi           = 300,
    out_dir       = ".",
    drive_dir     = NULL,
    rclone_remote = NULL,
    rclone_args   = "--progress",
    remove_tmp    = FALSE
) {
  if (is.null(filename)) stop("filename must be specified")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  path <- file.path(out_dir, filename)
  if (!is.null(plot_obj)) ggsave(path, plot = plot_obj, width = width, height = height, dpi = dpi, bg = "white")
  if (!is.null(rclone_remote) && !is.null(drive_dir)) {
    remote_path <- sprintf("%s:%s", rclone_remote, drive_dir)
    system(sprintf("rclone copy %s %s %s", shQuote(path), shQuote(remote_path), rclone_args))
    message(sprintf("Uploaded: %s -> %s", filename, remote_path))
    if (remove_tmp) file.remove(path)
  }
  invisible(path)
}


# ── max_anchor_bar_plot ───────────────────────────────────────────────────────
# Per query group.by level: % (or count) of cells anchored.
# metric     "percent" | "count"
# percent_of "group_total" (% of the group's own cells) | "all_anchored" (share of all anchored cells)
# color_by   "avg_score" (mean max anchor score), group.by itself (legend shows
#            n anchored), or any other query meta column (stacked composition).

max_anchor_bar_plot <- function(
    anchors,
    group.by,
    metric       = c("percent", "count"),
    color_by     = "avg_score",
    percent_of   = c("group_total", "all_anchored"),
    query_subset = NULL,
    ref_subset   = NULL,
    bar.width    = 0.6,
    title        = NULL
) {
  metric     <- match.arg(metric)
  percent_of <- match.arg(percent_of)
  if (is.null(anchors$max_anchor_df)) stop("anchors$max_anchor_df not found -- run max_anchors() first")
  q_meta <- anchors$query$meta
  if (!group.by %in% colnames(q_meta)) stop(sprintf("'%s' not found in anchors$query$meta", group.by))
  max_df <- anchors$max_anchor_df

  if (is.logical(query_subset)) query_subset <- q_meta$query_id_surface[query_subset]
  if (is.logical(ref_subset))   ref_subset   <- anchors$ref$meta$ref_id_surface[ref_subset]
  q_meta_use <- q_meta
  if (!is.null(query_subset)) q_meta_use <- q_meta_use[q_meta_use$query_id_surface %in% query_subset, ]
  max_df_use <- max_df
  if (!is.null(query_subset)) max_df_use <- max_df_use[max_df_use$query_id_surface %in% query_subset, ]
  if (!is.null(ref_subset))   max_df_use <- max_df_use[max_df_use$ref_id_surface %in% ref_subset, ]

  anchored_group  <- q_meta_use[[group.by]][match(max_df_use$query_id_surface, q_meta_use$query_id_surface)]
  total_counts    <- table(q_meta_use[[group.by]])
  anchored_counts <- table(anchored_group)
  groups     <- names(total_counts)
  n_anchored <- as.integer(anchored_counts[groups]); n_anchored[is.na(n_anchored)] <- 0
  n_total    <- as.integer(total_counts[groups])
  avg_score  <- sapply(groups, function(g) mean(max_df_use$score[anchored_group == g], na.rm = TRUE))
  percent_vals <- if (percent_of == "group_total") 100 * n_anchored / n_total else 100 * n_anchored / sum(n_anchored)

  plot_df <- data.frame(group = groups, n_anchored = n_anchored, n_total = n_total,
                        percent = percent_vals, avg_score = avg_score)
  y_col   <- if (metric == "percent") "percent" else "n_anchored"
  y_label <- if (metric == "percent") {
    if (percent_of == "group_total") "% of cells anchored" else "% of all anchored cells"
  } else "# of anchors"
  plot_df$group_label <- sprintf("%s\n(%.1f%%)", plot_df$group, plot_df$percent)

  if (color_by != "avg_score" && !identical(color_by, group.by)) {
    if (!color_by %in% colnames(q_meta_use)) stop(sprintf("'%s' not found in anchors$query$meta", color_by))
    anchored_color <- q_meta_use[[color_by]][match(max_df_use$query_id_surface, q_meta_use$query_id_surface)]
    n_total_named  <- setNames(n_total, groups)
    stack_df <- as.data.frame(table(group = anchored_group, color_val = anchored_color), stringsAsFactors = FALSE)
    names(stack_df) <- c("group", "color_val", "n")
    stack_df <- stack_df[stack_df$group %in% groups, ]
    stack_df$percent <- if (percent_of == "group_total") 100 * stack_df$n / n_total_named[stack_df$group] else
      100 * stack_df$n / sum(n_anchored)
    stack_df$group_label <- factor(plot_df$group_label[match(stack_df$group, plot_df$group)], levels = plot_df$group_label)
    y_stack_col <- if (metric == "percent") "percent" else "n"
    p <- ggplot(stack_df, aes(x = group_label, y = .data[[y_stack_col]], fill = color_val)) +
      geom_col(width = bar.width) +
      scale_x_discrete(expand = expansion(add = 0.8)) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
      theme_minimal() +
      labs(x = group.by, y = y_label, fill = color_by) +
      theme(axis.text.x = element_text(angle = 0, hjust = 0.5))
    return(.add_title(p, title, sprintf("%s by %s, split by %s", y_label, group.by, color_by)))
  }

  p <- ggplot(plot_df, aes(x = group_label, y = .data[[y_col]]))
  p <- if (color_by == "avg_score") {
    p + geom_col(aes(fill = avg_score), width = bar.width) +
      scale_fill_viridis_c(name = "avg max\nanchor score", limits = c(0, 1))
  } else {
    p + geom_col(aes(fill = group), width = bar.width) + labs(fill = color_by) +
      scale_fill_discrete(labels = setNames(sprintf("%s (%d)", plot_df$group, plot_df$n_anchored), plot_df$group))
  }
  p <- p +
    scale_x_discrete(expand = expansion(add = 0.8)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
    theme_minimal() +
    labs(x = group.by, y = y_label) +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5))
  .add_title(p, title, sprintf("%s by %s", y_label, group.by))
}
