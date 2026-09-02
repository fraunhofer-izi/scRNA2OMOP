#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# split_all_cell_by_orig_ident.R
# Writes Seurat subsets as anndata 0.1.0 .h5ad with X as CSC (cells x features),
# so anndata.read_h5ad reads the matrix with the correct orientation.

# ---- setup: libraries + paths --------------------------------------------
library(Seurat)
library(rhdf5)
library(digest)

DATA_ROOT    <- Sys.getenv("VTI_DATA_ROOT", "./data")
ALL_CELL_RDS <- file.path(DATA_ROOT, "meta/integration/05_seurat_harmony_all_new.Rds")
OUT_DIR      <- file.path(DATA_ROOT, "h5ad/splits")
MANIFEST_OUT <- file.path(OUT_DIR, "_split_manifest.csv")
TIMEPOINT_ORDER <- c("LP", "Late", "Very Late")  # rank for stable sort

# ---- setup: small helpers -------------------------------------------------
stop_msg <- function(...) stop(paste0(...), call. = FALSE)

# sha256 of a file via digest::digest (namespace-qualified so it works even if
# library(digest) failed to attach -- digest 0.6.35 throws a "cannot be
# unloaded" warning on attach because progressr/future/htmltools import it).
# A hashing failure stops the export: an un-checked hash on an audit sidecar
# is worse than no hash, because downstream consumers trust it silently.
sha256 <- function(p) {
  if (!file.exists(p)) stop_msg("sha256: file not found: ", p)
  h <- digest::digest(p, algo = "sha256", file = TRUE)
  stopifnot("sha256 returned NA" = !is.na(h) && nzchar(h))
  h
}

# Wrapper for rhdf5::h5writeAttribute. The default call writes character scalars
# as length-1 fixed-byte arrays (e.g. array([b'anndata'], dtype='|S8')); anndata
# 0.12's get_spec treats that as IOSpec(encoding_type=['anndata', ...]) -- a
# Python list, which is unhashable and raises TypeError at registry lookup.
# anndata writes these as variable-length UTF-8 string scalars (`str`) for
# length-1 string values and `object` arrays of `str` for length>1. We mirror
# that: pass variableLengthString=TRUE for character input, asScalar=TRUE for
# length-1 input (leaving length>1 string arrays as object arrays, matching
# anndata's column-order encoding). Numeric input is left to rhdf5 defaults.
h5_attr <- function(value, name, handle) {
  if (is.character(value)) {
    if (length(value) == 1L) {
      rhdf5::h5writeAttribute(value, handle, name,
                              variableLengthString = TRUE, asScalar = TRUE)
    } else {
      rhdf5::h5writeAttribute(value, handle, name,
                              variableLengthString = TRUE)
    }
  } else {
    rhdf5::h5writeAttribute(value, handle, name)
  }
}

h5_dset <- function(value, name, handle) {
  rhdf5::h5writeDataset(value, handle, name)
}

# ---- setup: anndata 0.1.0 schema writer -----------------------------------
# Writes a Seurat subset to a valid anndata 0.1.0 .h5ad. X is CSC of shape
# (cells, features): anndata reads it as AnnData with X.shape == (n_obs, n_var)
# and obs/var indexes matching the on-disk groups. A dgCMatrix from Seurat is
# (features x cells); transpose it before dumping @x/@i/@p so the CSC buffers
# describe (cells x features), matching the declared shape and encoding.
#
# NOTE on column-order for empty dataframes: rhdf5 cannot write a zero-length
# string attribute (max(nchar(character(0))) -> -Inf -> H5Tset_size -> NA ->
# H5Awrite "Bad value"). anndata accepts a length-1 empty string "" for an
# empty column-order, so we use "" instead of character(0).
write_h5ad <- function(seu, file) {
  da <- "RNA"
  stopifnot(da %in% Assays(seu))
  if (inherits(seu[[da]], "Assay5") && length(Layers(seu[[da]])) > 1) seu <- JoinLayers(seu)
  x <- GetAssayData(seu, assay = da, layer = "counts")   # dgCMatrix (features x cells)
  xt <- Matrix::t(x)                                     # dgCMatrix (cells x features)
  cells <- colnames(x); feats <- rownames(x)
  stopifnot(identical(cells, rownames(seu@meta.data)))
  meta <- seu@meta.data[cells, , drop = FALSE]
  for (c in names(meta)) meta[[c]] <- as.character(meta[[c]])   # factors -> char
  if (file.exists(file)) file.remove(file)
  h5 <- H5Fcreate(file, flags = "H5F_ACC_TRUNC")
  h5_attr("anndata", "encoding-type",    h5)
  h5_attr("0.1.0",   "encoding-version", h5)
  # X group (CSC, cells x features)
  h5createGroup(h5, "X")
  h5_dset(as.numeric(xt@x),  "X/data",    h5)
  h5_dset(as.integer(xt@i), "X/indices", h5)
  h5_dset(as.integer(xt@p), "X/indptr",  h5)
  g <- H5Gopen(h5, "X")
  h5_attr("csc_matrix",        "encoding-type",    g)
  h5_attr("0.1.0",              "encoding-version", g)
  h5_attr(c(ncol(x), nrow(x)), "shape",            g)
  H5Gclose(g)
  # obs group (dataframe)
  h5createGroup(h5, "obs")
  h5_dset(cells, "obs/_index", h5)
  for (c in names(meta)) h5_dset(meta[[c]], file.path("obs", c), h5)
  g <- H5Gopen(h5, "obs")
  h5_attr("_index",    "_index",          g)
  h5_attr("dataframe", "encoding-type",    g)
  h5_attr("0.2.0",     "encoding-version", g)
  h5_attr(names(meta), "column-order",     g)
  H5Gclose(g)
  # var group (dataframe, no extra columns). Use "" not character(0): rhdf5
  # cannot write a zero-length string attribute (see NOTE above).
  h5createGroup(h5, "var")
  h5_dset(feats, "var/_index", h5)
  g <- H5Gopen(h5, "var")
  h5_attr("_index",     "_index",          g)
  h5_attr("dataframe",  "encoding-type",    g)
  h5_attr("0.2.0",      "encoding-version", g)
  h5_attr("",           "column-order",     g)
  H5Gclose(g)
  # Empty dict groups (anndata expects these)
  for (gn in c("layers","obsm","obsp","varm","varp","uns")) {
    h5createGroup(h5, gn)
    g <- H5Gopen(h5, gn)
    h5_attr("dict",  "encoding-type",    g)
    h5_attr("0.1.0", "encoding-version", g)
    H5Gclose(g)
  }
  H5Fclose(h5)
}

# ---- setup: load all-cell object -----------------------------------------
# Run this block ONCE. seu stays in the global env; do NOT re-run on re-exports.
stopifnot(file.exists(ALL_CELL_RDS))
seu <- readRDS(ALL_CELL_RDS)
cat("object:", ncol(seu), "cells x", nrow(seu), "features; DefaultAssay =", DefaultAssay(seu),
    ";", length(unique(seu$orig.ident)), "orig.ident\n")

# ---- setup: derive specimen map from object metadata --------------------
# Build orig.ident -> (PATIENT_ID, TIMEPOINT) directly from seu@meta.data:
# the all-cell object is the only authoritative source for which specimens
# exist. Run ONCE after loading seu; sm stays in the global env.
needed <- c("orig.ident", "PATIENT_ID", "TIMEPOINT")
miss <- setdiff(needed, colnames(seu@meta.data))
if (length(miss)) stop_msg("seu@meta.data missing cols: ", paste(miss, collapse = ", "))
sm <- unique(seu@meta.data[, needed, drop = FALSE])
sm$orig.ident <- as.character(sm$orig.ident)
sm$PATIENT_ID <- as.character(sm$PATIENT_ID)
sm$TIMEPOINT  <- as.character(sm$TIMEPOINT)
stopifnot(!anyDuplicated(sm$orig.ident))
dup <- sm$orig.ident[ave(sm$PATIENT_ID, sm$orig.ident, FUN = function(z) length(unique(z))) > 1]
if (length(dup)) stop_msg("orig.ident maps to >1 PATIENT_ID: ", paste(dup, collapse = ", "))
if (any(is.na(sm$PATIENT_ID) | sm$PATIENT_ID == "")) stop_msg("missing PATIENT_ID in specimen map")
cat("specimen map:", nrow(sm), "rows,", length(unique(sm$orig.ident)), "orig.ident\n")

# ---- setup: join RNA layers ONCE on the full object ----------------------
# The saved object has the RNA assay in 135 split Assay5 layers (one per
# orig.ident) because 05_integration.R never calls JoinLayers. Subsetting
# such an object and then calling JoinLayers per-subset triggers a Seurat
# 5.1.0 bug ("object 'res' not found") inside JoinLayers. The robust fix is
# to join layers ONCE on the full object here, so each subset() below gets
# a clean single-layer assay and the JoinLayers call inside write_h5ad
# becomes a no-op. Run this block ONCE, right after loading.
da <- DefaultAssay(seu)
if (inherits(seu[[da]], "Assay5") && length(Layers(seu[[da]])) > 1) {
  cat("joining", length(Layers(seu[[da]])), "RNA layers on full object...\n")
  seu <- JoinLayers(seu)
  cat("after join:", length(Layers(seu[[da]])), "layer(s)\n")
} else {
  cat("RNA already single-layer (", length(Layers(seu[[da]])), ")\n", sep = "")
}

# ---- setup: stable export order ------------------------------------------
# Re-runnable. ordered is just a character vector of orig.ident.
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
tp_rank <- function(tp) match(tp, TIMEPOINT_ORDER, nomatch = 9)
ordered <- sm$orig.ident[order(factor(sm$PATIENT_ID), sapply(sm$TIMEPOINT, tp_rank), sm$orig.ident)]
cat("export order:", length(ordered), "specimens\n")

# ---- setup: partition invariants -----------------------------------------
stopifnot(!anyDuplicated(colnames(seu)))

expected <- table(as.character(seu$orig.ident))
observed <- setNames(
  vapply(ordered, function(oi) sum(seu$orig.ident == oi), integer(1)),
  ordered
)

stopifnot(
  identical(sort(names(expected)), sort(names(observed))),
  all(expected[names(observed)] == observed),
  sum(observed) == ncol(seu)
)
cat("partition invariants OK:", length(observed), "partitions cover",
    sum(observed), "of", ncol(seu), "cells\n")

# ---- export loop: all 135 specimens ---------------------------------------
# Re-runnable WITHOUT reloading seu. Resets rows/failed/exported each run.
# To reexport only failed ones, replace `ordered` with `sapply(failed, [[, "orig.ident")`.
exported <- 0L; failed <- list(); rows <- list()
cat("exporting", length(ordered), "samples...\n")
for (oi in ordered) {
  sub <- subset(seu, subset = orig.ident == oi)
  sub$orig.ident <- factor(as.character(sub$orig.ident), levels = oi)
  for (c in names(sub@meta.data)) if (is.factor(sub@meta.data[[c]])) sub@meta.data[[c]] <- droplevels(sub@meta.data[[c]])
  mr  <- sm[sm$orig.ident == oi, ]
  out <- file.path(OUT_DIR, paste0(oi, ".h5ad"))
  tryCatch({
    write_h5ad(sub, out)
    exported <- exported + 1L
    rows[[length(rows) + 1]] <- data.frame(
      orig.ident = oi, PATIENT_ID = mr$PATIENT_ID, TIMEPOINT = mr$TIMEPOINT,
      n_cells = ncol(sub), file = out, sha256 = sha256(out),
      r_version = R.version.string, seurat_version = as.character(packageVersion("Seurat")),
      stringsAsFactors = FALSE)
    cat(sprintf("  [%3d/%d] %s  cells=%d  patient=%s  tp=%s\n",
                exported, length(ordered), oi, ncol(sub), mr$PATIENT_ID, mr$TIMEPOINT))
  }, error = function(e) {
    failed[[length(failed) + 1]] <<- list(orig.ident = oi, error = conditionMessage(e))
    cat(sprintf("  FAILED  %s : %s\n", oi, conditionMessage(e)))
  })
}

# ---- single specimen: reexport one orig.ident ----------------------------
# Re-runnable. Use to retry a single failed specimen without touching the rest.
# Set `oi` to the orig.ident you want, then run this block.
# oi <- "MXMERZ002A_13"
# sub <- subset(seu, subset = orig.ident == oi)
# sub$orig.ident <- factor(as.character(sub$orig.ident), levels = oi)
# for (c in names(sub@meta.data)) if (is.factor(sub@meta.data[[c]])) sub@meta.data[[c]] <- droplevels(sub@meta.data[[c]])
# mr  <- sm[sm$orig.ident == oi, ]
# out <- file.path(OUT_DIR, paste0(oi, ".h5ad"))
# write_h5ad(sub, out)
# cat(sprintf("  %s  cells=%d  patient=%s  tp=%s -> %s\n", oi, ncol(sub), mr$PATIENT_ID, mr$TIMEPOINT, out))

# ---- manifest: write _split_manifest.csv + summary -----------------------
# Re-runnable after the export loop. Uses whatever is in `rows` right now.
manifest <- do.call(rbind, rows)
write.csv(manifest, MANIFEST_OUT, row.names = FALSE)
cat("\n=== SPLIT COMPLETE: ", exported, "/", length(ordered), " exported, ",
    length(failed), " failed, ", sum(manifest$n_cells), " cells ===\n", sep = "")
cat("manifest:", MANIFEST_OUT, "\n")
if (length(failed)) {
  for (f in failed) cat("  ", f$orig.ident, ": ", f$error, "\n", sep = "")
  stop_msg("split completed with ", length(failed), " failure(s).")
}
if (exported != length(ordered)) stop_msg("expected ", length(ordered), ", exported ", exported)
