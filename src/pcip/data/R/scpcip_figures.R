# scPCIP figures: organ map (ggPlantmap), UMAPs, marker heatmap and violins (SCP).
# Rscript scpcip_figures.R <input dir> <dataset id> <organ map: root|leaf|leaf_grass|flower|seed|none> <output dir>
# <input dir> holds <id>.predictions.tsv, <id>.umap.tsv, <id>.markers.tsv and expr/ (data.mtx genes x cells of
# log-normalised expression for the marker member genes, genes.tsv, cells.tsv).
args <- commandArgs(trailingOnly = TRUE)
dir_in <- args[1]; dsid <- args[2]; organ <- args[3]; dir_out <- args[4]
dir.create(dir_out, showWarnings = FALSE, recursive = TRUE)
suppressMessages({library(ggplot2); library(dplyr); library(patchwork); library(Matrix)})
options(Seurat.object.assay.version = "v3")

palette <- c(
  epidermis = "#8DB65A", cortex = "#E6B85C", endodermis = "#B85C38", pericycle = "#9B7BB0",
  procambium = "#7FB3D5", phloem = "#D9829B", xylem = "#34618E", root_cap = "#CDBB8A",
  mesophyll = "#4E8A3E", vascular = "#34618E",
  pollen = "#EBCB52", stamen = "#E39A55", carpel = "#B5546F", perianth = "#E2A9C8", non_floral = "#B3B3B3",
  embryo = "#5E9E6E", endosperm = "#EED79A", seed_coat = "#93683F",
  dermal = "#8DB65A", ground = "#E6B85C")
conf_palette <- c(Q1 = "#E3E1DC", Q2 = "#A9C6DD", Q3 = "#5B8DB8", Q4 = "#244A73")
unassigned <- "#ECEAE4"

pred <- read.delim(file.path(dir_in, paste0(dsid, ".predictions.tsv")), row.names = 1, colClasses = c(cell_id = "character"))
umap <- read.delim(file.path(dir_in, paste0(dsid, ".umap.tsv")), row.names = 1, colClasses = c(cell_id = "character"))
mk <- read.delim(file.path(dir_in, paste0(dsid, ".markers.tsv")))

extra <- setdiff(unique(c(pred$predicted_cell_type, pred$tissue_system)), c(names(palette), ""))
palette <- c(palette, setNames(grDevices::hcl.colors(max(length(extra), 1), "Dark 3")[seq_along(extra)], extra))

save_plot <- function(p, name, w, h, dpi = 300) {
  ggsave(file.path(dir_out, paste0("scPCIP_", name, ".png")), p, width = w, height = h, dpi = dpi, bg = "white")
  ggsave(file.path(dir_out, paste0("scPCIP_", name, ".pdf")), p, width = w, height = h, bg = "white")
}

if (organ != "none") {
  suppressMessages(library(ggPlantmap))
  stat <- pred |> group_by(type = predicted_cell_type) |>
    summarise(n = n(), hi = mean(confidence_quartile >= 3), .groups = "drop") |>
    mutate(share = n / sum(n))
  maps <- switch(organ,
    root = list(
      ggPm.At.roottip.longitudinal |>
        mutate(type = case_when(Level1 == "Meristem" ~ NA_character_,
                                Level1 %in% c("Columella", "Lateral Root Cap") ~ "root_cap",
                                Level2 == "Vascular Tissue" ~ "procambium",
                                Level2 == "Pericycle" ~ "pericycle",
                                TRUE ~ tolower(Level1))),
      ggPm.At.roottip.crosssection |>
        mutate(type = recode(ROI.name, Epidermis = "epidermis", Columella = "root_cap", Cortex = "cortex",
                             Endodermis = "endodermis", Pericycle = "pericycle", Phloem = "phloem",
                             Procambium = "procambium", Xylem = "xylem"))),
    leaf = list(
      ggPm.At.leaf.crosssection |>
        mutate(type = recode(Cell, epidermis = "epidermis", Parenchima = "mesophyll", vascularbundle = "vascular"))),
    leaf_grass = list(
      ggPm.Os.leaf.crosssection |>
        mutate(type = case_when(ROI.name %in% c("Epidermis", "Bulliform", "Guard.Cell", "Trichome") ~ "epidermis",
                                ROI.name == "Mesophyll" ~ "mesophyll",
                                ROI.name == "Large.Parenchyma" ~ NA_character_,
                                TRUE ~ "vascular"))),
    flower = list(
      ggPm.At.flower.diagram |>
        mutate(type = recode(Organ, Sepal = "perianth", Petal = "perianth", Stamen = "stamen", Pistil = "carpel"))),
    seed = list(
      ggPm.At.seed.devseries |> filter(Stage == "Maturegreen") |>
        mutate(type = recode(Part, seedcoat = "seed_coat", embryo = "embryo", endosperm = "endosperm"))),
    stop("unknown organ map: ", organ))

  # pollen grains drawn inside each stamen
  pollen_dots <- function(m) {
    st <- m |> filter(type == "stamen") |> group_by(ROI.id) |>
      summarise(cx = mean(x), cy = mean(y), r = 0.18 * min(diff(range(x)), diff(range(y))), .groups = "drop")
    merge(st, expand.grid(dx = c(-1, 0, 1), dy = c(-1, 0, 1))) |> mutate(x = cx + dx * r, y = cy + dy * r)
  }
  map_plot <- function(m) {
    m <- m |> mutate(fill = ifelse(is.na(type) | !type %in% stat$type, "unassigned", type))
    p <- ggplot(m, aes(x, y, group = ROI.id)) +
      geom_polygon(aes(fill = fill), colour = "white", linewidth = 0.35) +
      scale_fill_manual(values = c(palette, unassigned = unassigned), guide = "none") +
      coord_equal() + theme_void()
    if (organ == "flower" && "pollen" %in% stat$type)
      p <- p + geom_point(data = pollen_dots(m), aes(x, y), inherit.aes = FALSE, shape = 21, size = 2.2,
                          fill = palette[["pollen"]], colour = "white", stroke = 0.3)
    p
  }
  panels <- lapply(maps, map_plot)

  leg <- stat |> arrange(desc(n)) |> mutate(row = nrow(stat) + 1 - seq_along(type), label = gsub("_", " ", type))
  hdr <- nrow(stat) + 1
  p_leg <- ggplot(leg, aes(y = row)) +
    geom_tile(aes(x = 0, fill = type), width = 0.7, height = 0.7) +
    geom_text(aes(x = 0.6, label = label), hjust = 0, size = 3.8) +
    geom_text(aes(x = 4.6, label = format(n, big.mark = ",")), hjust = 1, size = 3.6, colour = "grey20") +
    geom_text(aes(x = 6.2, label = scales::percent(share, 0.1)), hjust = 1, size = 3.6, colour = "grey20") +
    geom_text(aes(x = 7.9, label = scales::percent(hi, 1)), hjust = 1, size = 3.6, colour = "grey20") +
    annotate("text", x = c(0.6, 4.6, 6.2, 7.9), y = hdr, label = c("cell type", "cells", "share", "high conf."),
             hjust = c(0, 1, 1, 1), size = 3.3, fontface = "bold", colour = "grey30") +
    annotate("segment", x = -0.35, xend = 7.9, y = hdr - 0.5, yend = hdr - 0.5, colour = "grey75", linewidth = 0.3) +
    scale_fill_manual(values = palette, guide = "none") +
    scale_x_continuous(limits = c(-0.5, 8.1)) + scale_y_continuous(limits = c(min(0.3, hdr - 8.7), hdr + 0.6)) +
    theme_void()
  widths <- switch(organ, root = c(0.7, 1.2, 1.25), leaf = c(1.6, 1.1), leaf_grass = c(1.5, 1.1),
                   flower = c(1.2, 1.1), seed = c(0.9, 1.1))
  w <- switch(organ, root = 12, leaf = 11, leaf_grass = 11, flower = 10, seed = 9)
  save_plot(wrap_plots(c(panels, list(p_leg)), nrow = 1, widths = widths), "organ_map", w, 5, dpi = 200)
}

suppressMessages({library(Seurat); library(SCP)})
# GroupHeatmap draws a second legend for the feature groups that repeats the cell-type legend; drop it
group_heatmap <- SCP::GroupHeatmap
body(group_heatmap) <- parse(text = gsub('lgd[["Cluster"]] <- Legend(', 'lgd_feature_split <- Legend(',
                                         deparse(body(SCP::GroupHeatmap)), fixed = TRUE))

expr <- readMM(file.path(dir_in, "expr", "data.mtx"))
rownames(expr) <- readLines(file.path(dir_in, "expr", "genes.tsv"))
colnames(expr) <- readLines(file.path(dir_in, "expr", "cells.tsv"))
expr <- as(expr, "CsparseMatrix")
stopifnot(identical(colnames(expr), rownames(pred)))
srt <- CreateSeuratObject(expr, meta.data = pred)
srt <- SetAssayData(srt, layer = "data", new.data = expr)
emb <- as.matrix(umap[colnames(srt), c("umap1", "umap2")]); colnames(emb) <- c("UMAP_1", "UMAP_2")
srt[["umap"]] <- CreateDimReducObject(embeddings = emb, key = "UMAP_", assay = "RNA")
types <- names(sort(table(srt$predicted_cell_type), decreasing = TRUE))
srt$cell_type <- factor(gsub("_", " ", srt$predicted_cell_type), levels = gsub("_", " ", types))
pal <- palette; names(pal) <- gsub("_", " ", names(pal))
srt$confidence <- factor(paste0("Q", srt$confidence_quartile), levels = paste0("Q", 1:4))

p <- CellDimPlot(srt, group.by = "cell_type", reduction = "umap", palcolor = pal[levels(srt$cell_type)],
                 label = TRUE, label_insitu = TRUE, label_repel = TRUE, label.size = 4, pt.size = 1.2,
                 theme_use = "theme_blank", show_stat = FALSE, title = NULL)
save_plot(p, "umap_cell_type", 6.5, 5.5)

p <- CellDimPlot(srt, group.by = "confidence", reduction = "umap", palcolor = unname(conf_palette), pt.size = 1.0,
                 theme_use = "theme_blank", show_stat = FALSE, title = NULL)
save_plot(p, "umap_confidence", 6.2, 5.5)

if ("tissue_system" %in% colnames(srt@meta.data)) {
  srt$tissue <- factor(gsub("_", " ", srt$tissue_system))
  p <- CellDimPlot(srt, group.by = "tissue", reduction = "umap", palcolor = pal[levels(srt$tissue)],
                   label = TRUE, label_insitu = TRUE, pt.size = 1.0, theme_use = "theme_blank", show_stat = FALSE,
                   title = NULL)
  save_plot(p, "umap_tissue_system", 6.2, 5.5)
}

feats <- intersect(c("margin", "expected_accuracy", "state_meristem_score"), colnames(srt@meta.data))
feats <- feats[vapply(feats, function(f) any(is.finite(srt@meta.data[[f]])), logical(1))]
p <- FeatureDimPlot(srt, features = feats, reduction = "umap", pt.size = 0.8, theme_use = "theme_blank",
                    ncol = length(feats), palette = "Blues", bg_cutoff = -Inf)
save_plot(p, "umap_scores", 4.2 * length(feats), 4.2)

# per program orthogroup, the member gene with the highest mean expression in the type it marks
data <- GetAssayData(srt, layer = "data")
mk <- mk |> filter(as.character(present_in_query) %in% c("True", "TRUE"), log2fc_in_vs_rest > 0,
                   detect_frac_in_type > 0.05, cell_type %in% types) |>
  rowwise() |>
  mutate(gene = {
    g <- intersect(strsplit(query_member_genes, ",")[[1]], rownames(srt))
    cells <- colnames(srt)[srt$predicted_cell_type == cell_type]
    if (length(g) == 0) NA_character_ else g[which.max(Matrix::rowMeans(data[g, cells, drop = FALSE]))]
  }) |> ungroup() |>
  filter(!is.na(gene)) |>
  group_by(cell_type) |> arrange(desc(log2fc_in_vs_rest), .by_group = TRUE) |> slice_head(n = 4) |> ungroup() |>
  mutate(cell_type = factor(gsub("_", " ", cell_type), levels = levels(srt$cell_type))) |> arrange(cell_type) |>
  distinct(gene, .keep_all = TRUE)
write.table(data.frame(cell_type = mk$cell_type, orthogroup = mk$orthogroup_id, gene = mk$gene),
            file.path(dir_out, "scPCIP_marker_genes_plotted.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
if (nrow(mk) >= 2) {
  ht <- group_heatmap(srt, features = mk$gene, group.by = "cell_type", feature_split = mk$cell_type, slot = "data",
                      group_palcolor = list(pal[levels(srt$cell_type)]),
                      feature_split_palcolor = list(pal[levels(droplevels(mk$cell_type))]),
                      add_dot = TRUE, add_bg = TRUE, nlabel = 0, show_row_names = TRUE,
                      heatmap_palette = "YlOrRd", width = 4, height = 6, row_names_side = "left", row_title = " ", column_title = " ",
                      row_title_rot = 0)
  save_plot(ht$plot, "marker_heatmap", 8, 8)
  top <- mk$gene[!duplicated(mk$cell_type)]  # the top marker of each type
  vd <- data.frame(cell_type = rep(srt$cell_type, length(top)),
                   gene = factor(rep(top, each = ncol(srt)), levels = top),
                   value = as.vector(t(as.matrix(data[top, , drop = FALSE]))))
  p <- ggplot(vd, aes(cell_type, value, fill = cell_type)) +
    geom_violin(scale = "width", linewidth = 0.25, colour = "grey25", adjust = 1.2) +
    facet_grid(gene ~ ., scales = "free_y", switch = "y") +
    scale_fill_manual(values = pal[levels(srt$cell_type)], guide = "none") +
    scale_y_continuous(position = "right", n.breaks = 3) +
    labs(x = NULL, y = NULL) + theme_classic(base_size = 11) +
    theme(strip.placement = "outside", strip.background = element_blank(),
          strip.text.y.left = element_text(angle = 0, hjust = 1), panel.spacing = unit(2, "pt"),
          axis.text.x = element_text(angle = 45, hjust = 1), axis.text.y.right = element_text(size = 7))
  save_plot(p, "marker_violin", 0.8 * length(levels(srt$cell_type)) + 2.5, 0.6 * length(top) + 1.5)
}
