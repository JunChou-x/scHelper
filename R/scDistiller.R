# ==============================================================================
# scDistiller - 单细胞聚类寻优与精炼工具
# ==============================================================================

# 注意：在正式的 R 包代码中不需要写 library()，依赖将由 DESCRIPTION 文件管理。
# 已经在控制台通过 usethis::use_package() 注册了以下依赖：
# Seurat, dplyr, cluster, ggplot2, presto


#' 自动化聚类寻优 (OptimizeRes)
#'
#' 基于 PCA 降维结果，通过计算轮廓系数 (Silhouette Score) 自动寻找并推荐最佳聚类分辨率。
#' 自动绘制不同分辨率下的平均轮廓系数折线图，并将最佳分辨率记录在对象的 misc 槽中。
#'
#' @param object Seurat 对象。必须已完成 `RunPCA`。
#' @param res_range 数值向量。聚类分辨率的测试范围，默认 \code{seq(0.2, 1.2, 0.2)}。
#' @param n_cells_sample 整数。用于计算轮廓系数的最大降采样细胞数（防内存溢出），默认 5000。
#' @param penalty_threshold 整数。聚类数量过少时的惩罚阈值（鼓励中等粒度），默认 4。
#'
#' @return 包含 \code{optimal_clusters} 新列与 \code{@misc[["scDistiller_best_res"]]} 记录的 Seurat 对象。
#' @export
#'
#' @examples
#' \dontrun{
#' pbmc <- OptimizeRes(pbmc, res_range = seq(0.2, 1.2, 0.2))
#' }
OptimizeRes <- function(object,
                        res_range         = seq(0.2, 1.2, 0.2),
                        n_cells_sample    = 5000,
                        penalty_threshold = 4) {

  if (!"pca" %in% names(object@reductions)) {
    stop("Error: 未检测到 PCA 降维结果，请先运行 RunPCA()。")
  }

  n_pcs      <- min(20, ncol(Embeddings(object, reduction = "pca")))
  pca_embed  <- Embeddings(object, reduction = "pca")[, 1:n_pcs]

  if (nrow(pca_embed) > n_cells_sample) {
    set.seed(42)
    sample_idx       <- sample(seq_len(nrow(pca_embed)), n_cells_sample)
    pca_embed_sample <- pca_embed[sample_idx, ]
  } else {
    sample_idx       <- seq_len(nrow(pca_embed))
    pca_embed_sample <- pca_embed
  }

  pca_dist   <- dist(pca_embed_sample)
  sil_scores <- numeric(length(res_range))
  best_res   <- res_range[1]
  max_sil    <- -Inf

  for (i in seq_along(res_range)) {
    res <- res_range[i]

    cols_before <- colnames(object@meta.data)
    object      <- FindClusters(object, resolution = res, verbose = FALSE)
    new_cols    <- setdiff(colnames(object@meta.data), cols_before)

    if (length(new_cols) == 0) {
      warning(paste("Warning: 分辨率", res, "未生成新聚类列，已跳过。"))
      next
    }

    cluster_col     <- new_cols[length(new_cols)]
    clusters_sample <- as.integer(object[[cluster_col]][sample_idx, 1])
    n_clusters      <- length(unique(clusters_sample))

    if (n_clusters > 1) {
      sil         <- silhouette(clusters_sample, pca_dist)
      penalty     <- ifelse(n_clusters < penalty_threshold,
                            0.1 * (penalty_threshold - n_clusters), 0)
      sil_scores[i] <- mean(sil[, 3]) - penalty
    } else {
      sil_scores[i] <- 0
    }

    if (sil_scores[i] > max_sil) {
      max_sil  <- sil_scores[i]
      best_res <- res
    }
  }

  cols_snap  <- colnames(object@meta.data)
  object     <- FindClusters(object, resolution = best_res, verbose = FALSE)
  best_col   <- setdiff(colnames(object@meta.data), cols_snap)

  if (length(best_col) > 0) {
    best_col <- best_col[length(best_col)]
  } else {
    best_col <- "seurat_clusters"
  }

  object$optimal_clusters             <- object[[best_col]]
  object@misc[["scDistiller_best_res"]] <- best_res

  plot_data         <- data.frame(Resolution = res_range, Silhouette_Score = sil_scores)
  plot_data$is_best <- plot_data$Resolution == best_res

  message(">> 输出可视化绘图数据 (Resolution vs Silhouette_Score):")
  print(plot_data)

  p <- ggplot(plot_data, aes(x = Resolution, y = Silhouette_Score)) +
    geom_line(color = "#2c3e50", linewidth = 1) +
    geom_point(aes(color = is_best, size = is_best)) +
    scale_color_manual(values = c("FALSE" = "#95a5a6", "TRUE" = "#c0392b")) +
    scale_size_manual(values  = c("FALSE" = 2.5, "TRUE" = 4.5)) +
    guides(color = "none", size = "none") +
    geom_vline(xintercept = best_res, linetype = "dashed", color = "#bdc3c7") +
    theme_classic() +
    theme(
      text          = element_text(color = "black"),
      axis.text     = element_text(color = "black"),
      axis.title    = element_text(face = "bold"),
      plot.title    = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5)
    ) +
    labs(title    = "Silhouette Score Optimization",
         subtitle = paste("Recommended resolution:", best_res),
         x = "Resolution", y = "Mean Silhouette Score (Adjusted)")
  print(p)

  message(">> 推荐的最佳 Resolution 为: ", best_res)

  if (best_res == min(res_range)) {
    warning("提示：推荐分辨率为测试范围最小值，可能存在欠聚类现象，建议结合下游生物学特征综合判断。")
  }

  return(object)
}


#' 极速特征提取与热图 (PlotTopMarkers)
#'
#' 利用 presto 包极速执行 Wilcoxon 秩和检验，提取各聚类群的 Top N 标志基因，
#' 并自动绘制出版级标准化热图。
#'
#' @param object Seurat 对象。
#' @param cluster_col 字符。用于提取 Marker 的聚类列名，默认 \code{"optimal_clusters"}。
#' @param top_n 整数。每个 Cluster 提取并在热图中展示的 Top 基因数量，默认 5。
#'
#' @return 一个包含 Top Marker 基因统计信息的 tibble 数据框。
#' @export
#'
#' @examples
#' \dontrun{
#' markers <- PlotTopMarkers(pbmc, cluster_col = "optimal_clusters", top_n = 5)
#' }
PlotTopMarkers <- function(object,
                           cluster_col = "optimal_clusters",
                           top_n       = 5) {

  if (!cluster_col %in% colnames(object@meta.data)) {
    stop(sprintf("Error: '%s' 列不存在于 meta.data 中。", cluster_col))
  }

  Idents(object) <- cluster_col
  markers        <- wilcoxauc(object, group_by = cluster_col)

  top_markers <- markers %>%
    filter(logFC > 0.25, pct_in >= 25, padj < 0.05) %>%
    group_by(feature) %>%
    slice_max(order_by = auc, n = 1) %>%
    ungroup() %>%
    group_by(group) %>%
    slice_max(order_by = auc, n = top_n) %>%
    ungroup()

  message(">> 输出可视化绘图数据 (Top Markers 列表):")
  print(top_markers %>% select(group, feature, auc, logFC), n = Inf)

  features_to_plot <- top_markers$feature
  object           <- ScaleData(object, features = features_to_plot, verbose = FALSE)

  # 使用 suppressMessages 屏蔽色标强行替换时的 ggplot 提示
  p <- suppressMessages(
    DoHeatmap(object, features = features_to_plot, group.by = cluster_col) +
      scale_fill_gradientn(colors = c("#4575b4", "white", "#d73027")) +
      theme(
        axis.text.y = element_text(face = "italic", size = 8, color = "black"),
        axis.text.x = element_text(color = "black")
      )
  )
  print(p)

  return(top_markers)
}


#' 通用锚点过滤引擎 (PurifyCells)
#'
#' 为特定下游分析（如 CellChat / NicheNet）精准切取目标细胞生态位。
#' 基于给定的核心基因和白名单基因表达阈值，将数据集二元划分为保留群体与背景群体。
#'
#' @param object Seurat 对象。
#' @param core_markers 字符向量。核心目标群体表达的基因列表。
#' @param whitelist_markers 字符向量。需要豁免保留的旁系细胞群基因列表，默认为 \code{NULL}。
#' @param expr_cutoff 数值。判定的最低 log-normalized 表达量阈值，默认 0.1。
#' @param group_by 字符。用于绘制纯化比例柱状图的分组列名，默认 \code{"condition"}。
#' @param label_keep 字符。达标细胞的注释标签，默认 \code{"Kept"}。
#' @param label_exclude 字符。未达标杂细胞的注释标签，默认 \code{"Excluded"}。
#' @param verbose 逻辑值。是否在终端打印过滤统计信息，默认 \code{TRUE}。
#'
#' @return 包含 \code{distilled_annotation} 新列的 Seurat 对象。
#' @export
#'
#' @examples
#' \dontrun{
#' pbmc <- PurifyCells(pbmc, core_markers = c("CD3D"), whitelist_markers = c("CD14"))
#' }
PurifyCells <- function(object,
                        core_markers,
                        whitelist_markers = NULL,
                        expr_cutoff       = 0.1,
                        group_by          = "condition",
                        label_keep        = "Kept",
                        label_exclude     = "Excluded",
                        verbose           = TRUE) {

  current_assay <- DefaultAssay(object)
  expr_data     <- GetAssayData(object, assay = current_assay, layer = "data")

  target_genes <- unique(c(core_markers, whitelist_markers))
  valid_genes  <- intersect(target_genes, rownames(expr_data))

  if (length(valid_genes) == 0) {
    stop("Error: 提供的核心或白名单基因均不在当前表达矩阵中，请检查基因名称。")
  }

  sub_expr     <- expr_data[valid_genes, , drop = FALSE]
  is_expressed <- sub_expr > expr_cutoff

  valid_core <- intersect(core_markers, valid_genes)
  is_core    <- if (length(valid_core) > 0) {
    colSums(is_expressed[valid_core, , drop = FALSE]) > 0
  } else {
    rep(FALSE, ncol(object))
  }

  is_whitelist <- rep(FALSE, ncol(object))
  if (!is.null(whitelist_markers)) {
    valid_whitelist <- intersect(whitelist_markers, valid_genes)
    if (length(valid_whitelist) > 0) {
      is_whitelist <- colSums(is_expressed[valid_whitelist, , drop = FALSE]) > 0
    }
  }

  final_status           <- rep(label_exclude, ncol(object))
  names(final_status)    <- colnames(object)
  final_status[is_core | is_whitelist] <- label_keep
  object$distilled_annotation <- final_status

  if (verbose) {
    message(sprintf(
      ">> 过滤汇总:\n - 核心目标命中: %d\n - 白名单豁免命中: %d\n - 分类为 %s: %d",
      sum(is_core),
      sum(is_whitelist & !is_core),
      label_exclude,
      sum(final_status == label_exclude)
    ))
  }

  if (group_by %in% colnames(object@meta.data)) {
    plot_data <- as.data.frame(object@meta.data)

    message(sprintf(">> 输出可视化绘图数据 (以 %s 为分组的纯化比例):", group_by))
    print(table(plot_data[[group_by]], plot_data$distilled_annotation))

    fill_colors <- c("#E64B35", "#3C5488")
    names(fill_colors) <- c(label_exclude, label_keep)

    p <- ggplot(plot_data, aes(x = .data[[group_by]], fill = distilled_annotation)) +
      geom_bar(position = "fill", color = "black", linewidth = 0.3, width = 0.7) +
      scale_fill_manual(values = fill_colors) +
      scale_y_continuous(expand = c(0, 0)) +
      theme_classic() +
      theme(
        text         = element_text(color = "black"),
        axis.text    = element_text(color = "black"),
        axis.title   = element_text(face = "bold"),
        legend.title = element_text(face = "bold")
      ) +
      labs(y = "Cell Proportion", x = group_by, fill = "Annotation")
    print(p)
  } else {
    warning(sprintf("提示: 在 meta.data 中未找到 '%s' 列，已跳过分组比例图绘制。", group_by))
  }

  return(object)
}


#' 聚类稳定性检验 (CheckStability)
#'
#' 对推荐的聚类分辨率进行无放回 Bootstrap 抽样检验，计算其拓扑结构的稳定性分数。
#'
#' @param object Seurat 对象。
#' @param res 数值。待检验的聚类分辨率。若为 \code{NULL}，将自动从 misc 读取推荐值。
#' @param n_iter 整数。Bootstrap 重抽样迭代次数，默认 10。
#' @param sample_ratio 数值。每次抽样保留的细胞比例，默认 0.8。
#'
#' @return 包含迭代聚类计数、主导聚类数、稳定性分数和语义解读的列表。
#' @export
#'
#' @examples
#' \dontrun{
#' stability_res <- CheckStability(pbmc, n_iter = 10)
#' }
CheckStability <- function(object,
                           res          = NULL,
                           n_iter       = 10,
                           sample_ratio = 0.8) {

  if (is.null(res)) {
    res <- object@misc[["scDistiller_best_res"]]
    if (is.null(res)) stop("Error: 请提供 res 参数，或先运行 OptimizeRes()。")
    message(">> 自动使用推荐的分辨率: ", res)
  }

  cells    <- colnames(object)
  n_sample <- floor(length(cells) * sample_ratio)

  message(sprintf(">> 开始执行聚类稳定性 Bootstrap 检验 (迭代 %d 次，无放回抽样)...", n_iter))
  cluster_counts <- integer(n_iter)

  for (i in seq_len(n_iter)) {
    set.seed(42 + i)

    sub_cells <- sample(cells, n_sample, replace = FALSE)
    sub_obj   <- subset(object, cells = sub_cells)

    sub_obj <- FindNeighbors(sub_obj, dims = 1:min(20, ncol(Embeddings(sub_obj, "pca"))),
                             verbose = FALSE)

    # 运行 FindClusters 后，Seurat 默认会将新聚类设为 Active Identity
    sub_obj <- FindClusters(sub_obj, resolution = res, verbose = FALSE)

    # 直接读取 Idents 数量，规避 subset 继承 meta.data 导致的新旧列名比对失败
    cluster_counts[i] <- length(unique(Idents(sub_obj)))
  }

  cluster_counts <- cluster_counts[!is.na(cluster_counts)]

  mode_val        <- as.integer(names(sort(table(cluster_counts), decreasing = TRUE)[1]))
  stability_score <- sum(cluster_counts == mode_val) / length(cluster_counts)

  message(sprintf(">> 检验完成。主导 Cluster 数量: %d，稳定性分数: %.2f",
                  mode_val, stability_score))

  interpretation <- dplyr::case_when(
    stability_score >= 0.8 ~ "高度稳定 (Highly Stable)",
    stability_score >= 0.5 ~ "中度稳定 (Moderately Stable)，建议结合生物学验证",
    TRUE                   ~ "不稳定 (Unstable)，建议调整分辨率"
  )
  message(">> 解读: ", interpretation)

  return(list(
    counts          = cluster_counts,
    mode_clusters   = mode_val,
    stability_score = stability_score,
    interpretation  = interpretation
  ))
}
