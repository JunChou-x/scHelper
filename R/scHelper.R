# ==============================================================================
# scHelper - 单细胞全生命周期自动化助手
# 命名规范：sh 前缀 + 功能词，按分析流程顺序编号
#
# 标准工作流：
#   shQC() -> shDoublet() -> shHVG() -> shPC() -> shCellCycle() ->
#   shBatch() -> shCluster() -> shMarkers() -> shFilter() -> shStability()
# ==============================================================================

# ==============================================================================
# 第一部分：前处理自动化 (Preprocessing)
# ==============================================================================

#' Step 1: 自动化质控寻优 (shQC)
#'
#' 利用中位数绝对偏差 (MAD) 算法自动识别离群细胞，替代人工看小提琴图划定阈值的步骤。
#'
#' @param object     Seurat 对象。
#' @param mt_pattern 字符。线粒体基因识别正则表达式，默认 "^MT-"（人）；小鼠数据请传入 "^mt-"。
#' @param n_mads     数值。离群值判定的 MAD 倍数，默认 3。
#'
#' @return 含 `qc_pass`（"Pass"/"Fail"）列的 Seurat 对象。
#' @export
shQC <- function(object, mt_pattern = "^MT-", n_mads = 3) {
  if (!"percent.mt" %in% colnames(object@meta.data)) {
    object[["percent.mt"]] <- PercentageFeatureSet(object, pattern = mt_pattern)
  }
  meta <- object@meta.data
  features_to_check  <- c("nFeature_RNA", "nCount_RNA", "percent.mt")
  cutoffs <- list()
  
  for (feat in features_to_check) {
    if (!feat %in% colnames(meta)) next
    vals  <- meta[[feat]]
    med   <- median(vals, na.rm = TRUE)
    mad_v <- mad(vals, na.rm = TRUE)
    if (feat == "percent.mt") {
      upper <- med + n_mads * mad_v
      cutoffs[[feat]] <- c(lower = 0, upper = upper)
      meta[[paste0(feat, "_outlier")]] <- vals > upper
    } else {
      lower <- max(0, med - n_mads * mad_v)
      upper <- med + n_mads * mad_v
      cutoffs[[feat]] <- c(lower = lower, upper = upper)
      meta[[paste0(feat, "_outlier")]] <- vals < lower | vals > upper
    }
  }
  
  existing_cols <- paste0(features_to_check, "_outlier")
  existing_cols <- existing_cols[existing_cols %in% colnames(meta)]
  meta$qc_pass <- ifelse(!apply(meta[, existing_cols, drop = FALSE], 1, any), "Pass", "Fail")
  object$qc_pass <- meta$qc_pass
  
  message(">> [shQC] 推荐阈值 (基于 ", n_mads, " 倍 MAD):")
  print(round(do.call(rbind, cutoffs), 2))
  message(">> [shQC] 过滤汇总:")
  print(table(object$qc_pass))
  
  plot_data_long <- tidyr::pivot_longer(
    meta[, c(features_to_check, "qc_pass")],
    cols = tidyselect::all_of(features_to_check),
    names_to = "variable", values_to = "value"
  )
  cutoff_lines <- do.call(rbind, lapply(names(cutoffs), function(f) {
    data.frame(variable = f, lower = unname(cutoffs[[f]]["lower"]), upper = unname(cutoffs[[f]]["upper"]))
  }))
  cutoff_lines$lower[cutoff_lines$lower == 0] <- NA
  
  p <- ggplot(plot_data_long, aes(x = variable, y = value)) +
    geom_violin(aes(fill = qc_pass), color = NA, alpha = 0.6) +
    geom_jitter(aes(color = qc_pass), size = 0.1, alpha = 0.3, width = 0.2) +
    geom_hline(data = cutoff_lines, aes(yintercept = upper), color = "#E64B35", linetype = "dashed", linewidth = 0.8) +
    geom_hline(data = cutoff_lines, aes(yintercept = lower), color = "#E64B35", linetype = "dashed", linewidth = 0.8, na.rm = TRUE) +
    scale_fill_manual(values = c("Fail" = "#E64B35", "Pass" = "#3C5488")) +
    scale_color_manual(values = c("Fail" = "#E64B35", "Pass" = "#3C5488")) +
    facet_wrap(~variable, scales = "free") +
    theme_classic() +
    theme(strip.background = element_blank(), strip.text = element_text(face = "bold", size = 12)) +
    labs(title = "Automated QC Filtering (MAD-based)", x = NULL, y = "Value")
  print(p)
  
  return(object)
}


#' Step 2: 双细胞率自动估算与鉴定 (shDoublet)
#'
#' 自适应单/多样本模式，基于 10x 官方经验公式自动推算预期双细胞率并鉴定。
#' 支持跨平台（包含 Windows）的并行加速。
#'
#' @param object          Seurat 对象。
#' @param run_scDblFinder 逻辑值。是否运行 scDblFinder，默认 TRUE。
#' @param batch_col       字符。样本分组列名，用于多样本独立测算，默认 "orig.ident"。
#' @param ncores          整数。并行计算的 CPU 核心数，默认 1（单线程）。
#'
#' @return 含 `scDblFinder_class` 和 `scDblFinder_score` 的 Seurat 对象。
#' @export
shDoublet <- function(object, run_scDblFinder = TRUE, batch_col = "orig.ident", ncores = 1) {
  if (batch_col %in% colnames(object@meta.data)) {
    n_samples <- length(unique(object@meta.data[[batch_col]]))
  } else {
    n_samples <- 1
    warning(sprintf("[shDoublet] 找不到分组列 '%s'，强制以单样本模式运行。", batch_col))
  }
  
  if (n_samples > 1) {
    mean_cells <- ncol(object) / n_samples
    expected_rate <- 0.008 * (mean_cells / 1000)
    message(sprintf(">> [shDoublet] 多样本模式：%d 个样本，平均 %d 细胞/样本", n_samples, as.integer(mean_cells)))
  } else {
    expected_rate <- 0.008 * (ncol(object) / 1000)
    message(sprintf(">> [shDoublet] 单样本模式：%d 个细胞", ncol(object)))
  }
  
  message(sprintf(">> [shDoublet] 推算预期双细胞率: %.2f%%", expected_rate * 100))
  object@misc[["scHelper_doublet_rate"]] <- expected_rate
  
  if (!run_scDblFinder) return(object)
  
  if (!requireNamespace("scDblFinder", quietly = TRUE)) {
    warning("[shDoublet] 未安装 scDblFinder。请运行 BiocManager::install('scDblFinder')。")
    return(object)
  }
  
  # Windows 兼容的多核处理
  bp_param <- if (ncores > 1 && requireNamespace("BiocParallel", quietly = TRUE)) {
    if (.Platform$OS.type == "windows") {
      BiocParallel::SnowParam(workers = ncores)
    } else {
      BiocParallel::MulticoreParam(workers = ncores)
    }
  } else {
    BiocParallel::SerialParam()
  }
  
  sce <- as.SingleCellExperiment(object)
  sce <- if (n_samples > 1) {
    scDblFinder::scDblFinder(sce, samples = object@meta.data[[batch_col]], dbr = expected_rate, BPPARAM = bp_param)
  } else {
    scDblFinder::scDblFinder(sce, dbr = expected_rate, BPPARAM = bp_param)
  }
  
  object$scDblFinder_class <- sce$scDblFinder.class
  object$scDblFinder_score <- sce$scDblFinder.score
  
  message(">> [shDoublet] 鉴定完成:")
  print(table(object$scDblFinder_class))
  return(object)
}


#' Step 3: HVG 数量寻优 (shHVG)
#'
#' 利用累积标准化方差饱和曲线与几何弦距法自动选择最优高变基因数量。
#'
#' @param object       Seurat 对象。
#' @param max_features 整数。评估的最大特征数，默认 5000。
#'
#' @return 完成 FindVariableFeatures 并记录最优点至 `scHelper_optimal_hvg` 的 Seurat 对象。
#' @export
shHVG <- function(object, max_features = 5000) {
  message(">> [shHVG] 计算全局标准化方差...")
  object <- FindVariableFeatures(object, selection.method = "vst", nfeatures = max_features, verbose = FALSE)
  hvf_info <- HVFInfo(object, method = "vst")
  var_data <- sort(hvf_info$variance.standardized, decreasing = TRUE)
  
  x <- 500:max_features
  y <- cumsum(var_data)[x]
  A <- y[length(y)] - y[1]; B <- x[1] - x[length(x)]; C <- x[length(x)] * y[1] - x[1] * y[length(y)]
  distances <- abs(A * x + B * y + C) / sqrt(A^2 + B^2)
  optimal_hvg <- x[which.max(distances)]
  
  object <- FindVariableFeatures(object, selection.method = "vst", nfeatures = optimal_hvg, verbose = FALSE)
  object@misc[["scHelper_optimal_hvg"]] <- optimal_hvg
  
  message(sprintf(">> [shHVG] 最优 HVG 数量 (弦距法): %d", optimal_hvg))
  
  p <- ggplot(data.frame(N_Features = x, Cumulative_Variance = y), aes(x = N_Features, y = Cumulative_Variance)) +
    geom_line(color = "#3C5488", linewidth = 1) +
    geom_vline(xintercept = optimal_hvg, color = "#E64B35", linetype = "dashed", linewidth = 1) +
    annotate("text", x = optimal_hvg + 100, y = min(y), label = paste("Optimal HVG =", optimal_hvg), color = "#E64B35", hjust = 0, fontface = "bold") +
    theme_classic() +
    labs(title = "Optimal HVG Selection", x = "Number of Features", y = "Cumulative Standardized Variance")
  print(p)
  
  return(object)
}


#' Step 4: 主成分数量寻优 (shPC)
#'
#' 利用几何弦距法在主成分方差碎石图上自动定位数学拐点。
#'
#' @param object Seurat 对象。必须已完成 RunPCA。
#' @param max_pc 整数。评估的最大 PC 数，默认 50。
#'
#' @return 记录最优点至 `scHelper_optimal_pc` 的 Seurat 对象。
#' @export
shPC <- function(object, max_pc = 50) {
  if (!"pca" %in% names(object@reductions)) stop("[shPC] 未检测到 PCA 降维结果，请先执行 RunPCA()。")
  
  stdev <- Stdev(object, reduction = "pca")
  n_pcs <- min(length(stdev), max_pc)
  var_exp <- (stdev[1:n_pcs]^2) / sum(stdev^2) * 100
  
  x <- seq_len(n_pcs); y <- var_exp
  A <- y[n_pcs] - y[1]; B <- x[1] - x[n_pcs]; C <- x[n_pcs] * y[1] - x[1] * y[n_pcs]
  distances <- abs(A * x + B * y + C) / sqrt(A^2 + B^2)
  optimal_pc <- which.max(distances)
  
  object@misc[["scHelper_optimal_pc"]] <- optimal_pc
  
  message(sprintf(">> [shPC] 最优 PC 数量 (弦距法): %d", optimal_pc))
  
  p <- ggplot(data.frame(PC = x, Variance = var_exp), aes(x = PC, y = Variance)) +
    geom_line(color = "#3C5488", linewidth = 1) +
    geom_point(color = "#3C5488", size = 2) +
    geom_vline(xintercept = optimal_pc, color = "#E64B35", linetype = "dashed", linewidth = 1) +
    annotate("text", x = optimal_pc + 0.5, y = max(var_exp) * 0.9, label = paste("Optimal PC =", optimal_pc), color = "#E64B35", hjust = 0, fontface = "bold") +
    theme_classic() +
    labs(title = "Optimal PC Selection (Elbow)", x = "Principal Component", y = "Variance Explained (%)")
  print(p)
  
  return(object)
}


#' Step 5: 细胞周期效应评估 (shCellCycle)
#'
#' 定量评估细胞周期打分对 PC1 的解释度。
#'
#' @param object       Seurat 对象。必须已完成 ScaleData 和 RunPCA。
#' @param s_genes      字符向量。S 期特征基因列表。
#' @param g2m_genes    字符向量。G2M 期特征基因列表。
#' @param r2_threshold 数值。判定细胞周期效应是否显著的 R-squared 阈值，默认 0.15。
#'
#' @return 含细胞周期评分的 Seurat 对象，并打印建议。
#' @export
shCellCycle <- function(object, s_genes = cc.genes.updated.2019$s.genes, g2m_genes = cc.genes.updated.2019$g2m.genes, r2_threshold = 0.15) {
  if (!"pca" %in% names(object@reductions)) stop("[shCellCycle] 未检测到 PCA，请先运行 RunPCA()。")
  
  object <- CellCycleScoring(object, s.features = s_genes, g2m.features = g2m_genes, set.ident = FALSE)
  meta <- object@meta.data
  pc1 <- Embeddings(object, "pca")[, 1]
  cc_data <- data.frame(pc1 = pc1, S.Score = meta$S.Score, G2M.Score = meta$G2M.Score)
  
  r_sq <- summary(lm(pc1 ~ S.Score + G2M.Score, data = cc_data))$r.squared
  message(sprintf(">> [shCellCycle] 细胞周期对 PC1 的方差解释度 R-squared: %.4f", r_sq))
  
  p <- ggplot(data.frame(PC1 = pc1, Phase = meta$Phase), aes(x = Phase, y = PC1, fill = Phase)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.2, size = 0.5, alpha = 0.5, color = "darkgray") +
    scale_fill_manual(values = c("G1" = "#4DBBD5", "S" = "#00A087", "G2M" = "#E64B35")) +
    theme_classic() +
    labs(title = "Cell Cycle Effect on PC1", subtitle = sprintf("R-squared = %.4f", r_sq))
  print(p)
  
  if (r_sq > r2_threshold) {
    message(sprintf(">> [shCellCycle] [警告] 效应显著 (R-squared > %.2f)。建议在 ScaleData 中回归 S.Score 和 G2M.Score。", r2_threshold))
  } else {
    message(sprintf(">> [shCellCycle] [正常] 效应微弱 (R-squared <= %.2f)，无需回归。", r2_threshold))
  }
  return(object)
}


# ==============================================================================
# 第二部分：降维自动化 (Dimensionality Reduction)
# ==============================================================================

#' Step 6: 批次效应检测 (shBatch)
#'
#' 利用 LISI 分数定量评估批次混合程度，并输出智能提示词。
#'
#' @param object    Seurat 对象。
#' @param batch_col 字符。批次信息列名（必填）。
#' @param reduction 字符。用于评估的降维空间，默认 "pca"。
#'
#' @return 以 invisible 方式返回 LISI 得分数据框，并渲染分布图。
#' @export
shBatch <- function(object, batch_col, reduction = "pca") {
  if (!requireNamespace("lisi", quietly = TRUE)) stop("[shBatch] 请安装 lisi 包: devtools::install_github('immunogenomics/lisi')")
  if (!batch_col %in% colnames(object@meta.data)) stop("[shBatch] 找不到批次列: ", batch_col)
  
  embed <- Embeddings(object, reduction = reduction)
  meta <- object@meta.data
  n_batches <- length(unique(meta[[batch_col]]))
  n_cells <- ncol(object)
  
  message(">> [shBatch] 计算 Local Inverse Simpson's Index (LISI)...")
  lisi_res <- lisi::compute_lisi(embed, meta, c(batch_col))
  mean_lisi <- mean(lisi_res[[batch_col]], na.rm = TRUE)
  
  message(sprintf(">> [shBatch] 平均 LISI: %.2f (满分: %d，接近 1 表示批次分离严重)", mean_lisi, n_batches))
  
  p <- ggplot(lisi_res, aes(x = .data[[batch_col]])) +
    geom_density(fill = "#3C5488", alpha = 0.5) +
    geom_vline(xintercept = mean_lisi, color = "#E64B35", linetype = "dashed", linewidth = 1) +
    theme_classic() +
    labs(title = "Batch Effect Assessment (LISI)", subtitle = paste("Mean LISI =", round(mean_lisi, 2)), x = "LISI Score", y = "Density")
  print(p)
  
  if (mean_lisi < n_batches * 0.6) {
    message(">> [shBatch] [警告] 批次混合不良，建议执行批次校正。")
    if (n_cells > 50000) {
      message(">> [shBatch] [建议] 细胞数 > 50k，推荐 Seurat RPCA 或 Harmony 以兼顾速度与内存。")
    } else {
      message(">> [shBatch] [建议] 细胞数适中，推荐 Harmony 或 CCA 整合。")
    }
  } else {
    message(">> [shBatch] [正常] 批次混合良好，通常无需额外校正。")
  }
  
  return(invisible(lisi_res))
}


# ==============================================================================
# 第三部分：聚类自动化 (Clustering)
# ==============================================================================

#' Step 7: 聚类分辨率寻优 (shCluster)
#'
#' 遍历测试一系列聚类分辨率，基于轮廓系数选择最优解，并施加欠聚类惩罚。
#'
#' @param object            Seurat 对象。
#' @param res_range         数值向量。测试的分辨率范围，默认 seq(0.2, 1.2, 0.2)。
#' @param n_cells_sample    整数。计算轮廓系数的采样细胞数上限，防止大样本内存溢出，默认 5000。
#' @param penalty_threshold 整数。当亚群数量小于此值时，按比例扣减轮廓系数分数，默认 4。
#'
#' @return 含 `optimal_clusters` 列及记录至 `scHelper_best_res` 的 Seurat 对象。
#' @export
shCluster <- function(object, res_range = seq(0.2, 1.2, 0.2), n_cells_sample = 5000, penalty_threshold = 4) {
  # 避免 %||% 依赖 rlang，使用基础 R 写法
  optimal_pc <- object@misc[["scHelper_optimal_pc"]]
  if (is.null(optimal_pc)) optimal_pc <- 20
  
  pca_embed <- Embeddings(object, reduction = "pca")[, 1:min(optimal_pc, ncol(Embeddings(object, "pca")))]
  
  if (nrow(pca_embed) > n_cells_sample) {
    set.seed(42)
    sample_idx <- sample(seq_len(nrow(pca_embed)), n_cells_sample)
  } else {
    sample_idx <- seq_len(nrow(pca_embed))
  }
  
  pca_dist <- dist(pca_embed[sample_idx, ])
  sil_scores <- numeric(length(res_range))
  best_res <- res_range[1]
  max_sil <- -Inf
  
  for (i in seq_along(res_range)) {
    res <- res_range[i]
    object <- FindClusters(object, resolution = res, verbose = FALSE)
    clusters_sample <- as.integer(Idents(object)[sample_idx])
    n_clusters <- length(unique(clusters_sample))
    
    if (n_clusters > 1) {
      sil <- cluster::silhouette(clusters_sample, pca_dist)
      penalty <- ifelse(n_clusters < penalty_threshold, 0.1 * (penalty_threshold - n_clusters), 0)
      sil_scores[i] <- mean(sil[, 3]) - penalty
    } else {
      sil_scores[i] <- 0
    }
    
    if (sil_scores[i] > max_sil) {
      max_sil <- sil_scores[i]
      best_res <- res
    }
  }
  
  object <- FindClusters(object, resolution = best_res, verbose = FALSE)
  object$optimal_clusters <- as.character(Idents(object))
  object@misc[["scHelper_best_res"]] <- best_res
  
  plot_data <- data.frame(Resolution = res_range, Silhouette_Score = sil_scores)
  plot_data$is_best <- plot_data$Resolution == best_res
  
  message(sprintf(">> [shCluster] 最优分辨率: %.2f", best_res))
  
  p <- ggplot(plot_data, aes(x = Resolution, y = Silhouette_Score)) +
    geom_line(color = "steelblue", linewidth = 1) +
    geom_point(aes(color = is_best, size = is_best)) +
    scale_color_manual(values = c("FALSE" = "darkred", "TRUE" = "gold")) +
    scale_size_manual(values  = c("FALSE" = 3, "TRUE" = 5)) +
    guides(color = "none", size = "none") +
    geom_vline(xintercept = best_res, linetype = "dashed", color = "grey50") +
    theme_classic() +
    labs(title = "Silhouette Score Optimization", subtitle = paste("Recommended resolution:", best_res), x = "Resolution", y = "Adjusted Silhouette Score")
  print(p)
  
  if (best_res == min(res_range)) warning("[shCluster] 推荐分辨率为测试下限，可能存在欠聚类，建议结合 shMarkers() 判断。")
  
  return(object)
}


#' Step 8: 特征提取与热图 (shMarkers)
#'
#' 利用 presto 包执行极速 Wilcoxon 检验提取 Marker，兼容所有 R 版本的语法结构。
#'
#' @param object      Seurat 对象。
#' @param cluster_col 字符。分组列名，默认 "optimal_clusters"。
#' @param top_n       整数。每个簇提取的 Marker 数量，默认 5。
#'
#' @return 提取的 Top Markers 数据框（tibble），并在视窗打印热图。
#' @export
shMarkers <- function(object, cluster_col = "optimal_clusters", top_n = 5) {
  if (!requireNamespace("presto", quietly = TRUE)) stop("[shMarkers] 请安装 presto。")
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("[shMarkers] 请安装 dplyr。")
  
  Idents(object) <- cluster_col
  markers <- presto::wilcoxauc(object, group_by = cluster_col)
  
  # 无管道符连续赋值，兼容低版本 R
  res <- dplyr::filter(markers, logFC > 0.25, pct_in >= 25, padj < 0.05)
  res <- dplyr::group_by(res, feature)
  res <- dplyr::slice_max(res, order_by = auc, n = 1, with_ties = FALSE)
  res <- dplyr::ungroup(res)
  res <- dplyr::group_by(res, group)
  res <- dplyr::slice_max(res, order_by = auc, n = top_n, with_ties = FALSE)
  top_markers <- dplyr::ungroup(res)
  
  message(">> [shMarkers] Top Marker 列表:")
  print(dplyr::select(top_markers, group, feature, auc, logFC), n = Inf)
  
  features_to_plot <- top_markers$feature
  object <- ScaleData(object, features = features_to_plot, verbose = FALSE)
  
  # 取消 |> 管道符，改为基础打印
  p <- suppressMessages(
    DoHeatmap(object, features = features_to_plot, group.by = cluster_col) + 
      scale_fill_gradientn(colors = c("#4575b4", "white", "#d73027"))
  )
  print(p)
  
  return(top_markers)
}


#' Step 9: 目标亚群切取 (shFilter)
#'
#' 基于核心和白名单 Marker 基因列表，精准切取所需亚群并标记保留状态。
#'
#' @param object            Seurat 对象。
#' @param core_markers      字符向量。核心需要保留的亚群特征基因。
#' @param whitelist_markers 字符向量。额外需要保留的亚群特征基因（可选）。
#' @param expr_cutoff       数值。基因表达量判定阈值，默认 0.1。
#' @param group_by          字符。用于统计提纯比例的分组列，默认 "condition"。
#' @param label_keep        字符。保留簇标签名，默认 "Kept"。
#' @param label_exclude     字符。排除簇标签名，默认 "Excluded"。
#'
#' @return 新增 `distilled_annotation` 二分类标签列的 Seurat 对象。
#' @export
shFilter <- function(object, core_markers, whitelist_markers = NULL, expr_cutoff = 0.1, group_by = "condition", label_keep = "Kept", label_exclude = "Excluded") {
  expr_data <- GetAssayData(object, assay = DefaultAssay(object), layer = "data")
  valid_genes <- intersect(unique(c(core_markers, whitelist_markers)), rownames(expr_data))
  
  if (length(valid_genes) == 0) stop("[shFilter] 错误：提供的基因无一在表达矩阵中。")
  
  sub_expr <- expr_data[valid_genes, , drop = FALSE] > expr_cutoff
  valid_core <- intersect(core_markers, valid_genes)
  
  is_core <- if (length(valid_core) > 0) {
    colSums(sub_expr[valid_core, , drop = FALSE]) > 0
  } else {
    rep(FALSE, ncol(object))
  }
  
  is_wl <- if (!is.null(whitelist_markers)) {
    valid_whitelist <- intersect(whitelist_markers, valid_genes)
    if (length(valid_whitelist) > 0) {
      colSums(sub_expr[valid_whitelist, , drop = FALSE]) > 0
    } else {
      rep(FALSE, ncol(object))
    }
  } else {
    rep(FALSE, ncol(object))
  }
  
  object$distilled_annotation <- ifelse(is_core | is_wl, label_keep, label_exclude)
  
  message(sprintf(">> [shFilter] 过滤汇总:\n - %s: %d\n - %s: %d", label_keep, sum(object$distilled_annotation == label_keep), label_exclude, sum(object$distilled_annotation == label_exclude)))
  
  if (group_by %in% colnames(object@meta.data)) {
    message(sprintf(">> [shFilter] 分组比例 (按 %s):", group_by))
    print(table(object@meta.data[[group_by]], object$distilled_annotation))
  }
  
  return(object)
}


#' Step 10: 聚类稳定性检验 (shStability)
#'
#' 利用 Bootstrap 重采样模拟数据扰动，检验指定分辨率下的聚类数量稳定性。
#'
#' @param object       Seurat 对象。
#' @param res          数值。需检验的分辨率，默认读取 `shCluster` 存入的 `scHelper_best_res`。
#' @param n_iter       整数。重采样迭代次数，默认 10。
#' @param sample_ratio 数值。每次随机抽取的细胞比例，默认 0.8。
#'
#' @return 包含 counts (每次迭代的聚类数), mode_clusters (众数聚类数) 和 stability_score (稳定性得分) 的列表。
#' @export
shStability <- function(object, res = NULL, n_iter = 10, sample_ratio = 0.8) {
  if (is.null(res)) {
    res <- object@misc[["scHelper_best_res"]]
    if (is.null(res)) stop("[shStability] 请提供 res 参数，或先运行 shCluster()。")
  }
  
  optimal_pc <- object@misc[["scHelper_optimal_pc"]]
  if (is.null(optimal_pc)) optimal_pc <- 20
  
  cells <- colnames(object)
  n_sample <- floor(length(cells) * sample_ratio)
  cluster_counts <- integer(n_iter)
  
  message(sprintf(">> [shStability] Bootstrap 检验开始 (%d 次迭代, 抽取 %.0f%%)...", n_iter, sample_ratio * 100))
  
  for (i in seq_len(n_iter)) {
    set.seed(42 + i)
    sub_obj <- subset(object, cells = sample(cells, n_sample, replace = FALSE))
    sub_obj <- FindNeighbors(sub_obj, dims = 1:min(optimal_pc, ncol(Embeddings(sub_obj, "pca"))), verbose = FALSE)
    sub_obj <- FindClusters(sub_obj, resolution = res, verbose = FALSE)
    cluster_counts[i] <- length(unique(Idents(sub_obj)))
  }
  
  mode_val <- as.integer(names(sort(table(cluster_counts), decreasing = TRUE)[1]))
  stability_score <- sum(cluster_counts == mode_val) / length(cluster_counts)
  
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    # 如果没装 dplyr，用 ifelse
    interpretation <- ifelse(stability_score >= 0.8, "高度稳定 (Highly Stable)", 
                             ifelse(stability_score >= 0.5, "中度稳定 (Moderately Stable)", "不稳定，建议调整分辨率 (Unstable)"))
  } else {
    interpretation <- dplyr::case_when(
      stability_score >= 0.8 ~ "高度稳定 (Highly Stable)",
      stability_score >= 0.5 ~ "中度稳定，建议结合生物学验证 (Moderately Stable)",
      TRUE                   ~ "不稳定，建议调整分辨率 (Unstable)"
    )
  }
  
  message(sprintf(">> [shStability] 主导 Cluster 数: %d | 稳定性分数: %.2f", mode_val, stability_score))
  message(">> [shStability] 专家解读: ", interpretation)
  
  return(list(
    counts          = cluster_counts,
    mode_clusters   = mode_val,
    stability_score = stability_score,
    interpretation  = interpretation
  ))
}