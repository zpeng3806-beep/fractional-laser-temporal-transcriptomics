#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(limma)
  library(msigdbr)
})

set.seed(20260823)
script_arg <- commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))]
script_path <- sub("^--file=", "", script_arg)
root <- normalizePath(file.path(dirname(script_path), "..", ".."))
out_dir <- file.path(root, "PROJECT_C_FRACTIONAL_LASER", "03_results")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

read_geo_matrix <- function(path) {
  d <- fread(cmd = paste("gzip -cd", shQuote(path)), skip = "!series_matrix_table_begin", header = TRUE, showProgress = FALSE)
  ids <- d[[1]]
  x <- as.matrix(d[, -1])
  storage.mode(x) <- "double"
  rownames(x) <- ids
  x
}

map_and_collapse <- function(x, platform_path, log2_needed) {
  annotation <- fread(platform_path, select = c("ID", "Gene Symbol"), showProgress = FALSE)
  setnames(annotation, c("probe", "symbol"))
  annotation <- annotation[!is.na(symbol) & symbol != "---" & !grepl("///", symbol, fixed = TRUE)]
  annotation <- annotation[match(rownames(x), probe)]
  keep <- !is.na(annotation$symbol)
  x <- x[keep, , drop = FALSE]
  if (log2_needed) x <- log2(x + 1)
  x <- avereps(x, ID = annotation$symbol[keep])
  x
}

make_programs <- function() {
  h <- as.data.table(msigdbr(species = "Homo sapiens", collection = "H"))
  c2 <- as.data.table(msigdbr(species = "Homo sapiens", collection = "C2"))
  c5 <- as.data.table(msigdbr(species = "Homo sapiens", collection = "C5"))
  wanted <- c(
    acute_injury = "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
    inflammation = "HALLMARK_INFLAMMATORY_RESPONSE",
    epidermal_repair = "GOBP_POSITIVE_REGULATION_OF_EPITHELIAL_CELL_PROLIFERATION_INVOLVED_IN_WOUND_HEALING",
    ecm_remodeling = "REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION",
    collagen_organization = "GOBP_COLLAGEN_FIBRIL_ORGANIZATION",
    proliferation = "HALLMARK_E2F_TARGETS",
    matrix_maturation = "GOBP_EXTRACELLULAR_MATRIX_ASSEMBLY"
  )
  all <- rbindlist(list(h[, .(gs_name, gene_symbol)], c2[, .(gs_name, gene_symbol)], c5[, .(gs_name, gene_symbol)]))
  lapply(wanted, function(name) unique(all[gs_name == name, gene_symbol]))
}

score_programs <- function(x, programs) {
  z <- t(scale(t(x)))
  z[!is.finite(z)] <- NA_real_
  scores <- vapply(programs, function(genes) {
    hit <- intersect(genes, rownames(z))
    if (length(hit) < 5) return(rep(NA_real_, ncol(z)))
    colMeans(z[hit, , drop = FALSE], na.rm = TRUE)
  }, numeric(ncol(z)))
  scores <- t(scores)
  colnames(scores) <- colnames(x)
  scores
}

exact_signflip_p <- function(d) {
  d <- d[is.finite(d)]
  n <- length(d)
  observed <- abs(mean(d))
  if (n <= 20) {
    total <- 2^n
    means <- numeric(total)
    for (i in 0:(total - 1)) {
      signs <- ifelse(as.logical(intToBits(i)[seq_len(n)]), 1, -1)
      means[i + 1] <- mean(d * signs)
    }
    return(mean(abs(means) >= observed - 1e-15))
  }
  means <- replicate(99999, mean(d * sample(c(-1, 1), n, replace = TRUE)))
  (1 + sum(abs(means) >= observed)) / 100000
}

summarize_pairs <- function(scores, meta, dataset, day_levels) {
  rows <- list()
  for (program in rownames(scores)) {
    for (day in day_levels) {
      people <- intersect(meta$participant_id[meta$day == 0], meta$participant_id[meta$day == day])
      diffs <- vapply(people, function(person) {
        b <- meta$geo_accession[meta$participant_id == person & meta$day == 0][1]
        t <- meta$geo_accession[meta$participant_id == person & meta$day == day][1]
        scores[program, t] - scores[program, b]
      }, numeric(1))
      diffs <- diffs[is.finite(diffs)]
      n <- length(diffs)
      mean_diff <- mean(diffs)
      sd_diff <- sd(diffs)
      se <- sd_diff / sqrt(n)
      ci <- mean_diff + c(-1, 1) * qt(0.975, df = n - 1) * se
      dz <- mean_diff / sd_diff
      hedges_g <- dz * (1 - 3 / (4 * n - 5))
      full_sign <- sign(mean_diff)
      lopo <- vapply(seq_along(diffs), function(i) sign(mean(diffs[-i])) == full_sign, logical(1))
      rows[[length(rows) + 1]] <- data.table(
        dataset = dataset, program = program, day = day, n_pairs = n,
        mean_paired_score_change = mean_diff, ci95_low = ci[1], ci95_high = ci[2],
        hedges_g_paired = hedges_g, exact_signflip_p = exact_signflip_p(diffs),
        lopo_direction_fraction = mean(lopo), direction = ifelse(mean_diff > 0, "UP", "DOWN")
      )
    }
  }
  rbindlist(rows)
}

programs <- make_programs()
fwrite(data.table(program = names(programs), n_genes_definition = lengths(programs)), file.path(out_dir, "laser_program_definitions.tsv"), sep = "\t")
fwrite(rbindlist(lapply(names(programs), function(nm) data.table(program = nm, gene_symbol = programs[[nm]]))),
       file.path(out_dir, "laser_program_genes.tsv"), sep = "\t")

# AFL: baseline and single-treatment time course.
x_a <- read_geo_matrix(file.path(root, "PROJECT_C_FRACTIONAL_LASER", "01_data", "GSE168760_series_matrix.txt.gz"))
x_a <- map_and_collapse(x_a, file.path(root, "PROJECT_C_FRACTIONAL_LASER", "00_metadata", "GPL13667_platform.soft.tsv"), TRUE)
m_a <- fread(file.path(root, "PROJECT_C_FRACTIONAL_LASER", "00_metadata", "GSE168760_metadata.tsv"))
m_a[, day := as.integer(sub("day ", "", time))]
m_a[, treatment_number := as.integer(treatment_number)]
m_a <- m_a[(treatment_number == 0 & day == 0) | treatment_number == 1]
x_a <- x_a[, m_a$geo_accession, drop = FALSE]
s_a <- score_programs(x_a, programs)
r_a <- summarize_pairs(s_a, m_a, "GSE168760_AFL", c(1, 3, 7, 14, 21, 28))

# NAFL: arm only, baseline and clean first-treatment observations. Day 29 is
# explicitly excluded because it is also one day after the second session.
x_n <- read_geo_matrix(file.path(root, "PROJECT_C_FRACTIONAL_LASER", "01_data", "GSE206495_series_matrix.txt.gz"))
x_n <- map_and_collapse(x_n, file.path(root, "PROJECT_C_FRACTIONAL_LASER", "00_metadata", "GPL15207_platform.soft.tsv"), FALSE)
m_n <- fread(file.path(root, "PROJECT_C_FRACTIONAL_LASER", "00_metadata", "GSE206495_metadata.tsv"))
m_n[, day := fcase(time == "baseline", 0L,
                   grepl("^1 day after first", time), 1L,
                   grepl("^7 days after first", time), 7L,
                   grepl("^14 days after first", time), 14L,
                   default = NA_integer_)]
m_n <- m_n[anatomical_site == "arm skin" & !is.na(day)]
x_n <- x_n[, m_n$geo_accession, drop = FALSE]
s_n <- score_programs(x_n, programs)
r_n <- summarize_pairs(s_n, m_n, "GSE206495_NAFL", c(1, 7, 14))

effects <- rbind(r_a, r_n)
fwrite(effects, file.path(out_dir, "laser_program_temporal_effects.tsv"), sep = "\t")

# Cross-platform conservation at matched days.
cross <- merge(r_a[day %in% c(1, 7, 14), .(program, day, afl = mean_paired_score_change, afl_lopo = lopo_direction_fraction)],
               r_n[, .(program, day, nafl = mean_paired_score_change, nafl_lopo = lopo_direction_fraction)],
               by = c("program", "day"))
cross[, same_direction := sign(afl) == sign(nafl)]
fwrite(cross, file.path(out_dir, "laser_cross_platform_conservation.tsv"), sep = "\t")

summary <- cross[, .(
  spearman_rho = cor(afl, nafl, method = "spearman"),
  directional_concordance = mean(same_direction),
  stable_same_direction_programs = sum(same_direction & afl_lopo == 1 & nafl_lopo == 1)
), by = day]
fwrite(summary, file.path(out_dir, "laser_cross_platform_summary.tsv"), sep = "\t")

# Peak timing is descriptive; no half-life model is fitted to sparse/nonmonotone data.
peak <- effects[, .SD[which.max(abs(mean_paired_score_change))], by = .(dataset, program)]
fwrite(peak, file.path(out_dir, "laser_descriptive_peak_timing.tsv"), sep = "\t")

capture.output(sessionInfo(), file = file.path(root, "PROJECT_C_FRACTIONAL_LASER", "04_logs", "sessionInfo.txt"))
