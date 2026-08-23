#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(limma)
})

SEED <- 20260823L
N_BOOT <- 5000L
N_RANDOM <- 2000L
PRIMARY_DAYS <- c(1L, 7L, 14L)
set.seed(SEED)

script_arg <- commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))]
script_path <- sub("^--file=", "", script_arg)
PROJECT <- normalizePath(Sys.getenv(
  "FRACTIONAL_LASER_PROJECT_ROOT",
  unset = file.path(dirname(script_path), "..")
), mustWork = FALSE)
META <- file.path(PROJECT, "metadata")
DATA <- Sys.getenv("FRACTIONAL_LASER_DATA_DIR", unset = file.path(PROJECT, "data"))
RESULTS <- file.path(PROJECT, "reproduction_output", "phase1b")
FIGURES <- file.path(PROJECT, "reproduction_output", "phase1b_figures")
LOGS <- file.path(PROJECT, "reproduction_output", "logs")
dir.create(RESULTS, recursive = TRUE, showWarnings = FALSE)
dir.create(FIGURES, recursive = TRUE, showWarnings = FALSE)
dir.create(LOGS, recursive = TRUE, showWarnings = FALSE)

read_geo_matrix <- function(path) {
  d <- fread(cmd = paste("gzip -cd", shQuote(path)), skip = "!series_matrix_table_begin",
             header = TRUE, showProgress = FALSE)
  ids <- d[[1L]]
  x <- as.matrix(d[, -1L])
  storage.mode(x) <- "double"
  rownames(x) <- ids
  x
}

build_mapping <- function(x, platform_path, platform, log2_needed) {
  ann <- fread(platform_path, select = c("ID", "Gene Symbol", "Entrez Gene"), showProgress = FALSE)
  setnames(ann, c("probe_id", "gene_symbol_raw", "entrez_raw"))
  ann <- ann[match(rownames(x), probe_id)]
  ann[, platform := platform]
  ann[, symbol_unique := !is.na(gene_symbol_raw) & gene_symbol_raw != "---" &
        !grepl("///", gene_symbol_raw, fixed = TRUE) & !grepl(",", gene_symbol_raw, fixed = TRUE)]
  ann[, entrez_unique := !is.na(entrez_raw) & entrez_raw != "---" &
        !grepl("///", entrez_raw, fixed = TRUE)]
  ann[, gene_symbol := fifelse(symbol_unique, trimws(gene_symbol_raw), NA_character_)]
  ann[, entrez_id := fifelse(entrez_unique, trimws(entrez_raw), NA_character_)]
  conflicts <- ann[symbol_unique & entrez_unique,
                   .(n_entrez_per_symbol = uniqueN(entrez_id)), by = gene_symbol]
  ann <- conflicts[ann, on = "gene_symbol"]
  ann[is.na(n_entrez_per_symbol), n_entrez_per_symbol := 0L]
  ann[, mapping_status := fcase(
    !symbol_unique, "MISSING_OR_AMBIGUOUS_SYMBOL",
    !entrez_unique, "MISSING_OR_AMBIGUOUS_ENTREZ",
    n_entrez_per_symbol > 1L, "SYMBOL_ENTREZ_CONFLICT",
    default = "RELIABLE_UNIQUE_SYMBOL_ENTREZ"
  )]
  ann[, used_reliable := mapping_status == "RELIABLE_UNIQUE_SYMBOL_ENTREZ"]

  x_use <- x[ann$used_reliable, , drop = FALSE]
  if (log2_needed) x_use <- log2(x_use + 1)
  ann_use <- ann[used_reliable == TRUE]
  primary <- avereps(x_use, ID = ann_use$gene_symbol)

  probe_iqr <- apply(x_use, 1L, IQR, na.rm = TRUE)
  choice <- data.table(probe_id = ann_use$probe_id, gene_symbol = ann_use$gene_symbol,
                       entrez_id = ann_use$entrez_id, probe_iqr = probe_iqr,
                       row_index = seq_len(nrow(x_use)))
  setorder(choice, gene_symbol, -probe_iqr, probe_id)
  choice <- choice[, .SD[1L], by = gene_symbol]
  alternative <- x_use[choice$row_index, , drop = FALSE]
  rownames(alternative) <- choice$gene_symbol
  ann[, selected_highest_iqr_probe := probe_id %in% choice$probe_id]
  list(primary = primary, alternative = alternative, audit = ann,
       gene_map = unique(ann_use[, .(gene_symbol, entrez_id)]))
}

z_by_gene <- function(x) {
  z <- t(scale(t(x)))
  z[!is.finite(z)] <- NA_real_
  z
}

score_programs <- function(z, programs) {
  ans <- vapply(programs, function(genes) {
    hit <- intersect(genes, rownames(z))
    if (length(hit) < 2L) return(rep(NA_real_, ncol(z)))
    colMeans(z[hit, , drop = FALSE], na.rm = TRUE)
  }, numeric(ncol(z)))
  ans <- t(ans)
  colnames(ans) <- colnames(z)
  ans
}

paired_hedges <- function(d) {
  d <- d[is.finite(d)]
  n <- length(d)
  if (n < 3L || sd(d) == 0) return(NA_real_)
  (mean(d) / sd(d)) * (1 - 3 / (4 * n - 5))
}

exact_signflip_p <- function(d) {
  d <- d[is.finite(d)]
  n <- length(d)
  observed <- abs(mean(d))
  if (n > 20L) stop("Exact sign-flip implementation is bounded to n <= 20")
  means <- numeric(2^n)
  for (i in 0:(2^n - 1L)) {
    signs <- ifelse(as.logical(intToBits(i)[seq_len(n)]), 1, -1)
    means[i + 1L] <- mean(d * signs)
  }
  mean(abs(means) >= observed - 1e-15)
}

make_score_long <- function(scores, meta, dataset) {
  d <- as.data.table(t(scores), keep.rownames = "geo_accession")
  d <- melt(d, id.vars = "geo_accession", variable.name = "program", value.name = "program_score")
  d <- meta[d, on = "geo_accession"]
  d[, dataset := dataset]
  d[]
}

fit_program_models <- function(score_long, dataset, primary_days = PRIMARY_DAYS) {
  effects <- list()
  diffs_out <- list()
  assumptions <- list()
  for (pr in unique(score_long$program)) {
    d <- score_long[program == pr & is.finite(program_score)]
    d[, participant_factor := factor(participant_id)]
    d[, day_factor := relevel(factor(day), ref = "0")]
    fit <- lm(program_score ~ day_factor + participant_factor, data = d)
    residuals <- residuals(fit)
    shapiro_p <- if (length(residuals) >= 3L && length(residuals) <= 5000L) shapiro.test(residuals)$p.value else NA_real_
    absdev <- abs(residuals - ave(residuals, d$day, FUN = median))
    bf_p <- tryCatch(anova(lm(absdev ~ factor(d$day)))[1L, "Pr(>F)"], error = function(e) NA_real_)
    assumptions[[length(assumptions) + 1L]] <- data.table(
      dataset = dataset, program = pr, residual_n = length(residuals),
      residual_shapiro_p = shapiro_p, brown_forsythe_p = bf_p,
      robust_companion = "exact_signflip_on_patient_differences"
    )
    for (dy in sort(unique(d$day[d$day != 0L]))) {
      coef_name <- paste0("day_factor", dy)
      ct <- summary(fit)$coefficients
      if (!coef_name %in% rownames(ct)) next
      people <- intersect(d[day == 0L, participant_id], d[day == dy, participant_id])
      pd <- rbindlist(lapply(people, function(pid) {
        base <- d[participant_id == pid & day == 0L, program_score][1L]
        treated <- d[participant_id == pid & day == dy, program_score][1L]
        data.table(dataset = dataset, program = pr, day = dy, participant_id = pid,
                   baseline_score = base, treated_score = treated, difference = treated - base)
      }))
      diff <- pd$difference
      lopo <- vapply(seq_along(diff), function(i) sign(mean(diff[-i])) == sign(mean(diff)), logical(1))
      effects[[length(effects) + 1L]] <- data.table(
        dataset = dataset, program = pr, day = dy, n_participants = length(diff),
        model_effect = ct[coef_name, "Estimate"], model_se = ct[coef_name, "Std. Error"],
        model_ci95_low = ct[coef_name, "Estimate"] - qt(.975, df.residual(fit)) * ct[coef_name, "Std. Error"],
        model_ci95_high = ct[coef_name, "Estimate"] + qt(.975, df.residual(fit)) * ct[coef_name, "Std. Error"],
        model_t = ct[coef_name, "t value"], model_df = df.residual(fit), model_p = ct[coef_name, "Pr(>|t|)"],
        paired_mean_change = mean(diff), paired_sd = sd(diff), paired_hedges_g = paired_hedges(diff),
        exact_signflip_p = exact_signflip_p(diff),
        lopo_direction_fraction = mean(lopo), direction = ifelse(mean(diff) > 0, "UP", "DOWN"),
        primary_timepoint = dy %in% primary_days
      )
      diffs_out[[length(diffs_out) + 1L]] <- pd
    }
  }
  e <- rbindlist(effects)
  e[, model_BH_q := p.adjust(model_p, method = "BH"), by = dataset]
  list(effects = e, differences = rbindlist(diffs_out), assumptions = rbindlist(assumptions))
}

effect_from_diff_matrix <- function(mat) {
  apply(mat, 2L, paired_hedges)
}

diff_matrix <- function(differences, dataset_value, day_value, programs) {
  d <- differences[dataset == dataset_value & day == day_value]
  w <- dcast(d, participant_id ~ program, value.var = "difference")
  ans <- as.matrix(w[, ..programs])
  rownames(ans) <- w$participant_id
  ans
}

spearman_safe <- function(a, b) {
  if (sum(is.finite(a) & is.finite(b)) < 3L) return(NA_real_)
  suppressWarnings(cor(a, b, method = "spearman", use = "complete.obs"))
}

all_permutations <- function(x) {
  if (length(x) == 1L) return(matrix(x, nrow = 1L))
  do.call(rbind, lapply(seq_along(x), function(i) {
    rest <- all_permutations(x[-i])
    cbind(x[i], rest)
  }))
}

cross_day_test <- function(mat_a, mat_n, day, n_boot = N_BOOT) {
  eff_a <- effect_from_diff_matrix(mat_a)
  eff_n <- effect_from_diff_matrix(mat_n)
  rho <- spearman_safe(eff_a, eff_n)
  raw_a <- colMeans(mat_a)
  raw_n <- colMeans(mat_n)
  raw_rho <- spearman_safe(raw_a, raw_n)
  boot <- numeric(n_boot)
  for (b in seq_len(n_boot)) {
    ba <- mat_a[sample(seq_len(nrow(mat_a)), nrow(mat_a), replace = TRUE), , drop = FALSE]
    bn <- mat_n[sample(seq_len(nrow(mat_n)), nrow(mat_n), replace = TRUE), , drop = FALSE]
    boot[b] <- spearman_safe(effect_from_diff_matrix(ba), effect_from_diff_matrix(bn))
  }
  boot <- boot[is.finite(boot)]
  perms <- all_permutations(seq_along(eff_n))
  perm_rho <- apply(perms, 1L, function(idx) spearman_safe(eff_a, eff_n[idx]))
  p_perm <- mean(abs(perm_rho) >= abs(rho) - 1e-15)
  data.table(day = day, standardized_spearman_rho = rho,
             bootstrap_ci95_low = unname(quantile(boot, .025)),
             bootstrap_ci95_high = unname(quantile(boot, .975)),
             bootstrap_replicates_valid = length(boot), exact_label_permutation_p = p_perm,
             exact_label_permutations = nrow(perms), raw_mean_change_spearman_rho = raw_rho,
             directional_concordance = mean(sign(eff_a) == sign(eff_n)))
}

quick_diff_matrices <- function(scores, meta, days = PRIMARY_DAYS) {
  lapply(days, function(dy) {
    people <- intersect(meta[day == 0L, participant_id], meta[day == dy, participant_id])
    ans <- do.call(rbind, lapply(people, function(pid) {
      b <- meta[participant_id == pid & day == 0L, geo_accession][1L]
      t <- meta[participant_id == pid & day == dy, geo_accession][1L]
      scores[, t] - scores[, b]
    }))
    if (is.null(dim(ans))) ans <- matrix(ans, ncol = nrow(scores))
    colnames(ans) <- rownames(scores)
    rownames(ans) <- people
    ans
  })
}

program_variant_effect <- function(z, meta, genes, program, dataset, day_value = 14L) {
  scores <- matrix(colMeans(z[genes, , drop = FALSE], na.rm = TRUE), nrow = 1L,
                   dimnames = list(program, colnames(z)))
  d <- quick_diff_matrices(scores, meta, day_value)[[1L]][, 1L]
  n <- length(d)
  se <- sd(d) / sqrt(n)
  ci <- mean(d) + c(-1, 1) * qt(.975, n - 1L) * se
  lopo <- vapply(seq_along(d), function(i) sign(mean(d[-i])) == sign(mean(d)), logical(1))
  data.table(dataset = dataset, program = program, day = day_value, n_participants = n,
             paired_mean_change = mean(d), paired_se = se, paired_ci95_low = ci[1L],
             paired_ci95_high = ci[2L], paired_t_p = 2 * pt(-abs(mean(d) / se), df = n - 1L),
             paired_hedges_g = paired_hedges(d), lopo_direction_fraction = mean(lopo),
             direction = ifelse(mean(d) > 0, "UP", "DOWN"))
}

# ---- Load, annotate, and define analysis samples ----
x_a_raw <- read_geo_matrix(file.path(DATA, "GSE168760_series_matrix.txt.gz"))
x_n_raw <- read_geo_matrix(file.path(DATA, "GSE206495_series_matrix.txt.gz"))
m_a <- fread(file.path(META, "GSE168760_metadata.tsv"))
m_a[, day := as.integer(sub("day ", "", time))]
m_a[, treatment_number := as.integer(treatment_number)]
m_a <- m_a[(treatment_number == 0L & day == 0L) | treatment_number == 1L]
setorder(m_a, participant_id, day)
m_n <- fread(file.path(META, "GSE206495_metadata.tsv"))
m_n[, day := fcase(time == "baseline", 0L,
                   grepl("^1 day after first", time), 1L,
                   grepl("^7 days after first", time), 7L,
                   grepl("^14 days after first", time), 14L,
                   default = NA_integer_)]
m_n <- m_n[anatomical_site == "arm skin" & !is.na(day)]
setorder(m_n, participant_id, day)
stopifnot(uniqueN(m_a$participant_id) == 14L, uniqueN(m_n$participant_id) == 17L)
stopifnot(!anyDuplicated(m_a[, .(participant_id, day)]), !anyDuplicated(m_n[, .(participant_id, day)]))

map_a <- build_mapping(x_a_raw[, m_a$geo_accession, drop = FALSE],
                       file.path(DATA, "GPL13667_platform.soft.tsv"), "GPL13667", TRUE)
map_n <- build_mapping(x_n_raw[, m_n$geo_accession, drop = FALSE],
                       file.path(DATA, "GPL15207_platform.soft.tsv"), "GPL15207", FALSE)
shared <- merge(map_a$gene_map, map_n$gene_map, by = "gene_symbol", suffixes = c("_a", "_n"))
shared[, entrez_agreement := entrez_id_a == entrez_id_n]
shared <- shared[entrez_agreement == TRUE]
setorder(shared, gene_symbol)
fwrite(shared, file.path(RESULTS, "shared_gene_universe.tsv"), sep = "\t")
audit <- rbindlist(list(map_a$audit, map_n$audit), fill = TRUE)
audit[, used_shared_universe := used_reliable & gene_symbol %in% shared$gene_symbol]
fwrite(audit, file.path(RESULTS, "platform_mapping_audit.tsv"), sep = "\t")

shared_genes <- shared$gene_symbol
x_a <- map_a$primary[shared_genes, , drop = FALSE]
x_n <- map_n$primary[shared_genes, , drop = FALSE]
x_a_alt <- map_a$alternative[shared_genes, , drop = FALSE]
x_n_alt <- map_n$alternative[shared_genes, , drop = FALSE]
z_a <- z_by_gene(x_a)
z_n <- z_by_gene(x_n)
z_a_alt <- z_by_gene(x_a_alt)
z_n_alt <- z_by_gene(x_n_alt)

# ---- Frozen Phase 1A programs ----
frozen <- fread(file.path(PROJECT, "program_definitions", "laser_program_genes.tsv"))
program_ids <- data.table(
  program = c("acute_injury", "inflammation", "epidermal_repair", "ecm_remodeling",
              "collagen_organization", "proliferation", "matrix_maturation"),
  source_database = c("MSigDB Hallmark", "MSigDB Hallmark", "MSigDB GO Biological Process",
                      "MSigDB Reactome", "MSigDB GO Biological Process", "MSigDB Hallmark",
                      "MSigDB GO Biological Process"),
  source_ID = c("HALLMARK_TNFA_SIGNALING_VIA_NFKB", "HALLMARK_INFLAMMATORY_RESPONSE",
                "GOBP_POSITIVE_REGULATION_OF_EPITHELIAL_CELL_PROLIFERATION_INVOLVED_IN_WOUND_HEALING",
                "REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION", "GOBP_COLLAGEN_FIBRIL_ORGANIZATION",
                "HALLMARK_E2F_TARGETS", "GOBP_EXTRACELLULAR_MATRIX_ASSEMBLY")
)
programs <- split(frozen$gene_symbol, frozen$program)
programs <- programs[program_ids$program]
programs_shared <- lapply(programs, intersect, y = shared_genes)
constitution <- program_ids[, .(program_name = program, source_database, source_ID)]
constitution[, original_gene_count := lengths(programs)]
constitution[, shared_universe_gene_count := lengths(programs_shared)]
constitution[, definition_status := "FROZEN_FROM_PHASE1A_NO_ADDITIONS_OR_DELETIONS"]
fwrite(constitution, file.path(RESULTS, "PROGRAM_CONSTITUTION.tsv"), sep = "\t")

# ---- Primary models and participant-level differences ----
s_a <- score_programs(z_a, programs_shared)
s_n <- score_programs(z_n, programs_shared)
sl_a <- make_score_long(s_a, m_a, "GSE168760_AFL")
sl_n <- make_score_long(s_n, m_n, "GSE206495_NAFL_FOREARM")
fwrite(rbind(sl_a, sl_n), file.path(RESULTS, "participant_program_scores.tsv"), sep = "\t")
fit_a <- fit_program_models(sl_a, "GSE168760_AFL")
fit_n <- fit_program_models(sl_n, "GSE206495_NAFL_FOREARM")
effects <- rbind(fit_a$effects, fit_n$effects)
differences <- rbind(fit_a$differences, fit_n$differences)
assumptions <- rbind(fit_a$assumptions, fit_n$assumptions)
fwrite(effects, file.path(RESULTS, "participant_blocked_program_effects.tsv"), sep = "\t")
fwrite(differences, file.path(RESULTS, "participant_program_differences.tsv"), sep = "\t")
fwrite(assumptions, file.path(RESULTS, "model_assumption_checks.tsv"), sep = "\t")

program_names <- program_ids$program
mats <- list()
cross_tests <- list()
cross_effects <- list()
for (dy in PRIMARY_DAYS) {
  ma <- diff_matrix(differences, "GSE168760_AFL", dy, program_names)
  mn <- diff_matrix(differences, "GSE206495_NAFL_FOREARM", dy, program_names)
  mats[[paste0("a", dy)]] <- ma
  mats[[paste0("n", dy)]] <- mn
  cross_tests[[length(cross_tests) + 1L]] <- cross_day_test(ma, mn, dy)
  cross_effects[[length(cross_effects) + 1L]] <- data.table(
    day = dy, program = program_names, afl_hedges_g = effect_from_diff_matrix(ma),
    nafl_hedges_g = effect_from_diff_matrix(mn), afl_mean_change = colMeans(ma),
    nafl_mean_change = colMeans(mn), same_direction = sign(effect_from_diff_matrix(ma)) == sign(effect_from_diff_matrix(mn))
  )
}
cross_summary <- rbindlist(cross_tests)
cross_effect <- rbindlist(cross_effects)
fwrite(cross_summary, file.path(RESULTS, "cross_cohort_conservation.tsv"), sep = "\t")
fwrite(cross_effect, file.path(RESULTS, "cross_cohort_program_effect_vectors.tsv"), sep = "\t")

# ---- Independent Phase 1A legacy recalculation (not copied from old result files) ----
legacy_map <- function(x, platform_path, log2_needed) {
  a <- fread(platform_path, select = c("ID", "Gene Symbol"), showProgress = FALSE)
  setnames(a, c("probe", "symbol"))
  a <- a[!is.na(symbol) & symbol != "---" & !grepl("///", symbol, fixed = TRUE)]
  a <- a[match(rownames(x), probe)]
  keep <- !is.na(a$symbol)
  xx <- x[keep, , drop = FALSE]
  if (log2_needed) xx <- log2(xx + 1)
  avereps(xx, ID = a$symbol[keep])
}
legacy_a <- legacy_map(x_a_raw[, m_a$geo_accession, drop = FALSE], file.path(DATA, "GPL13667_platform.soft.tsv"), TRUE)
legacy_n <- legacy_map(x_n_raw[, m_n$geo_accession, drop = FALSE], file.path(DATA, "GPL15207_platform.soft.tsv"), FALSE)
legacy_sa <- score_programs(z_by_gene(legacy_a), programs)
legacy_sn <- score_programs(z_by_gene(legacy_n), programs)
legacy_da <- fit_program_models(make_score_long(legacy_sa, m_a, "A"), "A")$differences
legacy_dn <- fit_program_models(make_score_long(legacy_sn, m_n, "N"), "N")$differences
expected <- c(`1` = 0.892857142857143, `7` = 0.857142857142857, `14` = 0.785714285714286)
legacy_audit <- rbindlist(lapply(PRIMARY_DAYS, function(dy) {
  ma <- diff_matrix(legacy_da, "A", dy, program_names)
  mn <- diff_matrix(legacy_dn, "N", dy, program_names)
  rho <- spearman_safe(colMeans(ma), colMeans(mn))
  data.table(day = dy, phase1a_reported_rho = expected[as.character(dy)], independently_recomputed_rho = rho,
             absolute_difference = abs(rho - expected[as.character(dy)]),
             method = "cohort-specific frozen-program raw mean changes; fresh calculation")
}))
fwrite(legacy_audit, file.path(RESULTS, "phase1a_recalculation_audit.tsv"), sep = "\t")

# ---- LOPO: re-estimate program vectors and cross-cohort correlations ----
lopo_program <- list()
lopo_cross <- list()
for (source in c("GSE168760_AFL", "GSE206495_NAFL_FOREARM")) {
  participants <- if (source == "GSE168760_AFL") rownames(mats$a1) else rownames(mats$n1)
  for (pid in participants) {
    for (dy in PRIMARY_DAYS) {
      ma <- mats[[paste0("a", dy)]]
      mn <- mats[[paste0("n", dy)]]
      if (source == "GSE168760_AFL") ma <- ma[rownames(ma) != pid, , drop = FALSE]
      if (source == "GSE206495_NAFL_FOREARM") mn <- mn[rownames(mn) != pid, , drop = FALSE]
      ea <- effect_from_diff_matrix(ma)
      en <- effect_from_diff_matrix(mn)
      lopo_cross[[length(lopo_cross) + 1L]] <- data.table(
        deleted_from = source, deleted_participant = pid, day = dy,
        standardized_spearman_rho = spearman_safe(ea, en),
        directional_concordance = mean(sign(ea) == sign(en)),
        all_program_effects_reestimated = TRUE
      )
      lopo_program[[length(lopo_program) + 1L]] <- data.table(
        deleted_from = source, deleted_participant = pid, day = dy, program = program_names,
        afl_hedges_g = ea, nafl_hedges_g = en,
        same_direction = sign(ea) == sign(en)
      )
    }
  }
}
lopo_cross <- rbindlist(lopo_cross)
lopo_program <- rbindlist(lopo_program)
fwrite(lopo_cross, file.path(RESULTS, "LOPO_temporal_conservation.tsv"), sep = "\t")
fwrite(lopo_program, file.path(RESULTS, "LOPO_program_effects.tsv"), sep = "\t")

# ---- Matched-random-program null on shared gene universe ----
perm_path <- file.path(RESULTS, "program_conservation_permutation.tsv")
if (file.exists(perm_path) && nrow(fread(perm_path, showProgress = FALSE)) == N_RANDOM + 1L) {
  perm_table <- fread(perm_path, showProgress = FALSE)
  random_null <- perm_table[is_real == FALSE]
  real_row <- perm_table[is_real == TRUE]
} else {
  random_null <- vector("list", N_RANDOM)
  for (iter in seq_len(N_RANDOM)) {
    random_programs <- lapply(lengths(programs_shared), function(k) sample(shared_genes, k, replace = FALSE))
    names(random_programs) <- program_names
    rsa <- score_programs(z_a, random_programs)
    rsn <- score_programs(z_n, random_programs)
    rda <- quick_diff_matrices(rsa, m_a, PRIMARY_DAYS)
    rdn <- quick_diff_matrices(rsn, m_n, PRIMARY_DAYS)
    rhos <- vapply(seq_along(PRIMARY_DAYS), function(i) {
      spearman_safe(effect_from_diff_matrix(rda[[i]]), effect_from_diff_matrix(rdn[[i]]))
    }, numeric(1))
    random_null[[iter]] <- data.table(iteration = iter, is_real = FALSE,
                                      day1_rho = rhos[1L], day7_rho = rhos[2L], day14_rho = rhos[3L],
                                      mean_rho = mean(rhos), min_rho = min(rhos))
  }
  random_null <- rbindlist(random_null)
  real_rhos <- cross_summary$standardized_spearman_rho
  real_row <- data.table(iteration = 0L, is_real = TRUE, day1_rho = real_rhos[1L], day7_rho = real_rhos[2L],
                         day14_rho = real_rhos[3L], mean_rho = mean(real_rhos), min_rho = min(real_rhos))
  perm_table <- rbind(real_row, random_null)
  fwrite(perm_table, perm_path, sep = "\t")
}
null_summary <- data.table(
  statistic = c("mean_rho_across_days", "minimum_rho_across_days", "day1_rho", "day7_rho", "day14_rho"),
  real_value = c(real_row$mean_rho, real_row$min_rho, real_row$day1_rho, real_row$day7_rho, real_row$day14_rho),
  null_mean = c(mean(random_null$mean_rho), mean(random_null$min_rho), mean(random_null$day1_rho),
                mean(random_null$day7_rho), mean(random_null$day14_rho)),
  null_q95 = c(quantile(random_null$mean_rho, .95), quantile(random_null$min_rho, .95),
               quantile(random_null$day1_rho, .95), quantile(random_null$day7_rho, .95),
               quantile(random_null$day14_rho, .95)),
  empirical_one_sided_p = c((1 + sum(random_null$mean_rho >= real_row$mean_rho)) / (N_RANDOM + 1),
                            (1 + sum(random_null$min_rho >= real_row$min_rho)) / (N_RANDOM + 1),
                            (1 + sum(random_null$day1_rho >= real_row$day1_rho)) / (N_RANDOM + 1),
                            (1 + sum(random_null$day7_rho >= real_row$day7_rho)) / (N_RANDOM + 1),
                            (1 + sum(random_null$day14_rho >= real_row$day14_rho)) / (N_RANDOM + 1)),
  permutations = N_RANDOM
)
fwrite(null_summary, file.path(RESULTS, "program_conservation_null_summary.tsv"), sep = "\t")

# ---- Collagen and proliferation sensitivities ----
collagen_genes <- programs_shared$collagen_organization
collagen_variants <- list(
  full = collagen_genes,
  minus_COL1A1 = setdiff(collagen_genes, "COL1A1"),
  minus_COL1A2 = setdiff(collagen_genes, "COL1A2"),
  minus_COL1A1_COL1A2 = setdiff(collagen_genes, c("COL1A1", "COL1A2")),
  minus_all_COL_family = collagen_genes[!grepl("^COL[0-9]", collagen_genes)]
)
collagen_sens <- rbindlist(lapply(names(collagen_variants), function(v) {
  ga <- collagen_variants[[v]]
  ea <- program_variant_effect(z_a, m_a, ga, "collagen_organization", "GSE168760_AFL", 14L)
  en <- program_variant_effect(z_n, m_n, ga, "collagen_organization", "GSE206495_NAFL_FOREARM", 14L)
  rbind(data.table(variant = v, genes_used = length(ga), ea), data.table(variant = v, genes_used = length(ga), en), fill = TRUE)
}))
fwrite(collagen_sens, file.path(RESULTS, "collagen_gene_family_sensitivity.tsv"), sep = "\t")

collagen_loo <- rbindlist(lapply(collagen_genes, function(g) {
  keep <- setdiff(collagen_genes, g)
  ea <- program_variant_effect(z_a, m_a, keep, "collagen_organization", "GSE168760_AFL", 14L)
  en <- program_variant_effect(z_n, m_n, keep, "collagen_organization", "GSE206495_NAFL_FOREARM", 14L)
  rbind(data.table(removed_gene = g, ea), data.table(removed_gene = g, en), fill = TRUE)
}))
fwrite(collagen_loo, file.path(RESULTS, "collagen_leave_one_gene_out.tsv"), sep = "\t")

gene_day14_effect <- function(z, meta, genes, dataset) {
  out <- rbindlist(lapply(genes, function(g) {
    sc <- matrix(z[g, ], nrow = 1L, dimnames = list(g, colnames(z)))
    d <- quick_diff_matrices(sc, meta, 14L)[[1L]][, 1L]
    data.table(dataset = dataset, gene_symbol = g, mean_gene_change = mean(d))
  }))
  out[, absolute_contribution_fraction := abs(mean_gene_change) / sum(abs(mean_gene_change))]
  out
}
collagen_influence <- rbind(gene_day14_effect(z_a, m_a, collagen_genes, "GSE168760_AFL"),
                            gene_day14_effect(z_n, m_n, collagen_genes, "GSE206495_NAFL_FOREARM"))
fwrite(collagen_influence, file.path(RESULTS, "collagen_gene_influence.tsv"), sep = "\t")

prolif_genes <- programs_shared$proliferation
prolif_loo <- rbindlist(lapply(prolif_genes, function(g) {
  keep <- setdiff(prolif_genes, g)
  ea <- program_variant_effect(z_a, m_a, keep, "proliferation", "GSE168760_AFL", 14L)
  en <- program_variant_effect(z_n, m_n, keep, "proliferation", "GSE206495_NAFL_FOREARM", 14L)
  rbind(data.table(removed_gene = g, ea), data.table(removed_gene = g, en), fill = TRUE)
}))
fwrite(prolif_loo, file.path(RESULTS, "proliferation_leave_one_gene_out.tsv"), sep = "\t")

program_halves <- split(sort(prolif_genes), rep(c("alphabetic_half_1", "alphabetic_half_2"), length.out = length(prolif_genes)))
prolif_halves <- rbindlist(lapply(names(program_halves), function(v) {
  keep <- program_halves[[v]]
  ea <- program_variant_effect(z_a, m_a, keep, "proliferation", "GSE168760_AFL", 14L)
  en <- program_variant_effect(z_n, m_n, keep, "proliferation", "GSE206495_NAFL_FOREARM", 14L)
  rbind(data.table(variant = v, genes_used = length(keep), ea), data.table(variant = v, genes_used = length(keep), en), fill = TRUE)
}))
fwrite(prolif_halves, file.path(RESULTS, "proliferation_gene_set_sensitivity.tsv"), sep = "\t")

s_a_alt <- score_programs(z_a_alt, programs_shared)
s_n_alt <- score_programs(z_n_alt, programs_shared)
alt_a <- fit_program_models(make_score_long(s_a_alt, m_a, "GSE168760_AFL"), "GSE168760_AFL")$effects
alt_n <- fit_program_models(make_score_long(s_n_alt, m_n, "GSE206495_NAFL_FOREARM"), "GSE206495_NAFL_FOREARM")$effects
mapping_sens <- rbind(alt_a, alt_n)[day %in% PRIMARY_DAYS]
mapping_sens[, mapping_method := "highest_IQR_single_probe_per_gene"]
primary_sens <- effects[day %in% PRIMARY_DAYS]
primary_sens[, mapping_method := "mean_all_reliable_probes_per_gene"]
mapping_sens <- rbind(primary_sens, mapping_sens, fill = TRUE)
fwrite(mapping_sens, file.path(RESULTS, "platform_mapping_sensitivity.tsv"), sep = "\t")

# Conservation must not depend on the secondary proliferation pattern.
leave_program_out <- rbindlist(lapply(program_names, function(drop_program) {
  cross_effect[program != drop_program, .(
    standardized_spearman_rho = spearman_safe(afl_hedges_g, nafl_hedges_g),
    directional_concordance = mean(sign(afl_hedges_g) == sign(nafl_hedges_g))
  ), by = day][, omitted_program := drop_program]
}))
setcolorder(leave_program_out, c("omitted_program", "day", "standardized_spearman_rho", "directional_concordance"))
fwrite(leave_program_out, file.path(RESULTS, "cross_cohort_leave_program_out.tsv"), sep = "\t")
proliferation_independence <- lopo_program[program != "proliferation", .(
  standardized_spearman_rho = spearman_safe(afl_hedges_g, nafl_hedges_g),
  directional_concordance = mean(sign(afl_hedges_g) == sign(nafl_hedges_g))
), by = .(deleted_from, deleted_participant, day)]
fwrite(proliferation_independence, file.path(RESULTS, "proliferation_independence_LOPO.tsv"), sep = "\t")

# ---- Deterministic program classification ----
class_rows <- list()
for (pr in program_names) {
  ce <- cross_effect[program == pr]
  ea <- effects[dataset == "GSE168760_AFL" & program == pr & day %in% PRIMARY_DAYS]
  en <- effects[dataset == "GSE206495_NAFL_FOREARM" & program == pr & day %in% PRIMARY_DAYS]
  setorder(ea, day); setorder(en, day); setorder(ce, day)
  same_n <- sum(ce$same_direction)
  peak_a <- ea$day[which.max(abs(ea$paired_hedges_g))]
  peak_n <- en$day[which.max(abs(en$paired_hedges_g))]
  opposite_strong <- any(sign(ea$paired_hedges_g) != sign(en$paired_hedges_g) &
                           ea$model_ci95_low * ea$model_ci95_high > 0 & en$model_ci95_low * en$model_ci95_high > 0 &
                           ea$lopo_direction_fraction == 1 & en$lopo_direction_fraction == 1)
  if (opposite_strong) {
    classification <- "DIVERGENT"
  } else if (same_n == 3L && peak_a == peak_n &&
             ea$lopo_direction_fraction[ea$day == peak_a] == 1 & en$lopo_direction_fraction[en$day == peak_n] == 1) {
    classification <- "CONSERVED"
  } else if (same_n >= 2L && peak_a == peak_n) {
    classification <- "PARTIALLY_CONSERVED"
  } else {
    classification <- "UNRESOLVED"
  }
  class_rows[[length(class_rows) + 1L]] <- data.table(
    program = pr, classification = classification, same_direction_days = same_n,
    afl_peak_day = peak_a, nafl_peak_day = peak_n,
    classification_rule = "DIVERGENT requires opposite significant stable effects; CONSERVED requires 3/3 same directions, same peak and stable peak; PARTIAL requires >=2/3 and same peak"
  )
}
classification <- rbindlist(class_rows)
fwrite(classification, file.path(RESULTS, "program_temporal_classification.tsv"), sep = "\t")

# ---- Exploratory, non-submission figures (raw/derived data retained in TSVs) ----
cols <- c("GSE168760_AFL" = "#0072B2", "GSE206495_NAFL_FOREARM" = "#D55E00")
png(file.path(FIGURES, "01_exploratory_temporal_trajectories.png"), width = 1800, height = 1500, res = 180)
par(mfrow = c(3, 3), mar = c(4, 4, 2.2, 1), oma = c(0, 0, 2, 0))
for (pr in program_names) {
  da <- effects[dataset == "GSE168760_AFL" & program == pr]
  dn <- effects[dataset == "GSE206495_NAFL_FOREARM" & program == pr]
  da <- rbind(data.table(day = 0L, paired_hedges_g = 0), da[, .(day, paired_hedges_g)])
  dn <- rbind(data.table(day = 0L, paired_hedges_g = 0), dn[, .(day, paired_hedges_g)])
  ylim <- range(c(da$paired_hedges_g, dn$paired_hedges_g), finite = TRUE)
  plot(da$day, da$paired_hedges_g, type = "b", pch = 16, col = cols[1], lty = 1,
       xlab = "Day", ylab = "Paired Hedges g", main = pr, ylim = ylim)
  lines(dn$day, dn$paired_hedges_g, type = "b", pch = 17, col = cols[2], lty = 2)
  abline(h = 0, col = "grey70", lty = 3)
}
plot.new(); legend("center", legend = c("GSE168760 AFL", "GSE206495 NAFL forearm"),
                   col = cols, pch = c(16, 17), lty = c(1, 2), bty = "n")
mtext("Exploratory participant-blocked program trajectories", outer = TRUE, cex = 1.2)
dev.off()

png(file.path(FIGURES, "02_exploratory_cross_cohort_effects.png"), width = 1800, height = 650, res = 180)
par(mfrow = c(1, 3), mar = c(4.5, 4.5, 3, 1))
for (dy in PRIMARY_DAYS) {
  d <- cross_effect[day == dy]
  max_abs <- max(abs(c(d$afl_hedges_g, d$nafl_hedges_g))) * 1.22
  lim <- c(-max_abs, max_abs)
  plot(d$afl_hedges_g, d$nafl_hedges_g, pch = 16, col = "#009E73", xlim = lim, ylim = lim,
       xlab = "AFL paired Hedges g", ylab = "NAFL paired Hedges g",
       main = paste0("Day ", dy, "; rho=", sprintf("%.2f", cross_summary[day == dy, standardized_spearman_rho])))
  abline(h = 0, v = 0, col = "grey70", lty = 3); abline(0, 1, col = "grey50", lty = 2)
  text(d$afl_hedges_g, d$nafl_hedges_g, labels = d$program, pos = 3, cex = .58, xpd = NA)
}
dev.off()

png(file.path(FIGURES, "03_exploratory_collagen_participant_distribution.png"), width = 1400, height = 650, res = 180)
par(mfrow = c(1, 2), mar = c(5, 4.5, 3, 1))
for (ds in c("GSE168760_AFL", "GSE206495_NAFL_FOREARM")) {
  d <- differences[dataset == ds & program == "collagen_organization" & day == 14L]
  plot(seq_len(nrow(d)), d$difference, pch = 16, col = if (ds == "GSE168760_AFL") cols[1] else cols[2],
       xlab = "Participant (ordered)", ylab = "Day14 - baseline score", main = ds, xaxt = "n")
  axis(1, seq_len(nrow(d)), labels = d$participant_id, las = 2, cex.axis = .7)
  abline(h = 0, col = "grey60", lty = 2); abline(h = mean(d$difference), col = "black", lwd = 2)
}
dev.off()

png(file.path(FIGURES, "04_exploratory_random_program_null.png"), width = 1100, height = 750, res = 180)
hist(random_null$mean_rho, breaks = 40, col = "grey80", border = "white",
     xlab = "Mean Spearman rho across day 1/7/14", main = "Gene-count-matched random-program null")
abline(v = real_row$mean_rho, col = "#D55E00", lwd = 3)
legend("topleft", legend = paste0("Real = ", sprintf("%.3f", real_row$mean_rho)), col = "#D55E00", lwd = 3, bty = "n")
dev.off()

# ---- QC summary and reproducibility ----
qc <- data.table(
  check = c("A_participants", "N_participants", "shared_reliable_gene_universe", "programs_frozen",
            "primary_days", "bootstrap_replicates", "random_program_permutations", "raw_data_modified",
            "GSE315794_retested", "formal_figure_created"),
  observed = c(uniqueN(m_a$participant_id), uniqueN(m_n$participant_id), length(shared_genes), length(programs_shared),
               paste(PRIMARY_DAYS, collapse = ","), N_BOOT, N_RANDOM, "NO", "NO", "NO"),
  pass = c(TRUE, TRUE, length(shared_genes) > 10000, length(programs_shared) == 7L, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE)
)
fwrite(qc, file.path(RESULTS, "PHASE1B_QC.tsv"), sep = "\t")
capture.output(sessionInfo(), file = file.path(LOGS, "phase1b_sessionInfo.txt"))
