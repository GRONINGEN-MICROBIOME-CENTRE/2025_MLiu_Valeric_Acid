suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(emmeans)
})

## -----------------------------------------------------------------------------
## File: 05_lm_intx_0908_emm_unified.R
## Goal: Refit the Cytokines × VA concentration models for each target gene,
##       compute estimated marginal means (EMMs), and save plots/tables that
##       share the same color palette, title, and subtitle formatting as the
##       inflame-score plots shown in the reference code snippet.
## -----------------------------------------------------------------------------

project_dir <- "~/OneDrive - UMCG/2025_09_02/R_Works/Valeric_acid"
data_dir <- file.path(project_dir, "results/Orga_bulkrna")
code_dir <- file.path(project_dir, "codes/VA_0810")
setwd(code_dir)

counts_path <- file.path(data_dir, "normalized_counts_clean_logcpm.csv")
metadata_path <- file.path(data_dir, "metadata_clean.csv")
gene_list_path <- file.path(project_dir, "data/Orga_bulk/VA_related_genes_for_organoids_0810.xlsx")

stopifnot(file.exists(counts_path), file.exists(metadata_path), file.exists(gene_list_path))

counts <- read.csv(counts_path, header = TRUE, row.names = 1, check.names = FALSE)
rownames(counts) <- gsub("\\..*", "", rownames(counts))
metadata <- read.csv(metadata_path, header = TRUE, row.names = 1, check.names = FALSE)
stopifnot(all(colnames(counts) %in% rownames(metadata)))

gene_sheet <- read_excel(gene_list_path)
valid_genes <- gene_sheet$Genes[gene_sheet$Ensembl_ID %in% rownames(counts)]
counts_subset <- counts[gene_sheet$Ensembl_ID[gene_sheet$Genes %in% valid_genes], , drop = FALSE]
rownames(counts_subset) <- gene_sheet$Genes[gene_sheet$Genes %in% valid_genes]

expr_long <- counts_subset %>%
  rownames_to_column("Gene") %>%
  pivot_longer(-Gene, names_to = "Sample", values_to = "Expression") %>%
  left_join(metadata %>% rownames_to_column("Sample"), by = "Sample") %>%
  mutate(
    Cytokines = factor(Cytokines, levels = c("No", "Yes")),
    Pretreat = factor(Pretreat, levels = c("No", "Yes")),
    Colon_organoids_cell_line = factor(Colon_organoids_cell_line),
    Timepoint = factor(Timepoint),
    Conc_num = as.numeric(as.character(Concentration)),
    Conc_group = factor(Concentration)
  )

va_levels <- c(0, 0.3, 1, 3, 10)
va_colors <- c("No" = "#4C72B0", "Yes" = "#DD8452")

output_root <- file.path(project_dir, "results/Orga_bulkrna/VA_0810/06_lm_intx/04_emm_consistent")
plot_dir <- file.path(output_root, "plots")
emm_dir <- file.path(output_root, "emm_tables")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(emm_dir, recursive = TRUE, showWarnings = FALSE)

safe_name <- function(x) gsub("[^A-Za-z0-9_-]", "_", x)

format_p <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 0.001) return("p<0.001")
  signif(p, 3)
}

fit_gene <- function(data_gene, gene_name) {
  if (nrow(data_gene) < 5) {
    warning(sprintf("Skipping %s: insufficient observations (%s rows).", gene_name, nrow(data_gene)))
    return(NULL)
  }
  
  model_formula <- Expression ~ Cytokines * Conc_num + Pretreat + Colon_organoids_cell_line + Timepoint
  fit <- lm(model_formula, data = data_gene)
  coef_table <- summary(fit)$coefficients
  
  get_p <- function(term) {
    if (term %in% rownames(coef_table)) signif(coef_table[term, 4], 3) else NA
  }
  
  p_ck <- get_p("CytokinesYes")
  p_va <- get_p("Conc_num")
  p_intx <- get_p("CytokinesYes:Conc_num")
  
  emm <- emmeans(fit, ~ Cytokines * Conc_num, at = list(Conc_num = va_levels))
  emm_df <- as.data.frame(emm) %>%
    mutate(
      Gene = gene_name,
      Concentration_mM = as.numeric(as.character(Conc_num)),
      Cytokines = factor(Cytokines, levels = c("No", "Yes"))
    )
  
  plot_df <- emm_df %>%
    mutate(Conc_factor = factor(Concentration_mM, levels = va_levels))
  
  plot_title <- paste(gene_name, "Estimated Marginal Means")
  plot_subtitle <- paste0(
    "P(CK)=", format_p(p_ck),
    "; P(VA)=", format_p(p_va),
    "; P(Intx)=", format_p(p_intx)
  )
  
  emm_plot <- ggplot(plot_df, aes(x = Conc_factor, y = emmean, group = Cytokines, color = Cytokines)) +
    geom_point(size = 3) +
    geom_line(linewidth = 1) +
    geom_errorbar(aes(ymin = emmean - SE, ymax = emmean + SE), width = 0.2) +
    scale_color_manual(values = va_colors, drop = FALSE) +
    labs(
      title = plot_title,
      subtitle = plot_subtitle,
      x = "VA Concentration (mM)",
      y = "Adjusted Expression",
      color = "Cytokines"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5),
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 12),
      legend.position = "right"
    )
  
  out_name <- safe_name(gene_name)
  ggsave(file.path(plot_dir, paste0(out_name, "_EMM_plot.pdf")), emm_plot, width = 6.5, height = 4.5)
  write_csv(emm_df, file.path(emm_dir, paste0(out_name, "_EMM_values.csv")))
  
  list(
    summary = tibble(
      Gene = gene_name,
      p_ck = p_ck,
      p_va = p_va,
      p_intx = p_intx,
      model = deparse(model_formula)
    ),
    emm = emm_df
  )
}

run_all_genes <- function(expr_df, genes) {
  summary_list <- list()
  emm_list <- list()
  for (gene in genes) {
    message("Processing: ", gene)
    gene_df <- expr_df %>% filter(Gene == gene)
    res <- tryCatch(
      fit_gene(gene_df, gene),
      error = function(e) {
        warning(sprintf("Failed for %s: %s", gene, e$message))
        NULL
      }
    )
    if (!is.null(res)) {
      summary_list[[gene]] <- res$summary
      emm_list[[gene]] <- res$emm
    }
  }
  summary_tbl <- if (length(summary_list)) bind_rows(summary_list) else tibble()
  emm_tbl <- if (length(emm_list)) bind_rows(emm_list) else tibble()
  
  list(
    summary = summary_tbl,
    emm = emm_tbl
  )
}

lm_results <- run_all_genes(expr_long, valid_genes)
lm_summary <- lm_results$summary
emm_table <- lm_results$emm

summary_path <- file.path(output_root, "VA_interaction_lm_summary_consistent.csv")
write_csv(lm_summary, summary_path)

emm_path <- file.path(output_root, "VA_emm_estimates_consistent.csv")
write_csv(emm_table, emm_path)

message("Finished. Summary table written to: ", summary_path)
message("Combined EMM estimates written to: ", emm_path)
