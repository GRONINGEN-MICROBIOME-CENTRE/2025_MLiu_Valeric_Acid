
library(readxl)
library(tidyverse)
library(rstatix)
library(viridis) 
library(ggpubr)
library(magrittr)

library(emmeans)
library(ggplot2)
library(dplyr)
library(ggrepel)

setwd('/Users/liumoting/R_Works/Valeric_acid/codes/')
countsrdyInf_organoids <- read.csv("../results/Inflame_score/countsrdyInf_organoids.csv", header = TRUE, row.names = 1)
df <- countsrdyInf_organoids


score_vars <- c("bMIS.IBD", "bMIS.CD", "bMIS.UC", "oxScore2007", "oxKeNf", "Barrier.score.Cld", "Barrier.score.noCld")


df$Cytokines <- factor(df$Cytokines, levels = c("No", "Yes"))
df$Pretreat <- factor(df$Pretreat, levels = c("No", "Yes"))
df$Colon_organoids_cell_line <- factor(df$Colon_organoids_cell_line)
df$Timepoint <- factor(df$Timepoint)
df$Concentration <- as.numeric(df$Concentration)


run_lm_and_emm_for_scores <- function(df, score_list, output_dir = ".", 
                                      va_levels = c(0, 0.3, 1, 3, 10),
                                      colors = c("No" = "#1f77b4", "Yes" = "#d62728")) {
  results_all <- list()
  
  for (score in score_list) {
    message("Processing: ", score)
    
    tryCatch({
      
      formula <- as.formula(paste(score, "~ Cytokines * Concentration + Pretreat + Colon_organoids_cell_line + Timepoint"))
      fit <- lm(formula, data = df)
      
      
      coef_table <- summary(fit)$coefficients
      get_p <- function(term) {
        if (term %in% rownames(coef_table)) signif(coef_table[term, 4], 3) else NA
      }
      p_ck <- get_p("CytokinesYes")
      p_va <- get_p("Concentration")
      p_intx <- get_p("CytokinesYes:Concentration")
      
      
      emm <- emmeans(fit, ~ Cytokines * Concentration, at = list(Concentration = va_levels))
      emm_df <- as.data.frame(emm)
      
      
      p <- ggplot(emm_df, aes(x = as.factor(Concentration), y = emmean, 
                              group = Cytokines, color = Cytokines)) +
        geom_point(size = 3) +
        geom_line(linewidth = 1) +
        geom_errorbar(aes(ymin = emmean - SE, ymax = emmean + SE), width = 0.2) +
        scale_color_manual(values = colors) +
        labs(
          title = paste(score, "Estimated Marginal Means"),
          subtitle = paste0("P(CK)=", p_ck, "; P(VA)=", p_va, "; P(Intx)=", p_intx),
          x = "VA Concentration (mM)",
          y = "Adjusted ", 
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
      
      
      file_out <- file.path(output_dir, paste0(score, "_EMM_plot.pdf"))
      ggsave(file_out, plot = p, width = 6.5, height = 4.5)
      
      
      results_all[[score]] <- data.frame(
        Score = score,
        p_ck = p_ck,
        p_va = p_va,
        p_intx = p_intx,
        stringsAsFactors = FALSE
      )
      
    }, error = function(e) {
      warning("Failed: ", score, " - ", e$message)
    })
  }
  
  
  df_pvals <- bind_rows(results_all)
  write.csv(df_pvals, file = file.path(output_dir, "LM_MainEffects_Pvalues.csv"), row.names = FALSE)
  
  return(df_pvals)
}

my_colors <- c("No" = "#4C72B0", "Yes" = "#DD8452") 
run_lm_and_emm_for_scores(df, score_list = score_vars, output_dir = "../results/Inflame_score/emms/plots/", colors = my_colors)


em_contrasts <- contrast(em, method = "revpairwise", by = "Cytokines", adjust = "fdr")
df_contrasts <- as.data.frame(em_contrasts)
df_contrasts$Score <- score
results_list[[score]] <- df_contrasts
df_all_stats <- bind_rows(results_list)
write.csv(df_all_stats, "../results/Inflame_score/emms/All_Stats_Contrasts_ByCytokines.csv", row.names = FALSE)
