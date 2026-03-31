
###############################################################
# Figure 3B: LINEAR REGRESSIONS with FIXED NPC LABEL POSITIONS
# R² and italic(p) placed LOW RIGHT inside each panel
###############################################################

library(readxl)
library(tidyverse)
library(grid)
library(gridExtra)

# set working directory to the directory the script is saved in
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


###############################################################
# 1. READ DATA
###############################################################

data <- read_xlsx("Biological Data file Tilstra et al.xlsx", sheet = "bio + env data", na = "NA")

# add % unit to the relevant colnames

colnames(data)[3:4]<-c("C host tissue (%)","N host tissue (%)")

# define variables which will be used for the plots

biological_vars <- c(
  "C host tissue (%)",
  "N host tissue (%)",
  "C:N ratio host tissue"
)

environmental_var <- "PAR"

###############################################################
# 2. MAKE RESULTS TABLE
###############################################################

results <- data.frame(
  Biological_Variable = character(),
  R2 = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

plots <- list()

###############################################################
# 3. LOOP THROUGH VARIABLES
###############################################################

for (biological_var in biological_vars) {
  
  summary_data <- data %>%
    group_by(.data[[environmental_var]]) %>%
    summarise(
      mean_value = mean(.data[[biological_var]], na.rm = TRUE),
      sd_value   = sd(.data[[biological_var]],   na.rm = TRUE),
      .groups = "drop"
    )
  
  lm_data <- data %>%
    select(
      bio = all_of(biological_var),
      env = all_of(environmental_var)
    ) %>%
    filter(!is.na(bio), !is.na(env))
  
  lm_fit <- lm(bio ~ env, data = lm_data)
  lm_sum <- summary(lm_fit)
  
  r2_val <- lm_sum$r.squared
  p_val  <- lm_sum$coefficients[2, 4]
  
  results <- rbind(
    results,
    data.frame(
      Biological_Variable = biological_var,
      R2 = round(r2_val, 3),
      p_value = signif(p_val, 3)
    )
  )
  
  ###############################################################
  # 4. BASE PLOT
  ###############################################################
  
  p <- ggplot(summary_data,
              aes(x = .data[[environmental_var]], y = mean_value)) +
    geom_errorbar(aes(ymin = mean_value - sd_value,
                      ymax = mean_value + sd_value),
                  width = 0.1, color = "sandybrown", linewidth = 1) +
    geom_point(size = 4, color = "sandybrown") +
    geom_smooth(method = "lm", se = TRUE,
                color = "sandybrown", fill = "lightgoldenrod3") +
    labs(x = NULL, y = biological_var) +
    theme_minimal() +
    theme(
      axis.text  = element_text(color = "black", size = 16),
      axis.title = element_text(size = 16),
      axis.ticks = element_line(color = "black"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
    )
  
  ###############################################################
  # 5. FIXED NPC (NORMALIZED) ANNOTATIONS
  ###############################################################
  
  # Build a grob for R² (normal)
  r2_grob <- textGrob(
    label = paste0("R² = ", round(r2_val, 3)),
    x = unit(0.90, "npc"),   # 90% from left
    y = unit(0.15, "npc"),   # 20% from bottom
    hjust = 1,
    gp = gpar(fontsize = 12)
  )
  
  # Build a grob for italic(p)
  p_grob <- textGrob(
    label = bquote(italic(p) == .(signif(p_val, 3))),
    x = unit(0.90, "npc"),
    y = unit(0.10, "npc"),
    hjust = 1,
    gp = gpar(fontsize = 12)
  )
  
  # Add grobs to the plot
  p <- p +
    annotation_custom(r2_grob) +
    annotation_custom(p_grob)
  
  plots[[biological_var]] <- p
}

###############################################################
# 6. PRINT RESULTS
###############################################################

cat("\n================ LINEAR REGRESSION SUMMARY ================\n")
print(results)

###############################################################
# 7. SAVE TO PDF WITH SHARED X-LABEL
###############################################################

xlab_grob <- textGrob(
  label = expression("PAR ("*mu*"mol photons m"^{-2}*" s"^{-1}*")"),
  gp = gpar(fontsize = 16)
)

fig3b<-grid.arrange(grobs = plots, ncol = 3, bottom = xlab_grob)

ggsave(
  "Figure_3B.pdf",
  fig3b,
  width =12,
  height = 6,
  dpi = 300
)

# In image processing software outside of R, the subplot letter were added and the plots were combined with the other plots of Figure 3.
