###############################################################
# Linear regression plots 
###############################################################

library(readxl)
library(ggplot2)
library(gridExtra)

# set working directory to the directory the script is saved in
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

###############################################################
# 1. Read data
###############################################################

cn_df <- read_xlsx("Biological Data file Tilstra et al.xlsx", sheet = "bio + env data", na = "NA")

print(colnames(cn_df))   # Just to verify the names

# subset data frame to the 4 columns we need

cn_df<-cn_df[,c("C host tissue","N host tissue","C algal symbiont","N algal symbiont")]

###############################################################
# 2. Helper function to build one regression plot
###############################################################

make_reg_plot <- function(cn_df, xcol, ycol,
                          x_label, y_label,
                          point_fill, point_shape,
                          line_color, line_fill) {
  
  # Fit linear model
  lm_fit <- lm(cn_df[[ycol]] ~ cn_df[[xcol]])
  sm     <- summary(lm_fit)
  
  r2_val <- sm$r.squared
  p_val  <- sm$coefficients[2, 4]
  
  # Compute annotation positions (near the top of the panel)
  x_range <- range(cn_df[[xcol]], na.rm = TRUE)
  y_range <- range(cn_df[[ycol]], na.rm = TRUE)
  
  x_pos    <- x_range[1] + 0.90 * diff(x_range)  # 70% from left
  y_max    <- y_range[2]
  y_span   <- diff(y_range)
  y_pos_r2 <- y_range[1] + 0.14  * y_span   # 25% from bottom
  y_pos_p  <- y_range[1] + 0.07 * y_span   # 12% from bottom
  
  ggplot(cn_df, aes(x = .data[[xcol]], y = .data[[ycol]])) +
    geom_point(color = "black", fill = point_fill,
               shape = point_shape, size = 3, stroke = 1) +
    geom_smooth(method = "lm", se = TRUE,
                color = line_color, fill = line_fill) +
    labs(
      x = x_label,
      y = y_label
    ) +
    # R² (normal, using R^2 syntax)
    annotate("text",
             x = x_pos, y = y_pos_r2,
             label = paste0("R^2 == ", round(r2_val, 3)),
             parse = TRUE,
             hjust = 1, vjust = 1,
             size = 5) +
    # p (italic p, numeric value normal)
    annotate("text",
             x = x_pos, y = y_pos_p,
             label = paste0("italic(p) == ", signif(p_val, 3)),
             parse = TRUE,
             hjust = 1, vjust = 1,
             size = 5) +
    theme_minimal() +
    theme(
      axis.text  = element_text(color = "black", size = 16),
      axis.title = element_text(size = 16),
      axis.ticks = element_line(color = "black"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 2)
    )
}

###############################################################
# 3. Build both plots using Excel names
###############################################################

plot_host <- make_reg_plot(
  cn_df          = cn_df,
  xcol        = "N host tissue",
  ycol        = "C host tissue",
  x_label     = "N host tissue (%)",
  y_label     = "C host tissue (%)",
  point_fill  = "sandybrown",
  point_shape = 21,
  line_color  = "sandybrown",
  line_fill   = "lightgoldenrod3"
)

plot_sym <- make_reg_plot(
  cn_df          = cn_df,
  xcol        = "N algal symbiont",
  ycol        = "C algal symbiont",
  x_label     = "N algal symbiont (%)",
  y_label     = "C algal symbiont (%)",
  point_fill  = "darkolivegreen4",
  point_shape = 24,
  line_color  = "darkolivegreen4",
  line_fill   = "mediumseagreen"
)

###############################################################
# 4. Print regression summaries to console
###############################################################

cat("\nHost model: C host tissue ~ N host tissue\n")
print(summary(lm(`C host tissue` ~ `N host tissue`, data = cn_df)))

cat("\nAlgal symbiont model: C algal symbiont ~ N algal symbiont\n")
print(summary(lm(`C algal symbiont` ~ `N algal symbiont`, data = cn_df)))

###############################################################
# 5. Show plots and save to Pdf
###############################################################

print(plot_host)
print(plot_sym)

fig4b<-grid.arrange(plot_host, plot_sym, ncol = 2)

ggsave(
  "Figure_4B.pdf",
  fig4b,
  width =8,
  height = 6,
  dpi = 300
)