###############################################################
# Combined Pearson/Spearman correlation heatmap with Arial
# - Data: Correlations_3Envs.xlsx
# - Env vars: last 3 columns (e.g. Temperature, PAR, DIN)
# - Bio vars: all other columns (in Excel order)
# - Row-wise method:
#     * Pearson for bio vars with Best != "none-normal"
#     * Spearman for bio vars with Best == "none-normal"
# - Pearson bio vars transformed: raw / log(x+1) / log(-x+1)
# - Spearman bio vars & env vars: raw
# - Bonferroni correction: ON
# - Font: Arial via showtext
###############################################################

########################
# 0. Font setup (Arial)
########################

library(showtext)

# Adjust this path accordingly

arial_path <- "C:/Windows/Fonts/Arial.ttf"

if (!file.exists(arial_path)) {
  stop("Arial.ttf not found at: ", arial_path)
}

font_add(family = "Arial", regular = arial_path)
showtext_auto()

########################
# 1. Packages & data
########################

library(readxl)
library(dplyr)
library(ggplot2)

# set working directory to the directory the script is saved in

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Read Excel file (must be in working directory)
df <- read_xlsx("Biological and Environmental Correlations input data.xlsx")

# Last 3 columns are environmental variables
env_vars <- tail(colnames(df), 3)

# Biological variables are all others, in Excel order
bio_vars <- colnames(df)[!(colnames(df) %in% env_vars)]

cat("Environmental variables:\n")
print(env_vars)
cat("\nBiological variables (Excel order):\n")
print(bio_vars)

bio_data <- df[, bio_vars, drop = FALSE]
env_data <- df[, env_vars, drop = FALSE]

###############################################################
# 2. Normality check + best transformation for each BIO variable
###############################################################

check_normality <- function(x) {
  z <- x[!is.na(x)]
  if (length(z) < 3) return(c(NA, NA, NA, "none-normal"))
  
  # raw
  raw_p <- tryCatch(shapiro.test(z)$p.value, error = function(e) NA)
  
  # log(x+1)
  if (all(z > -1)) {
    log_pos   <- log(z + 1)
    log_pos_p <- tryCatch(shapiro.test(log_pos)$p.value, error = function(e) NA)
  } else {
    log_pos_p <- NA
  }
  
  # log(-x+1)
  if (all(z <= 0)) {
    log_neg   <- log(-z + 1)
    log_neg_p <- tryCatch(shapiro.test(log_neg)$p.value, error = function(e) NA)
  } else {
    log_neg_p <- NA
  }
  
  best <- "none-normal"
  if (!is.na(raw_p)     && raw_p     > 0.05) best <- "raw"
  if (!is.na(log_pos_p) && log_pos_p > 0.05) best <- "log(x+1)"
  if (!is.na(log_neg_p) && log_neg_p > 0.05) best <- "log(-x+1)"
  
  c(raw_p, log_pos_p, log_neg_p, best)
}

norm_res <- data.frame(
  Variable = bio_vars,
  Raw_p    = NA_real_,
  LogPos_p = NA_real_,
  LogNeg_p = NA_real_,
  Best     = NA_character_
)

for (i in seq_along(bio_vars)) {
  r <- check_normality(bio_data[[bio_vars[i]]])
  norm_res$Raw_p[i]    <- as.numeric(r[1])
  norm_res$LogPos_p[i] <- as.numeric(r[2])
  norm_res$LogNeg_p[i] <- as.numeric(r[3])
  norm_res$Best[i]     <- r[4]
}

cat("\nNormality summary:\n")
print(norm_res)

pearson_vars  <- norm_res$Variable[norm_res$Best != "none-normal"]
spearman_vars <- norm_res$Variable[norm_res$Best == "none-normal"]

cat("\nPearson (normal) BIO vars:\n")
print(pearson_vars)
cat("\nSpearman (non-normal) BIO vars:\n")
print(spearman_vars)

###############################################################
# 3. Build mixed data matrix (transform Pearson BIO vars)
###############################################################

all_vars   <- c(bio_vars, env_vars)
mixed_data <- df[, all_vars, drop = FALSE]

# Apply transformations for Pearson biological variables
for (v in pearson_vars) {
  best <- norm_res$Best[norm_res$Variable == v]
  x    <- mixed_data[[v]]
  
  if (best == "log(x+1)") {
    mixed_data[[v]] <- log(x + 1)
  } else if (best == "log(-x+1)") {
    mixed_data[[v]] <- log(-x + 1)
  } else {
    mixed_data[[v]] <- x  # raw
  }
}

# Spearman BIO vars and ENV vars remain raw in mixed_data

###############################################################
# 4. Correlation matrix (row-wise: Pearson or Spearman)
###############################################################

row_vars <- bio_vars     # rows = all BIO in Excel order
col_vars <- all_vars     # columns = all BIO + ENV

n_rows <- length(row_vars)
n_cols <- length(col_vars)

R_mat <- matrix(NA_real_, n_rows, n_cols,
                dimnames = list(row_vars, col_vars))
P_mat <- matrix(NA_real_, n_rows, n_cols,
                dimnames = list(row_vars, col_vars))

for (i in seq_len(n_rows)) {
  row_name      <- row_vars[i]
  row_is_normal <- row_name %in% pearson_vars
  method        <- if (row_is_normal) "pearson" else "spearman"
  
  x <- mixed_data[[row_name]]
  
  for (j in seq_len(n_cols)) {
    col_name <- col_vars[j]
    y <- mixed_data[[col_name]]
    ok <- complete.cases(x, y)
    if (sum(ok) < 3) next
    
    ct <- cor.test(x[ok], y[ok], method = method)
    R_mat[i, j] <- unname(ct$estimate)
    P_mat[i, j] <- ct$p.value
  }
}

# will give warnings that for Spearman, exact p-values cannot be computed with ties. That's ok.

###############################################################
# 5. Bonferroni correction (ON)
###############################################################

use_bonferroni <- TRUE

P_for_plot <- P_mat

if (use_bonferroni) {
  p_flat <- as.vector(P_mat)
  p_adj  <- p.adjust(p_flat, method = "bonferroni")
  P_adj  <- matrix(p_adj, nrow = nrow(P_mat), ncol = ncol(P_mat),
                   dimnames = dimnames(P_mat))
  P_for_plot <- P_adj
}

###############################################################
# 6. Optional: r (p) matrix to inspect in console
###############################################################

r_p_matrix <- matrix(
  paste0(
    sprintf("%.2f", R_mat),
    " (",
    sprintf("%.3f", P_mat),
    ")"
  ),
  nrow = nrow(R_mat),
  ncol = ncol(R_mat),
  dimnames = list(rownames(R_mat), colnames(R_mat))
)

cat("\nCombined r (p) matrix (row-wise Pearson/Spearman, p raw):\n")
print(r_p_matrix, quote = FALSE)

###############################################################
# 7. Long-format data for ggplot
###############################################################

cor_df <- as.data.frame(as.table(R_mat))
colnames(cor_df) <- c("Biological", "Other", "Correlation")

p_df <- as.data.frame(as.table(P_for_plot))
colnames(p_df) <- c("Biological", "Other", "p_value")

cor_df <- merge(cor_df, p_df, by = c("Biological", "Other"))

cor_df$Annotation <- with(
  cor_df,
  ifelse(!is.na(p_value) & p_value < 0.001, "***",
         ifelse(!is.na(p_value) & p_value < 0.01,  "**",
                ifelse(!is.na(p_value) & p_value < 0.05,  "*", "")))
)

# Remove asterisks on the bio–bio diagonal (same variable)
same_pair <- as.character(cor_df$Biological) == as.character(cor_df$Other)
cor_df$Annotation[same_pair] <- ""

cor_df$Size <- abs(cor_df$Correlation)

cor_df$Biological <- factor(cor_df$Biological,
                            levels = rev(bio_vars))
cor_df$Other      <- factor(cor_df$Other,
                            levels = col_vars)

###############################################
# Custom axis labels with superscripts & italics
###############################################

label_map <- c(
  "C host tissue"           = "C host tissue",
  "N host tissue"           = "N host tissue",
  "C:N ratio host tissue"   = "C:N ratio host tissue",
  "C algal symbiont"        = "C algal symbiont",
  "N algal symbiont"        = "N algal symbiont",
  "C:N ratio algal symbiont"= "C:N ratio algal symbiont",
  
  # δ13C and δ15N superscripts
  "δ13C host tissue"        = expression(delta^13*C~"host tissue"),
  "δ15N host tissue"        = expression(delta^15*N~"host tissue"),
  
  "δ13C algal symbiont"     = expression(delta^13*C~"algal symbiont"),
  "δ15N algal symbiont"     = expression(delta^15*N~"algal symbiont"),
  
  "δ15N host-algal symbiont" = expression(delta^15*N~"host–algal symbiont"),
  "δ13C host-algal symbiont" = expression(delta^13*C~"host–algal symbiont"),
  
  # Other biological parameters
  "Algal symbiont density"   = "Algal symbiont density",
  "Mitotic index"           = "Mitotic index",
  
  # italic a in Chlorophyll a
  "Chlorophyll a"           = expression(italic("Chlorophyll")~italic("a")),
  
  # italic nifH but not (Ct values)
  "nifH (Ct values)"        = expression(italic("nifH")~"(Ct values)"),
  
  # Environmental variables
  "Temperature"             = "Temperature",
  "PAR"                     = "PAR",
  "DIN"                     = "DIN"
)

###############################################################
# 8. Supplemtary Figuire S5: Heatmap (Arial) + bio–bio & bio–env rectangles + Spearman rows
###############################################################

color_palette <- c(
  "#67001f","#b2182b","#d6604d","#f4a582","#fddbc7",
  "#f7f7f7",
  "#d1e5f0","#92c5de","#4393c3","#2166ac","#053061"
)

bio_count  <- length(bio_vars)
env_count  <- length(env_vars)
total_cols <- bio_count + env_count
row_count  <- length(bio_vars)

# Bio–Bio block
xmin1 <- 0.5
xmax1 <- bio_count + 0.5
ymin1 <- 0.5
ymax1 <- row_count + 0.5

# Bio–Env block
xmin2 <- bio_count + 0.5
xmax2 <- total_cols + 0.5
ymin2 <- 0.5
ymax2 <- row_count + 0.5

# Spearman rows: one red rectangle per row across full width
# spearman_vars comes from the normality test results
spearman_idx     <- match(spearman_vars, bio_vars)       # positions in original order (1..n)
spearman_idx_rev <- row_count - spearman_idx + 1         # reverse, because y-axis is reversed

spearman_rects <- data.frame(
  xmin = 0.5,
  xmax = total_cols + 0.5,
  ymin = spearman_idx_rev - 0.5,
  ymax = spearman_idx_rev + 0.5
)

p_mix <- ggplot(cor_df,
                aes(x = Other, y = Biological,
                    size = Size, color = Correlation)) +
  geom_tile(color = "grey", size = 0.5, alpha = 0) +
  geom_point(shape = 15) +
  geom_text(aes(label = Annotation),
            size = 5,
            color = "black",
            position = position_nudge(y = 0.10)) +
  # Black rectangle: bio–bio
  annotate("rect", xmin = xmin1, xmax = xmax1,
           ymin = ymin1, ymax = ymax1,
           fill = NA, color = "black", size = 1.2) +
  # Black rectangle: bio–env
  annotate("rect", xmin = xmin2, xmax = xmax2,
           ymin = ymin2, ymax = ymax2,
           fill = NA, color = "black", size = 1.2) +
  # Red rectangles: one per Spearman row
  geom_rect(data = spearman_rects,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            inherit.aes = FALSE,
            fill = NA, color = "red", size = 1.2) +
  scale_size_continuous(range = c(2, 10), guide = "none") +
  scale_color_gradientn(colors = color_palette, limits = c(-1, 1),
                        name = "Correlation (r / \u03c1)") +
  theme_minimal(base_family = "Arial") +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5, size = 9),
    axis.text.y = element_text(size = 9),
    panel.grid.major = element_line(color = "white", size = 1),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    legend.key.height = unit(1, "cm"),
    legend.key.width  = unit(0.3, "cm"),
    legend.title = element_text(size = 10),
    legend.text  = element_text(size = 9)
  ) +
  labs(x = NULL, y = NULL) +
  coord_fixed() +
  scale_x_discrete(
    position = "top",
    labels = label_map[levels(cor_df$Other)]
  ) +
  scale_y_discrete(
    labels = label_map[levels(cor_df$Biological)]
  )

# Check it out
print(p_mix)

ggsave("Supplementary_Fig_S5.pdf",
       plot = p_mix, width = 10, height = 7.5)

# Italics were added to Chlorophyll a and nifH (Ct values) in image processing software (setting labels as italics did not work)

###############################################################
# 9. Figure 2E – subset heatmap:
# Algal symbiont density, Mitotic index, Chlorophyll a
# vs. remaining environmental parameters (Temperature, PAR, DIN)
# Using existing R_mat and P_for_plot from the full analysis
###############################################################

# 1. Define the subsets
bio_subset <- c("Algal symbiont density",
                "Mitotic index",
                "Chlorophyll a")

env_subset <- c("Temperature", "PAR", "DIN")  # or env_vars if that equals these

# 2. Extract subset matrices from the full results
R_2e <- R_mat[bio_subset, env_subset, drop = FALSE]
P_2e <- P_for_plot[bio_subset, env_subset, drop = FALSE]

# Optional: inspect r (p) for Fig. 2E
r_p_matrix_2e <- matrix(
  paste0(
    sprintf("%.2f", R_2e),
    " (",
    sprintf("%.3f", P_2e),
    ")"
  ),
  nrow = nrow(R_2e),
  ncol = ncol(R_2e),
  dimnames = list(rownames(R_2e), colnames(R_2e))
)
cat("\nFig. 2E r (p) matrix (subset):\n")
print(r_p_matrix_2e, quote = FALSE)

# 3. Long-format data frame for ggplot
cor_df_2e <- as.data.frame(as.table(R_2e))
colnames(cor_df_2e) <- c("Biological", "Environmental", "Correlation")

p_df_2e <- as.data.frame(as.table(P_2e))
colnames(p_df_2e) <- c("Biological", "Environmental", "p_value")

cor_df_2e <- merge(cor_df_2e, p_df_2e, by = c("Biological", "Environmental"))

# Asterisks for significance (Bonferroni already applied in P_for_plot)
cor_df_2e$Annotation <- with(
  cor_df_2e,
  ifelse(!is.na(p_value) & p_value < 0.001, "***",
         ifelse(!is.na(p_value) & p_value < 0.01,  "**",
                ifelse(!is.na(p_value) & p_value < 0.05,  "*", "")))
)

cor_df_2e$Size <- abs(cor_df_2e$Correlation)

# Keep row order as in bio_subset (reversed on y-axis, like main heatmap)
cor_df_2e$Biological    <- factor(cor_df_2e$Biological,
                                  levels = rev(bio_subset))
cor_df_2e$Environmental <- factor(cor_df_2e$Environmental,
                                  levels = env_subset)

# 4. Plot Figure 2E

color_palette <- c(
  "#67001f","#b2182b","#d6604d","#f4a582","#fddbc7",
  "#f7f7f7",
  "#d1e5f0","#92c5de","#4393c3","#2166ac","#053061"
)

# Black frame around the 3x3 grid
xmin <- 0.5
xmax <- length(env_subset) + 0.5
ymin <- 0.5
ymax <- length(bio_subset) + 0.5

p_2e <- ggplot(cor_df_2e,
               aes(x = Environmental, y = Biological,
                   size = Size, color = Correlation)) +
  geom_tile(color = "grey", size = 0.5, alpha = 0) +
  geom_point(shape = 15) +
  geom_text(aes(label = Annotation),
            size = 5,
            color = "black",
            position = position_nudge(y = 0.10)) +
    scale_size_continuous(range = c(2, 10), guide = "none") +
  scale_color_gradientn(colors = color_palette, limits = c(-1, 1),
                        name = "Correlation") +
  theme_minimal(base_family = "Arial") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 0, vjust = 0.5, size = 9),
    axis.text.y = element_text(size = 9),
    panel.grid.major = element_line(color = "white", size = 1),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    legend.key.height = unit(0.6, "cm"),
    legend.key.width  = unit(0.3, "cm"),
    legend.title = element_text(size = 10),
    legend.text  = element_text(size = 9)
  ) +
  labs(x = NULL, y = NULL) +
  coord_fixed() +
  scale_x_discrete(position = "top")

# Check it out
print(p_2e)

# 5. Save Fig. 2E to file
ggsave("Figure_2E.pdf",
       plot = p_2e, width = 4, height = 3.2
      )

# Italics were added to Chlorophyll a in image processing software (setting labels as italics did not work)