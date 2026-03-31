###############################################################
# ENVIRONMENTAL AUTOCORRELATION ANALYSIS
# - Normality tests
# - log10-transform where appropriate
# - Pearson env × env correlation matrix
# - Heatmap in the same color palette as previous figures
# - List of high-correlation pairs (|r| ≥ 0.7)
# - Dendrogram for clustering variables
###############################################################

library(readxl)
library(ggplot2)

# set working directory to the directory the script is saved in

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

###############################################################
# STEP 1: Load data (minus ecologically irrelevant parameters)
###############################################################

env_data <- read_xlsx("Environmental Data file Tilstra et al.xlsx", sheet = "season_mean_sem_env")

# Keep only ecologically relevant / varying environmental variables and exclude standard error columns
env_vars <- c(
  "Temperature",
  "DO",
  "PAR",
  "Nitrate",
  "DIN",
  "DIN:DIP",
  "DOC",
  "DON",
  "DOC:DON"
)

env <- env_data[, env_vars]

# Check structure
print(colnames(env))
print(env)

###############################################################
# 2. Shapiro–Wilk normality tests (raw environmental data)
###############################################################

shapiro_results_env <- data.frame(
  variable = env_vars,
  p_value = NA_real_
)

for (i in seq_along(env_vars)) {
  v <- env_vars[i]
  x <- env[[v]]
  x_no_na <- x[!is.na(x)]
  
  if (length(x_no_na) >= 3) {
    shapiro_results_env$p_value[i] <- shapiro.test(x_no_na)$p.value
  }
}

print(shapiro_results_env)

###############################################################
# 3. Pearson correlation matrix (env × env)
###############################################################

# Correlation matrix (Pearson r)
R_env <- cor(env, use = "pairwise.complete.obs", method = "pearson")

# Matrix for p-values
p_env <- matrix(NA_real_, nrow = ncol(env), ncol = ncol(env))
rownames(p_env) <- colnames(env)
colnames(p_env) <- colnames(env)

n <- nrow(env)

for (i in seq_len(ncol(env))) {
  for (j in seq_len(ncol(env))) {
    if (i == j) {
      p_env[i, j] <- NA
    } else {
      x <- env[[i]]
      y <- env[[j]]
      ok <- complete.cases(x, y)
      if (sum(ok) >= 3) {
        p_env[i, j] <- cor.test(x[ok], y[ok], method = "pearson")$p.value
      } else {
        p_env[i, j] <- NA
      }
    }
  }
}

# Convert to long format for plotting
cor_env_df <- as.data.frame(as.table(R_env))
colnames(cor_env_df) <- c("Var1", "Var2", "Correlation")

p_env_df <- as.data.frame(as.table(p_env))
colnames(p_env_df) <- c("Var1", "Var2", "p_value")

cor_env_df <- merge(cor_env_df, p_env_df, by = c("Var1", "Var2"))

# Significance stars (note: n = 4 -> use with caution in interpretation)
cor_env_df$Annotation <- with(
  cor_env_df,
  ifelse(!is.na(p_value) & p_value < 0.001, "***",
         ifelse(!is.na(p_value) & p_value < 0.01, "**",
                ifelse(!is.na(p_value) & p_value < 0.05, "*", "")))
)

# Absolute r for point size
cor_env_df$Size <- abs(cor_env_df$Correlation)

# Factor ordering
cor_env_df$Var1 <- factor(cor_env_df$Var1, levels = env_vars)
cor_env_df$Var2 <- factor(cor_env_df$Var2, levels = env_vars)

###############################################################
# 4. Heatmap of environmental autocorrelation (Pearson)
###############################################################

color_palette <- c(
  "#67001f", "#b2182b", "#d6604d", "#f4a582", "#fddbc7",
  "#f7f7f7",
  "#d1e5f0", "#92c5de", "#4393c3", "#2166ac", "#053061"
)

p_env <- ggplot(cor_env_df,
                aes(x = Var1, y = Var2,
                    size = Size, color = Correlation)) +
  geom_tile(color = "grey", size = 0.5, alpha = 0) +
  geom_point(shape = 15) +
  geom_text(aes(label = Annotation),
            size = 5, color = "white",
            position = position_nudge(y = 0.1)) +
  
  # Remove the SIZE legend
  scale_size_continuous(range = c(2, 10), guide = "none") +
  
  # KEEP the COLOR legend (this is your correlation scale)
  scale_color_gradientn(colors = color_palette,
                        limits = c(-1, 1),
                        name = "Correlation") +
  
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    panel.grid.major = element_line(color = "white", size = 1),
    panel.grid.minor = element_blank(),
    
    # Keep only the color legend
    legend.position = "right",
    legend.key.size = unit(1.2, "cm"),
    legend.key.width = unit(0.3, "cm"),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10)
  ) +
  
  labs(title = NULL, x = NULL, y = NULL) +
  coord_fixed(ratio = 1)

print(p_env)

ggsave("Figure_S3A.pdf",plot = p_env, width = 6, height = 6,dpi=300)

# Legend title was moved further up and subplot letter (A) was added in image processing software outside of R.

###############################################################
# 5. Strongly correlated pairs |rho| ≥ 0.7
###############################################################

threshold <- 0.7
R_abs <- abs(R_env)

idx <- which(R_abs >= threshold & row(R_abs) < col(R_abs), arr.ind = TRUE)

high_pairs <- data.frame(
  Var1 = rownames(R_env)[idx[, "row"]],
  Var2 = colnames(R_env)[idx[, "col"]],
  r = R_env[idx]
)

high_pairs <- high_pairs[order(-abs(high_pairs$r)), ]

print(high_pairs)

###############################################################
# 6. Dendrogram of environmental variables
#    Distance = 1 - |Pearson rho|
###############################################################

dist_env <- as.dist(1 - abs(R_env))

hc_env <- hclust(dist_env, method = "average")

plot(hc_env,xlab = "", sub = "")

abline(h = 0.3, col = "red", lty = 2)  # optional guideline

# save to file

pdf("Figure_S3B.pdf", width = 8, height = 6)
plot(hc_env,xlab = "", sub = "",main=NULL)
abline(h = 0.3, col = "red", lty = 2)
dev.off()

# Subplot letter (B) was added in image processing software outside of R.

