
###################

### CAP Figure 3C ###

###################

library(readxl)
library(tidyverse)
library(vegan)
library(mice)
library(RColorBrewer)

# Set working directory to the directory the script is saved in
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Read raw data
data <- read_xlsx("Biological Data file Tilstra et al.xlsx", sheet = "bio + env data", na = "NA")

# The data also contains environmental variables and some biological variables which will not be included in the subsequent analysis
# Select the desired variables

data <- data %>% select(-c(13:14,19:31))
str(data)

# Use the package 'mice' to impute missing values (make 20 imputed data sets)

imputed_data <- mice(data, m = 20, maxit = 50, method = "pmm", seed = 123)

# getting a warning with 2 logged events, check them out

imputed_data$loggedEvents # just means that the first 2 columns (Season and Coral ID) had nothing to impute, mice leaves them unchanged

# Extract random data set of the 20 imputed ones and check it out

set.seed(123)
number<-sample(1:20, 1) # pick a random number between 1 and 20
completed_data <- complete(imputed_data, action=number) # pick the random data set
summary(completed_data)
str(completed_data)

# Add a constant to make all variables positive for Bray-Curtis dissimilarity matrix

completed_data[, -c(1, 2)] <- completed_data[, -c(1, 2)] + abs(min(completed_data[, -c(1, 2)])) + 1

# Calculate Bray-Curtis dissimilarity matrix

bray_matrix <- vegdist(completed_data[, -1:-2], method = "bray")

# Perform CAP analysis using the distance matrix

cap_result <- capscale(bray_matrix ~ completed_data$Season, permutations = 999)

# Get the percentages that explain the variation for each axis (found under Accumulated constrained eigenvalues Importance of components: in the row 'Proportion explained')

summary(cap_result)

# Proportion explained for CAP1 is 0.7635 = 76.35% and for CAP2 0.2042 = 20.42%
# This will be added later on to the CAP plot

# Fit vectors for the biological variables onto the ordination

cap_envfit <- envfit(cap_result, completed_data[, 3:16], permutations = 999)

## Preparation for plots

# Make season an ordered factor and set levels for legend later on

completed_data$Season <- factor(completed_data$Season, levels = c("Spring", "Summer", "Fall", "Winter"))
season_levels <- levels(completed_data$Season)

# Choose a pastel color palette from RColorBrewer with 4 colors (one for each season)

season_colors <- brewer.pal(4, "Pastel1")

# Define point symbols for each season

pch_values <- c(15, 16, 17, 18)

# Map colors and symbols to each sample

dot_colors <- season_colors[completed_data$Season]
dot_pch <- pch_values[completed_data$Season]

# Get scores for the plots

site_scores <- scores(cap_result, display = "sites")

##

## Make CAP plot with ellipses and a plot with the vectors next to it for Figure 3C and save as pdf

pdf("Figure_3C.pdf", width = 16, height = 8)

par(mfrow = c(1, 2), oma = c(2, 0, 0, 0)) # Two panels side by side

# Left plot: CAP + points + ellipses

par(mar = c(3,5,2,0.5)) # margin for left plot
plot(
  site_scores[, 1], site_scores[, 2],
  type = "n",
  xlab = "",
  ylab = "CAP2 (20.42%)",
  xlim = c(-0.65, 0.65),
  ylim = c(-0.7, 0.7),
  cex.lab = 1.5,
  cex.axis = 1.2
)

# Add points
points(
  cap_result,
  display = "sites",
  col = dot_colors,
  pch = dot_pch,
  cex = 1.5
)

# Add ellipses
ordiellipse(
  cap_result,
  display = "sites",
  groups = completed_data$Season,
  draw = "polygon",
  border = season_colors,
  col = season_colors,
  lwd = 2,
  kind = "sd",
  conf = 0.95,
  label = FALSE
)

# Add zero axes
abline(h = 0, v = 0, col = "darkgray", lty = "longdash")

# Add legend
legend(
  "bottomleft",
  legend = season_levels,
  col = season_colors,
  pch = pch_values,
  title = "Season",
  title.font=2,
  pt.cex = 2,
  cex = 1.2
)

# Right plot: vectors of biological variables
par(mar = c(3,0.5,2,1)) # margin for right plot
plot(
  site_scores[, 1], site_scores[, 2],
  type = "n",
  xlab = "",
  ylab = "",
  xlim = c(-0.65, 0.65),
  ylim = c(-0.7, 0.7),
  cex.lab = 1.5,
  cex.axis = 1.2,
  xaxt = "s",
  yaxt = "n"   # suppress default y-axis
)

# Add y-axis tick marks but no numbers
axis(2, labels = FALSE)

# Add zero axes
abline(h = 0, v = 0, col = "darkgray", lty = "longdash")

# Extract envfit vector coordinates
vec_scores <- scores(cap_envfit, display = "vectors")

# Keep only significant vectors with p.max = 0.001
sig <- cap_envfit$vectors$pvals <= 0.001
vec_scores_sig <- vec_scores[sig, , drop = FALSE]

# Adjust one vector name
rownames(vec_scores_sig)[rownames(vec_scores_sig) == "C:N ratio host tissue"]<-"C:N host tissue"

# Draw arrows
arrows(
  x0 = 0, y0 = 0,
  x1 = vec_scores_sig[, 1]*0.6, y1 = vec_scores_sig[, 2]*0.6, # scale the vectors to make them shorter
  length = 0.08,
  col = "black",
  lwd = 1.2
)

# Get vector label positions and scale them the same way as the vectors
label_x <- vec_scores_sig[, 1] * 0.6
label_y <- vec_scores_sig[, 2] * 0.6

# Manually adjust the label position where necessary

label_y["C algal symbiont"]  <- label_y["C algal symbiont"] - 0.05

label_y["Chlorophyll a"]  <- label_y["Chlorophyll a"] + 0.05

label_x["N host tissue"]  <- label_x["N host tissue"] - 0.05
label_y["N host tissue"]  <- label_y["N host tissue"] - 0.05

label_y["C host tissue"] <- label_y["C host tissue"] + 0.05

label_x["C:N host tissue"] <- label_x["C:N host tissue"] + 0.01
label_y["C:N host tissue"] <- label_y["C:N host tissue"] - 0.1

# Add all labels
text(
  x = label_x,
  y = label_y,
  labels = rownames(vec_scores_sig),
  cex = 1.2
)

# Add p.max note
usr <- par("usr")
text(
  x = usr[2] - 0.02 * diff(usr[1:2]),
  y = usr[4] - 0.02 * diff(usr[3:4]),
  labels = "p.max = 0.001",
  adj = c(1, 1),
  cex = 1.5
)

# Shared x-axis label across both panels
mtext("CAP1 (76.35%)", side = 1, outer = TRUE, line = 0, cex = 1.5)

dev.off()

# In image processing software outside of R, the subplot letter C was added and the plots were combined with the other plots of Figure 3.

## ANOSIM statistics ##

# Extract the CAP scores of the individual biological variables and print them
envmarker_scores <- cap_envfit$vectors
print(envmarker_scores)

# Perform ANOSIM and print results
anosim_result <- anosim(bray_matrix, completed_data$Season, permutations = 999)
print(anosim_result)
