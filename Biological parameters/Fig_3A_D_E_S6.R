library(readxl)
library(tidyverse)
library(patchwork)

# set working directory to the directory the script is saved in
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Read the data from the Excel file 

element_iso <- read_xlsx("Biological Data file Tilstra et al.xlsx", sheet = "elemental_and_isotopic_data", na = "NA")

# Change column names and reorder Compartment factor levels

colnames(element_iso)[1:8] <- c("Season","Coral ID","Compartment","d15N", "N", "d13C", "C", "C_N")

element_iso$Compartment <- factor(
  element_iso$Compartment,
  levels = c("Host tissue", "Algal symbiont")
)

# C:N ratio plot (Fig. 3A)

c_n<-ggplot(element_iso, aes(x = Season, y = C_N, fill = Compartment)) +
  geom_boxplot(outlier.shape = NA, lwd = 0.35) +
  ggforce::geom_sina(
    size = 1,
    alpha = 0.4,
    maxwidth = 0.75
  ) +
  scale_x_discrete(limits = c("Spring", "Summer", "Fall", "Winter")) +
  scale_y_continuous(limits = c(5, 25), breaks = seq(5, 25, 5)) +
  theme(
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 10),
    axis.text.y = element_text(size = 10)
  ) +
  labs(y = expression(C:N ~ " ratio")) +
  scale_fill_manual(values = c("sandybrown", "darkolivegreen4")) +
  stat_summary(
    fun = mean,
    geom = "point",
    shape = 20,
    size = 4,
    color = "black",
    position = position_dodge(0.75)
  ) +
  stat_summary(
    geom = 'text',
    size = 4,
    label = c("b", "y", "ab", "y", "a", "y", "b", "y"), # labels need to be set in the alphabetic order of the Seasons (1. Fall: Host, Symbio; 2. Spring: Host, Symbio;...)
    fun = max,
    vjust = -1,
    position = position_dodge(0.75)
  )

# Delta 13C plot (Fig. 3D)

delta13C<-ggplot(element_iso, aes(x = Season, y = d13C, fill = Compartment)) +
  geom_boxplot(outlier.shape = NA, lwd = 0.35) +
  ggforce::geom_sina(
    size = 1,
    alpha = 0.4,
    maxwidth = 0.75
  ) +
  scale_x_discrete(limits = c("Spring", "Summer", "Fall", "Winter")) +
  scale_y_continuous(limits = c(-18, -13), breaks = seq(-18, -13, 1)) +
  theme(
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 10),
    axis.text.y = element_text(size = 10)
  ) +
  labs(y = expression(delta^13*C ~ " (%)")) +
  scale_fill_manual(values = c("sandybrown", "darkolivegreen4")) +
  stat_summary(
    fun = mean,
    geom = "point",
    shape = 20,
    size = 4,
    color = "black",
    position = position_dodge(0.75)
  ) +
  stat_summary(
    geom = 'text',
    size = 4,
    label = c("a", "y", "a", "yz", "a", "y", "a", "z"), # labels need to be set in the alphabetic order of the Seasons (1. Fall: Host, Symbio; 2. Spring: Host, Symbio;...)
    fun = max,
    vjust = -1,
    position = position_dodge(0.75)
  )

# Delta 15N plot (Fig. 3E)

delta15n<-ggplot(element_iso, aes(x = Season, y = d15N, fill = Compartment)) +
  geom_boxplot(outlier.shape = NA, lwd = 0.35) +
  ggforce::geom_sina(
    size = 1,
    alpha = 0.4,
    maxwidth = 0.75
  ) +
  scale_x_discrete(limits = c("Spring", "Summer", "Fall", "Winter")) +
  scale_y_continuous(limits = c(2, 9), breaks = seq(2, 9, 2)) +
  theme(
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA),
    axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
    axis.title.y = element_text(size = 10),
    axis.text.y = element_text(size = 10)
  ) +
  labs(y = expression(delta^15*N ~ " (%)")) +
  scale_fill_manual(values = c("sandybrown", "darkolivegreen4")) +
  stat_summary(
    fun = mean,
    geom = "point",
    shape = 20,
    size = 4,
    color = "black",
    position = position_dodge(0.75)
  ) +
  stat_summary(
    geom = 'text',
    size = 4,
    label = c("a", "y", "a", "y", "a", "y", "a", "y"), # labels need to be set in the alphabetic order of the Seasons (1. Fall: Host, Symbio; 2. Spring: Host, Symbio;...)
    fun = max,
    vjust = -1,
    position = position_dodge(0.75)
  )

# Combine the three plots

isotope_stacked <- c_n / delta13C / delta15n +
  plot_layout(ncol = 1, heights = c(1, 1, 1))

# Save figure

ggsave(
  "Figure_3A_D_E.pdf",
  isotope_stacked,
  width = 4,
  height = 8,
  dpi = 300
)

# These three plots were later combined with the regression plots (Fig. 3B) and the CAP plots (Fig. 3C) in image processing software. Subplot letters were added there too.

# make a separate plot just for the legend to extract it later on 

legend<-ggplot(element_iso, aes(x = Season, y = C_N, fill = Compartment)) +
  geom_boxplot(outlier.shape = NA, lwd = 0.35) +
  ggforce::geom_sina(
    size = 1,
    alpha = 0.4,
    maxwidth = 0.75
  ) +
  scale_x_discrete(limits = c("Spring", "Summer", "Fall", "Winter")) +
  scale_y_continuous(limits = c(5, 25), breaks = seq(5, 25, 5)) +
  theme(
    legend.title = element_blank(),
    legend.key=element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 10),
    axis.text.y = element_text(size = 10)
  ) +
  labs(y = expression(C:N ~ " ratio")) +
  scale_fill_manual(values = c("sandybrown", "darkolivegreen4")) +
  stat_summary(
    fun = mean,
    geom = "point",
    shape = 20,
    size = 4,
    color = "black",
    position = position_dodge(0.75)
  ) +
  stat_summary(
    geom = 'text',
    size = 4,
    label = c("b", "y", "ab", "y", "a", "y", "b", "y"), # labels need to be set in the alphabetic order of the Seasons (1. Fall: Host, Symbio; 2. Spring: Host, Symbio;...)
    fun = max,
    vjust = -1,
    position = position_dodge(0.75)
  )

ggsave(
  "Figure_3A_D_E_vertical_legend.pdf",
  legend,
  width = 4,
  height = 4,
  dpi = 300
)

# This legend was then manually extracted in image processing software outside of R and added to Figure 3A of the previously created stacked Figure 3A,D,E

# Supplementary Figure S6

# C plot (Fig. S6A)

c_plot<-ggplot(element_iso, aes(x = Season, y = C, fill = Compartment)) +
  geom_boxplot(outlier.shape = NA, lwd = 0.35) +
  ggforce::geom_sina(
    size = 1,
    alpha = 0.4,
    maxwidth = 0.75
  ) +
  scale_x_discrete(limits = c("Spring", "Summer", "Fall", "Winter")) +
  scale_y_continuous(limits=c(20,70),breaks=seq(20,70,10)) +
  theme(
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
    axis.title.y = element_text(size = 16),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_blank()
  ) +
  labs(y = expression(C ~ " (%)")) +
  scale_fill_manual(values = c("sandybrown", "darkolivegreen4")) +
  stat_summary(
    fun = mean,
    geom = "point",
    shape = 20,
    size = 4,
    color = "black",
    position = position_dodge(0.75)
  ) +
  stat_summary(
    geom = 'text',
    size = 4,
    label = c("bc", "y", "ab", "yz", "a", "y", "c", "z"), # labels need to be set in the alphabetic order of the Seasons (1. Fall: Host, Symbio; 2. Spring: Host, Symbio;...)
    fun = max,
    vjust = -1,
    position = position_dodge(0.75)
  )

# N plot (Fig. S6B)

n_plot<-ggplot(element_iso, aes(x = Season, y = N, fill = Compartment)) +
  geom_boxplot(outlier.shape = NA, lwd = 0.35) +
  ggforce::geom_sina(
    size = 1,
    alpha = 0.4,
    maxwidth = 0.75
  ) +
  scale_x_discrete(limits = c("Spring", "Summer", "Fall", "Winter")) +
  scale_y_continuous(limits=c(2,8),breaks=seq(2,8,2)) +
  theme(
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
    axis.title.y = element_text(size = 16),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_blank()
  ) +
  labs(y = expression(N ~ " (%)")) +
  scale_fill_manual(values = c("sandybrown", "darkolivegreen4")) +
  stat_summary(
    fun = mean,
    geom = "point",
    shape = 20,
    size = 4,
    color = "black",
    position = position_dodge(0.75)
  ) +
  stat_summary(
    geom = 'text',
    size = 4,
    label = c("b","y","ab","yz","a","yz","b","z"), # labels need to be set in the alphabetic order of the Seasons (1. Fall: Host, Symbio; 2. Spring: Host, Symbio;...)
    fun = max,
    vjust = -1,
    position = position_dodge(0.75)
  )

# Combine the two plots

elements_stacked <- c_plot / n_plot +
  plot_layout(ncol = 2, widths = c(1, 1))

# Save figure

ggsave(
  "Figure_S6.pdf",
  elements_stacked,
  width =10,
  height = 4,
  dpi = 300
)

# In image processing software outside of R, subplot letters and the previously created legend were added.
