library(readxl)
library(tidyverse)

# set working directory to the directory the script is saved in
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Read the data
df <- read_excel("raw data/Algal symbiont data + nifH.xlsx", sheet = "nifH")

# 16S Ct value averages are in the 6th column, rename it to "Ct"

df <- df %>%
  rename(Ct = 6) %>%
  filter(!is.na(Ct)) %>%  # remove empty value
  
  # Set desired order for plotting
  mutate(Season = factor(Season,
                         levels = c("Spring", "Summer", "Fall", "Winter")))

# Summarize (mean + standard error)
summary_df <- df %>%
  group_by(Season) %>%
  summarise(
    mean_Ct = mean(Ct),
    se_Ct = sd(Ct) / sqrt(n()),
    .groups = "drop"
  )

# Plot
plot_16s<-ggplot(summary_df, aes(x = Season, y = mean_Ct)) +
  
  geom_bar(stat = "identity",
           fill = "#1565C0",
           color = "black",
           width = 0.6) +
  
  geom_errorbar(aes(ymin = mean_Ct,
                    ymax = mean_Ct + se_Ct),
                width = 0,
                size = 0.8) +
  scale_y_continuous(
    limits = c(0, 30),
    breaks = seq(0, 30, 5)
  ) +
  labs(
    y = "16S rRNA\n(Ct values from qPCR)",
    x = NULL
  ) +
  
  theme_classic(base_size = 14) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    axis.line = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave(plot_16s,file="Figure_S9.pdf",width=4,height = 5)
