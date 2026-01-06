# Load required libraries

library(tidyverse)
library(ggforce)
library(readxl)
library(patchwork)
library(ggtext)

# set working directory to the directory the script is saved in

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

### Figure(s) 1C-J ###

# Load the respective sheet of Excel file

env_data <- read_excel("Environmental Data file Tilstra et al.xlsx", sheet = "season_mean_sem_env")

# Reshape data from wide to long format

env_long <- env_data %>%
  pivot_longer(cols = -Season, names_to = "Parameter", values_to = "Value")

# Add a new column for standard errors of mean (SEM)

env_long <- env_long %>%
  mutate(Std_Error = lead(Value),   # lead() function shifts the values down by one row
         Std_Error = if_else(str_detect(lead(Parameter), "SEM"), Std_Error, NA_real_))

# Remove rows with SEM in the Parameter column

env_long <- env_long %>%
  filter(!str_detect(Parameter, "SEM"))

# Define which parameters to keep for this plot and their desired order

plot_order_1 <- c("Temperature", "PAR", "Nitrate", "DIN", "DIN:DIP", "DOC", "DON", "DOC:DON")

# Remove the environmental parameters not needed for this plot

env_long_1 <- env_long %>%
  filter(Parameter %in% plot_order_1)

# Convert Season to factor for correct ordering

env_long_1$Season <- factor(env_long_1$Season, levels = c("Spring", "Summer", "Fall", "Winter"))

# Reorder Parameter according to plot_order_1

env_long_1$Parameter <- factor(env_long_1$Parameter, levels = plot_order_1)

# Define the significance letters for each plot (based on statistics in different R script)

significance_letters <- list(
  c("a", "b", "c", "d"),
  c("a", "b", "ac", "c"),
  c("a", "a", "b", "a"),
  c("a", "a", "b", "a"),
  c("a", "ab", "b", "b"),
  c("ac", "b", "ab", "c"),
  c("a", "a", "b", "b"),
  c("a", "b", "c", "c")
)

# Match significance letters with the parameter order

env_long_1 <- env_long_1 %>% arrange(match(Parameter, plot_order_1))

# Repeat the significance letters for each parameter

env_long_1$significance <- rep(unlist(significance_letters), length.out = nrow(env_long_1))

# Make titles for each plot

facet_labels <- c(
  "Temperature" = "Water temperature (°C)",
  "PAR"         = "PAR (µmol m<sup>-2</sup> s<sup>-1</sup>)",
  "Nitrate"     = "Nitrate (µM)",
  "DIN"         = "DIN (µM)",
  "DIN:DIP"     = "DIN:DIP",
  "DOC"         = "DOC (µM)",
  "DON"         = "DON (µM)",
  "DOC:DON"     = "DOC:DON"
)

# Create plot

plot_1 <- env_long_1 %>%
  ggplot(aes(x = Season, y = Value, group = Parameter)) +
  geom_line(color = "deepskyblue3") +
  geom_point() +
  geom_errorbar(aes(ymin = Value - Std_Error, ymax = Value + Std_Error), width = 0.2, color = "black") +  
  geom_text(aes(label = significance), vjust = -0.5, size = 5, color = "black", position = position_nudge(y = -0.2, x = 0.3)) +  # Adjust position and aesthetics of significance letters
  facet_wrap(
    ~Parameter,
    scales = "free_y",
    ncol = 4,
    strip.position = "top",
    labeller = as_labeller(facet_labels)
  ) +
  labs(y = NULL, x = NULL) +
  guides(color = "none") + 
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    text = element_text(size = 16),
    strip.text = ggtext::element_markdown(size = 12),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  )

# Show the plot

print(plot_1)

# Save the plot

ggsave(
  "Figure_1C-J_blue_with_stats.pdf",
  plot_1,
  device = cairo_pdf,
  width = 9,
  height = 5
)

# significance letters and plot titles were adjusted and and panel letters added in image processing software outside of R


#### Supplementary Figure S1 ####

# Set order of parameters

plot_order_2 <- c(
  "Temperature", "PAR", "Salinity", "DO",
  "Nitrate", "Nitrite", "Ammonia",
  "DIN", "DON", "Phosphate",
  "DIN:DIP", "DOC", "DOC:DON"
)

env_long_2 <- env_long %>%
  filter(Parameter %in% plot_order_2) %>%
  mutate(
    Season    = factor(Season, levels = c("Spring","Summer","Fall","Winter")),
    Parameter = factor(Parameter, levels = plot_order_2)
  )

# Make Y-axis labels

y_labels_2 <- list(
  "Temperature" = "Water temperature (°C)",
  "PAR"         = expression(bold("PAR ("*mu*"mol m"^{-2}*" s"^{-1}*")")),
  "Salinity"    = "Salinity (PSU)",
  "DO"          = expression(bold("DO (mg L"^{-1}*")")),
  "Nitrate"     = "Nitrate (µM)",
  "Nitrite"     = "Nitrite (µM)",
  "Ammonia"     = "Ammonium (µM)",
  "DIN"         = "DIN (µM)",
  "DON"         = "DON (µM)",
  "Phosphate"   = "DIP (µM)",
  "DIN:DIP"     = "DIN:DIP",
  "DOC"         = "DOC (µM)",
  "DOC:DON"     = "DOC:DON"
)

# Define significance letters (based on statistics in the other R script)

significance_letters_2 <- list(
  c("a","b","c","d"),    # Temperature
  c("a","b","ac","c"),    # PAR
  c("a","bc","c","ab"),    # Salinity
  c("a","b","a","c"),   # DO
  c("a","a","b","a"),   # Nitrate
  c("a","a","a","a"),    # Nitrite
  c("a","a","a","a"),    # Ammonia
  c("a","a","b","a"),   # DIN
  c("a","a","b","b"),    # DON
  c("a","a","ab","b"),   # Phosphate
  c("a","ab","b","b"),    # DIN:DIP
  c("ac","b","ab","c"),  # DOC
  c("a","b","c","c")     # DOC:DON
)

env_long_2 <- env_long_2 %>%
  arrange(match(Parameter, plot_order_2), Season) %>%
  mutate(significance = unlist(significance_letters_2))

# Axis ticks

axis_specs <- tibble::tribble(
  ~Parameter,    ~breaks,                                 ~ymin,   ~ymax,
  "Temperature", 25:33,                                    25,      33,
  "PAR",         c(0, 200, 400, 600),                      0,       700,
  "Salinity",    38:42,                                    38,      42,
  "DO",          c(5.0, 5.5, 6.0, 6.5, 7.0),               5.0,     7.0,
  "Nitrate",     c(0.0, 0.2, 0.4, 0.6, 0.8, 1),            0.0,     1.1,
  "Nitrite",     c(0.00, 0.02, 0.04, 0.06, 0.08, 0.1),     0.00,    0.10,
  "Ammonia",     c(0.00, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30), 0.00, 0.30,
  "DIN",         c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4), 0.0,     1.4,
  "DON",         c(10, 11, 12, 13),                        9.5,      13,
  "Phosphate",   c(0.00, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30), 0.00, 0.30,
  "DIN:DIP",     c(0, seq(2, 16, 2)),                      0,       16,
  "DOC",         c(60, 65, 70, 75, 80),                    60,      80,
  "DOC:DON",     c(5.5, 6.0, 6.5, 7.0, 7.5, 8.0),          5.5,     8.0
)

# Make one square panel per parameter

plots_list <- env_long_2 %>%
  split(.$Parameter) %>%
  lapply(function(df_param) {
    
    par_name <- as.character(unique(df_param$Parameter))
    
    spec <- axis_specs %>%
      dplyr::filter(Parameter == par_name)
    
    brks <- spec$breaks[[1]]   
    ymin <- spec$ymin          
    ymax <- spec$ymax          
    ylab_val <- y_labels_2[[par_name]]
    
    ggplot(df_param, aes(x = Season, y = Value, group = 1)) +
      geom_errorbar(
        aes(ymin = Value - Std_Error,
            ymax = Value + Std_Error),
        width = 0.1,
        linewidth=0.25,
        color = "black"
      ) +
      geom_line(color = "dodgerblue1") +
      geom_point(size = 2.5, color = "dodgerblue1") +
      geom_text(
        aes(label = significance),
        vjust = -0.8,
        size  = 4,
        fontface = "bold"
      ) +
      scale_y_continuous(
        breaks = brks,
        limits = c(ymin, ymax),
        expand = c(0, 0)
      ) +
      scale_x_discrete(expand = expansion(mult = c(0.15, 0.15))) +
      labs(x = NULL, y = ylab_val) +
      theme_minimal(base_size = 10) +
      theme(
        axis.text.x  = element_text(angle = 45, hjust = 1,size=10),
        axis.text.y  = element_text(size = 10),
        axis.title.y = element_text(size = 10, face="bold"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
        panel.grid   = element_blank(),
        aspect.ratio = 1,
        plot.margin  = margin(5,5,5,5),
        panel.spacing = unit(1, "pt"), 
        axis.ticks   = element_line(color = "black"),
        axis.ticks.length = grid::unit(3, "pt")   # adjust or remove as you prefer
      )
  })

plot_2 <- wrap_plots(plots_list, ncol = 4)
plot_2

ggsave(
  "Supplementary_Fig_S1.pdf",
  plot_2,
  device = cairo_pdf,
  width = 12,
  height = 10
)

# Position of significance letters and plot titles were adjusted and panel letters added in image processing software outside of R

### Supplementary Figure S2 ###

# Load the respective sheet of Excel file

env_data_2 <- read_excel("Environmental Data file Tilstra et al.xlsx", sheet = "month_mean_sem_env")

# Reshape data from wide to long format

env_long_3 <- env_data_2 %>%
  pivot_longer(cols = -Month, names_to = "Parameter", values_to = "Value")

# Add a new column for standard errors of mean (SEM)

env_long_3 <- env_long_3 %>%
  mutate(Std_Error = lead(Value),   # lead() function shifts the values down by one row
         Std_Error = if_else(str_detect(lead(Parameter), "SEM"), Std_Error, NA_real_))

# Remove rows with SEM in the Parameter column

env_long_3 <- env_long_3 %>%
  filter(!str_detect(Parameter, "SEM"))

# Set order of parameters and transform month to proper date format

env_long_3 <- env_long_3 %>%
  filter(Parameter %in% plot_order_2) %>% # same order as in the previous plot S1
  mutate(
    Parameter = factor(Parameter, levels = plot_order_2),
    Month = as.Date(Month)
  )

# Axis ticks

axis_specs_2 <- tibble::tribble(
  ~Parameter,    ~breaks,                                 ~ymin,   ~ymax,
  "Temperature", seq(24, 34, 2),                           24,      34,
  "PAR",         c(0, 200, 400, 600),                      0,       700,
  "Salinity",    seq(39, 41, .5),                          39,      41,
  "DO",          c(4.5,5.0, 5.5, 6.0, 6.5, 7.0),           4.5,     7.0,
  "Nitrate",     c(0.0, 0.2, 0.4, 0.6, 0.8, 1,1.2),        0.0,     1.2,
  "Nitrite",     c(0.00, 0.02, 0.04, 0.06, 0.08, 0.1),     0.00,    0.10,
  "Ammonia",     c(0.00, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30), 0.00, 0.30,
  "DIN",         c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4,1.6), 0.0,     1.6,
  "DON",         seq(8, 13, 1),                        8,      13,
  "Phosphate",   c(0.00, 0.05, 0.10, 0.15, 0.20, 0.25), 0.00, 0.25,
  "DIN:DIP",     seq(0, 25, 5),                      0,       25,
  "DOC",         c(60, 65, 70, 75, 80),                    60,      80,
  "DOC:DON",     c(5.5, 6.0, 6.5, 7.0, 7.5, 8.0),          5.5,     8.0
)

plots_list_2 <- env_long_3 %>%
  split(.$Parameter) %>%
  lapply(function(df_param) {
    
    par_name <- as.character(unique(df_param$Parameter))
    
    spec <- axis_specs_2 %>%
      dplyr::filter(Parameter == par_name)
    
    brks <- spec$breaks[[1]]   
    ymin <- spec$ymin          
    ymax <- spec$ymax          
    ylab_val <- y_labels_2[[par_name]] # labels are the same as in S1 plot
    
    ggplot(df_param, aes(x = Month, y = Value, group = 1)) +
      geom_errorbar(
        aes(ymin = Value - Std_Error,
            ymax = Value + Std_Error),
        width = 10, # high width because it is expressed in time units
        linewidth=0.25,
        color = "black"
      ) +
      geom_line(color = "dodgerblue1") +
      geom_point(size = 2.5, color = "dodgerblue1") +
      scale_y_continuous(
        breaks = brks,
        limits = c(ymin, ymax),
        expand = c(0, 0)
      ) +
      scale_x_date(date_labels = "%b %y",date_breaks = "1 month", expand = expansion(mult = c(0.05, 0.05))) +
      labs(x = NULL, y = ylab_val) +
      theme_minimal(base_size = 10) +
      theme(
        axis.text.x  = element_text(angle = 45, hjust = 1,size=8),
        axis.text.y  = element_text(size = 8),
        axis.title.y = element_text(size = 8, face="bold"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
        panel.grid   = element_blank(),
        aspect.ratio = 0.5,
        plot.margin  = margin(5,5,5,5),
        panel.spacing = unit(1, "pt"), 
        axis.ticks   = element_line(color = "black"),
        axis.ticks.length = grid::unit(3, "pt")   # adjust or remove as you prefer
      )
  })

plot_3 <- wrap_plots(plots_list_2, ncol = 3)
plot_3

ggsave(
  "Supplementary_Fig_S2.pdf",
  plot_3,
  device = cairo_pdf,
  width = 15,
  height = 10
)

# Position of significance letters and plot titles were adjusted and panel letters added in image processing software outside of R

