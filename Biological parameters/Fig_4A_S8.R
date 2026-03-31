### SIBER and bootstrapping procedure to estimate isotopic niche overlap between host tissue and algal symbionts ###

library(readxl)
library(tidyverse)
library(gridExtra)
library(SIBER)
library(patchwork)
library(cowplot)
library(grid)

# set working directory to the directory the script is saved in
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Read data
iso_data <- read_xlsx("Biological Data file Tilstra et al.xlsx", sheet = "elemental_and_isotopic_data", na = "NA")

# Subset data to the isotopic values (keeping only samples where both 13C and 15N values are present) and compartment info and arrange as required by SIBER:
# iso1 = d13C
# iso2 = d15N
# group = Compartment (host tissue or algal symbiont)
# community = 1 (needs to sequential numbers, so we set it just as 1 for all samples here because we will not differentiate between seasons for the following analysis)

sib <- iso_data %>% 
  select(c(6, 4, 3)) %>% # change order at the same time as subsetting
  na.omit() %>%
  setNames(replace(names(.), c(1, 2, 3), c("iso1", "iso2", "group"))) %>%
  mutate(community = 1)

sib <- as.data.frame(sib) # has to be made into dataframe for next step

# Make SIBER object
sib.dat <- createSiberObject(sib)

## Calculate summary statistics for host tissue and algal symbionts based on the maximum likelihood estimates of the means and covariance matrices of each group: ##
# TA: total area (= total convex hull area) of the polygon enclosing all sample points in δ13C–δ15N space
# SEA: standard ellipse area (= area where ~40% of observations are expected to fall) = core isotopic niche area
# SEAc: standard ellipse area corrected for small sample size 

group.ML <- groupMetricsML(sib.dat)
print(group.ML)


## Get Bayesian estimate for SEA (SEAb) 

# fit the ellipses via JAGS using an Inverse Wishart prior on the covariance matrix Sigma, and a vague normal prior on the means

# options for running jags
parms <- list()
parms$n.iter <- 2 * 10^4   # number of iterations to run the model for
parms$n.burnin <- 1 * 10^3 # discard the first set of values
parms$n.thin <- 10         # thin the posterior by this many
parms$n.chains <- 2        # run this many chains

# define the priors
priors <- list()
priors$R <- 1 * diag(2)
priors$k <- 2
priors$tau.mu <- 1.0E-3

# Fitting is via the JAGS method
# You need to have JAGS installed on your machine for this
# https://sourceforge.net/projects/mcmc-jags/files/
# -> go to JAGS directory
# -> go to 4.x (or whatever is the latest available) directory
# -> go to Windows, MAC OS et.c directory based on what machine you have
# -> Download JAGS-4.3.2.exe (or whatever is the latest) and run it on your machine

ellipses.posterior <- siberMVN(sib.dat, parms, priors)

# Posterior estimates of the ellipses for each group can be used to calculate the SEAb for each group
SEA.B <- siberEllipses(ellipses.posterior)

### organize the data for plotting SEAb with 50% and 95% credible intervals and mode (peak of the distribution, i.e., the point in plot S8B)

# Get HDR (highest density region) intervals (50 and 95%) from each posterior distribution

SEA.B.credibles <- lapply(
  as.data.frame(SEA.B),
  function(x) hdrcde::hdr(x, prob = c(50, 95))$hdr
)

# make table for plotting
# hdrcde::hdr() does not guarantee rows are returned in 50/95 order.
# So assign CI labels by interval width: narrower = 50%, wider = 95%.
seaBcreds <- bind_rows(lapply(seq_along(SEA.B.credibles), function(i) {
  out <- as.data.frame(SEA.B.credibles[[i]])
  colnames(out) <- c("lower", "upper")
  out$width <- out$upper - out$lower
  out <- out %>% arrange(width)
  out$CI <- c(50, 95)
  out$group_id <- i
  out
}))

# Get posterior mode for each group

SEA.B.modes <- lapply(
  as.data.frame(SEA.B),
  function(x) hdrcde::hdr(x, prob = c(50, 95), all.modes = TRUE)$mode
)

SeaBmodes <- unlist(lapply(SEA.B.modes, function(x) x[1]))

# Format data

seaBcreds$mode <- SeaBmodes[seaBcreds$group_id] # add mode to each row
seaBcreds$Tissue <- factor(
  c("Host", "Host", "Algal symbiont", "Algal symbiont"), # add labels
  levels = c("Host", "Algal symbiont")
)

# Set theme for all plots below

newtheme <- theme_classic() + theme(text = element_text(size=11))+
  theme(axis.text.x = element_text(size=11,colour="black"), axis.text.y = element_text(size=11,colour="black"))+
  theme(plot.margin = unit(c(5.5,5.5,5.5,20), "pt"))+
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=1))+
  theme(panel.background = element_blank())+
  theme(axis.line = element_blank())+
  theme(legend.position = "none")


# plot SEAb estimate with 50 and 95% CI
SEAb.plot <- ggplot(seaBcreds, aes(x = Tissue, y = mode, color = Tissue, fill = Tissue)) +
  scale_color_manual(values = c("sandybrown", "darkolivegreen4")) +
  scale_fill_manual(values = c("sandybrown", "darkolivegreen4")) +
  geom_errorbar(
    data = subset(seaBcreds, CI == 95),
    aes(ymin = lower, ymax = upper),
    width = 0,
    size = 0.8,
    show.legend = FALSE
  ) +
  geom_errorbar(
    data = subset(seaBcreds, CI == 50),
    aes(ymin = lower, ymax = upper),
    width = 0,
    size = 2,
    show.legend = FALSE
  ) +
  geom_point(pch = 21, col = "black", size = 4) +
  ylab(expression('SEA'[b] ~ '(\u2030'^2 * ')')) +
  scale_y_continuous(breaks = c(0.5, 2.5, 4.5), limits = c(0.3, 4.8)) +
  xlab("") +
  newtheme

SEAb.plot


## SEAc overlap as a metric of coral trophic strategy ##

# Calculate the overlap between host and symbiont pooled across seasons following Conti-Jerpe et al. 2020 Sci Adv.
# Overlap is presented as the proportion of the overlapping area relative to the host tissue SEAc

# Set group labels used for overlap functions below
host <- "1.Host tissue"
sym  <- "1.Algal symbiont"

# SEAc ellipses are corrected for sample sizes and fitted using maximum likelihood (gives host and symbiont ellipse area plus overlap area)
sea.overlap <- as.data.frame(t(maxLikOverlap(host, sym, sib.dat, p.interval = NULL, n = 100))) # p.interval = NULL -> uses ellipse definition as used by Conti-Jerpe et al., rather than the 40% predictive ellipse setting
sea.overlap$Prop.overlap <- sea.overlap[3] / sea.overlap[1] # overlap area as proportion of hist tissue SEAc
colnames(sea.overlap) <- c("Mean.SEAc.Host", "Mean.SEAc.Sym", "Overlap", "Prop.overlap")

print(sea.overlap) # the first overlap is the absolute area, the second value is the overlap in percentage area of host tissue ellipse area


## Bootstrap procedure to estimate uncertainty

nested <- sib %>% # split data for both compartments
  group_by(group, community) %>%
  nest()

boots <- list()
dist <- NULL

n <- 10000   # bootstrapping: 10,000 times (use smaller values for test runs if needed)

set.seed(123)
for(i in 1:n){  
  cat(i, fill = TRUE)
  
  # randomly sample each source population with replacement for its own sample size
  samps <- nested %>%
    mutate(samp = map(data, ~ sample_n(.x, size = nrow(.x), replace = TRUE)))
  
  # unnest iterative sampled lists
  samps2 <- samps %>%
    select(-data) %>%
    unnest(samp) %>%
    select(iso1, iso2, group, community)
  samps2 <- as.data.frame(samps2)
  
  sib.boot <- createSiberObject(samps2)
  
  # calculate distance between centroids of host and symbiont
  cen <- samps2 %>%
    group_by(group) %>%
    summarize(
      d13C = mean(iso1),
      d15N = mean(iso2),
      .groups = "drop"
    )
  
  c.axis <- abs(cen$d13C[1] - cen$d13C[2])
  n.axis <- abs(cen$d15N[1] - cen$d15N[2])
  D <- sqrt((c.axis^2) + (n.axis^2))
  
  dist <- rbind(dist, D)
  
  # now calculate overlap for each bootstrap sample
  boot.overlap <- as.data.frame(t(maxLikOverlap(host, sym, sib.boot, p.interval = NULL, n = 100)))
  boot.overlap$Prop.overlap <- as.numeric(boot.overlap[3] / boot.overlap[1])
  
  # gather results into list to save memory
  boots[[length(boots) + 1]] <- boot.overlap
}

boots <- do.call("rbind", boots)
colnames(boots) <- c("Mean.SEAc.Host", "Mean.SEAc.Sym", "Overlap", "Prop.overlap")

colnames(dist) <- "Distance"
dist <- as.data.frame(dist)

# save the output if desired for re-runs
saveRDS(boots,"siber_boots.rds")
saveRDS(dist,"siber_dist.rds")

# you can also write .csv files if desired
# write.csv(boots, file = "Bootstrapped_overlaps.csv", row.names = FALSE)
# write.csv(dist, file = "Bootstrapped_distance_btwn_centroids.csv", row.names = FALSE)


# create mode function to calculate mode of each distribution
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}


boots.ci <- boots %>%
  summarize(
    Mean.prop = mean(Prop.overlap, na.rm = TRUE), # this is the estimated mean percent overlap between the host tissue and algal symbionts (11.23%)
    Median.prop = median(Prop.overlap, na.rm = TRUE),
    Mode = getmode(Prop.overlap),
    prop.95.upper = quantile(Prop.overlap, 0.975),
    prop.95.lower = quantile(Prop.overlap, 0.025),
    prop.75.upper = quantile(Prop.overlap, 0.875),
    prop.75.lower = quantile(Prop.overlap, 0.125)
  )

# do same for distance between centroids
dist.ci <- dist %>%
  summarize(
    Mean.dist = mean(Distance, na.rm = TRUE),
    Median.dist = median(Distance, na.rm = TRUE),
    Mode = getmode(Distance),
    prop.95.upper = quantile(Distance, 0.975),
    prop.95.lower = quantile(Distance, 0.025),
    prop.75.upper = quantile(Distance, 0.875),
    prop.75.lower = quantile(Distance, 0.125)
  )

#### Figures 4A and S8 ####

## Fig. 4A ##

# Make a biplot of the pooled data set (all isotopic data across all season) (Fig. 4A)
# The core niche area is 40% of data ellipse

niche.plot<-ggplot()+
  geom_point(data=sib,mapping=aes(x=iso1,y=iso2,shape=group,fill=group),col='black',alpha=0.75,size=2)+ 
  scale_color_manual(values=c("darkolivegreen4","sandybrown"))+
  scale_shape_manual(values=c(24,21))+ 
  scale_fill_manual(values=c("darkolivegreen4","sandybrown"))+
  guides(fill=guide_legend(override.aes=list(shape=21)))+ #force color into legend
  stat_ellipse(data=subset(sib,group=="Host tissue"),mapping=aes(iso1,iso2,color=group,fill=group),
               type="norm",geom="polygon",level=.4,alpha=.25,size=1.5,show.legend = FALSE)+
  stat_ellipse(data=subset(sib,group=="Algal symbiont"),mapping=aes(iso1,iso2,color=group,fill=group),
               type="norm",geom="polygon",level=.4,alpha=.25,size=1.5,show.legend = FALSE)+
  xlab(expression({delta}^13*C~'(\u2030)'))+
  ylab(expression({delta}^15*N~'(\u2030)'))+
  newtheme+
  theme(legend.direction = "vertical",legend.position =c(.22,.85),legend.background = element_blank(),
        legend.title = element_blank())
niche.plot

ggsave(niche.plot,file="Figure_4A.pdf",height=5,width=5)

# In image processing software,...
# ...the subplot letter A was added
# ...the overlap values calculated above during bootstrap procedure was added
# ...the legend position and group shapes were adjusted
# ...the plot was combined with Figure 4B

## Fig. S8 ##

# labels for trophic cutoffs
het <- grobTree(textGrob(
  "Heterotrophy", x = 0.02, y = .07, hjust = 0,
  gp = gpar(col = "black", fontsize = 10, fontface = "italic")
))
mix <- grobTree(textGrob(
  "Mixotrophy", x = .37, y = .07, hjust = 0,
  gp = gpar(col = "black", fontsize = 10, fontface = "italic")
))
auto <- grobTree(textGrob(
  "Autotrophy", x = .80, y = .07, hjust = 0,
  gp = gpar(col = "black", fontsize = 10, fontface = "italic")
))

# PLot for S8C

pmain <- ggplot(
  boots.ci,
  aes(x = Mean.prop, y = 1.5, xmin = prop.95.lower, xmax = prop.95.upper)
) +
  geom_errorbar(
    aes(xmin = prop.95.lower, xmax = prop.95.upper),
    orientation = "y",
    height = 0,
    size = 1.5,
    color = "sandybrown"
  ) +
  geom_errorbar(
    aes(xmin = prop.75.lower, xmax = prop.75.upper),
    orientation = "y",
    height = 0,
    size = 3.5,
    color = "sandybrown"
  ) +
  geom_point(size = 5, pch = 21, col = "black", fill = "sandybrown") +
  xlab(expression('Percent of overlap with Coral Host SEA'[c] ~ '(%)')) +
  scale_x_continuous(
    breaks = c(0, 0.1, 0.2, 0.4, 0.6, 0.8, 1),
    labels = c(0, 10, 20, 40, 60, 80, 100),
    limits = c(-0.02, 1)
  ) +
  geom_text(
    data = boots.ci,
    aes(x = Mean.prop, y = 1.5, label = round(Mean.prop * 100, 2)),
    vjust = -.8,
    hjust = .1
  ) +
  geom_vline(xintercept = 0.1, lty = 2) +
  geom_vline(xintercept = 0.7, lty = 2) +
  annotation_custom(het) +
  annotation_custom(mix) +
  annotation_custom(auto) +
  newtheme +
  theme(axis.text.y = element_blank()) +
  ylab("") +
  scale_y_discrete(breaks = NULL)

pmain

# marginal density of resampled overlap distribution
xdens <- axis_canvas(pmain, axis = "x") +
  geom_density(data = boots, aes(x = Prop.overlap), alpha = 0.8, size = .2, fill = "sandybrown")

xdens

overlap.plot <- insert_xaxis_grob(pmain, xdens, grid::unit(.5, "null"), position = "top")
ggdraw(overlap.plot)


### draw distance confidence intervals
dist.plot <- ggplot(dist.ci, aes(x = 1, y = Mean.dist)) +
  geom_errorbar(
    data = dist.ci,
    aes(x = 1, ymin = prop.95.lower, ymax = prop.95.upper),
    width = 0,
    size = 1,
    show.legend = FALSE,
    color = "sandybrown"
  ) +
  geom_errorbar(
    data = dist.ci,
    aes(x = 1, ymin = prop.75.lower, ymax = prop.75.upper),
    width = 0,
    size = 2,
    show.legend = FALSE,
    color = "sandybrown"
  ) +
  geom_point(
    pch = 21,
    col = "black",
    fill = "sandybrown",
    size = 4
  ) +
  ylab(expression('Distance between centroids' ~ '(\u2030)')) +
  newtheme +
  theme(axis.text.x = element_blank()) +
  xlab("") +
  scale_x_continuous(breaks = NULL)

dist.plot


### Combine all 4 figures created in this script ###

##custom plot layout
coral.layout <- "
AABB
CCCD
"

combined.fig <- niche.plot + SEAb.plot + ggdraw(overlap.plot) + dist.plot +
  plot_layout(design=coral.layout)

combined.fig

ggsave("Figure_S8.pdf", combined.fig, width = 12, height = 8)

# In image processing software,...
# ...the subplot letter were added
# ...the plots position was adjusted
# ...the legend position and group shapes were adjusted in S8A
