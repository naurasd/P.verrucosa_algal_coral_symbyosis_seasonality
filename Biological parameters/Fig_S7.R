library(readxl)
library(tidyverse)
library(ggforce)
library(patchwork)

# Set working directory to the directory the script is saved in
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Read raw data
data <- read_xlsx("Biological Data file Tilstra et al.xlsx", sheet = "elemental_and_isotopic_data", na = "NA")

# Subset data to the columns containing season info and delta d15N & delta d13C values

data<-data %>% select(c(11,13,14)) %>% filter(if_any(everything(), ~ !is.na(.)))

# rename columns
colnames(data)<-c("Season","delta_d15N","delta_d13C")

# Make plots

dd15n<-ggplot(data, aes(x = Season, y = delta_d15N)) +
  geom_boxplot(fill = "dodgerblue1",lwd=.35,outlier.shape = NA)+ 
  ggforce::geom_sina(
    ## draw bigger points
    size = 1,
    ## add some transparency
    alpha = .4,
    ## control range of the sina plot
    maxwidth = .5
  ) +
  scale_x_discrete(limits=c("Spring","Summer","Fall","Winter"))+
  scale_y_continuous(limits=c(-1,6),breaks=seq(-1,6,1))+
  theme(aspect.ratio=2,panel.grid.major = element_blank(), axis.title.x=element_blank(),panel.grid.minor = element_blank(), panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA),axis.text.x=element_text(size=10,angle=45,hjust=1),axis.title.y=element_text(size=16,margin = margin(t = 0, r = 3, b = 0, l = 0)),axis.text.y=element_text(size=10))+
  labs(y = expression(atop("",atop("Host tissue - algal symbiont",delta^15*N~("\211")))))+
  # labels need to be set in the alphabetic order of the Seasons
  stat_summary(geom = 'text', size=4,label = c("a","a","a","a"), fun = max, vjust = -1)+
  stat_summary(fun=mean, geom="point", shape=20, size=4, color="black") 

dd13c<-ggplot(data, aes(x = Season, y = delta_d13C)) +
  geom_boxplot(fill = "dodgerblue1",lwd=.35,outlier.shape = NA)+ 
  ggforce::geom_sina(
    ## draw bigger points
    size = 1,
    ## add some transparency
    alpha = .4,
    ## control range of the sina plot
    maxwidth = .5
  ) +
  scale_x_discrete(limits=c("Spring","Summer","Fall","Winter"))+
  scale_y_continuous(limits=c(-1.2,1.2),breaks=c(-1.2,-.8,-.4,0,.4,.8,1.2))+
  theme(aspect.ratio=2,panel.grid.major = element_blank(), axis.title.x=element_blank(),panel.grid.minor = element_blank(), panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA),axis.text.x=element_text(size=10,angle=45,hjust=1),axis.title.y=element_text(size=16,margin = margin(t = 0, r = 3, b = 0, l = 0)),axis.text.y=element_text(size=10))+
  labs(y = expression(atop("",atop("Host tissue - algal symbiont",delta^15*N~("\211")))))+
  # labels need to be set in the alphabetic order of the Seasons
  stat_summary(geom = 'text', size=4,label = c("ab","a","a","b"), fun = max, vjust = -1)+
  stat_summary(fun=mean, geom="point", shape=20, size=4, color="black")  

# Combine plots

combined<-dd13c + dd15n + plot_layout(ncol = 2)

ggsave("Figure_S7.pdf", combined, width = 5, height = 8)

# Subplot letters were added in image editing software outside of R
