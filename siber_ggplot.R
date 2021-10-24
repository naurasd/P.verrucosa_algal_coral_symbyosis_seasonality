library(devtools)
devtools::install_github("cmartin/ggConvexHull")
library(ggConvexHull)
library(SIBER)
library(tidyverse)

# import the data. 

setwd("C:/Users/Nauras/Dropbox/Nauras_Arjen/Analysis/R_scripts/SIBER")

siber.data <- read.csv("siber_data_plot.csv", header=T)

# make a copy of our data for use here in this example, and 
# set the columns group and community to be factor type using dplyr.
# Additionally rename the isotope data columns and drop the iso1 and iso2
# columns using the .keep option to keep only those that were not used to 
# create the new variables, i.e. keep only ones on the left of the "=". 
siber_data <- siber.data %>% mutate(group = factor(group), 
                                        Season = factor(Season),
                                        d13C = iso1, 
                                        d15N = iso2,
                                        .keep = "unused") 

# when plotting colors and shapes, we need to tell ggplot that these are to 
# be treated as categorical factor type data, and not numeric.

siber.plot <- ggplot(data = siber_data, 
                     aes(x = d13C, 
                         y = d15N)) + 
  geom_point(aes(color = group, shape = Season), size = 2) +
  scale_shape_manual(values=15:18,breaks=c("Spring","Summer","Fall","Winter"))+
  ylab(expression(paste(delta^{15}, "N (\u2030)"))) +
  theme_classic()+
  xlab(expression(paste(delta^{13}, "C (\u2030)"))) + 
  theme(aspect.ratio=.5,panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA),axis.text.x=element_text(size=10),legend.key.width= unit(.4, 'cm'),legend.key.height= unit(.4, 'cm'),axis.title.y=element_text(size=12,margin = margin(t = 0, r = 10, b = 0, l = 0)),axis.title.x=element_text(size=12,margin = margin(t = 10, r = 0, b = 0, l = 0)),axis.text.y=element_text(size=10))+
  scale_color_manual(name="",values = c("sandybrown","darkolivegreen4"))+
  scale_x_continuous(limits=c(-18,-13))+
  scale_y_continuous(limits=c(1.5,8.5),breaks=seq(2,8,2))+
  guides(color=guide_legend(order=1),shape=guide_legend(order=2))+
  geom_convexhull(alpha=0,aes(color = group),size=.2)+
  stat_ellipse(alpha=0,aes(group = interaction(group, community),color = group), 
               alpha = 0.25, 
               level = 0.4,
               type = "norm",
               geom = "polygon") 

ggsave(siber.plot,file="siber_plot.pdf",path="C:/Users/Nauras/Dropbox/Nauras_Arjen/Analysis/R_scripts/SIBER",width=7,height=5)








