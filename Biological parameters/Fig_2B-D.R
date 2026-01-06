library(readxl)
library(ggplot2)
library(ggforce)

# set working directory to the directory the script is saved in

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# read data

symbio<-read_xlsx("Biological Data file Tilstra et al.xlsx",sheet="symbionts_chla_MI_nifH_isotopic",na="NA")

# Change column names

colnames(symbio)[3:9]<-c("chla","mi","symbio","nifh","nifh_fold","delta_d15N","delta_d13C")

# Symbiont density

symbionts<-ggplot(symbio, aes(x = Season, y = symbio)) +
  geom_boxplot(fill = "darkolivegreen4",lwd=.35,outlier.shape = NA)+ 
  ggforce::geom_sina(
    ## draw bigger points
    size = 1,
    ## add some transparency
    alpha = .4,
    ## control range of the sina plot
    maxwidth = .5
  ) +
  scale_x_discrete(limits=c("Spring","Summer","Fall","Winter"))+
  scale_y_continuous(limits=c(0,1.6),breaks=seq(0,1.6,.2))+
  theme(aspect.ratio=2,panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA),axis.title.x=element_blank(),axis.text.x=element_text(size=10,angle=45,hjust=1),axis.title.y=element_text(size=16,margin = margin(t = 0, r = 3, b = 0, l = 0)),axis.text.y=element_text(size=10))+
  labs(y=expression(atop("",atop("Algal symbiont", paste("cell density (",x10^6, ~cm^-2,")")))))+
  # labels ned to be set in the alphabetic order of the Seasons
  stat_summary(geom = 'text', size=4,label = c("c","a","b","a"), fun = max, vjust = -1)+
  stat_summary(fun=mean, geom="point", shape=20, size=4, color="black") 

ggsave(
  "symbiont_density.pdf",
  symbionts,
  device = cairo_pdf,
  width = 5.8,
  height = 3.7
)

# Mitotic index

mi<-ggplot(symbio, aes(x = Season, y = mi)) +
  geom_boxplot(fill = "darkolivegreen4",lwd=.35)+ 
  ggforce::geom_sina(
    ## draw bigger points
    size = 1,
    ## add some transparency
    alpha = .4,
    ## control range of the sina plot
    maxwidth = .5
  ) +
  scale_x_discrete(limits=c("Spring","Summer","Fall","Winter"))+
  scale_y_continuous(limits=c(5,20),breaks=seq(5,20,5))+
  theme(aspect.ratio=2,panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA),axis.title.x=element_blank(),axis.text.x=element_text(size=10,angle=45,hjust=1),axis.title.y=element_text(size=16,margin = margin(t = 0, r = 3, b = 0, l = 0)),axis.text.y=element_text(size=10))+
  labs(y=expression(atop("",atop("","Mitotic index (%)"))))+
  # labels need to be set in the alphabetic order of the Seasons
  stat_summary(geom = 'text', size=4,label = c("a","ab","ab","b"), fun = max, vjust = -1)+
  stat_summary(fun=mean, geom="point", shape=20, size=4, color="black") 

ggsave(
  "mitotic_index.pdf",
  mi,
  device = cairo_pdf,
  width = 5.8,
  height = 3.7
)

#Chlorophyll a

chla<-ggplot(symbio, aes(x = Season, y = chla)) +
  geom_boxplot(fill = "darkolivegreen4",lwd=.35,outlier.shape = NA)+ 
  ggforce::geom_sina(
    ## draw bigger points
    size = 1,
    ## add some transparency
    alpha = .4,
    ## control range of the sina plot
    maxwidth = .5
  ) +
  scale_x_discrete(limits=c("Spring","Summer","Fall","Winter"))+
  scale_y_continuous(limits=c(0,9),breaks=seq(0,9,1.5))+
  theme(aspect.ratio=2,panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA),axis.title.x=element_blank(),axis.text.x=element_text(size=10,angle=45,hjust=1),axis.title.y=element_text(size=16,margin = margin(t = 0, r = 3, b = 0, l = 0)),axis.text.y=element_text(size=10))+
  labs(y=expression(atop("",atop("",paste("Chlorophyll",italic(" a "),"(",pg, ~cell^-1,")")))))+
  # labels ned to be set in the alphabetic order of the Seasons
  stat_summary(geom = 'text', size=4,label = c("a","a","a","b"), fun = max, vjust = -1)+
  stat_summary(fun=mean, geom="point", shape=20, size=4, color="black") 

ggsave(
  "chla.pdf",
  chla,
  device = cairo_pdf,
  width = 5.8,
  height = 3.7
)

# These three plots were later combined with the coral image and the correlation matrix (2E) in image processing software. Subplot letters were added there too.
