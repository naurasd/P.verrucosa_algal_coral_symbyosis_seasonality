# load required packages

library(readxl)
library(ggplot2)
library(ggforce)

# set working directory of input and output files 

setwd(...)

# read in the input data
# we used two files as input data, as the nifH data needed to be read in from pre-calculated means due to its visualization as a barplot rather than a boxplot
# all other data was read in as the bio_months object 
# see the input files R_bio_boxplots.xlsx (missing data should be stored as blank cells) and R_bio_means.xlsx in the repository

bio_months<-read_xlsx("R_bio_boxplots.xlsx")
nifh_mean<-read_xlsx("R_bio_means.xlsx")

# generate boxplots (except for nifH, for which a barplot is generated)
# labels (letters) in stat_summary define seasons with significant differences as determined by PERMANOVA 
  
mi<-ggplot(bio_months, aes(x = Season, y = MI)) +
  geom_boxplot(fill = "dodgerblue1",lwd=.35)+ 
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
  theme(aspect.ratio=2,panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA),axis.text.x=element_text(size=10,angle=45,hjust=1),axis.title.y=element_text(size=16,margin = margin(t = 0, r = 3, b = 0, l = 0)),axis.text.y=element_text(size=10))+
  labs(y=expression(atop("",atop("","Mitotic index (%)"))))+
  # labels ned to be set in the alphabetic order of the Seasons
  stat_summary(geom = 'text', size=4,label = c("a","a","a","b"), fun = max, vjust = -1)+
  stat_summary(fun=mean, geom="point", shape=20, size=4, color="black") 

ggsave(mi,file="mitotic.pdf")

zoox<-ggplot(bio_months, aes(x = Season, y = Zoox)) +
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
  scale_y_continuous(limits=c(0,1.4),breaks=seq(0,1.4,.2))+
  theme(aspect.ratio=2,panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA),axis.text.x=element_text(size=10,angle=45,hjust=1),axis.title.y=element_text(size=16,margin = margin(t = 0, r = 3, b = 0, l = 0)),axis.text.y=element_text(size=10))+
  labs(y=expression(atop("",atop("Zooxanthellae", paste("cell density (",x10^6, ~cm^-2,")")))))+
  # labels ned to be set in the alphabetic order of the Seasons
  stat_summary(geom = 'text', size=4,label = c("a","ac","b","c"), fun = max, vjust = -1)+
  stat_summary(fun=mean, geom="point", shape=20, size=4, color="black") 

ggsave(zoox,file="zoox.pdf")

chloro<-ggplot(bio_months, aes(x = Season, y = chla)) +
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
  scale_y_continuous(limits=c(0,19),breaks=seq(0,18,2))+
  theme(aspect.ratio=2,panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA),axis.text.x=element_text(size=10,angle=45,hjust=1),axis.title.y=element_text(size=16,margin = margin(t = 0, r = 3, b = 0, l = 0)),axis.text.y=element_text(size=10))+
  labs(y=expression(atop("",atop("",paste("Chlorophyll",italic(" a "),"(",pg, ~cell^-1,")")))))+
  # labels ned to be set in the alphabetic order of the Seasons
  stat_summary(geom = 'text', size=4,label = c("b","a","a","b"), fun = max, vjust = -1)+
  stat_summary(fun=mean, geom="point", shape=20, size=4, color="black") 

ggsave(chloro,file="chla.pdf")

nifH<- ggplot(nifh_mean, aes(x=Season, y=nifh)) +
  geom_bar(stat="identity",color="black", fill="dodgerblue1",width=0.5,position=position_dodge(),lwd=.35)+
  geom_linerange(aes(ymin=nifh, ymax=nifh+nifh_se),colour="black",lwd=.35)+
  scale_x_discrete(limits=c("Spring","Summer","Fall","Winter"))+
  theme(aspect.ratio=2,panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA),axis.text.x=element_text(size=10,angle=45,hjust=1),axis.title.y=element_text(size=16,margin = margin(t = 0, r = 3, b = 0, l = 0)),axis.text.y=element_text(size=10))+
  labs(y=expression(atop("",atop(paste("Relative", italic(" nifH "), "abundance"),"(fold change)"))))+
  scale_y_log10(limits=c(1,100),breaks=c(1,10,100))+
  geom_text(aes(label=nifh_diff,y=((nifh+nifh_se)*1.5)),size=4)

ggsave(nifH,file="nifH.pdf")

d15n<-ggplot(bio_months, aes(x = Season, y = d15N,fill=part)) +
  geom_boxplot(outlier.shape=NA,lwd=.35) +
  ggforce::geom_sina(
    ## draw bigger points
    size = 1,
    ## add some transparency
    alpha = .4,
    ## control range of the sina plot
    maxwidth = .75
  ) +
  scale_x_discrete(limits=c("Spring","Summer","Fall","Winter"))+
  scale_y_continuous(limits=c(2,9),breaks=seq(2,9,2))+
  theme(legend.position = "none",panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA),axis.text.x=element_text(size=10,angle=45,hjust=1),axis.title.y=element_text(size=16,margin = margin(t = 0, r = 3, b = 0, l = 0)),axis.text.y=element_text(size=10))+
  labs(y = expression(atop("",atop("",delta^15*N~("\211")))))+
  scale_fill_manual(values=c("sandybrown","darkolivegreen4"))+
  stat_summary(fun=mean, geom="point", shape=20, size=4, color="black",position=position_dodge(.75))+
  # labels ned to be set in the alphabetic order of the Seasons (1. Fall: Host, Zoox; 2. Spring: Host, Zoox;...)
  stat_summary(geom = 'text', size=4,label = c("a","yz","a","yz","ab","y","b","z"), fun = max, vjust = -1,position=position_dodge(.75))

ggsave(d15n,file="d15n.pdf")  

d13c<-ggplot(bio_months, aes(x = Season, y = d13C,fill=part)) +
  geom_boxplot(outlier.shape=NA,lwd=.35) +
  ggforce::geom_sina(
    ## draw bigger points
    size = 1,
    ## add some transparency
    alpha = .4,
    ## control range of the sina plot
    maxwidth = .75
  ) +
  scale_x_discrete(limits=c("Spring","Summer","Fall","Winter"))+
  scale_y_continuous(limits=c(-17,-13),breaks=seq(-17,-13,1))+
  theme(legend.position = "none",panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA),axis.text.x=element_text(size=10,angle=45,hjust=1),axis.title.y=element_text(size=16,margin = margin(t = 0, r = 3, b = 0, l = 0)),axis.text.y=element_text(size=10))+
  labs(y = expression(atop("",atop("",delta^13*C~("\211")))))+
  scale_fill_manual(values=c("sandybrown","darkolivegreen4"))+
  stat_summary(fun=mean, geom="point", shape=20, size=4, color="black",position=position_dodge(.75))+
  # labels ned to be set in the alphabetic order of the Seasons (1. Fall: Host, Zoox; 2. Spring: Host, Zoox;...)
  stat_summary(geom = 'text', size=4,label = c("ab","y","a","yz","b","y","a","z"), fun = max, vjust = -1,position=position_dodge(.75))

ggsave(d13c,file="d13c.pdf")  

dd15n<-ggplot(bio_months, aes(x = Season, y = delta_d15N)) +
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
  theme(aspect.ratio=2,panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA),axis.text.x=element_text(size=10,angle=45,hjust=1),axis.title.y=element_text(size=16,margin = margin(t = 0, r = 3, b = 0, l = 0)),axis.text.y=element_text(size=10))+
  labs(y = expression(atop("",atop("Host tissue - Symbiodiniaceae",delta^15*N~("\211")))))+
  # labels ned to be set in the alphabetic order of the Seasons
  stat_summary(geom = 'text', size=4,label = c("ab","a","ab","b"), fun = max, vjust = -1)+
  stat_summary(fun=mean, geom="point", shape=20, size=4, color="black") 

ggsave(dd15n,file="delta_d15n.pdf")

dd13c<-ggplot(bio_months, aes(x = Season, y = delta_d13C)) +
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
  theme(aspect.ratio=2,panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA),axis.text.x=element_text(size=10,angle=45,hjust=1),axis.title.y=element_text(size=16,margin = margin(t = 0, r = 3, b = 0, l = 0)),axis.text.y=element_text(size=10))+
  labs(y = expression(atop("",atop("Host tissue - Symbiodiniaceae",delta^13*C~("\211")))))+
  # labels ned to be set in the alphabetic order of the Seasons
  stat_summary(geom = 'text', size=4,label = c("abc","a","b","c"), fun = max, vjust = -1)+
  stat_summary(fun=mean, geom="point", shape=20, size=4, color="black") 

ggsave(dd13c,file="delta_d13c.pdf")

cn<-ggplot(bio_months, aes(x = Season, y = CN,fill=part)) +
  geom_boxplot(outlier.shape=NA,lwd=.35) +
  ggforce::geom_sina(
    ## draw bigger points
    size = 1,
    ## add some transparency
    alpha = .4,
    ## control range of the sina plot
    maxwidth = .75
  ) +
  scale_x_discrete(limits=c("Spring","Summer","Fall","Winter"))+
  scale_y_continuous(limits=c(5,25),breaks=seq(5,25,5))+
  theme(legend.position = "none",panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA),axis.text.x=element_text(size=10,angle=45,hjust=1),axis.title.y=element_text(size=16,margin = margin(t = 0, r = 3, b = 0, l = 0)),axis.text.y=element_text(size=10))+
  labs(y = expression(atop("",atop("",C[org]:N~ratio))))+
  scale_fill_manual(values=c("sandybrown","darkolivegreen4"))+
  stat_summary(fun=mean, geom="point", shape=20, size=4, color="black",position=position_dodge(.75))+
  # labels ned to be set in the alphabetic order of the Seasons (1. Fall: Host, Zoox; 2. Spring: Host, Zoox;...)
  stat_summary(geom = 'text', size=4,label = c("a","yz","a","y","b","yz","a","z"), fun = max, vjust = -1,position=position_dodge(.75))

ggsave(cn,file="CN.pdf")  

carbon<-ggplot(bio_months, aes(x = Season, y = carbon,fill=part)) +
  geom_boxplot(outlier.shape=NA,lwd=.35) +
  ggforce::geom_sina(
    ## draw bigger points
    size = 1,
    ## add some transparency
    alpha = .4,
    ## control range of the sina plot
    maxwidth = .75
  ) +
  scale_x_discrete(limits=c("Spring","Summer","Fall","Winter"))+
  scale_y_continuous(limits=c(20,70),breaks=seq(20,70,10))+
  theme(legend.position = "none",panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA),axis.text.x=element_text(size=10,angle=45,hjust=1),axis.title.y=element_text(size=16,margin = margin(t = 0, r = 3, b = 0, l = 0)),axis.text.y=element_text(size=10))+
  labs(y = expression(atop("",atop("",C[org]~("%")))))+
  scale_fill_manual(values=c("sandybrown","darkolivegreen4"))+
  stat_summary(fun=mean, geom="point", shape=20, size=4, color="black",position=position_dodge(.75))+
  # labels ned to be set in the alphabetic order of the Seasons (1. Fall: Host, Zoox; 2. Spring: Host, Zoox;...)
  stat_summary(geom = 'text', size=4,label = c("bc","yz","ab","y","a","y","c","z"), fun = max, vjust = -1,position=position_dodge(.75))

ggsave(carbon,file="C.pdf")  

nitro<-ggplot(bio_months, aes(x = Season, y = nitrogen,fill=part)) +
  geom_boxplot(outlier.shape=NA,lwd=.35) +
  ggforce::geom_sina(
    ## draw bigger points
    size = 1,
    ## add some transparency
    alpha = .4,
    ## control range of the sina plot
    maxwidth = .75
  ) +
  scale_x_discrete(limits=c("Spring","Summer","Fall","Winter"))+
  scale_y_continuous(limits=c(2,8),breaks=seq(2,8,2))+
  theme(legend.position = "none",panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA),axis.text.x=element_text(size=10,angle=45,hjust=1),axis.title.y=element_text(size=16,margin = margin(t = 0, r = 3, b = 0, l = 0)),axis.text.y=element_text(size=10))+
  labs(y = expression(atop("",atop("",N~("%")))))+
  scale_fill_manual(values=c("sandybrown","darkolivegreen4"))+
  stat_summary(fun=mean, geom="point", shape=20, size=4, color="black",position=position_dodge(.75))+
  # labels ned to be set in the alphabetic order of the Seasons (1. Fall: Host, Zoox; 2. Spring: Host, Zoox;...)
  stat_summary(geom = 'text', size=4,label = c("a","y","a","yz","b","yz","a","z"), fun = max, vjust = -1,position=position_dodge(.75))

ggsave(nitro,file="N.pdf")

# in case you want a legend for the plots, you may generate one of the plots again (here the nitrogen plot has been chosen) with a legend that can be used outside of R later on and added to an overall plot
# two versions are generated here, one with a horizontal and one with a vertical legend

# vertical legend

legend_vertical<-ggplot(bio_months, aes(x = Season, y = nitrogen,fill=part)) +
  geom_boxplot(outlier.shape=NA,lwd=.35) +
  ggforce::geom_sina(
    ## draw bigger points
    size = 1,
    ## add some transparency
    alpha = .4,
    ## control range of the sina plot
    maxwidth = .75
  ) +
  scale_x_discrete(limits=c("Spring","Summer","Fall","Winter"))+
  scale_y_continuous(limits=c(2,8),breaks=seq(2,8,2))+
  theme(legend.title = element_blank(),legend.key=element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA),axis.text.x=element_text(size=10,angle=45,hjust=1),axis.title.y=element_text(size=16,margin = margin(t = 0, r = 3, b = 0, l = 0)),axis.text.y=element_text(size=10))+
  labs(y = expression(atop("",atop("",N~("%")))))+
  scale_fill_manual(values=c("sandybrown","darkolivegreen4"))+
  stat_summary(fun=mean, geom="point", shape=20, size=4, color="black",position=position_dodge(.75))+
  # labels ned to be set in the alphabetic order of the Seasons (1. Fall: Host, Zoox; 2. Spring: Host, Zoox;...)
  stat_summary(geom = 'text', size=4,label = c("a","y","a","yz","b","yz","a","z"), fun = max, vjust = -1,position=position_dodge(.75))

ggsave(legend_vertical,file="legend_vertical.pdf")

# horizontal legend

legend_horizontal<-ggplot(bio_months, aes(x = Season, y = nitrogen,fill=part)) +
  geom_boxplot(outlier.shape=NA,lwd=.35) +
  ggforce::geom_sina(
    ## draw bigger points
    size = 1,
    ## add some transparency
    alpha = .4,
    ## control range of the sina plot
    maxwidth = .75
  ) +
  scale_x_discrete(limits=c("Spring","Summer","Fall","Winter"))+
  scale_y_continuous(limits=c(2,8),breaks=seq(2,8,2))+
  theme(legend.title = element_blank(),legend.key=element_blank(),legend.position="bottom",panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA),axis.text.x=element_text(size=10,angle=45,hjust=1),axis.title.y=element_text(size=16,margin = margin(t = 0, r = 3, b = 0, l = 0)),axis.text.y=element_text(size=10))+
  labs(y = expression(atop("",atop("",N~("%")))))+
  scale_fill_manual(values=c("sandybrown","darkolivegreen4"))+
  stat_summary(fun=mean, geom="point", shape=20, size=4, color="black",position=position_dodge(.75))+
  # labels ned to be set in the alphabetic order of the Seasons (1. Fall: Host, Zoox; 2. Spring: Host, Zoox;...)
  stat_summary(geom = 'text', size=4,label = c("a","y","a","yz","b","yz","a","z"), fun = max, vjust = -1,position=position_dodge(.75))

ggsave(legend_horizontal,file="legend_horizontal.pdf")  

