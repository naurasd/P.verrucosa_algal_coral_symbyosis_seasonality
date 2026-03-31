library(readxl)
library(ggplot2)
library(ggforce)

# set working directory to the directory the script is saved in
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Create plot

nifh <- read_xlsx("Biological Data file Tilstra et al.xlsx", sheet = "season_mean_sem_bio", na = "NA")

# Remove unit row and Change column names for nifH

nifh<-nifh[-1,]
colnames(nifh)[8:10]<-c("nifh","nifh_se","nifh_diff")

# Make plot

nifh_plot<-ggplot(nifh, aes(x=Season, y=nifh)) +
  geom_bar(stat="identity",color="black", fill="dodgerblue1",width=0.5,position=position_dodge(),lwd=.35)+
  geom_linerange(aes(ymin=nifh, ymax=nifh+nifh_se),colour="black",lwd=.35)+
  scale_x_discrete(limits=c("Spring","Summer","Fall","Winter"))+
  theme(aspect.ratio=2,panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA),axis.text.x=element_text(size=15,angle=45,hjust=1),axis.title.y=element_text(size=25,margin = margin(t = 0, r = 3, b = 0, l = 0)),axis.text.y=element_text(size=15))+
  labs(x=NULL,y=expression(atop("",atop(paste("Relative", italic(" nifH "), "abundance"),"(fold change)"))))+
  scale_y_log10(limits=c(1,100),breaks=c(1,10,100))+
  geom_text(aes(label=nifh_diff,y=((nifh+nifh_se)*1.5)),size=6)

ggsave("Figure_S4.pdf",plot = nifh_plot, width = 4, height = 6,dpi=300)

