## SIBER Script following the scripts by:

## Conti-Jerpe, I.E., Thompson, P.D., Wong, C.W.M., Oliveira, N.L., Duprey, N.N., Moynihan, M.A., BAker, D.M. 2020. 
## Trophic strategy and bleaching resistance in reef-building corals. 
## Science Advances, 6: eaaz5443. 
## DOI: https://doi.org/10.1126/sciadv.aaz5443

## which in turn is based on scripts and theory by the following works: 

## Jackson, A.L., Parnell, A.C., Inger R., & Bearhop, S. 2011. Comparing isotopic niche 
## widths among and within communities: SIBER - Stable Isotope Bayesian Ellipses in R. 
## Journal of Animal Ecology, 80: 595-602. 
## DOI: https://doi.org/10.1111/j.1365-2656.2011.01806.x 

## Turner, T.F., Collyer, M.L. & Krabbenhoft, T.J. 2010. A general hypothesis-testing framework 
## for stable isotope ratios in ecological studies. 
## Ecology, 91: 2227-2233. 
## DOI: https://doi.org/10.1890/09-1454.1

## Used in ...


### Calculate the overlap of SEAcs ### 

# install and load the SIBER package

install.packages("SIBER")
library(SIBER)
set.seed(1)

# read input data 

# headings of columns in .csv must be "iso1" (d13C), "iso2" (d15N), "group" (host or symbionts), and "community" (set to "1" for all of our entries, as we are only dealing with a single overall community)
# all headings & columns must be present
# only include coral fragments where values exist for both for groups (host and symbiont) and for both carbon and nitrogen stable isotopes
# check the file siber_data_total.csv in the repository

# set working directory of input files

setwd(...)

isotopes <- read.csv("siber_data_total.csv", header=T)

siber.data <- createSiberObject(isotopes)

# The first ellipse is referenced using a character string representation where 
# in "x.y", "x" is the community, and "y" is the group within that community.

ellipse1 <- "1.Host tissue" 

ellipse2 <- "1.Symbiodiniaceae"

# Overlap metric (units are per mil^2)
# The overlap of the maximum likelihood fitted standard ellipses are estimated as follows

sea.overlap <- maxLikOverlap(ellipse1, ellipse2, siber.data, 
                             p.interval = NULL, n = 100)
sea.overlap

# overlap as a proportion of group 1 area (host)

prop.sea.over.1 <- sea.overlap[3] / (sea.overlap[1])
prop.sea.over.1


### Caculate distance of host and symbiont SEAc centroids and other measures for subsequent permutation procedure ###

# Read in the the input data
# this should be a .csv file with 3 columns: "group" (host and symbionts), "x" (d13C), and "y" (d15N)
# only include coral fragments where values exist for both for groups (host and symbiont) and for both carbon and nitrogen stable isotopes
# check the file iso_niche_data.csv in the repository

iso<-read.csv("iso_niche_data.csv",header=T)

Y<-as.matrix(iso[,2:3])

# Designate groups
group<-as.factor(iso[,1])
gp<-length(levels(group)) # number of groups
n.comp<-(gp^2-gp)/2 # number of possible comparisons
rownames(Y)<-group 

lm.gp<-lm(Y~group,x=T,model=T) # for estimating group means
#This outputs the centroid of the Host tissue group and distance in x and y of the Symbiodiniaceae group from it 

res.gp<-resid(lm.gp) # residuals of groups from group means
yhat.gp<-predict(lm.gp) # predicted values

lm.gp.red<-lm(Y~1) # this is the model reduced by the group factor: only estimates an overall mean

# Read in the source code developed by Turner et al. (see above) for the following functions
# You can find this file in the repository and it should be placed in your working directory of input files

source('Turner.et.al.ecology.source.r')

# DISPERSION MEASURES
ex1<-ds.prep(res.gp,group) # see source file
ex1.ds<-disp.stat(ex1) # see source file
ex1.ds

# GROUP MEANS
gp.m<-group.means(Y,group) # finds the group means for the raw data
gp.m

# CONTRASTS (distance of group SEAc centroids (in per mil))
# This is the value reported as distance between host and symbiont SEAcs in the paper
mean.dif<-as.vector(dist(as.matrix(gp.m)))
mean.dif

# Calculate...
#mdc = mean centroid distance
#mnn = mean nearest neighbor
#ecc = eccentricity
#...needed for eprmutation

disp.dif<-ds.dif(ex1.ds) # Dimension names not given by output

disp.dif 
disp.dif$mdc
disp.dif$mnn
disp.dif$ecc


### PERMUTATION PROCEDURE to generate p-value ###

# Describe a function that needs
# x = the full linear model
# x2 = the reduced linear model
# y = the raw data
# g = group list
# p = the number of permutations to use (keep in mind that
# the observed values count as 1 random permutation)

# SETUP
permute<-function(x,x2,y,g,p){
  p.table<-NULL
  
  yhat<-predict(x)
  res<-resid(x)
  line<-nrow(res) # this defines the number of permutable objects
  yhat.2<-predict(x2)
  res.2<-resid(x2)
  
  mean.dif<-as.vector(dist(as.matrix(group.means(Y,g))))
  
  
  for(i in 1:p){
    
    #For dispersion tests
    
    line.rand<-sample(line,replace=FALSE) # This creates a line from 1 to the number
    # of permutable objects, then randomizes the order
    res.temp<-cbind(line.rand,res) # attaches the random order to the ordered matrix
    z<-(order(line.rand)) # randomizes matrix as it reorders line
    res.temp2<-as.matrix(res.temp[z,])# removes randomized line
    res.p<-res.temp2[,-1] # Rows of residuals are now randomized
    y.rand<-yhat+res.p # random values created
    
    # The resampling procedure above is the same as randomizing original values
    # (Since the reduced model only contains the overall mean)
    # But using residuals makes it applicable to multi-factor models
    
    lm.r<-lm(y.rand~g,x=T,model=T) # new linear model
    r<-resid(lm.r)
    yhat<-predict(lm.r)
    ex1.r<-ds.prep(r,g) #thanks to the functions, prep takes one step
    ex.r.ds<-disp.stat(ex1.r) # thanks to the function, dispersion stats require one step
    disp.dif.r<-ds.dif(ex.r.ds) # thanks to function, contrasts require one step
    
    
    # For means tests
    
    # uses different linear model
    # the null hypothesis that means are equal
    # means that an intercept model (i.e., defines only the overall mean)
    # is as viable as a group means model.
    # Thus, residuals are calculated from the intercept model (see above)
    # and random means are created despite no mean differences define by the model
    # this creates random distibutions of outcomes under the null hypothesis
    
    res.2.temp<-cbind(line.rand,res.2)
    z<-(order(line.rand))
    res.2.temp2<-as.matrix(res.2.temp[z,])
    res.2.p<-res.temp2[,-1] # Rows of residuals are now randomized
    y.rand.2<-yhat.2+res.2.p
    
    
    gm.r<-group.means(y.rand.2,g)
    md.r<-as.vector(dist(as.matrix(gm.r)))
    
    result<-c(i,md.r,disp.dif.r$mdc,disp.dif.r$mnn,disp.dif.r$ecc) # bind all results together
    p.table<-rbind(p.table,result) # add them to a table, row by row
    
  }
  
  head<-NULL # create a header
  line1<-as.vector(c(0,mean.dif,disp.dif$mdc,disp.dif$mnn,disp.dif$ecc))
  
  # The following is a bunch of code simply for generating column names in the output
  
  cn<-length(mean.dif) # cn = column name
  
  test.list<-NULL
  if (cn>1) for(i in 1:cn){
    l1<-rep(i,cn)
    l2<-array(1:cn)
    l12<-cbind(l1,l2)
    test.list<-rbind(test.list,l12)
  }
  
  test.list2<-NULL
  if (cn>1) for(j in 1:nrow(test.list)){
    t<-test.list[j,]
    if(t[2]>t[1]) test.list2<-rbind(test.list2,t)
  }
  
  test.list3<-NULL
  if (cn>1) for(k in 1:nrow(test.list2)){
    t<-test.list2[k,]
    lab<-paste(t[1],t[2],sep="--")
    test.list3<-rbind(test.list3,lab)}
  
  if (cn==1) test.list3<-c("1--2")
  
  lab2<-c(rep("MD",cn),rep("MDC",cn),rep("MNN",cn),rep("ECC",cn))
  test.list4<-paste(lab2,test.list3,sep=".")
  
  head<-c("iteration",test.list4)
  p.table<-rbind(head,line1,p.table)
  
}


test1<-permute(lm.gp,lm.gp.red,Y,group,999) # run the permutation test


# calculate P-value by using a function that ranks observed values as follows:

h<-test1[1,];h<-h[-1]
t<-test1[-1,];t<-t[,-1]
f<-function(t){
  r<-rank(t)
  p<-r[1]
  c<-length(t)
  pv<-(c-p+1)/c}
p.value<-apply(t,2,f)
names(p.value)<-h

# Check the p-value for the distance between SEAc centroids of the two groups
# this is the first p-value mentioned in the following output (MD.1--2)
p.value


