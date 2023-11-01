#############################################################################################
source("https://neuroconductor.org/neurocLite.R")
neuro_install("neurobase", release = "stable")
library(neurobase)
library(neuRosim)
library(ggplot2)
library(mosaic)
library(Hmisc)

##############################################################################################
# Reading and visualising fMRI data with R/neuroconductor.
setwd("C:/Users/6000600/Documents/Fulbright-Research/New_Course/Chapter6- fMRI-GLM/Scripts")
data <- readnii("sub-10159_task-rest_bold_space-MNI152NLin2009cAsym_preproc.nii.gz")
#This data was obtained from the OpenfMRI database. Its accession number is ds000030
# Data: UCLA Consortium for Neuropsychiatric Phenomics
# https://legacy.openfmri.org/dataset/ds000030/

dim(data)
#This shows us that the first three dimensions are 65 x 77 x 49. These correspond to the x-, y- and z-coordinates of the brain. The fourth dimension is time, which means there are 152 time points measured.

check_nifti(data)
#This shows some information that is stored in the “header” of the file. For example, it tells us the pixel dimensions. The first three numbers tell us the size of the voxels in space: 3 x 3 x 4 cm. The last number tells us the dimension in time: a scan was taken every 2 seconds.

#Visualisation of nifti’s. We only visualise one timepoint. R will automatically show the first timepoint.

orthographic(data)


dim(data)

check_nifti(data)
#3 x 3 x 4 cm x 2 seconds

data[20,20,30,1]
orthographic(data,xyz=c(20,20,30))
neurobase::ortho2(data,add.orient = FALSE)

##############################################################################################
#Visualization
#ts
data[20,20,30,]
plot(data[20,20,30,],type='l')

ortho2(data, y = data > quantile(data, 0.9))

ortho2(data, zlim = quantile(data, probs = c(0, 0.999)))


ortho2(data, y = data > 1000)

double_ortho(data,y = data > quantile(data, 0.9) ) 


image(data, z = 45)

image(data, z = 45, plot.type = "single", plane = "sagittal")


##############################################################################################
vals = c(data)
class(vals)
plot(density(vals))

#subset of voxels
timeseries= data.frame(t(data[20,20,c(18:37),]))
colnames(timeseries)=paste0("V",1:20)
head(timeseries)
str(timeseries)

##############################################################################################
# EV (Design matrix)
# Assumption: The resting state data is coming from a task experiment with a random design matrix. 
# Ref: Cluster failure: Why fMRI inferences for spatial extent have inflated false-positive rates
# Anders Eklund et al 2016; https://doi.org/10.1073/pnas.1602413113
# Assume we are designing an fMRI experiment, with a blocked design (10 seconds on / 10 seconds 
#off). We can use neuRosim to create our design.
totaltime <- 152*2
onsets <- seq(1,totaltime,40)
dur <- 20
s <- specifydesign(totaltime=totaltime,onsets=list(onsets),durations=dur,
                   accuracy=1,effectsize=1,TR=2)
plot(s,type='l')

#The signal that we expect looks like this:
s <- specifydesign(totaltime=totaltime,onsets=list(onsets),durations=dur,
                   accuracy=1,effectsize=1,TR=2,conv="double-gamma")
plot(s,type='l')


##############################################################################################
# GLM: Voxel level glm
# Is a particular voxel significantly related to the task?
plot(zscore(timeseries$V2),type="l")
lines(s,col=2)

timeseries$design <- s

model1 <- lm("zscore(V1) ~ design",data=timeseries)
summary(model1)


model11 <- lm("zscore(V11) ~ design",data=timeseries)
summary(model11)


##############################################################################################
# Loop
pvals <- c()
tvals <- c()

for (i in 1:(dim(timeseries)[2]-1)) {
  model <- lm(zscore(timeseries[,i]) ~ design,data=timeseries)
  pvals[i] <- summary(model)$coefficients[2,4]
  tvals[i] <- summary(model)$coefficients[2,3]
}

sum(pvals<0.05)/length(pvals)

sum(pvals<0.05/dim(timeseries)[2])/length(pvals)


stats=data.frame(cbind(as.double(tvals),pvals))
stats$sig=ifelse(stats$pvals<0.05, "sig", "no-sig")

  
ggplot(stats, aes(x=tvals)) + 
  geom_histogram(aes(y=..density..), colour="black", fill="white")+
  geom_density(alpha=.2, fill="#FF6666") 
# Color by groups
ggplot(stats, aes(x=tvals, color=sig, fill=sig)) + 
  geom_histogram(aes(y=..density..), alpha=0.5, 
                 position="identity")+
  geom_density(alpha=.2) 


hist(tvals,freq=FALSE,ylim=c(0,0.8))
x <- seq(-6,6,length=1000)
y <- dt(x,152)
lines(x,y)

#################################################################################
#Practice
#What happens when you choose a more random design (onsets not evenly spaced)?
#  Do you observe the same pattern in another?
