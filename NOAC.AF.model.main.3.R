# This research was funded by the National Institute for Health Research (NIHR) Health Technology Assessment programme project number 11/92/17 and NIHR Senior Investigator award NF-SI-0611-10168.
# Support for R in CEA provided by MRC Hubs for Trials Methodology Research ConDuCT-II (Collaboration and innovation in Difficult and Complex randomised controlled Trials In Invasive procedures) hub
# Howard Thom 31-October-2018. Bristol Medical School: Population Health Sciences. Bristol University, UK. howard.thom@bristol.ac.uk
# Details of model described in two publications:
# Sterne JA, Bodalia PN, Bryden PA, Davies PA, Lopez-Lopez JA, Okoli GN, et al. Oral anticoagulants for primary prevention, treatment and secondary prevention of venous thromboembolic disease, and for prevention of stroke in atrial fibrillation: systematic review, network meta-analysis and cost-effectiveness analysis. Health Technol Assess. 2017; 21(9):1-386.
# Lopez-Lopez JA, Sterne J, Thom H, Higgins J, Hingorani A, Okoli GN, et al. Oral anticoagulants for prevention of stroke in atrial fibrillation: systematic review, network meta-analysis, and cost effectiveness analysis. BMJ. 2017; 359(J5058).


# NOAC AF Cost-effectiveness model
# Main script to run the model
# Howard Thom 14-03-2015
# Updated for Wei Fang and Zhenru Wang to use noac.net.benefit() 29-05-2017

# Optimized to calculate cohort vectors for treatments in parallel (increases speed by factor of 3.4)
# Updated to sample n.cycle transition matrices for each age (very slow and memory intensive!)

# Updated to allow apixaban as reference treatment (set apixref=TRUE in the noac.net.benefit() function)


# Set to your baseline directory
baseline.directory<-"C:/Users/Howard/Documents/Bristol/R for CEA/DOACs model for GitHub"

library(expm)

setwd(paste(baseline.directory,"/code",sep=""))

n.samples<-1000
n.cycles<-120
initial.age<-70

# CE thresholds at which to generate the NB
lambdas<-c(1:50)*1000 

# Run requisite scripts and load functions
source("utility.functions.1.R")
source("generate.transition.matrix.7.R")
source("generate.transition.matrix.7.apixref.R")
source("NOAC.AF.net.benefit.3.R")

# Load outputs from WinBUGS for log hazard ratios, baseline log hazards, and no treamtent hazard ratios
bugs.loghr<-read.csv(file=paste(baseline.directory,"/data/bugs.loghr.csv",sep=""))
bugs.baseline<-read.csv(file=paste(baseline.directory,"/data/bugs.baseline.csv",sep=""))
hr.no.treatment<-read.csv(file=paste(baseline.directory,"/data/hr.no.treatment.csv",sep=""))
# Remove the empty first columns
bugs.loghr<-bugs.loghr[,-1]
bugs.baseline<-bugs.baseline[,-1]
hr.no.treatment<-hr.no.treatment[,-1]

# Set (fixed) NOAC drug costs. Allowed this to assist with threshold analysis
doac.costs<-c(200.42, 200.44, 200.44, 191.63)
names(doac.costs)<-treatment.names[2:5]

# Generate all model inputs
age.independent.samples<-age.independent.generate.probabilities(n.samples, doac.costs=doac.costs)

# Need the following for the threshold analysis
treatment.costs<-age.independent.samples$treatment.costs

# Generate the net benefits (And INB and CEAC)
model.outputs<-noac.net.benefit(n.samples=n.samples,n.cycles=n.cycles,initial.age=initial.age,lambdas=lambdas,age.independent.samples=age.independent.samples,
					apixref=FALSE,prob.treatment.samples=FALSE)

# Check if the parameterisation has been done correctly
#model.outputs.apixref<-noac.net.benefit(n.samples=n.samples,n.cycles=n.cycles,initial.age=initial.age,lambdas=lambdas,age.independent.samples=age.independent.samples,apixref=TRUE)
#model.outputs.standard<-noac.net.benefit(n.samples=n.samples,n.cycles=n.cycles,initial.age=initial.age,lambdas=lambdas,age.independent.samples=age.independent.samples,apixref=FALSE)


######################################################################################################
## Calculate results summaries #######################################################################

total.costs<-model.outputs$total.costs
total.qalys<-model.outputs$total.qalys
NB<-model.outputs$NB
INB<-model.outputs$INB

# Outputs for threshold analysis
prob.on.treatment<-model.outputs$prob.on.treatment
mean.years.treatment<-model.outputs$mean.years.treatment
prob.on.treatment.discounted<-model.outputs$prob.on.treatment.discounted
mean.years.treatment.discounted<-model.outputs$mean.years.treatment.discounted


# Export the results of the model to an rda file for later analysis 
#save(total.costs,total.qalys,cohort.vector,n.samples,sampled.probabilities,age.independent.samples,event.utilities,file=paste("model.run.",n.samples,".rda",sep=""))
# If you want to re-analyse previous data
#load(file=paste("model.run.",n.samples,".rda",sep=""))



# Calculate expected net benefit at lambda=20000

CEAC<-matrix(NA,n.treatments,length(lambdas))
	
# Results matrix presents incremental costs, incremental qalys, and INB
results.matrix<-matrix(" - (-, -)",nrow=6,ncol=(n.treatments-1))
colnames(results.matrix)<-treatment.names[1:(n.treatments-1)]
rownames(results.matrix)<-c("Costs","QALYs","Incremental Costs","Incremental QALYs","Incremental Net Benefit £20,000","Incremental Net Benefit £30,000")
incremental.nb2<-incremental.nb<-incremental.qalys<-incremental.costs<-matrix(NA,nrow=(n.treatments-1),ncol=3)
absolute.costs<-absolute.qalys<-matrix(NA,nrow=n.treatments,ncol=3)

rownames(incremental.nb)<-rownames(incremental.qalys)<-rownames(incremental.costs)<-treatment.names[-n.treatments]
# Evaluate at lambda=£20000
i.lambda<-20
# And at lambda=3000
i.lambda2<-30
for(i in 1:(n.treatments-1))
{
	absolute.costs[i,]<-c(mean(total.costs[,i]),quantile(total.costs[,i],probs=c(0.025,0.975)))
	absolute.qalys[i,]<-c(mean(total.qalys[,i]),quantile(total.qalys[,i],probs=c(0.025,0.975)))
	results.matrix[1,i]<-paste(format(absolute.costs[i,1],digits=4)," (",format(absolute.costs[i,2],digits=4),", ",format(absolute.costs[i,3],digits=4),") ",sep="")
	results.matrix[2,i]<-paste(format(absolute.qalys[i,1],digits=4)," (",format(absolute.qalys[i,2],digits=4),", ",format(absolute.qalys[i,3],digits=4),") ",sep="")

	if(i>1){
		incremental.costs[i,]<-c(mean(total.costs[,i]-total.costs[,1]),quantile(total.costs[,i]-total.costs[,1],probs=c(0.025,0.975)))
		incremental.qalys[i,]<-c(mean(total.qalys[,i]-total.qalys[,1]),quantile(total.qalys[,i]-total.qalys[,1],probs=c(0.025,0.975)))
		incremental.nb[i,]<-c(mean( lambdas[i.lambda]*(total.qalys[,i]-total.qalys[,1])-(total.costs[,i]-total.costs[,1])),quantile( lambdas[i.lambda]*(total.qalys[,i]-total.qalys[,1])-(total.costs[,i]-total.costs[,1]),probs=c(0.025,0.975)))
		incremental.nb2[i,]<-c(mean( lambdas[i.lambda2]*(total.qalys[,i]-total.qalys[,1])-(total.costs[,i]-total.costs[,1])),quantile( lambdas[i.lambda2]*(total.qalys[,i]-total.qalys[,1])-(total.costs[,i]-total.costs[,1]),probs=c(0.025,0.975)))
		results.matrix[3,i]<-paste(format(incremental.costs[i,1],digits=4)," (",format(incremental.costs[i,2],digits=4),", ",format(incremental.costs[i,3],digits=4),") ",sep="")
		results.matrix[4,i]<-paste(format(incremental.qalys[i,1],digits=4)," (",format(incremental.qalys[i,2],digits=4),", ",format(incremental.qalys[i,3],digits=4),") ",sep="")
		results.matrix[5,i]<-paste(format(incremental.nb[i,1],digits=4)," (",format(incremental.nb[i,2],digits=4),", ",format(incremental.nb[i,3],digits=4),") ",sep="")
		results.matrix[6,i]<-paste(format(incremental.nb2[i,1],digits=4)," (",format(incremental.nb2[i,2],digits=4),", ",format(incremental.nb2[i,3],digits=4),") ",sep="")
	}
}

# Threshold analysis for dabigatran
base.cost<-treatment.costs[1,3] # All the same as no randomness
#new.cost<-321.1632 # CE compared to Warfarin?
# New inb would be (solve this for different conditions):
#threshold.inb[3]<-threshold.inb[3]-(new.cost-base.cost)*mean.years.dabigatran

# Expected incremental net benefit
threshold.inb<-incremental.nb[,1]
# When is other NOAC cost-effective compared to warfarin (inb=0) (how high can price be?)
ce.warfarin.cost<-base.cost+threshold.inb[]/mean.years.treatment.discounted
# When is other NOACs inb higher than apixaban inb (inb[3]=inb[2]) (how low does price have to be?)
ce.apixaban.cost<-base.cost-(threshold.inb[2]-threshold.inb[])/mean.years.treatment.discounted



#write.csv(file=paste(baseline.directory,"/results/model.results.apixref.csv",sep=""),results.matrix)

# Summarise the probability of each type of event for the different treatments

#well.prob.summary<-matrix(0,nrow=n.treatments,ncol=n.events)
#colnames(well.prob.summary)<-event.names
#rownames(well.prob.summary)<-paste(treatment.names,"Well")
#for(i.treatment in 1:n.treatments)well.prob.summary[i.treatment,]<-colMeans(sampled.probabilities$probability.matrix[[1]][,(i.treatment-1)*16+1,])
#stroke.prob.summary<-matrix(0,nrow=n.treatments,ncol=n.events)
#colnames(stroke.prob.summary)<-event.names
#rownames(stroke.prob.summary)<-paste(treatment.names,"Stroke")
#for(i.treatment in 1:n.treatments)stroke.prob.summary[i.treatment,]<-colMeans(sampled.probabilities$probability.matrix[,(i.treatment-1)*16+5,])
#bleed.prob.summary<-matrix(0,nrow=n.treatments,ncol=n.events)
#colnames(bleed.prob.summary)<-event.names
#rownames(bleed.prob.summary)<-paste(treatment.names,"Bleed")
#for(i.treatment in 1:n.treatments)bleed.prob.summary[i.treatment,]<-colMeans(sampled.probabilities$probability.matrix[,(i.treatment-1)*16+2,])

#prob.summary<-rbind(well.prob.summary,stroke.prob.summary,bleed.prob.summary)

#write.csv(file=paste(baseline.directory,"/results/prob.summary.csv",sep=""),prob.summary)


# Calculate the EVPI for a range of lambda=20000
# Expected max NB - max Expected NB
EVPI<-rep(NA,length(lambdas))
EVPI.pop<-rep(NA,length(lambdas))
NB.max<-matrix(NA,nrow=n.samples,ncol=length(lambdas))
for(i.lambda in 1:length(lambdas))
{

	NB.which.max<-apply(NB[,,i.lambda],c(1),which.max)
	for(i.sample in 1:n.samples)NB.max[i.sample,i.lambda]<-NB[i.sample,NB.which.max[i.sample],i.lambda]
	EVPI[i.lambda]<-mean(NB.max[,i.lambda])-mean(NB[,which.max(colMeans(NB[,,i.lambda],na.rm=TRUE)),i.lambda])
	EVPI.pop[i.lambda]<-5000*sum(EVPI[i.lambda]*((1/1.035)^c(1:n.cycles)))
}

#jpeg(file=paste(baseline.directory,"/results/model.evpi.jpg",sep=""))
plot(c(0,0),xlim=c(min(lambdas),max(lambdas)),ylim=c(0,max(EVPI)),col=0,xlab="Willngess-to-pay (£)",ylab="Per-person EVPI (£)")
lines(lambdas,EVPI)
#dev.off()

# Calculate the probability that each treatment is most cost-effective and then
# plot the cost-effectiveness acceptability curve
for(i.treatment in 2:(n.treatments-1))
{
for(i.lambda in 1:length(lambdas))
{
	CEAC[i.treatment,i.lambda]<-mean(apply(INB[,,i.lambda],c(1),which.max)==(i.treatment-1))
}
}
# Calculate the probability that Warfarin (Coumarin) is the most cost-effective

for(i.lambda in 1:length(lambdas))
{
# Probability all treatments have negative INB
CEAC[1,i.lambda]<-1-mean(rowSums(INB[,,i.lambda]>0,na.rm=TRUE)>0)
CEAC[2:(n.treatments-1),i.lambda]<-CEAC[2:(n.treatments-1),i.lambda]*mean(rowSums(INB[,,i.lambda]>0,na.rm=TRUE)>0)
}

jpeg(file=paste(baseline.directory,"/results/model.ceac.jpg",sep=""))

plot(c(0,0),col=0,xlim=c(0,max(lambdas)),ylim=c(0,1.3),xlab="Willingness-to-pay (£)",ylab="Probability most cost-effective",yaxt="n")
axis(side=2,at=c(0,0.2,0.4,0.6,0.8,1.0))
for(i.treatment in 1:(n.treatments-1))
{
	lines(x=lambdas,y=CEAC[i.treatment,],lty=i.treatment,col=i.treatment,lwd=2)
}
legend("topright",c("Warfarin (INR 2-3)",treatment.names[2:(n.treatments-1)]),lty=1:(n.treatments-1),col=1:(n.treatments-1),lwd=2)

dev.off()

# Plot the Cost-effectiveness Acceptability Frontier Curve

# Which treatment has highest Expected NB at each threshold
expected.NB.which.max<-apply((apply(NB[,,],c(2,3),mean)),c(2),which.max)
CEAF<-rep(NA, length(lambdas))
for(i in 1:length(lambdas))
{
	CEAF[i]<-CEAC[expected.NB.which.max[i],i]
}

jpeg(file=paste(baseline.directory,"/results/model.ceaf.jpg",sep=""))

plot(c(0,0),col=0,xlim=c(0,max(lambdas)),ylim=c(0,1),xlab="Willingness-to-pay (£)",ylab="Probability most cost-effective")
for(i in 2:length(lambdas))
{
	lines(x=lambdas[(i-1):i],y=CEAF[(i-1):i],lty=expected.NB.which.max[i],col=expected.NB.which.max[i],lwd=2)
}
legend("topright",c("Warfarin (INR 2-3)",treatment.names[2:(n.treatments-1)]),lty=1:(n.treatments-1),col=1:(n.treatments-1),lwd=2)

dev.off()

# Plot the cost-effectiveness plane

all.incremental.costs<-total.costs-total.costs[,1]
all.incremental.qalys<-total.qalys-total.qalys[,1]

jpeg(file=paste(baseline.directory,"/results/model.ce.plane.jpg",sep=""))
scaling=1.75 # To ensure legend doesn't overlap with points
options(scipen=5)
plot(c(0,0),col=0,ylim=range(all.incremental.costs)*c(1,scaling),xlim=range(all.incremental.qalys)*c(1,2),ylab="Incremental costs (£)",xlab="Incremental QALYs")
# Include X and Y axes
lines(c(0,0),range(all.incremental.costs)*c(1,scaling))
lines(range(all.incremental.qalys)*c(1,2),c(0,0))
for(i in 2:(n.treatments-1))
{
	points(all.incremental.qalys[,i],all.incremental.costs[,i],col=i,pch=i) # Use pch=19 for simple dots
}
legend("topleft",c(treatment.names[2:(n.treatments-1)]),pch=2:(n.treatments-1),col=2:(n.treatments-1))
dev.off()
