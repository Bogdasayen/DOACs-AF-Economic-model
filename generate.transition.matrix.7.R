# This research was funded by the National Institute for Health Research (NIHR) Health Technology Assessment programme project number 11/92/17 and NIHR Senior Investigator award NF-SI-0611-10168.
# Support for R in CEA provided by MRC Hubs for Trials Methodology Research ConDuCT-II (Collaboration and innovation in Difficult and Complex randomised controlled Trials In Invasive procedures) hub
# Howard Thom 31-October-2018. Bristol Medical School: Population Health Sciences. Bristol University, UK. howard.thom@bristol.ac.uk
# Details of model described in two publications:
# Sterne JA, Bodalia PN, Bryden PA, Davies PA, Lopez-Lopez JA, Okoli GN, et al. Oral anticoagulants for primary prevention, treatment and secondary prevention of venous thromboembolic disease, and for prevention of stroke in atrial fibrillation: systematic review, network meta-analysis and cost-effectiveness analysis. Health Technol Assess. 2017; 21(9):1-386.
# Lopez-Lopez JA, Sterne J, Thom H, Higgins J, Hingorani A, Okoli GN, et al. Oral anticoagulants for prevention of stroke in atrial fibrillation: systematic review, network meta-analysis, and cost effectiveness analysis. BMJ. 2017; 359(J5058).



# NOAC AF Cost-effectiveness model
# Code to build a transition matrix for the AF Economic model

# Updated from version 6 so that state costs and utilities, "Kind age factor" and the hr.event.history
# are all sampled in the age.independent.samples() function

# Version 7 (among other fixes) allows apixaban to be reference and sets DOAC drug costs externally for threshold analysis


# Changed from generate.transition.matrix.4.R to include betrixaban and edoxaban

# There are 17 health states in the model: AF well, Stroke, Bleed, MI, ICH, S+B,
# S+MI, S+ICH, B+MI, B+ICH, MI+ICH, S+B+MI, S+B+ICH, S+MI+ICH, B+MI+ICH, 
# S+B+MI+ICH, Dead
# The 16 non-death health states are repeated for the five treatments (Rivaroxaban, 
# dabigatran, apixaban, warfarin, betrixaban 40mg od, edoxaban 60mg od, no treatment), giving a total of 81 states.

# At any stage, a patient can experience a stroke, clinically relevant bleed, MI, ICH, SE,
# TIA. The last two are transient events.


load(file=paste(baseline.directory,"/data/bugs.objects.rda",sep=""))

# Load the effect of prior events on future Ischemic strokes, TIA/SE, ICH and bleeds
hr.future.stroke<-read.csv(file=paste(baseline.directory,"/data/hr.event.history.Ischemic.stroke.csv",sep=""))
hr.future.tiase<-read.csv(file=paste(baseline.directory,"/data/hr.event.history.TIA.SE.csv",sep=""))
hr.future.ich<-read.csv(file=paste(baseline.directory,"/data/hr.event.history.ICH.csv",sep=""))
hr.future.bleed<-read.csv(file=paste(baseline.directory,"/data/hr.event.history.Bleed.csv",sep=""))
# Only limited evidence for effect on death
hr.future.death<-read.csv(file=paste(baseline.directory,"/data/hr.event.history.Death.csv",sep=""))
# Name the rows for easy access
rownames(hr.future.death)<-rownames(hr.future.bleed)<-rownames(hr.future.ich)<-rownames(hr.future.tiase)<-rownames(hr.future.stroke)<-hr.future.stroke[,1]

lifetables<-read.csv(file=paste(baseline.directory,"/data/EnglandWales-LifeTables-2011-13.csv",sep=""))

# Utility functions
source("generate.hr.death.age.1.R")
source("next.state.name.R")


# Define the event names and codes for ease of indexing
event.names<-c("MI","Ischemic stroke","Death (all causes)","Transient ischemic attack (TIA)","Clinically relevant bleeding","SE","ICH","Stay")
event.codes<-which(is.element(names(outcome.codes),event.names))
# Fixed event codes
event.codes<-c(5,1,6,7,17,15,16)
n.events<-length(event.names)

event.state.codes<-c(" M "," S "," Dead ",""," B ",""," I ","")
event.codes.alpha<-c(" B "," I "," M "," S ")


# Construct a list of state names for the rows and columns of the matrix
nondeath.health.states<-c("Well"," B "," I "," M "," S "," B + I "," B + M "," B + S ",
" I + M "," I + S "," M + S "," B + I + M "," B + I + S "," B + M + S ",
" I + M + S "," B + I + M + S ")
treatment.names<-c("Coumarin (INR 2-3)","Apixaban (5mg bd)","Dabigatran (150mg bd)","Edoxaban (60mg od)","Rivaroxaban (20mg od)","No treatment")
n.treatments<-length(treatment.names)
n.health.states<-length(nondeath.health.states)
n.states<-n.health.states*n.treatments+1
state.names<-rep("",n.states)
for(i in 1:n.treatments){
	for(j in 1:n.health.states){
		state.names[(i-1)*n.health.states+j]<-paste(treatment.names[i],nondeath.health.states[j])
	}
}
state.names[n.states]<-"Dead"

# Define number of indices to move if a patient switches treatment (non-ICH and non-MI)#
# NOAC - > Warfarin
# Warfarin - > No treatment
treatment.switch.indices<-c((n.states-16-1),-16*c(1:(n.treatments-2)),0)
# If ICH, NOAC/Warfarin - > No treatment
ich.treatment.switch.indices<-c(n.states-1-16*c(1:(n.treatments-1)),0)
# If MI, switch to Warfarin if on dabigatran, no switch for all other treatments
mi.treatment.switch.indices<-rep(0,n.treatments)
mi.treatment.switch.indices[treatment.names=="Dabigatran (150mg bd)"]<--32
names(treatment.switch.indices)<-names(ich.treatment.switch.indices)<-names(mi.treatment.switch.indices)<-treatment.names



# Load the life-tables from 2011-13 for England and Wales
# Probability of death is averaged over males and females

hr.death.age<-generate.hr.death.age(base.age=initial.age,lifetables=lifetables)

# Function to generate hazards and other probabilities that are age independent
age.independent.generate.probabilities<-function(n.samples, doac.costs=doac.costs)
{
	# bugs.object.baseline$sims.array[i.sample,1,j.outcome] is the baseline
	# log hazard for the j.outcome event
	# bugs.object.fixed$sims.array[i.sample,1,(k.treat-2)*17+j.outcome] is
	# the log hazard ratio of treatment k.treat relative to Coumarin (INR 2-3)

	# Mean costs for now
	event.costs<-t(matrix(c(2956,10844,0,1064,1751.5,2373,24234,0),ncol=n.samples,nrow=n.events))
	# SE acute costs
	event.costs[,6]<-runif(n.samples,1186.5,3559.5)
	# TIA acute costs
	event.costs[,4]<-runif(n.samples,532,1596)
	# Bleed acute costs
	event.costs[,5]<-runif(n.samples,875.75,2627.25)
	# MI acute costs
	event.costs[,1]<-runif(n.samples,2415.24,7245.72)
	# Ischemic stroke acute costs
	S.acute.cost.mean=11626
	S.acute.cost.SD=16868
	S.acute.cost.SE=S.acute.cost.SD/sqrt(162)
	event.costs[,2]<-rnorm(n.samples,mean=S.acute.cost.mean,sd=S.acute.cost.SE)
	# ICH acute costs
	I.acute.cost.mean=11453
	I.acute.cost.SD=13815
	I.acute.cost.SE=I.acute.cost.SD/sqrt(17)
	event.costs[,7]<-rnorm(n.samples,mean=I.acute.cost.mean,sd=I.acute.cost.SE)

	colnames(event.costs)<-event.names

	# Probability patient will switch given they have each type of event
	event.switch.probs<-matrix(0,nrow=n.samples,ncol=length(event.names))
	colnames(event.switch.probs)<-event.names
	# Only switch after MI if on Dabigatran
	event.switch.probs[,"MI"]<-rep(1,n.samples)
	is.disc.params<-list("alpha"=0.1,"beta"=0.9) #beta.parameters(0.3,0.7) # (0.5,0.25)
	event.switch.probs[,"Ischemic stroke"]<-rbeta(n.samples,is.disc.params$alpha,is.disc.params$beta)
	b.disc.params<-list("alpha"=0.3,"beta"=0.7) #beta.parameters(0.3,0.7) #(0.5,0.4)
	event.switch.probs[,"Clinically relevant bleeding"]<-rbeta(n.samples,b.disc.params$alpha,b.disc.params$beta)
	tia.disc.params<-list("alpha"=0.1,"beta"=0.9) #beta.parameters(0.1,0.9) #(0.05,0.1)
	event.switch.probs[,"Transient ischemic attack (TIA)"]<-rbeta(n.samples,tia.disc.params$alpha,tia.disc.params$beta)
	se.disc.params<-list("alpha"=0.1,"beta"=0.9) #beta.parameters(0.1,0.9) #(0.05,0.1)
	event.switch.probs[,"SE"]<-rbeta(n.samples,se.disc.params$alpha,se.disc.params$beta)
	event.switch.probs[,"ICH"]<-rep(1,n.samples)
	event.switch.probs[,"Stay"]<-rep(0,n.samples)
	event.switch.probs[,"Death (all causes)"]<-rep(0,n.samples)

	# Hazard ratios for effect of previous events on future events
	hr.event.history<-array(1,dim=c(n.samples,n.events,n.events-1))
	# Effect of prior bleeds on future events
	hr.event.history[,event.state.codes==" B ",event.names[-n.events]=="MI"]<-1 # No evidence
	hr.event.history[,event.state.codes==" B ",event.names[-n.events]=="Ischemic stroke"]<-exp(rnorm(n=n.samples,mean=hr.future.stroke["Bleed","Logmean"],sd=hr.future.stroke["Bleed","logSD"]))
	# Had to asssume effect on death same as that of stroke
	hr.event.history[,event.state.codes==" B ",event.names[-n.events]=="Death (all causes)"]<-exp(rnorm(n=n.samples,mean=hr.future.death["Stroke","Logmean"],sd=hr.future.death["Stroke","logSD"]))
	hr.event.history[,event.state.codes==" B ",event.names[-n.events]=="Transient ischemic attach (TIA)"]<-exp(rnorm(n=n.samples,mean=hr.future.tiase["Bleed","Logmean"],sd=hr.future.tiase["Bleed","logSD"]))
	hr.event.history[,event.state.codes==" B ",event.names[-n.events]=="Clinically relevant bleeding"]<-exp(rnorm(n=n.samples,mean=hr.future.bleed["Bleed","Logmean"],sd=hr.future.bleed["Bleed","logSD"]))
	hr.event.history[,event.state.codes==" B ",event.names[-n.events]=="SE"]<-exp(rnorm(n=n.samples,mean=hr.future.tiase["Bleed","Logmean"],sd=hr.future.tiase["Bleed","logSD"]))
	hr.event.history[,event.state.codes==" B ",event.names[-n.events]=="ICH"]<-exp(rnorm(n=n.samples,mean=hr.future.ich["Bleed","Logmean"],sd=hr.future.ich["Bleed","logSD"]))
	# Effect of prior Intracranial hemorrhage of future events
	hr.event.history[,event.state.codes==" I ",event.names[-n.events]=="MI"]<-1
	hr.event.history[,event.state.codes==" I ",event.names[-n.events]=="Ischemic stroke"]<-exp(rnorm(n=n.samples,mean=hr.future.stroke["ICH","Logmean"],sd=hr.future.stroke["ICH","logSD"]))
	# Had to assume effect on death the same as that of stroke
	hr.event.history[,event.state.codes==" I ",event.names[-n.events]=="Death (all causes)"]<-exp(rnorm(n=n.samples,mean=hr.future.death["Stroke","Logmean"],sd=hr.future.death["Stroke","logSD"]))
	hr.event.history[,event.state.codes==" I ",event.names[-n.events]=="Transient ischemic attach (TIA)"]<-exp(rnorm(n=n.samples,mean=hr.future.tiase["ICH","Logmean"],sd=hr.future.tiase["ICH","logSD"]))
	hr.event.history[,event.state.codes==" I ",event.names[-n.events]=="Clinically relevant bleeding"]<-exp(rnorm(n=n.samples,mean=hr.future.bleed["ICH","Logmean"],sd=hr.future.bleed["ICH","logSD"]))
	hr.event.history[,event.state.codes==" I ",event.names[-n.events]=="SE"]<-exp(rnorm(n=n.samples,mean=hr.future.tiase["ICH","Logmean"],sd=hr.future.tiase["ICH","logSD"]))
	hr.event.history[,event.state.codes==" I ",event.names[-n.events]=="ICH"]<-exp(rnorm(n=n.samples,mean=hr.future.ich["ICH","Logmean"],sd=hr.future.ich["ICH","logSD"]))
	# Effect of prior MI on future events
	hr.event.history[,event.state.codes==" M ",event.names[-n.events]=="MI"]<-1
	hr.event.history[,event.state.codes==" M ",event.names[-n.events]=="Ischemic stroke"]<-exp(rnorm(n=n.samples,mean=hr.future.stroke["MI","Logmean"],sd=hr.future.stroke["MI","logSD"]))
	hr.event.history[,event.state.codes==" M ",event.names[-n.events]=="Death (all causes)"]<-exp(rnorm(n=n.samples,mean=hr.future.death["MI","Logmean"],sd=hr.future.death["MI","logSD"]))
	hr.event.history[,event.state.codes==" M ",event.names[-n.events]=="Transient ischemic attach (TIA)"]<-exp(rnorm(n=n.samples,mean=hr.future.tiase["MI","Logmean"],sd=hr.future.tiase["MI","logSD"]))
	hr.event.history[,event.state.codes==" M ",event.names[-n.events]=="Clinically relevant bleeding"]<-exp(rnorm(n=n.samples,mean=hr.future.bleed["MI","Logmean"],sd=hr.future.bleed["MI","logSD"]))
	hr.event.history[,event.state.codes==" M ",event.names[-n.events]=="SE"]<-exp(rnorm(n=n.samples,mean=hr.future.tiase["MI","Logmean"],sd=hr.future.tiase["MI","logSD"]))
	hr.event.history[,event.state.codes==" M ",event.names[-n.events]=="ICH"]<-exp(rnorm(n=n.samples,mean=hr.future.ich["MI","Logmean"],sd=hr.future.ich["MI","logSD"]))
	# Effect of prior strokes on future events
	hr.event.history[,event.state.codes==" S ",event.names[-n.events]=="MI"]<-1
	hr.event.history[,event.state.codes==" S ",event.names[-n.events]=="Ischemic stroke"]<-exp(rnorm(n=n.samples,mean=hr.future.stroke["Stroke","Logmean"],sd=hr.future.stroke["Stroke","logSD"]))
	hr.event.history[,event.state.codes==" S ",event.names[-n.events]=="Death (all causes)"]<-exp(rnorm(n=n.samples,mean=hr.future.death["Stroke","Logmean"],sd=hr.future.death["Stroke","logSD"]))
	hr.event.history[,event.state.codes==" S ",event.names[-n.events]=="Transient ischemic attach (TIA)"]<-exp(rnorm(n=n.samples,mean=hr.future.tiase["Stroke","Logmean"],sd=hr.future.tiase["Stroke","logSD"]))
	hr.event.history[,event.state.codes==" S ",event.names[-n.events]=="Clinically relevant bleeding"]<-exp(rnorm(n=n.samples,mean=hr.future.bleed["Stroke","Logmean"],sd=hr.future.bleed["Stroke","logSD"]))
	hr.event.history[,event.state.codes==" S ",event.names[-n.events]=="SE"]<-exp(rnorm(n=n.samples,mean=hr.future.tiase["Stroke","Logmean"],sd=hr.future.tiase["Stroke","logSD"]))
	hr.event.history[,event.state.codes==" S ",event.names[-n.events]=="ICH"]<-exp(rnorm(n=n.samples,mean=hr.future.ich["Stroke","Logmean"],sd=hr.future.ich["Stroke","logSD"]))	
	

	# Proportional utility decrements from Kind et al. 1999
	# Beta distributions for each age range estimated by Pete Bryden
	# Use weighted average of males (60%) and females (40%) to represent AF pop
	# Ratio of each category to 70 year olds is the proportional decrement/increment

	# Matrix containing absolute utilities
	kind.age.utility<-matrix(NA,nrow=n.samples,ncol=8)
	colnames(kind.age.utility)<-c(35,45,55,65,75,85,95,105)

	kind.age.utility[,"35"]<-0.6*rbeta(n.samples,656.7,64.95)+0.4*rbeta(n.samples,1006.6,99.5)
	kind.age.utility[,"45"]<-0.6*rbeta(n.samples,341.41,65.03)+0.4*rbeta(n.samples,544.1,96.02)
	kind.age.utility[,"55"]<-0.6*rbeta(n.samples,330.43,93.2)+0.4*rbeta(n.samples,526.59,123.52)
	kind.age.utility[,"65"]<-0.6*rbeta(n.samples,388.47,109.57)+0.4*rbeta(n.samples,551.74,155.62)
	kind.age.utility[,"75"]<-0.6*rbeta(n.samples,191.17,63.72)+0.4*rbeta(n.samples,406.37,165.98)

	# Assume absolute utilities for ages 90, 100 and 110 are the same as 80 (not true)
	kind.age.utility[,"85"]<-kind.age.utility[,"95"]<-kind.age.utility[,"105"]<-kind.age.utility[,"75"]

	age.utility.factor<-kind.age.utility/kind.age.utility[,4]


	# Treatment.costs costs (all but Warfarin costs are fixed)
	treatment.costs<-matrix(0,nrow=n.samples,ncol=n.treatments)
	# Uniform distribution on Warfarin costs for now
	treatment.costs[,1]<-runif(n.samples,52.57,157.70)
	treatment.costs[,2]<-doac.costs["Apixaban (5mg bd)"] #200.42
	treatment.costs[,3]<-doac.costs["Dabigatran (150mg bd)"] #200.44
	treatment.costs[,4]<-doac.costs["Edoxaban (60mg od)"] #200.44
	treatment.costs[,5]<-doac.costs["Rivaroxaban (20mg od)"] #191.63
	names(treatment.costs)<-treatment.names

	# Health state costs (divided by four to go from annual to quarterly costs)
	# Stroke
	S.cost.mean=3613
	S.cost.SD=4235
	S.cost.SE=S.cost.SD/sqrt(136)
	S.cost<-rnorm(n.samples,mean=S.cost.mean,sd=S.cost.SE)
	# ICH (Assume it is similar to stroke)
	I.cost.mean=3613
	I.cost.SD=4235
	I.cost.SE=I.cost.SD/sqrt(136)
	I.cost<-rnorm(n.samples,mean=I.cost.mean,sd=I.cost.SE)
	# MI adds only an instant cost, so this post-state has 0 management cost
	M.cost<-rep(0,n.samples)

	# Health state utilities (from sources identified in Bayer Table 49)
	# These are combined proportionally
	# All utilities are later divided by 4 to make them 3-monthly
	AF.utility<-rnorm(n.samples,0.779,0.0045)
	S.utility<-rnorm(n.samples,0.69,0.025)/AF.utility
	M.utility<-rnorm(n.samples,0.718,0.0163)/AF.utility
	I.utility<-rbeta(n.samples,3.941,1.385)/AF.utility
	# Need evidence for post bleed utility (for now assume same as Stroke)
	B.utility<-S.utility

	# Cap the utilities at 1
	AF.utility[AF.utility>1]<-1
	S.utility[S.utility>1]<-1
	M.utility[M.utility>1]<-1
	I.utility[I.utility>1]<-1
	B.utility[B.utility>1]<-1


	state.utilities<-state.costs<-matrix(0,nrow=n.samples,ncol=length(state.names))
	colnames(state.utilities)<-colnames(state.costs)<-state.names
	for(i.treatment in 1:n.treatments){
		for(i.health.state in 1:n.health.states){
			state.costs[,(i.treatment-1)*n.health.states+i.health.state]<-treatment.costs[,i.treatment]
		}
	}

	# Use a for loop to add effects of previous events (should use something faster)
	for(i.sample in 1:n.samples)
	{
		# History of 1 event
		# Add MI costs 
		state.costs[i.sample,gregexpr(" M ",state.names)!=-1 & gregexpr(" I ",state.names)==-1 & gregexpr(" S ",state.names)==-1]<-state.costs[i.sample,gregexpr(" M ",state.names)!=-1 & gregexpr(" I ",state.names)==-1 & gregexpr(" S ",state.names)==-1]+M.cost[i.sample]
		# Add ICH costs
		state.costs[i.sample,gregexpr(" M ",state.names)==-1 & gregexpr(" I ",state.names)!=-1 & gregexpr(" S ",state.names)==-1]<-state.costs[i.sample,gregexpr(" M ",state.names)==-1 & gregexpr(" I ",state.names)!=-1 & gregexpr(" S ",state.names)==-1]+I.cost[i.sample]
		# Add stroke costs 
		state.costs[i.sample,gregexpr(" M ",state.names)==-1 & gregexpr(" I ",state.names)==-1 & gregexpr(" S ",state.names)!=-1]<-state.costs[i.sample,gregexpr(" M ",state.names)==-1 & gregexpr(" I ",state.names)==-1 & gregexpr(" S ",state.names)!=-1]+S.cost[i.sample]
		# History of 2 events
		# Add bleed costs (excluding those that are both ICH)
		# Add max of stroke or ICH cost to states with history of both events but not Bleed
		state.costs[i.sample,gregexpr(" M ",state.names)==-1 & gregexpr(" I ",state.names)!=-1 & gregexpr(" S ",state.names)!=-1]<-state.costs[i.sample,gregexpr(" M ",state.names)==-1 & gregexpr(" I ",state.names)!=-1 & gregexpr(" S ",state.names)!=-1]+max(I.cost[i.sample],S.cost[i.sample])
		# Add max of MI or ICH cost to states with history of both events but not Stroke
		state.costs[i.sample,gregexpr(" M ",state.names)!=-1 & gregexpr(" I ",state.names)!=-1 & gregexpr(" S ",state.names)==-1]<-state.costs[i.sample,gregexpr(" M ",state.names)!=-1 & gregexpr(" I ",state.names)!=-1 & gregexpr(" S ",state.names)==-1]+max(M.cost[i.sample],I.cost[i.sample])
		# Add max of MI or stroke cost to states with history of both events but not ICH
		state.costs[i.sample,gregexpr(" M ",state.names)!=-1 & gregexpr(" I ",state.names)==-1 & gregexpr(" S ",state.names)!=-1]<-state.costs[i.sample,gregexpr(" M ",state.names)!=-1 & gregexpr(" I ",state.names)==-1 & gregexpr(" S ",state.names)!=-1]+max(M.cost[i.sample],S.cost[i.sample])
		# History of 3 events
		# Add max of MI, stroke, or ICH costs
		state.costs[i.sample,gregexpr(" M ",state.names)!=-1 & gregexpr(" I ",state.names)!=-1 & gregexpr(" S ",state.names)!=-1]<-state.costs[i.sample,gregexpr(" M ",state.names)!=-1 & gregexpr(" I ",state.names)!=-1 & gregexpr(" S ",state.names)!=-1]+max(M.cost[i.sample],I.cost[i.sample],S.cost[i.sample])

		# Utilities are multiplicative		
		state.utilities[i.sample,]<-AF.utility[i.sample]
		state.utilities[i.sample,gregexpr(" B ",state.names)!=-1]<-state.utilities[i.sample,gregexpr(" B ",state.names)!=-1]*B.utility[i.sample]
		state.utilities[i.sample,gregexpr(" I ",state.names)!=-1]<-state.utilities[i.sample,gregexpr(" I ",state.names)!=-1]*I.utility[i.sample]
		state.utilities[i.sample,gregexpr(" M ",state.names)!=-1]<-state.utilities[i.sample,gregexpr(" M ",state.names)!=-1]*M.utility[i.sample]
		state.utilities[i.sample,gregexpr(" S ",state.names)!=-1]<-state.utilities[i.sample,gregexpr(" S ",state.names)!=-1]*S.utility[i.sample]
	}


	# Costs and utilities for the death state should be zero
	state.utilities[,"Dead"]<-0
	state.costs[,"Dead"]<-0

	# Disutilities
	# Old dummy disutilities
	#event.utilities=t(matrix(c(-0.1,-0.1,0,-0.01,-0.03,-0.01,-0.1,0),ncol=n.samples,nrow=n.events))
	event.utilities<-matrix(0,nrow=n.samples,ncol=n.events)
	colnames(event.utilities)<-event.names
	event.utilities[,"MI"]<-rnorm(n.samples,0.683,0.0156)-AF.utility
	event.utilities[,"Ischemic stroke"]<-runif(n.samples,min=1.5*(-0.59),max=0.5*(-0.59))
	event.utilities[,"Transient ischemic attack (TIA)"]<-runif(n.samples,min=1.5*(-0.131),max=1.5*(-0.131))
	event.utilities[,"Clinically relevant bleeding"]<--rnorm(n.samples,0.03,0.001531)
	event.utilities[,"SE"]<-runif(n.samples,min=1.5*(-0.131),max=1.5*(-0.131))
	event.utilities[,"ICH"]<-rnorm(n.samples,0.6,0.064)-AF.utility

	# Finally, scale the state utilities to be 3-monthly
	state.utilities<-state.utilities/4

	
	return(list("event.costs"=event.costs,"event.utilities"=event.utilities,
			"event.switch.probs"=event.switch.probs,"hr.event.history"=hr.event.history,
			"bugs.loghr"=bugs.loghr[1:n.samples,],"bugs.baseline"=bugs.baseline[1:n.samples,],
			"bugs.loghr.apixref"=bugs.loghr.apixref[1:n.samples,],"bugs.baseline.apixref"=bugs.baseline.apixref[1:n.samples,],
			"hr.no.treatment"=hr.no.treatment[1:n.samples,],
			"age.utility.factor"=age.utility.factor,
			"state.costs"=state.costs,"state.utilities"=state.utilities,
			"treatment.costs"=treatment.costs))
}

# Function to generate the transition matrix, costs, and utilities
# Transform log hazards from ith sample into probability.matrix, a matrix of
# event probabilities for each state
# Transform the probability.matrix into event utilities and event costs for
# each state
# Transform the probability.matrix into a n.states x n.states transition.matrix
generate.probabilities<-function(n.samples,ages=c(70:100),event.costs,event.utilities,event.switch.probs,hr.event.history,
						bugs.loghr,bugs.baseline,hr.no.treatment)
{
	# Will return age specific parameters (using common hazards, costs and utilities but adjusting for mortality)
	transition.matrix<-transient.costs<-transient.utilities<-probability.matrix<-list()
	for(age in ages)
	{
	# Probability of each event that patient can have from each state
	probability.matrix[[age]]<-array(0,dim=c(n.samples,n.states,n.events))
	for(i.treatment in 1:n.treatments)
	{
	for(i.health.state in 1:n.health.states)
	{
		# Want to define
		# probability.matrix[(i.treatment-1)*n.health.states+i.health.state,]
		# Find treatment index for baseline hazards and hazard ratios from BUGS
		k.treat<-which(t.names==treatment.names[i.treatment])

		# Add the baseline hazards
		#hazards<-exp(bugs.object.baseline$sims.array[1:n.samples,1,event.codes])
		hazards<-exp(bugs.baseline[1:n.samples,])
		colnames(hazards)<-event.names[1:(n.events-1)]

		# Use life tables for hazard of death
		# Life tables give probability for male/female death within one year

		# Assume the baseline hazard of death comes from the trials (more representative of AF patients)
		# Adjust by age specific hazard ratios from lifetables
		hazards[,"Death (all causes)"]<-hazards[,"Death (all causes)"]*hr.death.age[age-initial.age+1]

		# If not on Warfarin (coumarin) or no treatment, apply hazard ratios treatment effects
		#if(treatment.names[i.treatment]!="Coumarin (INR 2-3)" & treatment.names[i.treatment]!="No treatment")hazards<-hazards*exp(bugs.object.fixed$sims.array[1:n.samples,1,(k.treat-2)*17+event.codes])
		if(treatment.names[i.treatment]!="Coumarin (INR 2-3)" & treatment.names[i.treatment]!="No treatment")hazards<-hazards*exp(bugs.loghr[1:n.samples,(i.treatment-1)*length(event.codes)+c(1:length(event.codes))])


		# If on no treatment, apply hazard ratios from Hart et al. 2007
		if(treatment.names[i.treatment]=="No treatment"){hazards<-hazards*hr.no.treatment}

		# Apply event history hazard ratios
		if(gregexpr(" B ",nondeath.health.states[i.health.state])[[1]][1]!=-1){
			hazards<-hazards*hr.event.history[,event.state.codes==" B ",]}	
		if(gregexpr(" I ",nondeath.health.states[i.health.state])[[1]][1]!=-1){
			hazards<-hazards*hr.event.history[,event.state.codes==" I ",]}
		if(gregexpr(" M ",nondeath.health.states[i.health.state])[[1]][1]!=-1){
			hazards<-hazards*hr.event.history[,event.state.codes==" M ",]}
		if(gregexpr(" S ",nondeath.health.states[i.health.state])[[1]][1]!=-1){
			hazards<-hazards*hr.event.history[,event.state.codes==" S ",]}

		# Competing risks hazards
		probs<-(1-exp(-0.25*rowSums(hazards)))*hazards/rowSums(hazards)
		# Add the probs plus the probability of staying in the same state (no event)
		probability.matrix[[age]][,(i.treatment-1)*n.health.states+i.health.state,]<-t(rbind(t(probs),1-rowSums(probs)))
	}
	}
	# Calculate the event utilities and costs.
	# This includes transient events such as TIA and SE
	transient.costs[[age]]<-transient.utilities[[age]]<-matrix(NA,nrow=n.samples,ncol=n.states)
	for(i.sample in 1:n.samples)
	{
		transient.utilities[[age]][i.sample,]<-probability.matrix[[age]][i.sample,,]%*%event.utilities[i.sample,]
		transient.costs[[age]][i.sample,]<-probability.matrix[[age]][i.sample,,]%*%event.costs[i.sample,]
	}
	
	# Change the probability.matrix into transition.matrix
	# Need rules to assign event probabilities to state probabilities
	# These are conditional on i.treatment and i.health.state
	transition.matrix[[age]]<-array(0,dim=c(n.samples,n.states,n.states))
	for(i.treatment in 1:n.treatments)
	{
	for(i.health.state in 1:n.health.states)
	{
		# The old state is this
		old.state.name<-paste(treatment.names[i.treatment],nondeath.health.states[i.health.state])
		i.state<-which(state.names==old.state.name)
		# Need the name and index of the state patients in this state go
		# to if they have each event.state.codes

		# Find the name and index (in the state.names vector) of the new state (a bit messy)
		new.state.indices<-new.state.name<-rep(NA,n.events)

		# State following Clinically relevant bleeding (B)
		new.state.name[event.state.codes==" B "]<-next.state.name(old.state.name," B ")
		new.state.indices[event.state.codes==" B "]<-which(next.state.name(old.state.name," B ")==state.names)
		# If no discontinuation/switching
		transition.matrix[[age]][,i.state,new.state.indices[event.state.codes==" B "]]<-transition.matrix[[age]][,i.state,new.state.indices[event.state.codes==" B "]]+probability.matrix[[age]][,i.state,event.state.codes==" B "]*(1-event.switch.probs[,event.state.codes==" B "])
		# If discontinuation/switching
		transition.matrix[[age]][,i.state,new.state.indices[event.state.codes==" B "]+treatment.switch.indices[i.treatment]]<-transition.matrix[[age]][,i.state,new.state.indices[event.state.codes==" B "]+treatment.switch.indices[i.treatment]]+probability.matrix[[age]][,i.state,event.state.codes==" B "]*event.switch.probs[,event.state.codes==" B "]
		
		# State following ICH (I)
		new.state.name[event.state.codes==" I "]<-next.state.name(old.state.name," I ")
		new.state.indices[event.state.codes==" I "]<-which(next.state.name(old.state.name," I ")==state.names)
		# If no discontinuation/switching
		transition.matrix[[age]][,i.state,new.state.indices[event.state.codes==" I "]]<-transition.matrix[[age]][,i.state,new.state.indices[event.state.codes==" I "]]+probability.matrix[[age]][,i.state,event.state.codes==" I "]*(1-event.switch.probs[,event.state.codes==" I "])
		# If discontinuation/switching
		transition.matrix[[age]][,i.state,new.state.indices[event.state.codes==" I "]+ich.treatment.switch.indices[i.treatment]]<-transition.matrix[[age]][,i.state,new.state.indices[event.state.codes==" I "]+ich.treatment.switch.indices[i.treatment]]+probability.matrix[[age]][,i.state,event.state.codes==" I "]*event.switch.probs[,event.state.codes==" I "]

		# State following MI (M)
		new.state.name[event.state.codes==" M "]<-next.state.name(old.state.name," M ")
		new.state.indices[event.state.codes==" M "]<-which(next.state.name(old.state.name," M ")==state.names)
		# If no discontinuation/switching
		transition.matrix[[age]][,i.state,new.state.indices[event.state.codes==" M "]]<-transition.matrix[[age]][,i.state,new.state.indices[event.state.codes==" M "]]+probability.matrix[[age]][,i.state,event.state.codes==" M "]*(1-event.switch.probs[,event.state.codes==" M "])
		# If discontinuation/switching
		transition.matrix[[age]][,i.state,new.state.indices[event.state.codes==" M "]+mi.treatment.switch.indices[i.treatment]]<-transition.matrix[[age]][,i.state,new.state.indices[event.state.codes==" M "]+mi.treatment.switch.indices[i.treatment]]+probability.matrix[[age]][,i.state,event.state.codes==" M "]*event.switch.probs[,event.state.codes==" M "]
		
		# State following Ischemic stroke (S)
		new.state.name[event.state.codes==" S "]<-next.state.name(old.state.name," S ")
		new.state.indices[event.state.codes==" S "]<-which(next.state.name(old.state.name," S ")==state.names)
		# If no discontinuation/switching
		transition.matrix[[age]][,i.state,new.state.indices[event.state.codes==" S "]]<-transition.matrix[[age]][,i.state,new.state.indices[event.state.codes==" S "]]+probability.matrix[[age]][,i.state,event.state.codes==" S "]*(1-event.switch.probs[,event.state.codes==" S "])
		# If discontinuation/switching
		transition.matrix[[age]][,i.state,new.state.indices[event.state.codes==" S "]+treatment.switch.indices[i.treatment]]<-transition.matrix[[age]][,i.state,new.state.indices[event.state.codes==" S "]+treatment.switch.indices[i.treatment]]+probability.matrix[[age]][,i.state,event.state.codes==" S "]*event.switch.probs[,event.state.codes==" S "]


		#transition.matrix[i.state,new.state.indices]<-probability.matrix[i.state,event.state.codes!="" & event.names!="Death (all causes)"]
		# Probability stay (always sum of "Stay" and transient states
		transition.matrix[[age]][,i.state,i.state]<-transition.matrix[[age]][,i.state,i.state]+ rowSums(probability.matrix[[age]][,i.state,event.state.codes=="" & event.names!="Death (all causes)"]*(1-event.switch.probs[,event.state.codes==""]))                  
		# If discontinuation/switching (sum of transient event switching probabilities and no event switching probability)
		transition.matrix[[age]][,i.state,i.state+treatment.switch.indices[i.treatment]]<-transition.matrix[[age]][,i.state,i.state+treatment.switch.indices[i.treatment]]+rowSums(probability.matrix[[age]][,i.state,event.state.codes=="" & event.names!="Death (all causes)"]*(event.switch.probs[,event.state.codes==""]))
		
		# Probability death (always just the probability of death)
		transition.matrix[[age]][,i.state,n.states]<-probability.matrix[[age]][,i.state,event.names=="Death (all causes)"]
	}
	}

	# Dead patients stay dead with probability 1
	transition.matrix[[age]][,n.states,n.states]<-1
	} # End loop over ages
	return(list("transition.matrix"=transition.matrix,"transient.utilities"=transient.utilities,"transient.costs"=transient.costs,"probability.matrix"=probability.matrix[[1]][1:n.samples,,]))
}

