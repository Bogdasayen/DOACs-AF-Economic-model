# This research was funded by the National Institute for Health Research (NIHR) Health Technology Assessment programme project number 11/92/17 and NIHR Senior Investigator award NF-SI-0611-10168.
# Support for R in CEA provided by MRC Hubs for Trials Methodology Research ConDuCT-II (Collaboration and innovation in Difficult and Complex randomised controlled Trials In Invasive procedures) hub
# Howard Thom 31-October-2018. Bristol Medical School: Population Health Sciences. Bristol University, UK. howard.thom@bristol.ac.uk
# Details of model described in two publications:
# Sterne JA, Bodalia PN, Bryden PA, Davies PA, Lopez-Lopez JA, Okoli GN, et al. Oral anticoagulants for primary prevention, treatment and secondary prevention of venous thromboembolic disease, and for prevention of stroke in atrial fibrillation: systematic review, network meta-analysis and cost-effectiveness analysis. Health Technol Assess. 2017; 21(9):1-386.
# Lopez-Lopez JA, Sterne J, Thom H, Higgins J, Hingorani A, Okoli GN, et al. Oral anticoagulants for prevention of stroke in atrial fibrillation: systematic review, network meta-analysis, and cost effectiveness analysis. BMJ. 2017; 359(J5058).


# Modified version of generate.probabilities so that apixaban is the reference

bugs.loghr.apixref<-read.csv(file=paste(baseline.directory,"/data/bugs.loghr.apixref.csv",sep=""))
bugs.baseline.apixref<-read.csv(file=paste(baseline.directory,"/data/bugs.baseline.apixref.csv",sep=""))
# Remove empty first columns
bugs.loghr.apixref<-bugs.loghr.apixref[,-1]
bugs.baseline.apixref<-bugs.baseline.apixref[,-1]


# Function to generate the transition matrix, costs, and utilities
# Transform log hazards from ith sample into probability.matrix, a matrix of
# event probabilities for each state
# Transform the probability.matrix into event utilities and event costs for
# each state
# Transform the probability.matrix into a n.states x n.states transition.matrix
generate.probabilities.apixref<-function(n.samples,ages=c(70:100),event.costs,event.utilities,event.switch.probs,hr.event.history,
					bugs.loghr.apixref,bugs.baseline.apixref,hr.no.treatment)
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
		hazards<-exp(bugs.baseline.apixref[1:n.samples,])
		colnames(hazards)<-event.names[1:(n.events-1)]

		# Add the baseline hazards for warfarin; needed for no treatment as it was a meta-analysis vs warfarin, not apixaban
		hazards.warfarin<-exp(bugs.baseline.apixref[1:n.samples,])
		colnames(hazards.warfarin)<-event.names[1:(n.events-1)]
		hazards.warfarin<-hazards.warfarin*exp(bugs.loghr.apixref[1:n.samples,c(1:length(event.codes))])


		# Use life tables for hazard of death
		# Life tables give probability for male/female death within one year

		# Assume the baseline hazard of death comes from the trials (more representative of AF patients)
		# Adjust by age specific hazard ratios from lifetables
		hazards[,"Death (all causes)"]<-hazards[,"Death (all causes)"]*hr.death.age[age-initial.age+1]
		# Also need to adjust warfarin hazards (by same amount)
		hazards.warfarin[,"Death (all causes)"]<-hazards.warfarin[,"Death (all causes)"]*hr.death.age[age-initial.age+1]

		# If not on Warfarin (coumarin) or no treatment, apply hazard ratios treatment effects
		#if(treatment.names[i.treatment]!="Coumarin (INR 2-3)" & treatment.names[i.treatment]!="No treatment")hazards<-hazards*exp(bugs.object.fixed$sims.array[1:n.samples,1,(k.treat-2)*17+event.codes])
		#if(treatment.names[i.treatment]!="Coumarin (INR 2-3)" & treatment.names[i.treatment]!="No treatment")hazards<-hazards*exp(bugs.loghr[1:n.samples,(i.treatment-1)*length(event.codes)+c(1:length(event.codes))])
		if(treatment.names[i.treatment]!="Apixaban (5mg bd)" & treatment.names[i.treatment]!="No treatment")hazards<-hazards*exp(bugs.loghr.apixref[1:n.samples,(i.treatment-1)*length(event.codes)+c(1:length(event.codes))])


		# If on no treatment, apply hazard ratios from Hart et al. 2007
		if(treatment.names[i.treatment]=="No treatment"){hazards<-hazards.warfarin*hr.no.treatment}

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
		# If no discontinuation/switching
		transition.matrix[[age]][,i.state,i.state]<-transition.matrix[[age]][,i.state,i.state]+rowSums(probability.matrix[[age]][,i.state,event.state.codes=="" & event.names!="Death (all causes)"])*(1-rowSums(event.switch.probs[,event.state.codes==""]))
		# If discontinuation/switching (sum of transient event switching probabilities and no event switching probability)
		transition.matrix[[age]][,i.state,i.state+treatment.switch.indices[i.treatment]]<-transition.matrix[[age]][,i.state,i.state+treatment.switch.indices[i.treatment]]+rowSums(probability.matrix[[age]][,i.state,event.state.codes=="" & event.names!="Death (all causes)"])*rowSums(event.switch.probs[,event.state.codes==""])

		# Probability death (always just the probability of death)
		transition.matrix[[age]][,i.state,n.states]<-probability.matrix[[age]][,i.state,event.names=="Death (all causes)"]
	}
	}

	# Dead patients stay dead with probability 1
	transition.matrix[[age]][,n.states,n.states]<-1
	} # End loop over ages
	return(list("transition.matrix"=transition.matrix,"transient.utilities"=transient.utilities,"transient.costs"=transient.costs,"probability.matrix"=probability.matrix[[1]][1:n.samples,,]))
}

