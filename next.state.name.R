# This research was funded by the National Institute for Health Research (NIHR), NIHR Senior Investigator award NF-SI-0611-10168.
# Support for R in CEA provided by MRC Hubs for Trials Methodology Research ConDuCT-II (Collaboration and innovation in Difficult and Complex randomised controlled Trials In Invasive procedures) hub
# Howard Thom 31-October-2018. Bristol Medical School: Population Health Sciences. Bristol University, UK. howard.thom@bristol.ac.uk
# Details of model described in two publications:
# Sterne JA, Bodalia PN, Bryden PA, Davies PA, Lopez-Lopez JA, Okoli GN, et al. Oral anticoagulants for primary prevention, treatment and secondary prevention of venous thromboembolic disease, and for prevention of stroke in atrial fibrillation: systematic review, network meta-analysis and cost-effectiveness analysis. Health Technol Assess. 2017; 21(9):1-386.
# Lopez-Lopez JA, Sterne J, Thom H, Higgins J, Hingorani A, Okoli GN, et al. Oral anticoagulants for prevention of stroke in atrial fibrillation: systematic review, network meta-analysis, and cost effectiveness analysis. BMJ. 2017; 359(J5058).


# Function to work out which state follows which event

# Laborious version
next.state.name<-function(old.state.name,event.code)
{
	# Extract the nondeath.health.state
	i.health.state<-n.health.states
	while(substr(old.state.name,nchar(old.state.name)-nchar(nondeath.health.states[i.health.state])+1,nchar(old.state.name))!=nondeath.health.states[i.health.state] & i.health.state>0){i.health.state<-i.health.state-1}
	# Extract the treatment name
	treatment.name<-substr(old.state.name,0,nchar(old.state.name)-nchar(nondeath.health.states[i.health.state]))
	# Extract the old health state
	health.state<-nondeath.health.states[i.health.state]
	if(event.code==" B ")
	{
		if(health.state=="Well"){
			return(paste(treatment.name," B ",sep=""))
		}
		if(health.state==" B "){
			return(paste(treatment.name," B ",sep=""))
		}
		if(health.state==" I "){
			return(paste(treatment.name," B + I ",sep=""))
		}
		if(health.state==" M "){
			return(paste(treatment.name," B + M ",sep=""))
		}
		if(health.state==" S "){
			return(paste(treatment.name," B + S ",sep=""))
		}
		if(health.state==" B + I "){
			return(paste(treatment.name," B + I ",sep=""))
		}
		if(health.state==" B + M "){
			return(paste(treatment.name," B + M ",sep=""))
		}
		if(health.state==" B + S "){
			return(paste(treatment.name," B + S ",sep=""))
		}
		if(health.state==" I + M "){
			return(paste(treatment.name," B + I + M ",sep=""))
		}
		if(health.state==" I + S "){
			return(paste(treatment.name," B + I + S ",sep=""))
		}
		if(health.state==" M + S "){
			return(paste(treatment.name," B + M + S ",sep=""))
		}
		if(health.state==" B + I + M "){
			return(paste(treatment.name, " B + I + M ",sep=""))
		}
		if(health.state==" B + I + S "){
			return(paste(treatment.name, " B + I + S ",sep=""))
		}
		if(health.state==" B + M + S "){			
			return(paste(treatment.name," B + M + S ",sep=""))
		}
		if(health.state==" I + M + S "){
			return(paste(treatment.name," B + I + M + S ",sep=""))
		}
		if(health.state==" B + I + M + S "){
			return(paste(treatment.name,health.state,sep=""))
		}
	}
	if(event.code==" I ")
	{
		if(health.state=="Well"){
			return(paste(treatment.name," I ",sep=""))
		}
		if(health.state==" B "){
			return(paste(treatment.name," B + I ",sep=""))
		}
		if(health.state==" I "){
			return(paste(treatment.name," I ",sep=""))
		}
		if(health.state==" M "){
			return(paste(treatment.name," I + M ",sep=""))
		}
		if(health.state==" S "){
			return(paste(treatment.name," I + S ",sep=""))
		}
		if(health.state==" B + I "){
			return(paste(treatment.name," B + I ",sep=""))
		}
		if(health.state==" B + M "){
			return(paste(treatment.name," B + I + M ",sep=""))
		}
		if(health.state==" B + S "){
			return(paste(treatment.name," B + I + S ",sep=""))
		}
		if(health.state==" I + M "){
			return(paste(treatment.name," I + M ",sep=""))
		}
		if(health.state==" I + S "){
			return(paste(treatment.name," I + S ",sep=""))
		}
		if(health.state==" M + S "){
			return(paste(treatment.name," I + M + S ",sep=""))
		}
		if(health.state==" B + I + M "){
			return(paste(treatment.name, " B + I + M ",sep=""))
		}
		if(health.state==" B + I + S "){
			return(paste(treatment.name, " B + I + S ",sep=""))
		}
		if(health.state==" B + M + S "){			
			return(paste(treatment.name," B + I + M + S ",sep=""))
		}
		if(health.state==" I + M + S "){
			return(paste(treatment.name," I + M + S ",sep=""))
		}
		if(health.state==" B + I + M + S "){
			return(paste(treatment.name,health.state,sep=""))
		}
	}
	if(event.code==" M ")
	{
		if(health.state=="Well"){
			return(paste(treatment.name," M ",sep=""))
		}
		if(health.state==" B "){
			return(paste(treatment.name," B + M ",sep=""))
		}
		if(health.state==" I "){
			return(paste(treatment.name," I + M ",sep=""))
		}
		if(health.state==" M "){
			return(paste(treatment.name," M ",sep=""))
		}
		if(health.state==" S "){
			return(paste(treatment.name," M + S ",sep=""))
		}
		if(health.state==" B + I "){
			return(paste(treatment.name," B + I + M ",sep=""))
		}
		if(health.state==" B + M "){
			return(paste(treatment.name," B + M ",sep=""))
		}
		if(health.state==" B + S "){
			return(paste(treatment.name," B + M + S ",sep=""))
		}
		if(health.state==" I + M "){
			return(paste(treatment.name," I + M ",sep=""))
		}
		if(health.state==" I + S "){
			return(paste(treatment.name," I + M + S ",sep=""))
		}
		if(health.state==" M + S "){
			return(paste(treatment.name," M + S ",sep=""))
		}
		if(health.state==" B + I + M "){
			return(paste(treatment.name, " B + I + M ",sep=""))
		}
		if(health.state==" B + I + S "){
			return(paste(treatment.name, " B + I + M + S ",sep=""))
		}
		if(health.state==" B + M + S "){			
			return(paste(treatment.name," B + M + S ",sep=""))
		}
		if(health.state==" I + M + S "){
			return(paste(treatment.name," I + M + S ",sep=""))
		}
		if(health.state==" B + I + M + S "){
			return(paste(treatment.name,health.state,sep=""))
		}
	}
	if(event.code==" S ")
	{
		if(health.state=="Well"){
			return(paste(treatment.name," S ",sep=""))
		}
		if(health.state==" B "){
			return(paste(treatment.name," B + S ",sep=""))
		}
		if(health.state==" I "){
			return(paste(treatment.name," I + S ",sep=""))
		}
		if(health.state==" M "){
			return(paste(treatment.name," M + S ",sep=""))
		}
		if(health.state==" S "){
			return(paste(treatment.name," S ",sep=""))
		}
		if(health.state==" B + I "){
			return(paste(treatment.name," B + I + S ",sep=""))
		}
		if(health.state==" B + M "){
			return(paste(treatment.name," B + M + S ",sep=""))
		}
		if(health.state==" B + S "){
			return(paste(treatment.name," B + S ",sep=""))
		}
		if(health.state==" I + M "){
			return(paste(treatment.name," I + M + S ",sep=""))
		}
		if(health.state==" I + S "){
			return(paste(treatment.name," I + S ",sep=""))
		}
		if(health.state==" M + S "){
			return(paste(treatment.name," M + S ",sep=""))
		}
		if(health.state==" B + I + M "){
			return(paste(treatment.name, " B + I + M + S ",sep=""))
		}
		if(health.state==" B + I + S "){
			return(paste(treatment.name, " B + I + S ",sep=""))
		}
		if(health.state==" B + M + S "){			
			return(paste(treatment.name," B + M + S ",sep=""))
		}
		if(health.state==" I + M + S "){
			return(paste(treatment.name," I + M + S ",sep=""))
		}
		if(health.state==" B + I + M + S "){
			return(paste(treatment.name,health.state,sep=""))
		}
	}
}


# Failed attempt at smart version
next.state.name.alt<-function(old.state.name,event.code)
{
	# If history already includes event, new state is old state
	if(gregexpr(event.code,old.state.name)[[1]][1]!=-1)
	{
		new.state.name<-old.state.name
	}else{
	# History does not include the event so it must be alphabetically inserted
	# Find alphabetically highest event code
	#i.max.code<-0
	#while(gregexpr(event.codes.alpha[i.max.code+1],old.state.name)[[1]][1]!=-1){i.max.code<-i.max.code+1}
	# Find the position (or -1) of each event
	event.code.index<-rep(-1, length(event.codes.alpha))
	for(i in 1:length(event.code.index))event.code.index[i]<-gregexpr(event.codes.alpha[i],old.state.name)[[1]][1]
	i.max.code<-which.max(event.code.index)
	# Find the location in the state name string at which to place new code
	# If any events have occured
	if(i.max.code>0){
		# Event is alphabetically higher than previous max (eg. "M">"B")
		if(i.max.code<which(event.codes.alpha==event.code))
		{
			char.position<-gregexpr(event.codes.alpha[i.max.code],old.state.name)[[1]][1]+nchar(event.codes.alpha[i.max.code])
			new.state.name<-paste(substr(old.state.name, 1, char.position-1), "+", event.code, substr(old.state.name, char.position, nchar(old.state.name)), sep = "")
		}else{
		# Event is alphabetically lower than previous max (eg. "M"<"S")
			char.position<-gregexpr(event.codes.alpha[i.max.code],old.state.name)[[1]][1]
	###########################################################
			new.state.name<-0
	#############################################################
		}
	}else{
	# If this is the first event put it at the end
		#char position<-nchar(old.state.name)
		new.state.name<-paste(old.state.name, event.code , sep = "")
	}

	}

	return(new.state.name)
}

