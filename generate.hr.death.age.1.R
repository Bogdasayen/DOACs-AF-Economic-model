# This research was funded by the National Institute for Health Research (NIHR), NIHR Senior Investigator award NF-SI-0611-10168.
# Support for R in CEA provided by MRC Hubs for Trials Methodology Research ConDuCT-II (Collaboration and innovation in Difficult and Complex randomised controlled Trials In Invasive procedures) hub
# Howard Thom 31-October-2018. Bristol Medical School: Population Health Sciences. Bristol University, UK. howard.thom@bristol.ac.uk
# Details of model described in two publications:
# Sterne JA, Bodalia PN, Bryden PA, Davies PA, Lopez-Lopez JA, Okoli GN, et al. Oral anticoagulants for primary prevention, treatment and secondary prevention of venous thromboembolic disease, and for prevention of stroke in atrial fibrillation: systematic review, network meta-analysis and cost-effectiveness analysis. Health Technol Assess. 2017; 21(9):1-386.
# Lopez-Lopez JA, Sterne J, Thom H, Higgins J, Hingorani A, Okoli GN, et al. Oral anticoagulants for prevention of stroke in atrial fibrillation: systematic review, network meta-analysis, and cost effectiveness analysis. BMJ. 2017; 359(J5058).


# Function to convert life-tables into hazard ratios for age


generate.hr.death.age<-function(base.age,lifetables)
{
		# Use life tables for hazard of death
		# Life tables give probability for male/female death within one year
		# Assume a 60/40 split male/female
		# Calculate a hazard from this using Prob=1-exp(-hazard)
		# So annual hazard=-log(1-Prob)
		hazard.death<-rep(NA,max(lifetables[,"Age"])-base.age+1)
		for(i in 1:(max(lifetables[,"Age"])-base.age+1)){
			hazard.death[i]<--log(1-((0.6*lifetables[lifetables[,"Age"]==(base.age+i-1),"male.qx"]+0.4*lifetables[lifetables[,"Age"]==(base.age+i-1),"female.qx"])))
		}
		hr.death<-hazard.death/hazard.death[1]
		names(hr.death)<-c(base.age:(max(lifetables[,"Age"])))
		return(hr.death)
}