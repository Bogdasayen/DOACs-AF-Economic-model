# This research was funded by the National Institute for Health Research (NIHR) Health Technology Assessment programme project number 11/92/17 and NIHR Senior Investigator award NF-SI-0611-10168.
# Support for R in CEA provided by MRC Hubs for Trials Methodology Research ConDuCT-II (Collaboration and innovation in Difficult and Complex randomised controlled Trials In Invasive procedures) hub
# Howard Thom 31-October-2018. Bristol Medical School: Population Health Sciences. Bristol University, UK. howard.thom@bristol.ac.uk
# Details of model described in two publications:
# Sterne JA, Bodalia PN, Bryden PA, Davies PA, Lopez-Lopez JA, Okoli GN, et al. Oral anticoagulants for primary prevention, treatment and secondary prevention of venous thromboembolic disease, and for prevention of stroke in atrial fibrillation: systematic review, network meta-analysis and cost-effectiveness analysis. Health Technol Assess. 2017; 21(9):1-386.
# Lopez-Lopez JA, Sterne J, Thom H, Higgins J, Hingorani A, Okoli GN, et al. Oral anticoagulants for prevention of stroke in atrial fibrillation: systematic review, network meta-analysis, and cost effectiveness analysis. BMJ. 2017; 359(J5058).


# NOAC AF Cost-effectiveness model
# Script with utility functions

gamma.parameters<-function(gamma.mean,gamma.sd)
{
	gamma.shape=(gamma.mean^2)/(gamma.sd^2)
	gamma.scale=(gamma.sd^2)/gamma.mean
	return(list("shape"=gamma.shape,"scale"=gamma.scale))
}

beta.parameters<-function(beta.mean,beta.sd)
{
	beta.var<-beta.sd^2
  alpha <- ((1 - beta.mean) / beta.var - 1 / beta.mean) * beta.mean ^ 2
  beta <- alpha * (1 / beta.mean - 1)
  return(list(alpha = alpha, beta = beta))
}
