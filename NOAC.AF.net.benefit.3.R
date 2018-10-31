# This research was funded by the National Institute for Health Research (NIHR), NIHR Senior Investigator award NF-SI-0611-10168.
# Support for R in CEA provided by MRC Hubs for Trials Methodology Research ConDuCT-II (Collaboration and innovation in Difficult and Complex randomised controlled Trials In Invasive procedures) hub
# Howard Thom 31-October-2018. Bristol Medical School: Population Health Sciences. Bristol University, UK. howard.thom@bristol.ac.uk
# Details of model described in two publications:
# Sterne JA, Bodalia PN, Bryden PA, Davies PA, Lopez-Lopez JA, Okoli GN, et al. Oral anticoagulants for primary prevention, treatment and secondary prevention of venous thromboembolic disease, and for prevention of stroke in atrial fibrillation: systematic review, network meta-analysis and cost-effectiveness analysis. Health Technol Assess. 2017; 21(9):1-386.
# Lopez-Lopez JA, Sterne J, Thom H, Higgins J, Hingorani A, Okoli GN, et al. Oral anticoagulants for prevention of stroke in atrial fibrillation: systematic review, network meta-analysis, and cost effectiveness analysis. BMJ. 2017; 359(J5058).


# Function to generate the costs, qalys, and net benefits of the NOACs model
# V3 (16-Oct-2017) includes option to use apixaban as the reference (apixref=TRUE)


noac.net.benefit<-function(n.samples,n.cycles,initial.age,lambdas,age.independent.samples,
			apixref=FALSE,prob.treatment.samples=FALSE)
{
	# Vectors of total costs and total qalys output from the model
	total.qalys<-total.costs<-matrix(0,nrow=n.samples,ncol=(n.treatments-1))
	colnames(total.costs)<-colnames(total.qalys)<-treatment.names[1:(n.treatments-1)]

	prob.on.treatment<-matrix(NA,ncol=n.treatments-1,nrow=n.cycles)
	if(prob.treatment.samples)prob.on.treatment.samples<-array(NA,dim=c(n.samples,n.cycles,n.treatments-1),dimnames=list(NULL,NULL,treatment.names[1:(n.treatments-1)]))
	colnames(prob.on.treatment)<-treatment.names[1:(n.treatments-1)]


	# Set up cohort vector
	cohort.vector<-array(0,dim=c(n.samples,(n.treatments-1),n.states))
	# Initialise the cohort vectors for each treatment
	for(i.treatment in 1:(n.treatments-1))cohort.vector[,i.treatment,(n.health.states)*(i.treatment-1)+1]<-1

	# Run model for n.cycles (eg. n.cycles=120 means years=30)
	age<-initial.age
	old.age<-age-1

	for(i.cycle in 1:n.cycles)
	{

		# Calculate mean time on each treatment(for price threshold analysis)
		for(i.treatment in 1:(n.treatments-1))
		{
			prob.on.treatment[i.cycle,i.treatment]<-sum(colMeans(cohort.vector[,i.treatment,(16*(i.treatment-1)+1):(16*(i.treatment-1)+16)]))
			if(prob.treatment.samples)prob.on.treatment.samples[,i.cycle,i.treatment]<-rowSums(cohort.vector[,i.treatment,(16*(i.treatment-1)+1):(16*(i.treatment-1)+16)])
		}

		# Only resample every four cycles
		if((floor(age)-old.age)==1)
		{
			if(apixref==FALSE){
			sampled.probabilities<-generate.probabilities(n.samples,ages=floor(age),event.costs=age.independent.samples$event.costs,
							event.utilities=age.independent.samples$event.utilities,event.switch.probs=age.independent.samples$event.switch.probs,
							hr.event.history=age.independent.samples$hr.event.history,bugs.loghr=age.independent.samples$bugs.loghr,
							bugs.baseline=age.independent.samples$bugs.baseline,hr.no.treatment=age.independent.samples$hr.no.treatment)
			# Some code used to make sure the reparameterisation has been done correctly
			# sampled.probabilities.standard<-sampled.probabilities
			# sampled.probabilities.standard$transition.matrix[[70]][100,1,]
			}
			if(apixref==TRUE){
			sampled.probabilities<-generate.probabilities.apixref(n.samples,ages=floor(age),event.costs=age.independent.samples$event.costs,
							event.utilities=age.independent.samples$event.utilities,event.switch.probs=age.independent.samples$event.switch.probs,
							hr.event.history=age.independent.samples$hr.event.history,bugs.baseline.apixref=age.independent.samples$bugs.baseline.apixref,
							bugs.loghr.apixref=age.independent.samples$bugs.loghr.apixref,hr.no.treatment=age.independent.samples$hr.no.treatment)
			# Some code used to make sure the reparameterisation has been done correctly
			# sampled.probabilities.apixref<-sampled.probabilities
			# sampled.probabilities.apixref$transition.matrix[[70]][100,1,]
			
			}
			old.age<-age
		}

		# Half cycle correction (subtract half of initial QALYs and costs)
		if(i.cycle==1)
		{
			for(i.sample in 1:n.samples)
			{
				total.costs[i.sample,]<--0.5*((1/1.035)^(1/4))*cohort.vector[i.sample,,]%*%(age.independent.samples$state.costs+sampled.probabilities$transient.costs[[floor(age)]])[i.sample,]
				total.qalys[i.sample,]<--age.independent.samples$age.utility.factor[1,toString(round(age/10)*10-5)]*0.5*((1/1.035)^(1/4))*cohort.vector[i.sample,,]%*%(age.independent.samples$state.utilities+sampled.probabilities$transient.utilities[[floor(age)]])[i.sample,]
			}
		}

		for(i.sample in 1:n.samples)
		{
			# Add costs and QALYs for current cohort vector
			total.costs[i.sample,]<-total.costs[i.sample,]+((1/1.035)^(i.cycle/4))*cohort.vector[i.sample,,]%*%(age.independent.samples$state.costs+sampled.probabilities$transient.costs[[floor(age)]])[i.sample,]
			total.qalys[i.sample,]<-total.qalys[i.sample,]+age.independent.samples$age.utility.factor[1,toString(round(age/10)*10-5)]*((1/1.035)^(i.cycle/4))*cohort.vector[i.sample,,]%*%(age.independent.samples$state.utilities+sampled.probabilities$transient.utilities[[floor(age)]])[i.sample,]
			# Matrix multiplication to update cohort vector
			cohort.vector[i.sample,,]<-cohort.vector[i.sample,,]%*%sampled.probabilities$transition.matrix[[floor(age)]][i.sample,,]		
		} # End sample loop

		age<-age+0.25

	} # End cycle loop

	# Half cycle correction (add half of final QALYs and costs)
	age<-age+0.25
	# Only resample every four cycles
	if((floor(age)-old.age)==1)
	{
		if(apixref==FALSE){
		sampled.probabilities<-generate.probabilities(n.samples,ages=floor(age),event.costs=age.independent.samples$event.costs,
						event.utilities=age.independent.samples$event.utilities,event.switch.probs=age.independent.samples$event.switch.probs,
						hr.event.history=age.independent.samples$hr.event.history,bugs.loghr=age.independent.samples$bugs.loghr,
						bugs.baseline=age.independent.samples$bugs.baseline,hr.no.treatment=age.independent.samples$hr.no.treatment)
		}
		if(apixref==TRUE){
		sampled.probabilities<-generate.probabilities.apixref(n.samples,ages=floor(age),event.costs=age.independent.samples$event.costs,
							event.utilities=age.independent.samples$event.utilities,event.switch.probs=age.independent.samples$event.switch.probs,
							hr.event.history=age.independent.samples$hr.event.history,bugs.baseline.apixref=age.independent.samples$bugs.baseline.apixref,
							bugs.loghr.apixref=age.independent.samples$bugs.loghr.apixref,hr.no.treatment=age.independent.samples$hr.no.treatment)
		}
		old.age<-age
	}
	for(i.sample in 1:n.samples)
	{
		# Sum the costs and QALYs over all n.cycles, discounting at 3.5% per year.
		total.costs[i.sample,]<-total.costs[i.sample,]+0.5*((1/1.035)^((n.cycles+1)/4))*cohort.vector[i.sample,,]%*%(age.independent.samples$state.costs+sampled.probabilities$transient.costs[[floor(age)]])[i.sample,]
		total.qalys[i.sample,]<-total.qalys[i.sample,]+age.independent.samples$age.utility.factor[1,toString(round(age/10)*10-5)]*0.5*((1/1.035)^((n.cycles+1)/4))*cohort.vector[i.sample,,]%*%(age.independent.samples$state.utilities+sampled.probabilities$transient.utilities[[floor(age)]])[i.sample,]
	}

	# Average time spent on treatment (if starting on it) used for threshold analysis
	mean.years.treatment<-colSums(prob.on.treatment)/4
	if(prob.treatment.samples)mean.years.treatment.samples<-apply(prob.on.treatment.samples,c(1,3),sum)/4

	# Need to discount at this stage so that later years contribute less
	prob.on.treatment.discounted<-prob.on.treatment*rep((1/1.035)^c(0:((n.cycles/4)-1)),each=4)
	mean.years.treatment.discounted<-colSums(prob.on.treatment.discounted)/4

	# Also need the discounted prob and years for each sample
	if(prob.treatment.samples){
		prob.on.treatment.samples.discounted<-prob.on.treatment.samples	
		for(i.treatment in 1:dim(prob.on.treatment.samples)[3])
		{
			prob.on.treatment.samples.discounted[,,i.treatment]<-t(t(prob.on.treatment.samples[,,i.treatment])*
					rep((1/1.035)^c(0:((n.cycles/4)-1)),each=4))
		}
		mean.years.treatment.samples.discounted<-apply(prob.on.treatment.samples.discounted,c(1,3),sum)/4
	}


	# Generate the net benefit and incremental net benefit for output
	INB<-array(NA, dim=c(n.samples,n.treatments-1,length(lambdas)))
	NB<-array(NA, dim=c(n.samples,n.treatments,length(lambdas)))
	# Calculate the incremental net benefit for each treatment
	for(i.treatment in 2:(n.treatments-1))
	{
	for(i.lambda in 1:length(lambdas))
	{
		INB[,i.treatment-1,i.lambda]<-lambdas[i.lambda]*(total.qalys[,i.treatment]-total.qalys[,1])-(total.costs[,i.treatment]-total.costs[,1])
	}
	}
	for(i.treatment in 1:(n.treatments-1))for(i.lambda in 1:length(lambdas)){NB[,i.treatment,i.lambda]<-lambdas[i.lambda]*total.qalys[,i.treatment]-total.costs[,i.treatment]}

	if(!prob.treatment.samples)mean.years.treatment.samples<-mean.years.treatment.samples.discounted<-NULL

	return(list("total.costs"=total.costs,"total.qalys"=total.qalys,"NB"=NB,"INB"=INB,
			"prob.on.treatment"=prob.on.treatment,
			"mean.years.treatment"=mean.years.treatment,
			"prob.on.treatment.discounted"=prob.on.treatment.discounted,
			"mean.years.treatment.discounted"=mean.years.treatment.discounted,
			"mean.years.treatment.samples"=mean.years.treatment.samples,
			"mean.years.treatment.samples.discounted"=mean.years.treatment.samples.discounted))

} # End net benefit function

