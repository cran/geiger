`negloglike` <-
function(theta, ds)
{
	#----------------------------------
	#-----       PARSE DATA       -----
	#----------------------------------
	y			<- ds$data				# TIP data
	tree		<- ds$tree			# Tree
	meserr		<- ds$meserr
	#----------------------------------
	#-----  DETERMINE OTHER INFO  -----
	#----------------------------------
	n			<- length(y)
	bounds		<- ds$bounds
	model		<- ds$model			
	if (sum(names(model)==c("lambda", "kappa", "delta", "ou", "eb"))!=5)
		stop("Error with \"model\" internal paramter")

	#----------------------------------
	#-----       DEFAULT FIT      -----
	#----------------------------------
	if (sum(model)==0) {
		mu		<-	theta[1]
		beta	<-  inv.logit(theta[2], min=bounds$beta[1], max=bounds$beta[2])
		means	<- 	rep(mu,n)
		vcv<-vcv.phylo(tree)
		diag(vcv)=diag(vcv)+meserr^2
		#---- RETURN NEGLOGLIKELIHOOD ---
		return( -dmvnorm(y,means,beta*vcv, log=TRUE))
	#----------------------------------
	#-----       LAMBDA ONLY      -----
	#----------------------------------
	} else if (model[1] & !(model[2] | model[3])){
		mu		<-	theta[1]
		beta	<-  inv.logit(theta[2], min=bounds$beta[1], max=bounds$beta[2])
		lambda	<- 	inv.logit(theta[3], min=bounds$lambda[1], max=bounds$lambda[2])
		vcv<-vcv.phylo(tree)
		#---- multiply off diagonals
		index			<-	matrix(TRUE, n,n)
		diag(index)		<- FALSE
		vcv[index] 	<- vcv[index]*lambda
		means		    <- rep(mu,n)
		diag(vcv)=diag(vcv)+meserr^2

		#---- RETURN NEGLOGLIKELIHOOD ---
		return( -dmvnorm(y,means,beta*vcv, log=TRUE))
	#----------------------------------
	#-----        KAPPA ONLY      -----
	#----------------------------------
	} else if (model[2] & !(model[1] | model[3])){
		beta.bound 	 <- bounds$beta
		kappa.bound  <- bounds$kappa
		mu		<-	theta[1]
		beta	<-  inv.logit(theta[2], min=bounds$beta[1], max=bounds$beta[2])
		kappa	<- 	inv.logit(theta[3], min=bounds$kappa[1], max=bounds$kappa[2])
		if (kappa==0)
			stop("kappa = 0, procedure haulted")
		#---- RAISE BRANCH LENGTHS BY KAPPA
		tree$edge.length<-tree$edge.length^kappa
		#---- DETERMINE VCV
		vcv <- vcv.phylo(tree)   
		means	<- 	rep(mu,n)
		diag(vcv)=diag(vcv)+meserr^2

		#---- RETURN NEGLOGLIKELIHOOD ---
		return( -dmvnorm(y,means,beta*vcv, log=TRUE))
	#----------------------------------
	#-----        DELTA ONLY      -----
	#----------------------------------	
	} else if (model[3] & !(model[1] | model[2])){
		mu		<-	theta[1]
		beta	<-  inv.logit(theta[2], min=bounds$beta[1], max=bounds$beta[2])
		delta	<-  inv.logit(theta[3], min=bounds$delta[1], max=bounds$delta[2])
		means	<- rep(mu,n)
		mb=max(branching.times(tree))
		tree<-deltaTree(tree, delta=delta)
		vcv<-vcv.phylo(tree)

		rescale <- mb/max(vcv)
		vcv	<-(vcv*rescale)
		diag(vcv)=diag(vcv)+meserr^2

		#---- RETURN NEGLOGLIKELIHOOD ---
		return( -dmvnorm(y,means,beta*vcv, log=TRUE))
	#----------------------------------
	#-----        ALPHA ONLY      -----
	#----------------------------------			
	} else if (model[4] & !(model[1] | model[2] | model[3])){
		mu		<-	theta[1]
		beta	<-  inv.logit(theta[2], min=bounds$beta[1], max=bounds$beta[2])
		alpha	<-  inv.logit(theta[3], min=bounds$alpha[1], max=bounds$alpha[2])

		#---- OU TRANSFORMATION
		tree<-ouTree(tree, alpha=alpha)
		vcv	<-vcv.phylo(tree)
		means		    <- rep(mu,n)
		diag(vcv)=diag(vcv)+meserr^2

		#---- RETURN NEGLOGLIKELIHOOD ---
		return( -dmvnorm(y,means,beta*vcv, log=TRUE))

	#----------------------------------
	#-----        EXPON ONLY      -----
	#----------------------------------	
	} else if(model[5]&!(model[1] | model[2] | model[3] | model[4])){
		mu		<-	theta[1]
		beta	<-  inv.logit(theta[2], min=bounds$beta[1], max=bounds$beta[2])
		endRate	<-  inv.logit(theta[3], min=bounds$endRate[1], max=bounds$endRate[2])
		means		    <- rep(mu,n)

		#---- OU TRANSFORMATION
		tree<-exponentialchangeTree(tree, endRate)
		vcv	<-vcv.phylo(tree)
		diag(vcv)=diag(vcv)+meserr^2

		#---- RETURN NEGLOGLIKELIHOOD ---
		return( -dmvnorm(y,means,beta*vcv, log=TRUE))	}else{
		stop("Parameters  \"lambda, \"kappa\" and \"delta\" can only be fit one at a time currently")
	}
}
