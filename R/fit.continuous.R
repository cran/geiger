`fit.continuous` <-
function(phy, data, data.names=NULL, lambda=FALSE, kappa=FALSE, delta=FALSE, ou=FALSE, eb=FALSE, bounds=NULL,  meserr=NULL)
{
	
	# sort is T because sub-functions assume data are in
	# this particular order
	
	td<-treedata(phy, data, data.names, sort=T)

	ntax=length(td$phy$tip.label)

	if(is.null(meserr)) {
		me=td$data
		me[]=0
		meserr=me	
	} else if(length(meserr)==1) {
		me=td$data
		me[]=meserr
		meserr=me
	} else if(is.vector(meserr)) {
		if(!is.null(names(meserr))) {
			o<-match(rownames(td$data), names(meserr))
			if(length(o)!=ntax) stop("meserr is missing some taxa from the tree")
			meserr<-as.matrix(meserr[o,])
		} else {
			if(length(meserr)!=ntax) stop("No taxon names in meserr, and the number of taxa does not match the tree")
			me<-td$data
			me[]=meserr
			meserr=me
		}
	} else {
		if(!is.null(rownames(meserr))) {
			o<-match(rownames(td$data), rownames(meserr))
			meserr=meserr[o,]
		} else {
			if(sum(dim(meserr)!=dim(td$data))!=0)
				stop("No taxon names in meserr, and the number of taxa does not match the tree")
			print("No names in meserr; assuming that taxa are in the same order as tree")	
		}
	}

	#--------------------------------
    #---    PREPARE DATA LIST     ---
    #--------------------------------
	ds			<- list()
   		ds$tree 		<- td$phy          # TIP data 
    #--------------------------------
    #--- SET MODEL SPECIFICATIONS ---
    #--------------------------------
    model<-c(lambda, kappa, delta, ou, eb)
	names(model)<- c("lambda", "kappa", "delta", "ou", "eb")
	if (sum(model) > 1)
		stop("Currently, lambda, kappa, delta, ou, and eb can only be fit one at a time")
    #-----------------------------
    #---  SET PARAMETER BOUNDS ---
    #-----------------------------
    #---- DEFAULT BOUNDS
    bounds.default			 <- matrix(c(0.00001, 20, 0.0000001,1, 0.000001, 1, 0.00001, 5, 0, 5, -10, 10), nrow=6, ncol=2, byrow=TRUE)
    rownames(bounds.default) <- c("beta", "lambda", "kappa", "delta", "alpha", "a");
    colnames(bounds.default) <- c("min", "max")

 	#---- USER DEFINED PARAMETER BOUNDS
 	if (is.null(bounds)) {
 		bounds <- bounds.default       # USE DEFAULTS
 	}else{
 		if (class(bounds)!="list"){
 			stop("Please specify user defined parameter bounds as a list()")
 		}else{
 			specified   <- !c(is.null(bounds$beta), is.null(bounds$lambda), 
 							  is.null(bounds$kappa), is.null(bounds$delta),  is.null(bounds$alpha), is.null(bounds$r)
 							  )
 			bounds.user <- matrix(c(bounds$beta, bounds$lambda, bounds$kappa, bounds$delta, bounds$alpha, bounds$endRate), 
 								  nrow=sum(specified), ncol=2, byrow=TRUE
 								  )
 			rownames(bounds.user) <- c("beta", "lambda", "kappa", "delta", "alpha", "endRate")[specified]
   	 		colnames(bounds.user) <- c("min", "max")
  
   	 		#----  SET FINAL SEARCH BOUNDS
 			bounds <- bounds.default
 			bounds[specified,] <- bounds.user     # Final Bounds
   		} # END if list
   	}  # END user bound if loop
   	#--------------------------------
    #---   APPEND MODEL SETTINGS  ---
    #--------------------------------
  	ds$bounds <- data.frame(t(bounds))
  	ds$model  <- model
  	#--------------------------------
    #---        FIT MODEL         ---
    #--------------------------------
    result<-list()
    for(i in 1:ncol(td$data)) {
    	ds$data=td$data[,i]
    	ds$meserr=meserr[,i]
  		result[[i]]<-fitContinuousModel(ds, print=print)
  	}
  	result
}


`fitContinuousModel` <-
function(ds, print=TRUE)
{
	bounds 	<- ds$bounds
	model 	<- ds$model
	np 		<- sum(model)
	n 		<- length(ds$data)
	#--- INITIALIZE RESULTS MATRIX --
    results <- numeric(2+np)
    names(results)<- c("mu", "beta", c("lambda", "kappa", "delta", "alpha", "endRate")[model])
	#----- MINIMIZE NEGATIVE LOG LIKELIHOOD
	
	beta.start<-var(ds$data)/max(branching.times(ds$tree))
	theta.start <-c(log(beta.start),c(log(0.5), log(0.5), log(0.5), log(0.1), 0.0001)[model])     # Starting point for profile search
	lower=log(bounds[1,][c(1, model)==1])
	upper=log(bounds[2,][c(1, model)==1])

	out         <- NULL
	
	y			<- ds$data				# TIP data
	tree		<- ds$tree			# Tree
	meserr		<- ds$meserr
	n			<- length(y)
	
	if (sum(names(model)==c("lambda", "kappa", "delta", "ou", "eb"))!=5)
		stop("Error with \"model\" internal paramter")


	#----------------------------------
	#-----       DEFAULT FIT      -----
	#----------------------------------
	if (sum(model)==0) {
		

		vcv<-vcv.phylo(tree)

		
		foo<-function(x) {
			vv<-exp(x)*vcv
			diag(vv)<-diag(vv)+meserr^2
			mu<-phylogMean(vv, y)
			mu<-rep(mu, n)
			-dmvnorm(y, mu, vv, log=T)
		}
		o<-optim(foo, p=theta.start, lower=lower, upper=upper, method="L")
			
		results<-list(lnl=-o$value, beta= exp(o$par))

	#----------------------------------
	#-----       LAMBDA ONLY      -----
	#----------------------------------
	} else if (model[1] & !(model[2] | model[3])){
		
				
		
		foo<-function(x) {


			vcv<-vcv.phylo(tree)

			index			<-	matrix(TRUE, n,n)
			diag(index)		<- FALSE
			vcv[index] 	<- vcv[index]*exp(x[2])
			
			vv<-exp(x[1])*vcv

			
			diag(vv)<-diag(vv)+meserr^2
			mu<-phylogMean(vv, y)
			mu<-rep(mu, n)
			
			-dmvnorm(y, mu, vv, log=T)
		}
		o<-optim(foo, p=theta.start, lower=lower, upper=upper, method="L")
			
		results<-list(lnl=-o$value, beta= exp(o$par[1]), lambda=exp(o$par[2]))

	#----------------------------------
	#-----        KAPPA ONLY      -----
	#----------------------------------
	} else if (model[2] & !(model[1] | model[3])){
		
				
		
		foo<-function(x) {

			t<-kappaTree(tree, kappa=exp(x[2]))
			vcv<-vcv.phylo(t)

			
			vv<-exp(x[1])*vcv

			
			diag(vv)<-diag(vv)+meserr^2
			mu<-phylogMean(vv, y)
			mu<-rep(mu, n)
			
			-dmvnorm(y, mu, vv, log=T)
		}
		o<-optim(foo, p=theta.start, lower=lower, upper=upper, method="L")
		
		results<-list(lnl=-o$value, beta= exp(o$par[1]), lambda=exp(o$par[2]))


	#----------------------------------
	#-----        DELTA ONLY      -----
	#----------------------------------	
	} else if (model[3] & !(model[1] | model[2])){
		
		foo<-function(x) {

			t<-deltaTree(tree, delta=exp(x[2]))
			vcv<-vcv.phylo(t)

			
			vv<-exp(x[1])*vcv

			
			diag(vv)<-diag(vv)+meserr^2
			mu<-phylogMean(vv, y)
			mu<-rep(mu, n)
			
			-dmvnorm(y, mu, vv, log=T)
		}
		o<-optim(foo, p=theta.start, lower=lower, upper=upper, method="L")
		
		results<-list(lnl=-o$value, beta= exp(o$par[1]), delta=exp(o$par[2]))	#----------------------------------
	#-----        ALPHA ONLY      -----
	#----------------------------------			
	} else if (model[4] & !(model[1] | model[2] | model[3])){
		
		
		
		foo<-function(x) {
			t<-ouTree(tree, exp(x[2]))

			vcv<-vcv.phylo(t)
			
			vv<-exp(x[1])*vcv
			diag(vv)<-diag(vv)+meserr^2
			
			mu<-phylogMean(vv, y)
			mu<-rep(mu, n)
			
			-dmvnorm(y, mu, vv, log=T)
		}
		o<-optim(foo, p=theta.start, lower=lower, upper=upper, method="L")
			
		results<-list(lnl=-o$value, beta= exp(o$par[1]), alpha=exp(o$par[2]))


	#----------------------------------
	#-----        EXPON ONLY      -----
	#----------------------------------	
	} else if(model[5]&!(model[1] | model[2] | model[3] | model[4])){

		
		foo<-function(x) {
			t<-exponentialchangeTree(tree, a=x[2])

			vcv<-vcv.phylo(t)
			
			vv<-exp(x[1])*vcv
			diag(vv)<-diag(vv)+meserr^2
			
			mu<-phylogMean(vv, y)
			mu<-rep(mu, n)
			
			-dmvnorm(y, mu, vv, log=T)
		}
		o<-optim(foo, p=theta.start, lower=lower, upper=upper, method="L")
			
		results<-list(lnl=-o$value, beta= exp(o$par[1]), a=o$par[2])	}else{
		stop("Parameters  \"lambda, \"kappa\" and \"delta\" can only be fit one at a time currently")
	}
	
	k=np+1
	results$aic<-2*k-2*results$lnl
	results$aicc<-2*k*(n-1)/(n-k-2)-2*results$lnl

	return(results) 

}

phylogMean<-function(phyvcv, data) 
{
	o<-rep(1, length(data))
	ci<-solve(phyvcv)
	
	m1<-solve(t(o) %*% ci %*% o)
	m2<-t(o) %*% ci %*% data
	
	return(m1 %*% m2)
	
	}