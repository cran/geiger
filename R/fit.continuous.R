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
    bounds.default			 <- matrix(c(0.00001, 20, 0,1, 0.000001, 1, 0.00001, 5, 0, 5, 0.00001, 100), nrow=6, ncol=2, byrow=TRUE)
    rownames(bounds.default) <- c("beta", "lambda", "kappa", "delta", "alpha", "endRate");
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
   	 		#----  NOTIFICATION
   	 		print("Warning: The following user defined parameter bounds have been set:", quote=FALSE)
   	 		print(bounds.user)
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
  		result[[i]]<-fit.continuous.model(ds, print=print)
  	}
  	result
}

`fit.continuous.model` <-
function(ds, print=TRUE)
{
	bounds 	<- ds$bounds
	model 	<- ds$model
	np 		<- sum(model)
	n 		<- length(ds$data)
	#--- INITIALIZE RESULTS MATRIX --
    results <- matrix(nrow=2+np, ncol=4)
    colnames(results)<- c("Estimates", "SE", "Low.CI", "Upper.CI")
    rownames(results)<- c("mu", "beta", c("lambda", "kappa", "delta", "alpha", "endRate")[model])
	#----- MINIMIZE NEGATIVE LOG LIKELIHOOD
	theta.start <-c(0, 0.1,c(0.5, 0.5, 0.5, 0.5, 0)[model])     # Starting point for profile search
	out         <- NULL
	out	<- nlm(negloglike, theta.start, hessian=TRUE, ds=ds)
	
	model.full <-c(TRUE, model)
	bounds <- t(bounds)
	#----------------------------------------
	#-----  POINT ESTIMATES & VARIANCES -----
	#----------------------------------------
	results[1,1]	<- 	out$estimate[1]			   #---   MU - No back transformation
	i <- 1
	while (i < sum(model) + 2){
		i <- i+1
		results[i,1] <-	inv.logit(out$estimate[i], 
								min=bounds[which(model.full)[i-1],1], 
								max=bounds[which(model.full)[i-1],2]
							)
	}
	#----- VARIANCE AND CI, REGULR PARAMETER SPACE
	par.var			<-	NULL
	inv.fish.info 	<- solve(out$hessian)     # inverse of fisher info = variance if untransformed
	fit.var 		<- diag(inv.fish.info)   
	par.var[1]  	<- fit.var[1]             # variance of mu untransformed
	i <- 1
	while (i < sum(model)+2){
		i <- i+1
		y <- out$estimate[i]
		#--- Delta Transformation of error term
		par.var[i]  <- fit.var[i]*(((exp(y)*(1+exp(y))-exp(2*y))/(1+exp(y))^2) * 
						(bounds[which(model.full)[i-1],2] - bounds[which(model.full)[i-1],1]))^2
	}	
	#-------- CONFIDENCE INTERVALS -----
	conf.level  <- 0.95
	crit.val    <- qnorm((1 + conf.level) / 2)
	for (i in 1:length(out$estimate))
	{
		results[i,2]   <- round(sqrt(par.var[i]), digits=4)
	   	results[i,3:4] <- results[i,1] + c(-1, 1) * crit.val *sqrt(par.var[i])
	}
	lnl=-out$minimum
	k=1+sum(model)
	return(list(estimates=results, lnl=-out$minimum, aic=2*k-2*lnl)) 
	
}

inv.logit<-function (x, min = 0, max = 1) 
{
    p <- exp(x)/(1 + exp(x))
    p <- ifelse(is.na(p) & !is.na(x), 1, p)
    p * (max - min) + min
}
