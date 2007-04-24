`fit.continuous` <-
function(tips, phy, lambda=FALSE, kappa=FALSE, delta=FALSE, alpha=FALSE, r=FALSE, bounds=NULL, print=TRUE, meserr=0)
{

	
    if (!is.vector(tips))
    	stop("Currently only one set of TIP values can be analysed")
    if (is.null(names(tips)))
    	stop("the variable \"tips\" must have names specified for elements, use \"names()\" function or see help")
    #--------------------------------
    #---    CALCULATIONS       ---
    #--------------------------------	
 	vcv 	    <- vcv.phylo(phy, mode="Brownian") # Variance - Covariance Matrix
 	if(nrow(vcv)!=length(tips))
 		stop("Object \"tips\" and object \"phy\" are not of the same length")
	#--------------------------------
    #---     SORT TIP DATA        ---
    #--------------------------------
    spp.dat <- names(tips)
    spp.vcv <- rownames(vcv)
    n <- length(tips)
    if (sum(spp.dat == spp.vcv) != length(spp.vcv)){
    	#warning("TIP data was sorted to match tree")
        index <- numeric()              
        for (i in 1:n){
            index[i] <- (1:n)[spp.vcv[i] == spp.dat]
        }
        if(sum(is.na(index)) > 0)
        	stop("Cannot match object \"tips\" to object \"phy\", species names do not match")
    	tips <- tips[index]
    }
    spp.dat <- names(tips)
    spp.vcv <- rownames(vcv)
    if (sum(spp.dat == spp.vcv) != length(spp.vcv))
    	stop("Cannot match object \"tips\" to object \"phy\"")
	#--------------------------------
    #---    PREPARE DATA LIST     ---
    #--------------------------------
	data			<- list()
   		data$obs 		<- tips          # TIP data 
    	data$spp.name	<- names(tips)	 # SPP Names of TIP data
    	data$tree 		<- phy
    	data$meserr		<- meserr		# Standard errors of the mean
    #--------------------------------
    #--- SET MODEL SPECIFICATIONS ---
    #--------------------------------
    model<-c(lambda, kappa, delta, alpha, r)
	names(model)<- c("lambda", "kappa", "delta", "alpha", "r")
	if (sum(model) > 1)
		stop("Currently, lambda, kappa, delta, and alpha can only be fit one at a time")
    #-----------------------------
    #---  SET PARAMETER BOUNDS ---
    #-----------------------------
    #---- DEFAULT BOUNDS
    bounds.default			 <- matrix(c(0.00001, 20, 0,1, 0.000001, 5, 0, 12, 0, 5, -100, 100), nrow=6, ncol=2, byrow=TRUE)
    rownames(bounds.default) <- c("beta", "lambda", "kappa", "delta", "alpha", "r");
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
 			bounds.user <- matrix(c(bounds$beta, bounds$lambda, bounds$kappa, bounds$delta, bounds$alpha, bounds$r), 
 								  nrow=sum(specified), ncol=2, byrow=TRUE
 								  )
 			rownames(bounds.user) <- c("beta", "lambda", "kappa", "delta", "alpha", "r")[specified]
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
  	data$design$bounds <- data.frame(t(bounds))
  	data$design$model  <- model
  	#--------------------------------
    #---        FIT MODEL         ---
    #--------------------------------
  	fit.continuous.model(data, print=print)
}

`fit.continuous.model` <-
function(data, print=TRUE)
{
	bounds 	<- data$design$bounds
	model 	<- data$design$model
	np 		<- sum(model)
	n 		<- length(data$obs)
	#--- INITIALIZE RESULTS MATRIX --
    results <- matrix(nrow=2+np, ncol=4)
    colnames(results)<- c("Estimates", "SE", "Low.CI", "Upper.CI")
    rownames(results)<- c("mu", "beta", c("lambda", "kappa", "delta", "alpha", "r")[model])
	#----- MINIMIZE NEGATIVE LOG LIKELIHOOD
	theta.start <-c(0, 0.1,c(0.5, 0.5, 0.5, 0.5, 0)[model])     # Starting point for profile search
	out         <- NULL
	out	<- nlm(negloglike, theta.start, hessian=TRUE, data=data)
	
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
	#----------------------------------------
	#--- OUTPUT RESULTS (PRINT OR RETURN) ---
	#----------------------------------------
	if (print==TRUE){
		print(paste("Maximum Likelihood:", round(-out$minimum, digits=3)), quote=FALSE)
		print(results)
	}else{
	 	return(list(estimates=results, lnl=-out$minimum)) 
	}
}


