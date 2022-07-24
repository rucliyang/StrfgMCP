############################################################
#####                      fSCAD2                      #####
##### Method and Algorithm for two continous responses #####
##### and functional covariates X(t):                  #####
############################################################



##################### L2 norm of a vector #####################
vecNorm <- function(v)
{
    vectorNorm <- sqrt(sum(v^2))
    return(vectorNorm)
}
###############################################################




################## Derivative of SCAD penalty #################
DSCAD <- function(u,lambda,a)
{
    if(u<=lambda) Dp <- lambda
    else if(u< (a*lambda)) Dp <- -(u-a*lambda)/(a-1)
    else Dp <- 0
    return(Dp)
}
###############################################################




################## SCAD penalty #################
PSCAD <- function(u,lambda,a)
{
    if(u<=lambda) Pp <- lambda*u
    else if(u<(a*lambda)) Pp <- -(u^2-2*a*lambda*u+lambda^2)/(2*(a-1))
    else Pp <- (a+1)*lambda^2/2
    return(Pp)
}
###############################################################




################### Derivative of MCP penalty #################
DMCP <- function(u,lambda,a)
{
    if(u<=a*lambda) Dp <- lambda - u/a
    else Dp <- 0
    return(Dp)
}
###############################################################



########################## MCP penalty ########################
PMCP <- function(u,lambda,a)
{
    if(u<=a*lambda) Dp <- lambda*u - u^2/(2*a)
    else Dp <- (a*u^2)/2
    return(Dp)
}
###############################################################






############ Compute weight matrix for W (single) #############
# This function helps to get nonzero part of Wj,
# where Wj is an (M+d) by (M+d) matrix and 
# w_uv is an 5 by 5 matrix with j <= u,v <= j+d 
# j=1,...,M no. of subintervals.
# i.e. for j=10 (10th subinterval), w_uv is nonezero at j=10 to 14.

computeW <- function(basisfd)
{
	L <- basisfd$nbasis
    rng <- getbasisrange(basisfd)
    breaks <- c(rng[1],basisfd$params,rng[2])
    M <- length(breaks) - 1 #number of subintervals of rangeval
    norder <- L-M+1
    
    # creates an array of size M, and each element is a norder x norder matrix (i.e. M list of 5x5 matrix)
	W <- array(0,dim=c(norder,norder,M))
    for (j in 1:M)
    {
    	# calculate inprod of two basis per subinterval
    	# i.e. jth subinterval has range from breaks[j] to breaks[j+1].
        temp <- inprod(basisfd,basisfd,rng=c(breaks[j],breaks[j+1]))
        # W collects nonzero parts per subinterval
        W[,,j] <- temp[j:(j+norder-1),j:(j+norder-1)]
    }
    return(W)
}
###############################################################




############################## PMSE ###########################
calcPMSE <- function(y, fittedY)
{
    pmseY <- mean((y-fittedY)^2)
	return(pmseY) 
}
###############################################################



############################ ajjMat ###########################
# Compute group level ajjMat from xfd or xHat,
# return ajjMat as an J-by-J matrix.
calcAjj <- function(J,K,xHat)
{
	JK=J*K
	ajjMat=matrix(0,J,J)
	evalpoint=seq(0,1,1/(JK-1))
	corMat=cor.fd(evalpoint, xHat) #JK by JK
	for(j1 in 1:J)
   	{	
    	j1Ind <- c(((j1-1)*K+1):(j1*K))
   		for(j2 in j1:J)
    		{
    			j2Ind <- c(((j2-1)*K+1):(j2*K))
    			# corM=corMat[j1Ind, j2Ind]
    			# ajjMat[j1, j2] <- mean(corM[corM>0.2])
    			ajjMat[j1, j2] <- mean(corMat[j1Ind, j2Ind])
    		}
	}
	ajjMat[which(ajjMat <=0.4)]=0
	ajjMat = ajjMat + t(ajjMat)
	diag(ajjMat)=1
	
	return(ajjMat)
}
###############################################################	


calcAjjKvec <- function(J,JK,Kvec,xHat)
{
	ajjMat=matrix(0,J,J)
	evalpoint=cumsum(Kvec)/JK - ((1/JK)/2)
	corMat=cor.fd(evalpoint, xHat) #JK by JK
	for(j1 in 1:J)
   	{	
    	j1Ind=j1
   		for(j2 in j1:J)
    		{
    			j2Ind <- j2
    			ajjMat[j1, j2] <- corMat[j1Ind, j2Ind]
    		}
	}
	ajjMat[which(ajjMat <=0.4)]=0
	ajjMat = ajjMat + t(ajjMat)
	diag(ajjMat)=1
	
	return(ajjMat)
}


RVcoef_xHat <- function(J,JK,Kvec,xHat)
{
  tt=seq(0,1,by=1/(JK-1))
  xEE=t(eval.fd(tt,xHat))
  ind=c(0,cumsum(Kvec))
  ajjGroupMat=matrix(0,J,J)
  for(j1 in 1:J)
  {	
    j1Ind <- c((ind[j1]+1):ind[j1+1])
    for(j2 in j1:J)
    {
      if(j2>j1){
        j2Ind <- c((ind[j2]+1):ind[j2+1])
        # ajjGroupMat[j1, j2] <- mean(covMat[j1Ind, j2Ind])
        # ajjMat[j1Ind, j2Ind] <- mean(covMat[j1Ind, j2Ind])
        ajjGroupMat[j1, j2] <- coeffRV(xEE[,j1Ind], xEE[,j2Ind])$rv
      }
    }
	}
	aaa=ajjGroupMat
	ajjGroupMat[which(ajjGroupMat <=0.2)]=0 #1e-3 for 50 groups
	ajjGroupMat=ajjGroupMat+t(ajjGroupMat)
	diag(ajjGroupMat)=1
	return(ajjGroupMat)
}

RVcoef_xGeno <- function(J,JK, Kcumsum,xGeno)
{
  ind=c(0, Kcumsum)
  ajjGroupMat=matrix(0,J,J)
  for(j1 in 1:J)
  {	
    print(j1)
    j1Ind <- c((ind[j1]+1):ind[j1+1])
    for(j2 in j1:J)
    {
      if(j2>j1){
        j2Ind <- c((ind[j2]+1):ind[j2+1])
        ajjGroupMat[j1, j2] <- coeffRV(xGeno[,j1Ind], xGeno[,j2Ind])$rv
      }
    }
  }
  aaa=ajjGroupMat
  # ajjGroupMat[which(ajjGroupMat <=0.2)]=0 #1e-3 for 50 groups
  ajjGroupMat=ajjGroupMat+t(ajjGroupMat)
  diag(ajjGroupMat)=1
  return(ajjGroupMat)
}



########################## Predict Y ###########################
predict.slos <- function(beta, testX)
{
	x <- testX
	yHat <- 0
	g <- inprod(x$basis, beta$basis)
	yHat <- yHat + t(x$coef) %*% g %*% beta$coef

	return(yHat)
}
###################################################################



############################ TP / FP ##############################
dTPFP <- function(Beta, Beta.i, dimP)
{	
	
	rng <- getbasisrange(Beta$basis)
	rng.a <- rng[1]
	rng.b <- rng[2]
	t <- seq.int(rng.a, rng.b, rng.b/(dimP-1)) 
	fullInd <- c(1: dimP)
	
	# estimated beta
	bEval <- t(eval.fd(t,Beta)) 
	bEvalZero <- (which(bEval==0))
	bEvalNonZero <- fullInd[-bEvalZero]
	
	# true beta
	bTru= Beta.i(t)
	bTruZero <- (which(bTru==0))
	bTruNonZero <- fullInd[-bTruZero]
	
	# TP
	TP = length(intersect(bEvalNonZero, bTruNonZero))/length(bTruNonZero)
	# FP
	FP = length(intersect(bEvalNonZero, bTruZero))/length(bTruZero)
	
	TPFP =list(TP=TP, FP=FP)
}


# This fn is modified for NULL snp inserted data
dTPFP_ins <- function(Beta, Beta.i)
{	
  
  rng <- getbasisrange(Beta$basis)
  rng.a <- rng[1]
  rng.b <- rng[2]
  t.bTru <- seq.int(rng.a, rng.b, rng.b/(5000-1)) 
  t.bEval <- seq.int(rng.a, rng.b, rng.b/(5050-1)) 
  fullInd <- c(1:5050)
  
  # estimated beta
  bEval <- t(eval.fd(t.bEval,Beta)) 
  bEvalZero <- (which(bEval==0))
  bEvalNonZero <- fullInd[-bEvalZero]
  
  # true beta
  bTru= Beta.i(t.bTru)
  bTruNew <- c(bTru[c(1:75)], rep(0,25), bTru[c(76:900)], rep(0,25), bTru[c(901:5000)])
  bTruZero <- (which(bTruNew==0))
  bTruNonZero <- fullInd[-bTruZero]
  
  # TP
  TP = length(intersect(bEvalNonZero, bTruNonZero))/length(bTruNonZero)
  # FP
  FP = length(intersect(bEvalNonZero, bTruZero))/length(bTruZero)
  
  TPFP =list(TP=TP, FP=FP)
}




######################### algorithm ###########################
### method estimation function
method.alg <- function(xfd, xHat, y, J, Kvec, Ajj, lambda1, lambda2, gamma, POpar, a, int, basisfd.beta, Jmat, Ws, max.iter, changeThres, cutoff, estX)
{
	# get some constants
	Ls <- basisfd.beta$nbasis
	JKs <- length(basisfd.beta$params) +1
	ds <- Ls - JKs
	n <- length(y)
	
	Js <- J
	Kvecs <- Kvec
	Tvecs <- Kvecs*(diff(getbasisrange(xfd$basis))/JKs)
	L2NNerVec <- sqrt(Kvecs/Tvecs) #sqrt(K/T_j)
	L2NNerj= L2NNerVec[1]
	# L2NNer <- sqrt(JKs/diff(getbasisrange(xfd$basis))) #sqrt(K/T)
	
	Vs <- NULL
	if(is.null(Vs))
	{
	  Vs <- eval.penalty(basisfd.beta,int2Lfd(2))
	}
    
	# Ws <- NULL
	# if(is.null(Ws))
  #  {
  #      Ws <- computeW(basisfd.beta)  
  #  }
    
	# compute matrix U
    if (estX) {xcoef <- xHat$coefs }else {xcoef <- xfd$coefs} # M+d by n
    Uvec <- t(xcoef)%*%Jmat # n by M+d, each row is an M+d dim vector
    # if (xcoef.mod) {Uvec <- Uvec + 0.001}
    
    if(int){Uvec <- cbind(matrix(1,n,1),Uvec)}
    U <- Uvec   

    #compute matrix V
    VV <- n*gamma*Vs # M+d by M+d
    if(int){VV <- bdiag(0, VV)}
    
    
    # Compute the initial estimate
    bTilde <- solve(t(U)%*%U + VV, t(U)%*%y)
    bHat0 <- bTilde
	# plot(bHat0, type='l')
    
	if(int){muHat <- bHat0[1]}
    
        
	est <- lqa.algCD.last.rPartial(y, Uvec, bHat0, int, VV, Ws, Js, Kvecs, Tvecs, JKs, Ls, ds, L2NNerj, Ajj, lambda1, lambda2, gamma, POpar, a, max.iter, changeThres, cutoff, estX)
	bHatNew <- est$bHat
	SS <- est$lqaS
    
    #compute projected matrix  
    bZero <- (abs(bHatNew) < 1e-2) 
    	bNonZero <- !bZero
    bHatNew[bZero] <- 0  #set those elements to 0 if TRUE
		
	UW <- U[,bNonZero]
    	Vnew <- VV[bNonZero,bNonZero]
    	Snew <- SS[bNonZero,bNonZero]
    	bnew <- solve(t(UW)%*%UW + Snew + Vnew, t(UW)%*%y)
    		
    	bBrave <- matrix(0,sum(Ls),1)
    bBrave[bNonZero,1] <- matrix(bnew,length(bnew),1)

    projMat <- UW %*% solve(t(UW)%*%UW + Snew + Vnew,t(UW))
    fittedY <- projMat %*% y
    
    betaRes <- list(beta=fd(coef=bBrave, basisobj=basisfd.beta), fittedY=fittedY, projMat=projMat)
    return(betaRes)
}
###############################################################




############################ lqa ##############################
### method penalization function for beta1 or beta2.  call this function for only one beta at a time.  penalization is performed seperated for beta1 and beta2.

lqa.algCD.last.rPartial <- function(y, Uvec, bHat0, int, VV, Ws, Js, Kvecs, Tvecs, JKs, Ls, ds, L2NNerj, Ajj, lambda1, lambda2, gamma, POpar, a, max.iter, changeThres, cutoff, estX)
{
	JKmax <- max(JKs)
    Lmax <- max(Ls)
    n <- length(y)
	POvec=Kvecs*a*lambda1/POpar #was b, outer penalty parameter
	
	kcumsum <- c(0,cumsum(Kvecs))
	
	betaNormj <- matrix(0,JKmax,1) # 80 by 1
	betaNormjk <- matrix(0,JKmax,1) # 80 by 1
	betaNormGroupJK <- matrix(0,JKmax,1) # 80 by 1
    bZeroMat <- matrix(FALSE,1,sum(Ls)) # 1 by 84
    bZeroj <- matrix(FALSE,1,sum(Ls))
    bZeroGroupj <- matrix(FALSE,1,sum(Ls))
	betaNorm <- rep(Inf,1)
	
	bZeroMat2 <- matrix(FALSE,1,sum(Ls))
	betaNormjk2 <- matrix(0,JKmax,1)
	# penalization iteration for beta1 and beta2
	record=NULL
	record.rss=NULL
	
	bHatLoop=bHat0
	bHat=bHat0
	r <- y- Uvec%*%bHat
	
	it <- 1
	while (it <= max.iter)
	{
		# Stop criterias
		# change of betaNorm must less than changeThres
		betaNormOld <- betaNorm
		betaNorm <- vecNorm(bHat)
		
		change <- max((betaNormOld-betaNorm)^2)
		if(betaNorm >1e4) stop("betaNorm is too large!", call.=FALSE) 
		if(betaNorm < 1) stop("betaNorm is too small!", call.=FALSE) 
    if(change < changeThres || betaNorm >1e2 || betaNorm < 1) break
		# cat('\n','betaNormOld =', betaNormOld, "; betaNorm =", betaNorm, "; change = ", change, ", at iteration", it)
		
		# Coordinate Descent
        update.bk=matrix(0,Ls,1)
        rss=matrix(0,Ls,1)
        num.pen.bk=matrix(0,Ls,1)
        utu.k=matrix(0,Ls,1)
        den.pen.bk=matrix(0,Ls,1)
        sumWk = matrix(0,Ls,1)
        sumSk = matrix(0,Ls,1)
        wkk=matrix(0,Ls,1)
        skk=matrix(0,Ls,1)

    		###
		# Compute W:
		W <- Ws
		lqaW <- matrix(0,Ls,Ls) # nbasis by nbasis
		PIjk <- matrix(0,JKs,1) 
		PIj <- matrix(0,Js,1)
		DPOcj <- matrix(0,Js,1)
		DPIrecord <- matrix(0,JKs,2)
		for(j in 1:J)
		{
			zeroIndj <- NULL
			index.j <- c((kcumsum[j]+1):(kcumsum[j+1]))
			#Compute Inner Penalty
			lqaWj <- matrix(0,Ls,Ls) # nbasis by nbasis
			for(k in index.j) 
			{
				index.Wjk <- c(k:(k+ds)) # 5-dim vec
				bjk <- bHat[index.Wjk] # 5-dim vec
				betaNormjk[k,1] <- sqrt(t(bjk) %*% W[,,k] %*% bjk)
				DPIcjk <- DMCP(betaNormjk[k,1]*L2NNerj, lambda1, a) #(8) and (10.1)
				DPIrecord[k,1] <- DPIcjk
				DPIrecord[k,2] <- betaNormjk[k,1]
				if(DPIcjk !=0) # power of penalty is nonzero
				{
					# if betaNormj is extremely small, set to zero.
					#if(betaNormjk[k,1] < changeThres) bZeroMat[index.Wjk] <- TRUE
					# o.w. apply shrinkage, i.e. eq (12) from ref. paper.
					lqaWj[index.Wjk, index.Wjk] <- lqaWj[index.Wjk, index.Wjk] + DPIcjk*(L2NNerj/betaNormjk[k,1])*W[,,k] # 5-by-5 matrix
				}
				# Get PIjk to compute outter penalty
				PIjk[k,] <- PMCP(betaNormjk[k,1]*L2NNerj, lambda1, a)
			}
				
			# Compute outter penalty
			PIj[j,1] <- sum(PIjk[index.j, ])
			DPOcj[j,1] <- DMCP(PIj[j,1]/Kvecs[j], lambda1, POvec[j]) #(8) and (10.1)
				
			lqaWj <- DPOcj[j,1]*lqaWj/2
			lqaW <- lqaW + lqaWj
		}
		bZeroVec <- bZeroMat
    		bNonZeroVec <- !bZeroVec
    			
    		if(int)
    		{
    			lqaW <- direct.sum(0,lqaW)
    			bNonZeroVec <- as.vector(c(TRUE,bNonZeroVec))
    		}
			
		# Compute L2-norm of beta_j (i.e. group level beta)
		betaNormGroupJ <- matrix(0,Js,1) 
		for(bj in 1:Js)
    		{
    			index.bj <- c((kcumsum[bj]+1):(kcumsum[bj+1]))
    			for(k in index.bj)
    			{
				bjk <- bHat[c(k:(k+ds))]
				betaNormGroupJK[k,1] <- t(bjk) %*% W[,,k] %*% bjk
    			}
    			betaNormGroupJ[bj,1] <- sqrt(sum(betaNormGroupJK[index.bj,]))
    		}
    			
    		# Compute lqaS matrix
    		Tvecs_group=Tvecs
    		Wjk <- Ws
    		lqaS <- matrix(0,Ls,Ls) 
    		for(j1 in 1:Js)
    		{	
    			index.j1 <- c((kcumsum[j1]+1):(kcumsum[j1+1]))
    			Tvecs_j1=Tvecs[j1] #Compute length of j1^th network
    			Dlapj <- 0
   			for(j2 in 1:Js)
    			{
    				index.j2 <- c((kcumsum[j2]+1):(kcumsum[j2+1]))
    				Tvecs_j2=Tvecs[j2] #Compute length of j1^th network
    				if(j2 != j1)
    				{
    					Dlap.jj=NULL
    					ajjValue <- Ajj[j1, j2] 
    					if(ajjValue != 0)
    					{
    						Dlap.jj <- abs(ajjValue)*(betaNormGroupJ[j1,1]/sqrt(abs(Tvecs_j1)) - betaNormGroupJ[j2,1]/sqrt(abs(Tvecs_j1)))
    					}else{Dlap.jj <- 0} 
    					Dlapj <- Dlapj + Dlap.jj
    				}
    			}
    			if(Dlapj != 0)
    			{
    				lqaS.Wjk <- matrix(0,Ls,Ls) 
    				for(k in index.j1)
    				{
    					index.jk <- c(k:(k+ds))
    					lqaS.Wjk[index.jk,index.jk] <- lqaS.Wjk[index.jk,index.jk] + Wjk[,,k]	
    				}
    				lqaS <- lqaS + (lambda2*Dlapj/(betaNormGroupJ[j1,1]*sqrt(abs(Tvecs_j1)))) *lqaS.Wjk
    			}else{lqaS <- lqaS}
   	 	}
		lqaS= lqaS/2
			
        for(ki in 1:Ls)
        {
        		rp <- r + Uvec[,ki]*bHat[ki,]
        	
        		Wp <- lqaW[ki,-ki,drop=F]
        		Sp <- lqaS[ki,-ki,drop=F]
        		Vp <- VV[ki,-ki,drop=F]
        		WSVp <- (n*Wp + n*Sp + Vp)%*%bHat[-ki,]
        		# Coordinate Descent 
        		update.bk[ki,] = (sum(Uvec[,ki]*rp) - WSVp) / (t(Uvec[,ki])%*% Uvec[,ki] + n*lqaW[ki,ki]+n*lqaS[ki,ki]+VV[ki,ki])
        		bHat[ki,] = update.bk[ki,]
			
			r <- rp - Uvec[,ki]*bHat[ki,]
        } #end of ki
        
		it <- it + 1
		#print(cbind(bHat, update.bk, b.it))
	} #end of it
	#plot(bHat,type="l")
	
	b.it=bHat
	bNew = bHat
	for(k2 in 1:JKs)
    {
    		ind.k2 <- c(k2:(k2+ds))
		bNewjk <- bNew[ind.k2]
		betaNormjk2[k2,1] <- t(bNewjk) %*% W[,,k2] %*% bNewjk
		if(betaNormjk2[k2,1] < cutoff) bZeroMat2[ind.k2] <- TRUE
   	}
    bZeroVec2 <- bZeroMat2
    bNonZeroVec2 <- !bZeroVec2
    	
    bHat <- matrix(0,JKs+ds,1)
    bHat[bNonZeroVec2,1]=bNew[bNonZeroVec2,1]
	bHatLoop=cbind(bHatLoop,bHat)
	
	result <- list(bHat=bHat, b.it=b.it, lqaS=lqaS)		
	return(result)
	#result <- return(bHat=bHat, betaNormjk2= betaNormjk2, bNew=bNew)		
}
###############################################################




########################### CV Tuning ########################
# Options of tuning: BIC, AIC, RIC, VALIDATION, CV
tuneProgram <- function(tuning="AIC/BIC/RIC", xfd, xHat, y, J, Kvec, Ajj, lambda1, lambda2, gamma, POpar, a, int, basisfd.beta, Jmat, Ws, max.iter, changeThres, cutoff, estX)
{
	n <- length(y)
	ss <- NULL
	for(i in 1:length(gamma))
	{
		for(j in 1:length(lambda1))
		{
			for(k in 1: length(lambda2))
			{
				if(tuning != 'CV')
          		{
					no.error <- TRUE
					tryCatch( {fit=method.alg(xfd, xHat, y, J, Kvec, Ajj, lambda1[j], lambda2[k], gamma[i], POpar, a, int, basisfd.beta, Jmat, Ws, max.iter, changeThres, cutoff, estX)}, error = function(e) {no.error <<- FALSE})
					if(no.error) {
						resid <- y-fit$fittedY
                		rss <- sum(resid^2)
                		df <- sum(diag(fit$projMat))

                    	sig2 <- rss/(n-df)
                    	errRIC <- (n-df)*log(sig2)+df*(log(n)-1)+4/(n-df-2)

                    	errAIC = n*log(rss/n)+2*df
                    	errBIC = n*log(rss/n) + log(n)*df
               		}else{
               			errRIC=NA
               			errAIC=NA
               			errBIC=NA
               		}
               		errCV=NA
                }else if(tuning == 'CV')
            	{
		        		Kfold <- 10
		        		idx <- sample(1:n,n)
		        		s <- ceiling(n/Kfold)
		       		err.cv <- matrix(0,1,Kfold)
		        		for(f in 1:Kfold)
	            		{
	            			test.index <- idx[(s*(f-1)+1):min(s*f,n)]
	                		train.index <- idx[-c((s*(f-1)+1):min(s*f,n))]
						
	                		no.error3 <- TRUE
	                		tryCatch( {fit.cv=method.alg(xfd[train.index], xHat[train.index], y[train.index,,drop=F], J, Kvec, Ajj, lambda1[j], lambda2[k], gamma[i], POpar, a, int, basisfd.beta, Jmat, Ws, max.iter, changeThres, cutoff, estX)}, error=function(e) {no.error2 <<- FALSE})
	                		if(no.error3){
	                		# yhat <- predict(fit.cv,xHat,test.index)
                    		# resid <- y[test.index] - yhat
                    		# err.cv[f] <- sum(resid^2)
						beta.cv <- fit.cv$beta
						yhat.cv <- predict.slos(beta.cv, xHat[test.index])
						resid <- y[test.index] - yhat.cv
						err.cv[f] <- sum(resid^2)
						}else{
						err.cv[f] <- NA	
						}
					}
	            		errCV <- sum(err.cv)	
	            		errRIC=NA
               		errAIC=NA
               		errBIC=NA
            	}
            	
	            ss.i <- cbind(gamma[i], lambda1[j], lambda2[k], errAIC, errBIC, errRIC) 
	            ss <- rbind(ss, ss.i)
			}
		}
	}
	colnames(ss) <- c('gamma', 'lambda1', 'lambda2', 'errAIC', 'errBIC', 'errRIC')
	row.has.na <- apply(ss, 1, function(x){any(is.na(x))})
	ss1LL <- ss[!row.has.na,]

	if(is.matrix(ss1LL))
	{
	  if(dim(ss1LL)[1]==0){
	    selMat = rep(NA,8)
	  }else{		
		selectAIC <- which(ss1LL[,'errAIC']==min(as.numeric(ss1LL[,'errAIC'])), arr.ind=T)
		selectBIC <- which(ss1LL[,'errBIC']==min(as.numeric(ss1LL[,'errBIC'])), arr.ind=T)
		selectRIC <- which(ss1LL[,'errRIC']==min(as.numeric(ss1LL[,'errRIC'])), arr.ind=T)
		
		selMat <- ss1LL[c(selectAIC[1], selectBIC[1], selectRIC[1]),]
		nam <- c('selectAIC','selectBIC','selectRIC' )
		selMat <- cbind(nam,selMat)
	  }
	}else if(is.vector(ss1LL)){
		selMat <- ss1LL
		nam='AIC/BIC/RIC'
		selMat <- cbind(nam,selMat)
	}

	res <- list(selMat=selMat, tuneRes=ss1LL)
}
###############################################################




################################# Estimate datasets ###########################
estimate.datasets <- function(nruns, rangeval, data, lambda1, lambda2, gamma, Kvec, J, POpar, a, max.iter, changeThres, cutoff, int, estX, Beta.i, MAF.text)
{
	Kvec <- Kvec
	J <- J
	subint <- 25 # 25 obs/subint
	Kobs <- Kvec*subint
	p <- sum(Kobs) #5000
	JK <- p/subint # no. of subintervals = M = JK

	nknots <- JK+1 #nknots=M+1
	norder <- 5 #norder=d+1
	d=norder-1
	nb <- nknots + norder - 2
	basisfd <- create.bspline.basis(rangeval= rangeval, nbasis=nb, 	norder=norder)
	basisfd.beta <- basisfd
	
	#compute Ws
	Ws <- NULL
	if(is.null(Ws)) Ws <- computeW(basisfd.beta) 
	
	t <- seq.int(0, 1, 1/(JK-1))
	
	coefMatAIC <- NULL
  	coefMatBIC <- NULL
  	coefMatRIC <- NULL
	evalInd <- matrix(0,nruns,18)
	for (irun in 1:nruns)
	{
		cat("g50 iteration = ", irun, "\n")
		
		### Get training data
		xfd <- data$train[[irun]]$xfd
		xTru <- data$train[[irun]]$xTru
		y <- data$train[[irun]]$y
		if(MAF.text=="MAF0.5"){
			xHat <- data$train[[irun]]$xHat1
			xGeno <- data$train[[irun]]$xGeno1
		}else if(MAF.text=="MAF0.4"){
			xHat <- data$train[[irun]]$xHat2
			xGeno <- data$train[[irun]]$xGeno2
		}
		
		# Compute Ajj
		Ajj <- calcAjjKvec(J,JK,Kvec,xHat)
		diag(Ajj)=0

		# Compute Jmat
		Jmat <- inprod(xfd$basis, basisfd.beta)

		### Get test data
		if(estX) 
		{
			if(MAF.text=="MAF0.5")
			{
				testX <- data$test[[irun]]$xHat1
			}else if(MAF.text=="MAF0.4"){
				testX <- data$test[[irun]]$xHat2
			}
		}else 
		{testX <- data$test[[irun]]$xfd}
		testY <- data$test[[irun]]$y
		
		
		#10-fold CV tuning
		tunedRes <- tuneProgram(tuning="AIC/BIC/RIC", xfd, xHat, y, J, Kvec, Ajj, lambda1, lambda2, gamma, POpar, a, int, basisfd.beta, Jmat, Ws, max.iter, changeThres, cutoff, estX)
		AICgam=as.numeric(tunedRes$selMat[1,'gamma'])
		AIClam1=as.numeric(tunedRes$selMat[1,'lambda1'])
		AIClam2=as.numeric(tunedRes$selMat[1,'lambda2'])
		
		BICgam=as.numeric(tunedRes$selMat[2,'gamma'])
		BIClam1=as.numeric(tunedRes$selMat[2,'lambda1'])
		BIClam2=as.numeric(tunedRes$selMat[2,'lambda2'])
		
		RICgam=as.numeric(tunedRes$selMat[3,'gamma'])
		RIClam1=as.numeric(tunedRes$selMat[3,'lambda1'])
		RIClam2=as.numeric(tunedRes$selMat[3,'lambda2'])
		
		
		if(is.na(AICgam) | is.na(BICgam) | is.na(RICgam)){
		  # beta1.evalfd[irun,] <- NA
		  # beta2.evalfd[irun,] <- NA
		  evalInd[irun,] <- cbind(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,NA, NA, NA, NA) 
		}else{
  		
  		tryCatch( {AICfit=method.alg(xfd, xHat, y, J, Kvec, Ajj, AIClam1, AIClam2, AICgam, POpar, a, int, basisfd.beta, Jmat, Ws, max.iter, changeThres, cutoff, estX)}, error=function(e) {cat("ERROR :",conditionMessage(e), "irun =",irun, "lambda1 =", lam1, ", lambda2 =", lam2, ", gamma =", gam, "\n")})
		AICbeta=AICfit$beta
		# Compute TP/FP
		AIC.tpfp.disc=dTPFP(AICbeta, Beta.i, dimP=5000)
		AIC.dTP=AIC.tpfp.disc$TP
		AIC.dFP=AIC.tpfp.disc$FP
  		# Predict y
  		AIC.predY <- predict.slos(AICbeta, testX)
  		AIC.pmseY <- calcPMSE(y, AIC.predY)
  		
  		
  		tryCatch( {BICfit=method.alg(xfd, xHat, y, J, Kvec, Ajj, BIClam1, BIClam2, BICgam, POpar, a, int, basisfd.beta, Jmat, Ws, max.iter, changeThres, cutoff, estX)}, error=function(e) {cat("ERROR :",conditionMessage(e), "irun =",irun, "lambda1 =", lam1, ", lambda2 =", lam2, ", gamma =", gam, "\n")})
		BICbeta=BICfit$beta
		# Compute TP/FP
		BIC.tpfp.disc=dTPFP(BICbeta, Beta.i, dimP=5000)
		BIC.dTP=BIC.tpfp.disc$TP
		BIC.dFP=BIC.tpfp.disc$FP
  		# Predict y
  		BIC.predY <- predict.slos(BICbeta, testX)
  		BIC.pmseY <- calcPMSE(y, BIC.predY)
  		
  		
  		tryCatch( {RICfit=method.alg(xfd, xHat, y, J, Kvec, Ajj, RIClam1, RIClam2, RICgam, POpar, a, int, basisfd.beta, Jmat, Ws, max.iter, changeThres, cutoff, estX)}, error=function(e) {cat("ERROR :",conditionMessage(e), "irun =",irun, "lambda1 =", lam1, ", lambda2 =", lam2, ", gamma =", gam, "\n")})
		RICbeta=RICfit$beta
		# Compute TP/FP
		RIC.tpfp.disc=dTPFP(RICbeta, Beta.i, dimP=5000)
		RIC.dTP=RIC.tpfp.disc$TP
		RIC.dFP=RIC.tpfp.disc$FP
  		# Predict y
  		RIC.predY <- predict.slos(RICbeta, testX)
  		RIC.pmseY <- calcPMSE(y, RIC.predY)
  		
  		
  		coefMatAIC <- cbind(coefMatAIC, AICbeta$coefs)
  		coefMatBIC <- cbind(coefMatBIC, BICbeta$coefs)
  		coefMatRIC <- cbind(coefMatRIC, RICbeta$coefs)
		
  		evalInd[irun,] <- cbind(AICgam, AIClam1, AIClam2, AIC.dTP, AIC.dFP, AIC.pmseY, BICgam, BIClam1, BIClam2, BIC.dTP, BIC.dFP, BIC.pmseY, RICgam, RIClam1, RIClam2, RIC.dTP, RIC.dFP, RIC.pmseY) 
  	}
	}
	
	AIC.beta.fd <- fd(coef= coefMatAIC, basisobj=basisfd.beta)
	BIC.beta.fd <- fd(coef= coefMatBIC, basisobj=basisfd.beta)
	RIC.beta.fd <- fd(coef= coefMatRIC, basisobj=basisfd.beta)

	colnames(evalInd) <- c('AICgam', 'AIClam1', 'AIClam2', 'AIC.dTP', 'AIC.dFP', 'AIC.pmseY', 'BICgam', 'BIClam1', 'BIClam2', 'BIC.dTP', 'BIC.dFP', 'BIC.pmseY', 'RICgam', 'RIClam1', 'RIClam2', 'RIC.dTP', 'RIC.dFP', 'RIC.pmseY') 
	eval = na.omit(evalInd)
	eval.ave <- round(apply(evalInd[,c('AIC.dTP', 'AIC.dFP', 'AIC.pmseY', 'BIC.dTP', 'BIC.dFP', 'BIC.pmseY', 'RIC.dTP', 'RIC.dFP', 'RIC.pmseY')],2,mean), digits=4)
	eval.sd <- round(apply(evalInd[,c('AIC.dTP', 'AIC.dFP', 'AIC.pmseY', 'BIC.dTP', 'BIC.dFP', 'BIC.pmseY', 'RIC.dTP', 'RIC.dFP', 'RIC.pmseY')],2,sd), digits=4)

	result <- list(evalInd=evalInd, eval.ave=eval.ave, eval.sd=eval.sd, Ajj=Ajj, AIC.beta.fd=AIC.beta.fd, BIC.beta.fd=BIC.beta.fd, RIC.beta.fd=RIC.beta.fd)
}
###########################################################################





