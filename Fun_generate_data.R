
####################### Calculate InnerProd #######################
InnerProd <- function(Betaf,basisfd,j) {
	# compute the <beta_j, B_j>, integral of beta_j and B_j. 
	# Betaf: beta function, i.e. Beta1, Beta2, Beta3, Beta4, predefined.
	# basis: basis function, jth column of basismatrix (eval.basis) (t rows and nbasis columns)
	
	rng <- getbasisrange(basisfd)
	knots <- basisfd$param
	knots <- c(rng[1], knots, rng[2])
	nbasis <- basisfd$nbasis
	norder <- nbasis - length(knots) + 2
	
	a <-rng[1]
	if(j-norder > 0) {a <- knots[j-norder+1]}
	
	b <- rng[2]
	if(j <= (nbasis-norder)) {b <- knots[j+1]}
	
	BFun <- function(t) {
		basismatrix <- eval.basis(t,basisfd) # 71 by 74 matrix, t rows and nbasis column
		basismatrix.j <- t(basismatrix[,j]) #get jth column of basismatrix
	return(basismatrix.j)
	}
	
	flog <- function(t) {Betaf(t)*BFun(t)}
	in.prod <- integrate(flog,a,b) 
	return(in.prod$value)
}
###########################################################

################### Generate Data with EstX ##################
### Estimation of genetic variant function                 ###
### Convert continous values to 0,1,2                      ###
### Then estimate it to hat(X_i(t))                        ###

EstX <- function(xTru, p, basisfd){
	nbasis <- basisfd$nbasis
	rng <- getbasisrange(basisfd)
	rng.a <- rng[1]
	rng.b <- rng[2]
	t <- seq.int(rng.a, rng.b, rng.b/(p-1)) 
	evalX <- xTru
	
	fullInd <- c(1:p)
		
	# Hardy-Weinberg equilibrium
	# My original way of cutting xTru into SNPs
	# set continous evalX to 0,1,2 in "geno"
	# geno<- evalX*0
	p0=0.5^2
	p2=(1-0.5)^2
	p1=2*0.5*(1-0.5)
	geno1 <- XtoG(evalX,p0,p1,p2)
	
	# Way 2
	p00=0.4^2
	p22=(1-0.4)^2
	p11=2*0.4*(1-0.4)
	geno2 <- XtoG(evalX,p00,p11,p22)
	
	
	nsample <- nrow(xTru)
   	nsnp <- ncol(xTru)
   
	betabasis <- basisfd
    genobasis <- basisfd
    
    Phi <- eval.basis(t, genobasis)
    
    U1 <- geno1 %*% Phi %*% ginv(t(Phi) %*% Phi)
    U2 <- geno2 %*% Phi %*% ginv(t(Phi) %*% Phi)
		
	#return(data)	
	results <- list(xHat1=fd(coef=t(U1), basisobj=basisfd), geno1=geno1, xHat2=fd(coef=t(U2), basisobj=basisfd), geno2=geno2)
	return(results)		
}

###############################################################



###################### Transform X to G ######################
XtoG <- function(evalX,p0,p1,p2){
	geno <- matrix(1, dim(evalX)[1], dim(evalX)[2])
	for (j in 1:ncol(evalX))
	{
		genoInd_0 <- evalX[,j] <= quantile(evalX[,j], p0)
		geno[genoInd_0,j]=0
		genoInd_2 <- evalX[,j] > quantile(evalX[,j], p0+p1)
		geno[genoInd_2,j]=2
	}
	return(geno)
}
###############################################################



# rmultnorm<-function(n,muvec,sigmat,tSD){
#   # n is the number of random vectors to be generated
#   # mu is the mean vector, sigmat is the variance-covariance matrix
#   # the function returns an n by p matrix with each row being a
#   # random vector from a multivariate normal with mean vector muvec
#   # and covariance matrix sigmat
#   if(length(muvec)==1){
#     temp<-rnorm(n,muvec,sqrt(sigmat))
#     amat=NA
#     return(temp)
#   }
#   else{
#     sigeigen<-eigen(sigmat)
#     amat<-sigeigen$vectors%*%diag(sqrt(sigeigen$values))
#     temp<-matrix(rnorm(n*length(muvec),0,tSD),ncol=n)
#     temp<-(amat%*%temp)+muvec
#     temp<-t(temp)
#     return(temp)
#   }
#   return(temp, amat)
# }


# eigenMat<-function(sigmat){
#   # sigmat is the variance-covariance matrix
#   # the function return the eigen matrix from sigmat
#     sigeigen<-eigen(sigmat)
#     emat<-sigeigen$vectors%*%diag(sqrt(sigeigen$values))
#   return(emat)
# }


#################### Calc eigen matrix with corr ####################
corrEigenMat_Kvec <- function(p, Kvec, Kobs, Gsig, network, corr_interG, corr_interG_s1s2, corr_intraG, slen){
	
	kseq=c(0, Kvec)
	J=length(Kvec)
	kcumsum=c(0, cumsum(Kobs))
	
	covar <- matrix(rep(0,p*p),nrow=p)
	
	for(nwi in 1:length(network))
	{
		nwi_corrG = NULL
		nwi_corrG = network[[nwi]] #get the list of component per network
		nw_corrInd=NULL
		for(nwij in nwi_corrG)
		{
			ind.i <- c((kcumsum[nwij]+1):(kcumsum[nwij+1]))
			nw_corrInd <- c(nw_corrInd, ind.i)
		}
		# btw group correlation
		covar[nw_corrInd, nw_corrInd]=corr_intraG
	}
	
	# corrG = unlist(network)
	Gsig = Gsig
	# nullG = setdiff(corrG,signalG)
	
	# Within each group correlation	
	for(j in 1:J)
	{
		indj=c((kcumsum[j]+1):kcumsum[j+1])
		covar[indj,indj]= corr_interG
	}
	
	for(s in Gsig)
	{
		s1=c((kcumsum[s+1]-slen+1): kcumsum[s+1])
		s2=c((kcumsum[s]+1):(kcumsum[s+1]-slen))
		covar[s1,s2]= corr_interG_s1s2
		covar[s2,s1]= corr_interG_s1s2
	}
	

	# Set 1 to diagonal items
	diag(covar)=1
  	
  	# Set correlated matrix
  	# diagcovar <- diag(1,p,p)
  	
  	sigeigen<-eigen(covar)
    emat<-sigeigen$vectors%*%diag(sqrt(sigeigen$values))
    
    # way 2
    # X <- mvrnorm(n, mu=rep(0,p), Sigma=covar)
  	
  	res <- list(covar=covar, emat=emat, network=network, corr_interG=corr_interG, corr_intraG= corr_intraG)
  	return (res)
}
####################################################################





#################### Generate Data (single data) ####################
DataYBetaDisc <- function(n, p, rangeval, Beta.i, estX, eigenR, snr){
  # sigR50/75 use nnknots=201
  # sigR25 use nnknots=301
  
  # nnknots <- 301 # for sigR25
  nnknots <- 201 # for sigR50/75
	nnbasis <- nnknots + 5 - 2
	nnorder <- 5
	basisfd <- create.bspline.basis(rangeval=rangeval, nbasis=nnbasis, norder=nnorder)
	
	
	# Generate the vector of beta, p-dim
	tBeta <- seq.int(0, 1, 1/(p-1))
  	
	emat <- eigenR$emat
	covar <- eigenR$covar
	temp<-matrix(rnorm(n*p,0,1),ncol=n)
	xTru <- t(emat%*%temp) # without mean

	# Convert data to fdobj
	xfd=Data2fd(t(xTru), argvals=seq(0,1,len=p), basisobj=basisfd, lambda=0)
	
	#inner product of <beta#_j, B_j>
	G1 <- matrix(0,nrow=nnbasis, ncol=1)
	for (j in 1:nnbasis){
		G1[j,1] <- InnerProd(Beta.i,basisfd,j)
	} 	
	
	aMat <- xfd$coef 
	yTru <- t(aMat)%*%G1
	
	# Signal to noise
	ep <- rnorm(n=n, mean=0,sd=1)
	ye <- ep*sd(yTru)/snr
	# Generate y
	y <- yTru + ye
	
	XX_xTru <- EstX(xTru, p, basisfd)
	xHat1 <- XX_xTru$xHat1
	xGeno1 <- XX_xTru$geno1
	xHat2 <- XX_xTru$xHat2
	xGeno2 <- XX_xTru$geno2
	
	# return(data)	
	data <- list(y=y, yTru=yTru, xfd=xfd, xTru=xTru, xHat1=xHat1, xGeno1=xGeno1, xHat2=xHat2, xGeno2=xGeno2, nnbasis=nnbasis)
	return(data)		
}
###############################################################




########################## Replicates #########################
### Generate datasets over replicates                       ###
generate.datasets <- function(nruns, n, nTest, p, rangeval, Beta.i, estX, eigenR, snr)
{
	train <- list()
	test <- list()
	
	i <- 1
	while (i <= nruns)
	{
		print(sprintf('running=%d',i))
		train[[i]] <- DataYBetaDisc(n, p, rangeval, Beta.i, estX, eigenR, snr)
        if(nTest > 0) test[[i]] <- DataYBetaDisc(n, p, rangeval, Beta.i, estX, eigenR, snr) 

		i <- i + 1	
	}
		
	dat <- list(train=train, test=test)
}
###############################################################






Beta50g5_5sigR50 <- function(t) {
  beta <- rep(0,length(t))
  for(i in 1:length(t))
  {
  	# signalG = c(6:10)
  	# 50% sparsity for group 6 only
    if(0.18 < t[i] && t[i] <= 0.2) {beta[i] <- -2*(1-t[i])*sin(50*pi*(t[i]+0.2))}
    else if(0.11 < t[i] && t[i] <= 0.12) {beta[i] <- 0.5*sin(80*pi*(t[i])+0.6)} 
    else if(0.12 < t[i] && t[i] <= 0.14) {beta[i] <- 0.5*sin(30*pi*(t[i])-1.2)}
    else if(0.14 < t[i] && t[i] <= 0.16) {beta[i] <- -0.5*sin(30*pi*(t[i]))}
    else if(0.16 < t[i] && t[i] <= 0.18) {beta[i] <- -0.5*sin(40*pi*(t[i])-0.6)}
    else{beta[i] <- 0}}
  return(beta)
}
# t <- seq.int(0, 1, 1/(5000-1)) 
# bTru= Beta50g5_5sigR50(t)
# plot(t, bTru, col='blue',ylim=c(-0.6,2.5))


Beta50g5_6sigR50 <- function(t) {
  beta <- rep(0,length(t))
  for(i in 1:length(t))
  {
  	# signalG = c(6:10)
  	# 50% sparsity for group 6 only
    if(0.18 < t[i] && t[i] <= 0.2) {beta[i] <- -3*(1-t[i])*sin(50*pi*(t[i]+0.2))}
    else if(0.11 < t[i] && t[i] <= 0.12) {beta[i] <- 0.4*sin(80*pi*(t[i])+0.6)} 
    else if(0.12 < t[i] && t[i] <= 0.14) {beta[i] <- 0.4*sin(30*pi*(t[i])-1.2)}
    else if(0.14 < t[i] && t[i] <= 0.16) {beta[i] <- -0.4*sin(30*pi*(t[i]))}
    else if(0.16 < t[i] && t[i] <= 0.18) {beta[i] <- -0.4*sin(40*pi*(t[i])-0.6)}
    else{beta[i] <- 0}}
  return(beta)
}
# t <- seq.int(0, 1, 1/(5000-1)) 
# bTru= Beta50g5_6sigR50(t)
# plot(t, bTru, col='blue',ylim=c(-0.6,2.5))


Beta50g5_7sigR50 <- function(t) {
  beta <- rep(0,length(t))
  for(i in 1:length(t))
  {
  	# signalG = c(6:10)
  	# 50% sparsity for group 6 only
    if(0.18 < t[i] && t[i] <= 0.2) {beta[i] <- -2*(1-t[i])*sin(50*pi*(t[i]+0.2))}
    else if(0.11 < t[i] && t[i] <= 0.12) {beta[i] <- 0.4*sin(80*pi*(t[i])+0.7)} 
    else if(0.12 < t[i] && t[i] <= 0.14) {beta[i] <- -0.4*sin(35*pi*(t[i])-0.1)}
    else if(0.14 < t[i] && t[i] <= 0.16) {beta[i] <- 0.6*sin(35*pi*(t[i])+0.6)}
    else if(0.16 < t[i] && t[i] <= 0.18) {beta[i] <- -0.6*sin(40*pi*(t[i])-0.6)}
    else{beta[i] <- 0}}
  return(beta)
}
# t <- seq.int(0, 1, 1/(5000-1)) 
# bTru= Beta50g5_7sigR50(t)
# plot(t, bTru, col='blue',ylim=c(-0.6,2.5))


Beta50g5_8sigR50 <- function(t) {
  beta <- rep(0,length(t))
  for(i in 1:length(t))
  {
  	# signalG = c(6:10)
  	# 50% sparsity for group 6 only
    if(0.18 < t[i] && t[i] <= 0.2) {beta[i] <- -2*(1-t[i])*sin(50*pi*(t[i]+0.2))}
    else if(0.11 < t[i] && t[i] <= 0.12) {beta[i] <- 0.2*sin(80*pi*(t[i])+0.7)} 
    else if(0.12 < t[i] && t[i] <= 0.14) {beta[i] <- -0.2*sin(35*pi*(t[i])-0.1)}
    else if(0.14 < t[i] && t[i] <= 0.16) {beta[i] <- 0.6*sin(31*pi*(t[i])+2.2)}
    else if(0.16 < t[i] && t[i] <= 0.18) {beta[i] <- 0.6*sin(34*pi*(t[i])-0.4)}
    else{beta[i] <- 0}}
  return(beta)
}
# t <- seq.int(0, 1, 1/(5000-1)) 
# bTru= Beta50g5_8sigR50(t)
# plot(t, bTru, col='blue',ylim=c(-0.6,2.5))



Beta50g5_9sigR50 <- function(t) {
  beta <- rep(0,length(t))
  for(i in 1:length(t))
  {
  	# signalG = c(6:10)
  	# 50% sparsity for group 6 only
    if(0.18 < t[i] && t[i] <= 0.2) {beta[i] <- -2*(1-t[i])*sin(50*pi*(t[i]+0.2))}
    else if(0.11 < t[i] && t[i] <= 0.12) {beta[i] <- 0.2*sin(80*pi*(t[i])+0.7)} 
    else if(0.12 < t[i] && t[i] <= 0.14) {beta[i] <- -0.2*sin(35*pi*(t[i])-0.1)}
    else if(0.14 < t[i] && t[i] <= 0.16) {beta[i] <- 0.4*sin(32*pi*(t[i])+1.7)}
    else if(0.16 < t[i] && t[i] <= 0.18) {beta[i] <- 0.4*sin(33*pi*(t[i])+0.2)}
    else{beta[i] <- 0}}
  return(beta)
}
# t <- seq.int(0, 1, 1/(5000-1)) 
# bTru= Beta50g5_9sigR50(t)
# plot(t, bTru, col='blue',ylim=c(-0.6,2.5))






