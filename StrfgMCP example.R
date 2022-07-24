
library(fda)
library(pracma)
library(mvtnorm)
library(MASS)
library(matrixcalc)
library(FactoMineR)

source("Fun_StrfgMCP.R")
source("Fun_generate_data.R")


			#################################################
			#########   Step 1: Generate data   #############
			
						######    n=150     ######
						######     J=50     ######
						######   Eigen946   ######

rangeval=c(0,1) #domain of the function
# Kvec=c(3,5,3,5,2,6,2,6,4,4,4,4,3,5,3,5,2,6,2,6, rep(4,30)) 
# vector of 50 groups of genes with varous group sizes
Kvec=rep(4,50) # vector of 50 groups of genes with equal size (4)
Kobs=Kvec*25 # set 25 snps per gene
J=length(Kvec) # no. of groups
p=sum(Kobs) # total no. of snp

network=list(nw1=c(1:5), nw2=c(6:10), nw3=c(11:15), nw4=c(16:20), nw5=c(21:25), nw6=c(26:30), nw7=c(31:35), nw8=c(36:40), nw9=c(41:45), nw10=c(46:50))
Gsig=6

n=150 #sample size
eigenR=NULL
system.time(eigenR <- corrEigenMat_Kvec(p, Kvec, Kobs, Gsig, network, corr_interG=0.9, corr_interG_s1s2=0.4, corr_intraG=0.6, slen=50))
# saveRDS(eigenR, file="G50_nw2_sigR50_eigenR946")
# eigenR=readRDS(file="G50_nw2_sigR50_eigenR946")
###
data=NULL
set.seed(666)
Beta.i=Beta50g5_5sigR50
data=generate.datasets(nruns=5, n=n, nTest=n, p=p, rangeval, Beta.i, estX=T, eigenR, snr=10)
# saveRDS(data, file="dataG50_Beta50g5_5sigR50_eigenR946_snr10_n150")

                      



			#################################################
			#########     Step 2: Estimation    #############

a=3
POpar=2
max.iter=500
changeThres <- 1e-4
cutoff <- 1e-4
estX=T
xcoef.mod=F
int=F
subint=25 # 25 obs/subint
JK=p/subint # no. of subintervals = M = JK

### (1) Estimate a single dataset without tuning ###
irun=1
xHat <- data$train[[irun]]$xHat1 #MAF0.5; xHat2 for MAF0.4;
xfd <- data$train[[irun]]$xfd
xTru <- data$train[[irun]]$xTru
xGeno <- data$train[[irun]]$xGeno1 #MAF0.5; xGeno2 for MAF0.4;
y <- data$train[[irun]]$y

nknots <- JK+1 #nknots=M+1
norder <- 5 #norder=d+1
d=norder-1
nb <- nknots + norder - 2
basisfd <- create.bspline.basis(rangeval= rangeval, nbasis=nb, norder=norder)
basisfd.beta <- basisfd
Jmat <- inprod(xfd$basis, basisfd.beta)
int <- F
estX <- T

Ws <- NULL
if(is.null(Ws)) Ws <- computeW(basisfd.beta)  
Ajj <- RVcoef_xHat(J,JK,Kvec,xHat)
diag(Ajj)=0

gamma=6e-14
lambda1=0.1
lambda2=0.001
fit=method.alg(xfd, xHat, y, J, Kvec, Ajj, lambda1, lambda2, gamma, POpar, a, int, basisfd.beta, Jmat, Ws, max.iter, changeThres, cutoff, estX)
beta=fit$beta
tpfp.d=dTPFP(beta, Beta.i, dimP=5000)
tpfp.d
plot(fit$beta, ylim=c(-2,2)) # estimated curve
t <- seq.int(0, 1, 1/(5000-1)) 
bTru=Beta.i(t) # true curve
lines(t, bTru, col='blue')


### (2) Estimate a single dataset with tunings  ###
gamma=c(5.5e-14,6e-14,6.5e-14) 
lambda1=c(0.06,0.08,0.1,0.12)
lambda2=seq(from=0.0001,to=0.02,by=0.005)
tunedRes <- tuneProgram(tuning="AIC/BIC/RIC", xfd, xHat, y, J, Kvec, Ajj, lambda1, lambda2, gamma, POpar, a, int, basisfd.beta, Jmat, Ws, max.iter, changeThres, cutoff, estX)

gam=as.numeric(tunedRes$selMat[2,'gamma'])
lam1=as.numeric(tunedRes$selMat[2,'lambda1'])
lam2=as.numeric(tunedRes$selMat[2,'lambda2'])
fit=method.alg(xfd, xHat, y, J, Kvec, Ajj, lam1, lam2, gam, POpar, a, int, basisfd.beta, Jmat, Ws, max.iter, changeThres, cutoff, estX)
beta=fit$beta
tpfp.d=dTPFP(beta, Beta.i, dimP=5000)
tpfp.d
plot(fit$beta, ylim=c(-2,2)) # estimated curve
t <- seq.int(0, 1, 1/(5000-1)) 
bTru=Beta.i(t) # true curve
lines(t, bTru, col='blue')



### (3) Averaged estimates over nruns replicates ###
res=NULL
res = estimate.datasets(nruns=5, rangeval, data, lambda1, lambda2, gamma, Kvec, J, POpar, a, max.iter, changeThres, cutoff, int, estX, Beta.i, MAF.text="MAF0.5")
res$eval.ave
res$eval.sd


