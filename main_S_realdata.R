set.seed(2)
source("functions_S.R")
usePackage("splines2")
usePackage("TruncatedNormal")
usePackage("mvtnorm")
usePackage("matrixStats")
source("exchangable_cov.R")
usePackage("foreach")

library(doParallel)
library(doRNG)
registerDoParallel(cores=11)
options(warn=0)

# Loading Data ---------------------------------------------------
load('biocard.rda')

K <- 11 # Number of biomarkers
# All non-age covariates
X <- as.matrix(cbind(intercept=1,df[,c("apoe","SEX","education")])) 
Y <- df[,1:K] # Biomarkers array
t <- df$ageori # Age in original scale
dfi <- 20 # DoF of Spline
qknot <- (1:(dfi-3))/(dfi-2) # Quantiles to determine knots
VIF <- 0.1 # Variance inflation factor for BETAKDE
group <- 1:11

library(ggplot2)

# Construct biomarker-specific design matrix
covar.list <- as.list(rep(NA,K))
knot.list <- as.list(rep(NA,K))
boundary.knot <- c(0,120)#range(t)
#remove <- 2# Removing the first .. and last .. basis
for(i in 1:K){
  # Calculate knot points
  t01 <- (t-boundary.knot[1])/(boundary.knot[2]-boundary.knot[1])
  t01 <- t01[dplyr::between(t01,0,1)]
  knot.list[[i]] <- betaKDE(t01,s=VIF,q=qknot)$quantile
  knot.list[[i]] <- (boundary.knot[2]-boundary.knot[1])*knot.list[[i]]+boundary.knot[1]
  B <- ibs(pmin(pmax(t,min(boundary.knot)),max(boundary.knot)), 
           knots=knot.list[[i]], Boundary.knots = boundary.knot, 
           degree=2, intercept=TRUE) # IBSpline Basis
  #B <- iSpline(t, knots=knot.list[[i]], Boundary.knots = boundary.knot, 
  #             degree=2, intercept=TRUE) # IBSpline Basis
  #B <- B/min(apply(B,2,max))
  B <- B[,(3):(ncol(B)-2)]
  covar.list[[i]] <- cbind(X,B)
}


# Create consecutive pseudo-IDs for each individual for easy coding
unique.IDs <- sort(unique(df$Study.ID))
df$ID <- match(df$Study.ID,unique.IDs)
# Pre-calculate longitudinal sample size 
# for each individual-biomarker combination
long_ss <- matrix(0,nrow=length(unique.IDs),ncol=K)
for(i in 1:length(unique.IDs))
  for(j in 1:K){
    long_ss[i,j] <- sum((df$ID==i)&(!is.na(Y[,j])))
  }

long_all_ss <- rep(0,length(unique.IDs))
for(i in seq_along(unique.IDs)){
  long_all_ss[i] <- sum(df$ID==i)
}

R <- 1e4 # Set Number of Iterations
Burnin <- R/2 # Set Number of Burn-ins

library(progress)
pb <- progress_bar$new(total = R-1)

# Set Priors -----------------------------------------------------
# Beta parameter: Coefficients for adjusting covariates
beta.prior <- list(mean=rep(0,ncol(X)),
                   variance=diag(rep(10000,ncol(X))),
                   #variance=10000,
                   precision=NULL)
beta.prior$precision <- solve(beta.prior$variance)
# Gamma parameter: Coefficients for splines
gamma.prior <- list(mean=rep(0,ncol(B)),
                    variance=NULL)
coef.prior <- list(mean=c(beta.prior$mean, gamma.prior$mean),
                   variance=NULL,
                   precision=NULL)
#coef.prior$precision <- solve(coef.prior$variance)


# Set initial guess ----------------------------------------------

# Fixed Effect of X & All-positive Spline Coefs
coefs <- array(0,c(ncol(covar.list[[1]]),ncol(Y),R)) 
nX <- ncol(X)
# Variance of Biomarkers
sigmays <- rep(0,R)
sigmaws <- rep(0,R)
pens <- array(0, c(2, R))
REs <- array(0, c(dim(long_ss),R))
offsets <- array(0, c(dim(long_ss),R))

coefs[1:nX,,1] <- 
  t(rtmvnorm(ncol(Y),mu=beta.prior$mean,
             sigma=beta.prior$variance))
coefs[(nX+1):ncol(covar.list[[1]]),,1] <- 
  t(rtmvnorm(ncol(Y),mu=gamma.prior$mean,
             sigma=penalty_Matrix(ncol(B),
                                  smooth.sigma = 1,flat.sigma = 1)$V,
             lb=rep(0,ncol(B))))
sigmays[1] <- 1/rgamma(1,shape=3,rate=0.5)
pens[,1] <- c(0.1,0.5)
sigmaws[1] <- 1/rgamma(1,shape=3,rate=0.5)

# Prior density for penalties
lpd <- function(s) 
  log(2)+dnorm(s[1], sd=1/50, log=TRUE) +
  log(2)+dnorm(s[2], sd=1/50, log=TRUE)
ls <- -2 # Log of Jump Standard Deviation
acc <- 0 # Accepted Proposals in one batch
lss <- ls # Sequence of LS for reference
w <- planck_taper(ncol(B), eps=0.1) # Window Function
#w <- rep(1,ncol(B))
#w <- NULL

M_coef <- lincon(nX+dfi-4,nX)
M_pen <- lincon(dfi-4,0)
verbose=FALSE
# Perform MCMC ----------------------------------------------------
for(i in 1:(R-1)){
  prec <- block_Matrix(beta.prior$precision,
                       penalty_Matrix(ncol(B),pens[1,i],pens[2,i],
                                      weight=w)$prec)
  sigmays[i+1] <- update_sigmay(covar.list,Y,as.matrix(REs[df$ID,,i],ncol=1),
  as.matrix(coefs[,,i],ncol=1),
  3,0.5)
  u <- update_coef_grouped(covar.list,nX,Y,as.matrix(REs[df$ID,,i],ncol=1),
                   sigmays[i+1],sigmaws[i],df$ID,
                   coef.prior$mean,
                   prec,M_coef,verbose,group,samples=1)
  coefs[,,i+1] <- aperm(u$res,c(2,3,1))
  REs[,,i+1] <- update_W(covar.list,Y,as.matrix(coefs[,,i+1],ncol=K),long_ss,
                         df$ID,sigmays[i+1],sigmaws[i])
  new_pens <- update_pens_grouped(gamma=as.matrix(coefs[(nX+1):ncol(covar.list[[1]]),,i+1],ncol=1),
                          mu=gamma.prior$mean,
                          lambda=pens[,i],
                          lpd=lpd,
                          ls=ls,
                          weight=w,
                          Ms=M_pen,
                          verbose=verbose,
                          group=group
  )
  pens[,i+1] <- new_pens$new
  acc <- acc + new_pens$acc_status
  if(i%%50==0){
    delta <- min(0.1, 1/sqrt(i/50))
    rate <- acc/50
    if(rate >= 0.234){
      ls <- ls + delta
    }
    else{
      ls <- ls - delta
    }
    acc <- 0
    lss <- c(lss, ls)
  }
  sigmaws[i+1] <- update_sigmaw(REs[,,i+1],3,0.5)
  pb$tick()
}
stopImplicitCluster()
save.image('biocard_result_group16nonzeros_ungrouped.RData')