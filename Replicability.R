#packages required
library(mclust)
library(mvtnorm)
library(mvnfast)
library(MixGHD)
library(mclustAddons)
library(fpc)
library(igraph)
library(compositions)

# SET HERE YOUR WORKING DIRECTORY
#setwd()

#Function required
source("GMM_uni.R")
source("Hellinger_mixture.R")

# ----------------------------
# Motivating example 1
# ----------------------------

# data
data(Baudry_etal_2010_JCGS_examples)

# GMM
set.seed(123)
mod = Mclust(ex4.1)

# clustHel
hel = Hellingher_algorithm(mod, 0.92)

# H matrix
hel$H

# merging components
hel$merg

# plot
par(mfrow=c(1,2))
z = mod$classification
plot(ex4.1, col = z, type = "n")
text(ex4.1[z==1,],label=1,col=1)
text(ex4.1[z==2,],label=2,col=2)
text(ex4.1[z==3,],label=3,col=3)
text(ex4.1[z==4,],label=4,col=4)
text(ex4.1[z==5,],label=5,col=5)
text(ex4.1[z==6,],label=6,col=6)

z = ifelse(z==1, 6, z)
z = ifelse(z==3, 4, z)
plot(ex4.1, col = z)

# ----------------------------
# Motivating example 2
# ----------------------------

# data
n = 150
Sigma = matrix(c(0.2,0.1,0.1,0.2),byrow=TRUE,nrow=2)
Mu = c(1,1)
set.seed(123)
x = rlnorm.rplus(n,log(Mu),Sigma)
x = matrix(x,dim(x)[1],dim(x)[2])

# GMM
mod = Mclust(x)

# clustHel
hel = Hellingher_algorithm(mod, 0.92)
hel$H
hel$merg

# new clas
classification = ifelse(mod$classification==1, 1, 1)

# plot
par(mfrow = c(1, 2), mar = c(5.1, 6.1, 4.1, 1.1))
plot(x, type = "n", col = mod$classification, xlab = "X1", ylab = "X2", cex.axis=1, cex.lab = 1)
text(x[mod$classification==1,],label=1,col=1)
text(x[mod$classification==2,],label=2,col=2)
plot(x, col = classification, xlab = "X1", ylab = "X2", cex.axis=1, cex.lab = 1)

# -----------------------------
# Motivating example 3
# -----------------------------
mu1 = c(1,1)
mu2 = c(3,1)
sigma1 = matrix(0.6,2,2)
diag(sigma1) = 1
sigma2 = matrix(-0.5,2,2)
diag(sigma2) = 1

set.seed(123)
x1 = rmvnorm(50, mu1, sigma1)
x2 = rmvnorm(75, mu2, sigma2)

data = rbind(x1,x2)
plot(data)

# GMM
mod = Mclust(data)

# PMLE
mod_pmle = GMMuni(data, 100, 2, 1, 1, model = "yes", mod, force_est = 0, den = FALSE)

# plot
xx <- seq(-1, 6, length.out = 100)  
yy <- seq(-2, 5, length.out = 100)  

grid <- expand.grid(x = xx, y = yy)
pdf <- rep(0,10000)
for(k in 1:2){
  pdf <- pdf + mod_pmle$PiMclust[k]*dmvnorm(grid, mod_pmle$MuMclust[,k], mod_pmle$SigmaMclust[,,k])
}
z <- matrix(pdf, 100, 100)

pdf <- rep(0,10000)
for(k in 1:2){
  pdf <- pdf + mod_pmle$pis[k]*dmvnorm(grid, mod_pmle$Mus[,k], mod_pmle$Sigmas[,,k])
}
z2 <- matrix(pdf, 100, 100)

par(mfrow=c(1,2))
plot(mod_pmle$x)
contour(xx, yy, z, main = "Contour Plot no-unimodality", xlab = "X1", ylab = "X2", add = TRUE)
plot(mod_pmle$x)
contour(xx, yy, z2, main = "Contour Plot unimodality", xlab = "X1", ylab = "X2", add = TRUE)

# -----------------------------
# clustHel Simulation
# -----------------------------

hellinger.dist.mod = function(mu1, mu2, sigma1, sigma2, w1, w2){
  hd = sqrt(w1+w2-(2*sqrt(w1*w2))*exp(-bhattacharyya.dist(mu1, mu2, sigma1, sigma2)))
  return(as.numeric(hd))
}

# dimensionalty
n = 100

# B
nrep = 100

# epsilon grid
n_epsilon = seq(0.50, 0.99, 0.01)
res = matrix(NA, nrep, length(n_epsilon))

hellinger.dist.mod = function(mu1, mu2, sigma1, sigma2, w1, w2){
  hd = sqrt(w1+w2-(2*sqrt(w1*w2))*exp(-bhattacharyya.dist(mu1, mu2, sigma1, sigma2)))
  return(as.numeric(hd))
}

try(for(sim in 1:nrep){
  set.seed(122 + sim)
  x1 = rmvnorm(n/2, c(2,1), diag(1, 2)) #data x1
  x2 = rmvnorm(n/2, c(4,5), matrix(c(1, 0.6, 0.6, 1), 2, 2)) #data x2
  x = rbind(x1,x2)
  
  mod = Mclust(x,3) # overfitted GMM
  G = mod$G
  H = matrix(0, G, G)
  
  for(i in 1:(G-1)){
    for(j in (i+1):G){
      mu1 = mod$parameters$mean[,i]
      mu2 = mod$parameters$mean[,j]
      sigma1 = mod$parameters$variance$sigma[,,i]
      sigma2 = mod$parameters$variance$sigma[,,j]
      w1 = mod$parameters$pro[i]/(mod$parameters$pro[i] + mod$parameters$pro[j])
      w2 = mod$parameters$pro[j]/(mod$parameters$pro[i] + mod$parameters$pro[j])
      H[i,j] = hellinger.dist.mod(mu1, mu2, sigma1, sigma2, w1, w2)
    }
  }
  
  for(ep in 1:length(n_epsilon)){
    # algorithm
    IH = 1*(H > 0 & H <= n_epsilon[ep])
    merg = matrix(0, G, G)
    V = 1:G
    E = as.vector(t(which(IH == 1, arr.ind = TRUE)))
    if(length(E)==0){
      val = apply(merg, 2, function(x) sum(x!=0))-1
      res[sim,ep] = G-sum(val[val>0])
    }else{
      H_ = H+t(H)
      Gr = graph(E, directed = FALSE)
      C = largest_cliques(Gr)
      l = 0
      dimC = length(C)
      dimCi = length(C[[1]])
      
      if(dimC == 1){
        l = l+1
        merg[1:dimCi,l] = C[[1]]
      }else{
        while(dimC > 1){
          l = l+1
          Gr = graph(E, directed = FALSE)
          C = largest_cliques(Gr)
          dimC = length(C)
          dimCi = length(C[[1]])
          
          if(dimC > 1){
            d = array(0, dim = c(dimCi, dimCi, dimC))
            for(i in 1:dimC){
              for(j in 1:(dimCi-1)){
                for(k in (j+1):dimCi){
                  d[j,k,i] = H_[C[[i]][j],C[[i]][k]]
                }
              }
            }
            dmax = apply(d, 3, max)
            dmin = which.min(dmax)
            merg[1:dimCi,l] = C[[dmin]]
            A = matrix(E, ncol = 2, byrow = TRUE)
            E = as.vector(t(A[c(which(matrixStats::rowProds(t(apply(A, 1, function(x) !as.vector(C[[dmin]]) %in% x))) == 1)),]))
            rm(dmin)
            if(length(E) == 0){
              break
            }
          }else{
            merg[1:dimCi,l] = C[[1]]
            break
          }
        } 
      }
    }
    val = apply(merg, 2, function(x) sum(x!=0))-1
    res[sim,ep] = G-sum(val[val>0])
  }
})

corren100t2g3 = rep(NA, 50)
for(i in 1:50){corren100t2g3[i] = sum(res[,i] == 2)/length(res[,i])}

n_epsilon[which.max(corren100t2g3)]

# -----------------------------
# Simulations GMMuni
# -----------------------------

# simulation setting G=2
nrep = 1 # B
sampleSize = 100 # sample size
G = 2
pi = c(0.5, 0.7) # mixing proportion
mu = matrix(c(0, 1, 2, 5), 2, 2) # mu
sigma = array(c(1, 0.7, 0.7, 1, 1, -0.6, -0.6, 1), dim = c(2, 2, 2)) #Sigma
delta = 1 #tuning parameter

simpar = expand.grid(sim =  1:nrep,
                     sampleSize = sampleSize, G = G, pi = pi, delta = delta)

seed = seq(133, 132+nrep, by = 1)
simpar = cbind(simpar, seed)

sim = list()

for(i in 1:nrow(simpar)){
  set.seed(simpar$seed[i])
  x1 = rmvnorm(simpar$sampleSize[i]*simpar$pi[i], mu[,1], sigma[,,1])
  x2 = rmvnorm(simpar$sampleSize[i]*(1-simpar$pi[i]), mu[,2], sigma[,,2])
  x = rbind(x1,x2)
  sim[[i]] = try(GMMuni(x, 100, simpar$G[i], simpar$delta[i], 1, model = "no", force_est = 0, den = FALSE))
}

# -----------------------------
# Simulation informative features, X3 noise, partial code provided by clustvarsel
# -----------------------------

# simulated data
n = 200      # sample size
pro = 0.5    # mixing proportion
mu1 = c(0,0) # mean vector for the first cluster
mu2 = c(3,3) # mean vector for the second cluster
sigma1 = matrix(c(1,0.5,0.5,1),2,2)       # covar matrix for the first cluster
sigma2 = matrix(c(1.5,-0.7,-0.7,1.5),2,2) # covar matrix for the second cluster
X = matrix(0, n, 3, dimnames = list(NULL, paste0("X", 1:3)))
set.seed(123) 
u = runif(n)
Class = ifelse(u < pro, 1, 2)
X[u < pro, 1:2]  = MASS::mvrnorm(sum(u < pro), mu = mu1, Sigma = sigma1)
X[u >= pro, 1:2] = MASS::mvrnorm(sum(u >= pro), mu = mu2, Sigma = sigma2)
X[, 3] = rnorm(n, mean = 1.5, sd = 2)

# overfitted GMM
mod = Mclust(X, G = 3)

# clustHel
hel = Hellingher_algorithm(mod, 0.94)

# merging component and classification
merg = hel$merg
mergin_component = X[mod$classification==merg[1,1] | mod$classification==merg[2,1],]
classification = ifelse(mod$classification==2, 2, 1)

# plot
d = dim(X)[2]
par(mfrow = c(d, d), mar = rep(0.2/2, 4), oma = rep(3, 4))
plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
text(1, 1, "X1", cex = 2)
plot(X[,1],X[,2], col = mod$classification, xaxt="n", yaxt="n", ann=FALSE)
axis(side = 3, las = 1, cex.axis=2, cex.lab = 2.5) 
plot(X[,1],X[,3], col = mod$classification, xaxt="n", yaxt="n", ann=FALSE)
axis(side = 4, las = 2, cex.axis=2, cex.lab = 2.) 
plot(X[,2],X[,1], col = mod$classification, xaxt="n", yaxt="n", ann=FALSE)
axis(side = 2, las = 2, cex.axis=2, cex.lab = 2.) 
plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
text(1, 1, "X2", cex = 2, cex.axis=2, cex.lab = 2.)
plot(X[,2],X[,3], col = mod$classification, xaxt="n", yaxt="n", ann=FALSE)
plot(X[,3],X[,1], col = mod$classification, xaxt="n", yaxt="n", ann=FALSE)
axis(side = 1, las = 1, cex.axis=2, cex.lab = 2.) 
plot(X[,3],X[,2], col = mod$classification, xaxt="n", yaxt="n", ann=FALSE)
plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
text(1, 1, "X3", cex = 2)

d = dim(X)[2]
par(mfrow = c(d, d), mar = rep(0.2/2, 4), oma = rep(3, 4))
plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
text(1, 1, "X1", cex = 2)
plot(X[,1],X[,2], col = classification, xaxt="n", yaxt="n", ann=FALSE)
axis(side = 3, las = 1, cex.axis=2, cex.lab = 2.) 
plot(X[,1],X[,3], col = classification, xaxt="n", yaxt="n", ann=FALSE)
axis(side = 4, las = 2, cex.axis=2, cex.lab = 2.) 
plot(X[,2],X[,1], col = classification, xaxt="n", yaxt="n", ann=FALSE)
axis(side = 2, las = 2, cex.axis=2, cex.lab = 2.) 
plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
text(1, 1, "X2", cex = 2, cex.axis=2, cex.lab = 2.)
plot(X[,2],X[,3], col = classification, xaxt="n", yaxt="n", ann=FALSE)
plot(X[,3],X[,1], col = classification, xaxt="n", yaxt="n", ann=FALSE)
axis(side = 1, las = 1, cex.axis=2, cex.lab = 2.) 
plot(X[,3],X[,2], col = classification, xaxt="n", yaxt="n", ann=FALSE)
plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
text(1, 1, "X3", cex = 2)

# save original model
mod_original = mod

# subset of GMM
mod$parameters$pro = mod$parameters$pro[-2]/sum(mod$parameters$pro[-2])
mod$parameters$mean = mod$parameters$mean[,-2]
mod$parameters$variance$sigma = mod$parameters$variance$sigma[,,-2]
mod$n = nrow(mergin_component)

# PMLE
mergin_component = matrix(unlist(mergin_component), mod$n, 3)
newmod = GMMuni(mergin_component, 100, 2, 0.5, 1, "yes", mod, force_est = 0, den = FALSE)

# -----------------------------
# Simulation informative features, X3 correlated, partial code provided by clustvarsel
# -----------------------------

n = 200      # sample size
pro = 0.5    # mixing proportion
mu1 = c(0,0) # mean vector for the first cluster
mu2 = c(3,3) # mean vector for the second cluster
sigma1 = matrix(c(1,0.5,0.5,1),2,2)       # covar matrix for the first cluster
sigma2 = matrix(c(1.5,-0.7,-0.7,1.5),2,2) # covar matrix for the second cluster
X = matrix(0, n, 3, dimnames = list(NULL, paste0("X", 1:3)))
set.seed(123) 
u = runif(n)
Class = ifelse(u < pro, 1, 2)
X[u < pro, 1:2]  = MASS::mvrnorm(sum(u < pro), mu = mu1, Sigma = sigma1)
X[u >= pro, 1:2] = MASS::mvrnorm(sum(u >= pro), mu = mu2, Sigma = sigma2)
X[, 3] = X[, 1] + rnorm(n)

# overfitted GMM
mod = Mclust(X, G = 3)

# clustHel
hel = Hellingher_algorithm(mod, 0.94)

# merging component and classification
merg = hel$merg
mergin_component = X[mod$classification==merg[1,1] | mod$classification==merg[2,1],]
classification = ifelse(mod$classification==1, 1, 2)

# plot
d = dim(X)[2]
par(mfrow = c(d, d), mar = rep(0.2/2, 4), oma = rep(3, 4))
plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
text(1, 1, "X1", cex = 2)
plot(X[,1],X[,2], col = mod$classification, xaxt="n", yaxt="n", ann=FALSE)
axis(side = 3, las = 1, cex.axis=2, cex.lab = 2) 
plot(X[,1],X[,3], col = mod$classification, xaxt="n", yaxt="n", ann=FALSE)
axis(side = 4, las = 2, cex.axis=2, cex.lab = 2) 
plot(X[,2],X[,1], col = mod$classification, xaxt="n", yaxt="n", ann=FALSE)
axis(side = 2, las = 2, cex.axis=2, cex.lab = 2) 
plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
text(1, 1, "X2", cex = 2)
plot(X[,2],X[,3], col = mod$classification, xaxt="n", yaxt="n", ann=FALSE)
plot(X[,3],X[,1], col = mod$classification, xaxt="n", yaxt="n", ann=FALSE)
axis(side = 1, las = 1, cex.axis=2, cex.lab = 2) 
plot(X[,3],X[,2], col = mod$classification, xaxt="n", yaxt="n", ann=FALSE)
plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
text(1, 1, "X3", cex = 2)

d = dim(X)[2]
par(mfrow = c(d, d), mar = rep(0.2/2, 4), oma = rep(3, 4))
plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
text(1, 1, "X1", cex = 2)
plot(X[,1],X[,2], col = classification, xaxt="n", yaxt="n", ann=FALSE)
axis(side = 3, las = 1, cex.axis=2, cex.lab = 2) 
plot(X[,1],X[,3], col = classification, xaxt="n", yaxt="n", ann=FALSE)
axis(side = 4, las = 2, cex.axis=2, cex.lab = 2) 
plot(X[,2],X[,1], col = classification, xaxt="n", yaxt="n", ann=FALSE)
axis(side = 2, las = 2, cex.axis=2, cex.lab = 2) 
plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
text(1, 1, "X2", cex = 2, cex.axis=2, cex.lab = 2)
plot(X[,2],X[,3], col = classification, xaxt="n", yaxt="n", ann=FALSE)
plot(X[,3],X[,1], col = classification, xaxt="n", yaxt="n", ann=FALSE)
axis(side = 1, las = 1, cex.axis=2, cex.lab = 2) 
plot(X[,3],X[,2], col = classification, xaxt="n", yaxt="n", ann=FALSE)
plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
text(1, 1, "X3", cex = 2)

# save original model
mod_original = mod

# subset of GMM
mod$parameters$pro = mod$parameters$pro[-1]/sum(mod$parameters$pro[-1])
mod$parameters$mean = mod$parameters$mean[,-1]
mod$parameters$variance$sigma = mod$parameters$variance$sigma[,,-1]
mod$n = nrow(mergin_component)

# PMLE
mergin_component = matrix(unlist(mergin_component), mod$n, 3)
newmod = GMMuni(mergin_component, 100, 2, 0.1, 1, "yes", mod, force_est = 0, den = FALSE)

# -----------------------------
# Simulation informative features, X3 noise,  X4 correlated, partial code provided by clustvarsel
# -----------------------------

n = 200      # sample size
pro = 0.5    # mixing proportion
mu1 = c(0,0) # mean vector for the first cluster
mu2 = c(3,3) # mean vector for the second cluster
sigma1 = matrix(c(1,0.5,0.5,1),2,2)       # covar matrix for the first cluster
sigma2 = matrix(c(1.5,-0.7,-0.7,1.5),2,2) # covar matrix for the second cluster
X = matrix(0, n, 4, dimnames = list(NULL, paste0("X", 1:4)))
set.seed(123) 
u = runif(n)
Class = ifelse(u < pro, 1, 2)
X[u < pro, 1:2]  = MASS::mvrnorm(sum(u < pro), mu = mu1, Sigma = sigma1)
X[u >= pro, 1:2] = MASS::mvrnorm(sum(u >= pro), mu = mu2, Sigma = sigma2)
X[, 3] = X[, 1] + rnorm(n)
X[, 4] = rnorm(n, mean = 1.5, sd = 2)

# GMM
mod = Mclust(X, G = 3)

# clustHel
hel = Hellingher_algorithm(mod, 0.96)

# merging component and classification
merg = hel$merg
mergin_component = X[mod$classification==merg[1,1] | mod$classification==merg[2,1],]
classification = ifelse(mod$classification==1, 1, 2)

# plot
d = dim(X)[2]
par(mfrow = c(d, d), mar = rep(0.2/2, 4), oma = rep(3, 4))
plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
text(1, 1, "X1", cex = 2)
plot(X[,1],X[,2], col = mod$classification, xaxt="n", yaxt="n", ann=FALSE)
axis(side = 3, las = 1, cex.axis=2, cex.lab = 2) 
plot(X[,1],X[,3], col = mod$classification, xaxt="n", yaxt="n", ann=FALSE)
axis(side = 3, las = 1, cex.axis=2, cex.lab = 2) 
plot(X[,1],X[,4], col = mod$classification, xaxt="n", yaxt="n", ann=FALSE)
axis(side = 4, las = 2, cex.axis=2, cex.lab = 2) 
plot(X[,2],X[,1], col = mod$classification, xaxt="n", yaxt="n", ann=FALSE)
axis(side = 2, las = 2, cex.axis=2, cex.lab = 2) 
plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
text(1, 1, "X2", cex = 2)
plot(X[,2],X[,3], col = mod$classification, xaxt="n", yaxt="n", ann=FALSE)
plot(X[,2],X[,4], col = mod$classification, xaxt="n", yaxt="n", ann=FALSE)
plot(X[,3],X[,1], col = mod$classification, xaxt="n", yaxt="n", ann=FALSE)
axis(side = 2, las = 2, cex.axis=2, cex.lab = 2) 
plot(X[,3],X[,2], col = mod$classification, xaxt="n", yaxt="n", ann=FALSE)
plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
text(1, 1, "X3", cex = 2)
plot(X[,3],X[,4], col = mod$classification, xaxt="n", yaxt="n", ann=FALSE)
plot(X[,4],X[,1], col = mod$classification, xaxt="n", yaxt="n", ann=FALSE)
axis(side = 1, las = 1, cex.axis=2, cex.lab = 2) 
plot(X[,4],X[,2], col = mod$classification, xaxt="n", yaxt="n", ann=FALSE)
plot(X[,4],X[,3], col = mod$classification, xaxt="n", yaxt="n", ann=FALSE)
plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
text(1, 1, "X4", cex = 2)

d = dim(X)[2]
par(mfrow = c(d, d), mar = rep(0.2/2, 4), oma = rep(3, 4))
plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
text(1, 1, "X1", cex = 2)
plot(X[,1],X[,2], col = classification, xaxt="n", yaxt="n", ann=FALSE)
axis(side = 3, las = 1, cex.axis=2, cex.lab = 2) 
plot(X[,1],X[,3], col = classification, xaxt="n", yaxt="n", ann=FALSE)
axis(side = 3, las = 1,cex.axis=2, cex.lab = 2) 
plot(X[,1],X[,4], col = classification, xaxt="n", yaxt="n", ann=FALSE)
axis(side = 4, las = 2, cex.axis=2, cex.lab = 2) 
plot(X[,2],X[,1], col = classification, xaxt="n", yaxt="n", ann=FALSE)
axis(side = 2, las = 2, cex.axis=2, cex.lab = 2) 
plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
text(1, 1, "X2", cex = 2)
plot(X[,2],X[,3], col = classification, xaxt="n", yaxt="n", ann=FALSE)
plot(X[,2],X[,4], col = classification, xaxt="n", yaxt="n", ann=FALSE)
plot(X[,3],X[,1], col = classification, xaxt="n", yaxt="n", ann=FALSE)
axis(side = 2, las = 2, cex.axis=2, cex.lab = 2) 
plot(X[,3],X[,2], col = classification, xaxt="n", yaxt="n", ann=FALSE)
plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
text(1, 1, "X3", cex = 2)
plot(X[,3],X[,4], col = classification, xaxt="n", yaxt="n", ann=FALSE)
plot(X[,4],X[,1], col = classification, xaxt="n", yaxt="n", ann=FALSE)
axis(side = 1, las = 1, cex.axis=2, cex.lab = 2) 
plot(X[,4],X[,2], col = classification, xaxt="n", yaxt="n", ann=FALSE)
plot(X[,4],X[,3], col = classification, xaxt="n", yaxt="n", ann=FALSE)
plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
text(1, 1, "X4", cex = 2)

# save original model
mod_original = mod

# subset of GMM
mod$parameters$pro = mod$parameters$pro[-1]/sum(mod$parameters$pro[-1])
mod$parameters$mean = mod$parameters$mean[,-1]
mod$parameters$variance$sigma = mod$parameters$variance$sigma[,,-1]
mod$n = nrow(mergin_component)

# PMLE
mergin_component = matrix(unlist(mergin_component), mod$n, 4)
newmod = GMMuni(mergin_component, 100, 2, 0.5, 1, "yes", mod, force_est = 0, den = FALSE)

# -----------------------------
# Application 1
# Bankruptcy dataset
# -----------------------------

# dataset
data("bankruptcy")

set.seed(123)

# model estimation
mod = Mclust(bankruptcy[,2:3])
mod_icl = mclustICL(bankruptcy[,2:3])
out_comb = clustCombi(mod)
clustCombi = clustCombiOptim(out_comb)

# Hellingher algorithm
hel = Hellingher_algorithm(mod, 0.92)
merg = hel$merg
mergin_component = bankruptcy[mod$classification==merg[1,1] | mod$classification==merg[2,1], 2:3]

# plot
classification = ifelse(mod$classification==3, 3, 2)
par(mfrow=c(1,3))
par(mar = c(5,6,4,1))
plot(bankruptcy[,2], bankruptcy[,3], col = mod$classification, pch = 19, xlab = "RE", ylab = "EBIT", cex.lab = 2.5, cex.axis = 2.5, cex=2.5)
plot(bankruptcy[,2], bankruptcy[,3], col = classification, pch = 19, xlab = "RE", ylab = "EBIT", cex.lab = 2.5, cex.axis = 2.5, cex=2.5)
plot(bankruptcy[,2], bankruptcy[,3], col = ifelse(bankruptcy[,1]==0,2,3), pch = 19, xlab = "RE", ylab = "EBIT", cex.lab = 2.5, cex.axis = 2.5, cex=2.5)

#Comparison clustHell vs clustCombi
table(clustCombi$cluster.combi, classification)

# save original model
mod_original = mod

# subset of GMM
mod$parameters$pro = mod$parameters$pro[-3]/sum(mod$parameters$pro[-3])
mod$parameters$mean = mod$parameters$mean[,-3]
mod$parameters$variance$sigma = mod$parameters$variance$sigma[,,-3]
mod$n = nrow(mergin_component)

# data to be merged
mergin_component = matrix(unlist(mergin_component), mod$n, 2)

# PMLE
newmod = GMMuni(mergin_component, 100, 2, 0.5, 1, "yes", mod, force_est = 0, den = FALSE)

# -----------------------------
# Application 2
# Thyroid
# -----------------------------

# dataset
data("thyroid")

set.seed(123)

# model estimation
mod = Mclust(thyroid[,c(2,6)])
mod_icl = mclustICL(thyroid[,c(2,6)])
out_comb = clustCombi(mod)
clustCombi = clustCombiOptim(out_comb)

# Hellingher algorithm
hel = Hellingher_algorithm(mod, 0.92)
merg = hel$merg
mergin_component = thyroid[mod$classification==merg[1,1] | mod$classification==merg[2,1], c(2,6)]

# plot
classification = ifelse(mod$classification == 1, 3, mod$classification)
par(mfrow=c(1,3))
par(mar = c(5,6,4,1))
plot(thyroid[,2], thyroid[,6], col = mod$classification, pch = 19, xlab = "RT3U", ylab = "DTSH", cex.lab = 2.5, cex.axis = 2.5, cex=2.5)
plot(thyroid[,2], thyroid[,6], col = classification, pch = 19, xlab = "RT3U", ylab = "DTSH", cex.lab = 2.5, cex.axis = 2.5, cex=2.5)
plot(thyroid[,2], thyroid[,6], col = ifelse(thyroid[,1] == "Normal", 3, ifelse(thyroid[,1] == "Hyper", 4, 2)), pch = 19, xlab = "RT3U", ylab = "DTSH", cex.lab = 2.5, cex.axis = 2.5, cex=2.5)

#Comparison clustHell vs clustCombi
table(clustCombi$cluster.combi, classification)

# save original model
mod_original = mod

# subset of GMM
mod$parameters$pro = mod$parameters$pro[-c(2,4)]/sum(mod$parameters$pro[-c(2,4)])
mod$parameters$mean = mod$parameters$mean[,-c(2,4)]
mod$parameters$variance$sigma = mod$parameters$variance$sigma[,,-c(2,4)]
mod$n = nrow(mergin_component)

# data to be merged
mergin_component = matrix(unlist(mergin_component), mod$n, 2)

# PMLE
newmod = GMMuni(mergin_component, 100, 2, 0.5, 1, "yes", mod, force_est = 0, den = FALSE)
