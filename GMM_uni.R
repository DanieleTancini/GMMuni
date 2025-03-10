#Function for GMMuni

GMMuni =   function(x, R, G, gamma, set, model, model_est, force_est, den){
  
  #Initial estimates
  if(model == "yes"){
    mod = model_est
  }else{
    mod = Mclust(x, G,  verbose = FALSE)
  }
  
  mode = GaussianMixtureMEM(x, mod$parameters$pro, mod$parameters$mean, mod$parameters$variance$sigma, control = list(denoise = den))
  if(force_est == 1) mode$nmodes = 2
  #Check if unimodal
  if(mode$nmodes == 1){
    print("unimodal")
  }else{
    if(set == 1){
      m = t(mode$modes)%*%(exp(mode$logdens)/sum(exp(mode$logdens)))
    }else{
      m = mod$parameters$mean%*%mod$parameters$pro
    }
    
    #Initial settings
    n = mod$n
    d = dim(x)[2]
    q = d*G
    w = q+G*((d*(d+1)/2))+1
    logp = rep(0,R)
    
    #Starting values
    Mus = MuMclust = mod$parameters$mean
    Sigmas = SigmaMclust = mod$parameters$variance$sigma
    pis = PiMclust = mod$parameters$pro
    
    #Cholesky for numerical stability
    Chol = c()
    for(i in 1:G){
      chole = c(diag(chol(mod$parameters$variance$sigma[,,i])),chol(mod$parameters$variance$sigma[,,i])[upper.tri(chol(mod$parameters$variance$sigma[,,i]))])  
      Chol = c(Chol,chole)
    }
    
    #Vector of starting values
    par1 = c(Mus, Chol, pis)
    
    #Initial settings 2 
    loglo = 0
    p2 = 0
    gr2 = matrix(0, d, G)
    mx = apply(x, 1, function(x) m-x)
    norm_vec = rep(0, n)
    
    #GMM-uni
    for(i in 1:n){
      lo = 0
      for(k in 1:G){
        pr = pis[k]*mvnfast::dmvn(t(x[i,]), Mus[, k], Sigmas[, , k])
        gr2[, k] = pr*chol2inv(chol(Sigmas[, , k])) %*% (Mus[, k] - x[i,])
        lo = lo + pr
      }
      gradient2 = rowSums(gr2)
      norm_vec[i] = sqrt(sum((mx[,i])^2))
      norm_grad2 = sqrt(sum(gradient2^2))
      cos2 = (gradient2 %*% (mx[,i])) / (norm_vec[i] * norm_grad2)
      p2 = p2 + (cos2-1)*gamma
      loglo = loglo + log(lo)
    }
    logp[1] = loglo + p2
    
    comp = matrix(0, n, G)
    
    for(r in 1:(R-1)){
      #E-step
      for(k in 1:G){
        comp[,k] =  pis[k] * mvnfast::dmvn(x, Mus[,k], Sigmas[,,k])
      }
      
      z = comp/rowSums(comp)
      L = array(0,dim = c(d,d,G))
      Sigma = array(0,dim = c(d,d,G))
      
      #M-step
      loglikp = function(par, x, z) {
        Mu = matrix(par[1:q], d, G)
        L = array(0,dim = c(d,d,G))
        for(l in 1:G){
          L[,,l] = diag(par[(q+1+((d*(d+1)/2)*(l-1))):(q+d+((d*(d+1)/2)*(l-1)))])
          L[,,l][upper.tri(L[,,l])] = par[(q + d*l + 1 + (d*(d+1)/2 -d)*(l-1)) : (q + d*l + (d*(d+1)/2 -d)*l)]
          Sigma[,,l] = t(L[,,l])%*%L[,,l]
        }
        pi = exp(par[w:(w+(G-1))])/sum(exp(par[w:(w+(G-1))])) 
        lc = 0
        p = 0
        gr = matrix(0, d, G)
        for(i in 1:n){
          for(j in 1:G){
            pr = pi[j] * mvnfast::dmvn(t(x[i,]), Mu[, j], Sigma[, , j])
            gr[, j] = pr * chol2inv(chol(Sigma[, , j])) %*% (Mu[, j] - x[i,])
            lc = lc + z[i,j]*log(pr)
          }
          gradient = rowSums(gr)
          norm_grad = sqrt(sum(gradient^2))
          cos = (gradient %*% (mx[,i])) / (norm_vec[i] * norm_grad)
          p = p + (cos-1)*gamma
        }
        lp = lc + p
        return(lp)
      }
      
      result = optim(par1, loglikp, x=x, z=z, method = "BFGS", control=list(fnscale=-1))
      npar = par1 = result$par
      
      Mus = matrix(npar[1:q], d, G)
      Ls = array(0,dim = c(d,d,G))
      for(l in 1:G){
        Ls[,,l] = diag(npar[(q+1+((d*(d+1)/2)*(l-1))):(q+d+((d*(d+1)/2)*(l-1)))])
        Ls[,,l][upper.tri(Ls[,,l])] = npar[(q + d*l + 1 + (d*(d+1)/2 -d)*(l-1)) : (q + d*l + (d*(d+1)/2 -d)*l)]
        Sigmas[,,l] = t(Ls[,,l])%*%Ls[,,l]
      }
      pis = exp(npar[w:(w+(G-1))])/sum(exp(npar[w:(w+(G-1))])) 
      
      #Penalized log-lik evaluation
      loglo = 0
      p2 = 0
      gr2 = matrix(0, d, G)
      for(i in 1:n){
        lo = 0
        for(k in 1:G){
          pr = pis[k] * mvnfast::dmvn(t(x[i,]), Mus[, k], Sigmas[, , k])
          gr2[, k] = pr*chol2inv(chol(Sigmas[, , k])) %*% (Mus[, k] - x[i,])
          lo = lo + pr
        }
        gradient2 = rowSums(gr2)
        norm_grad2 = sqrt(sum(gradient2^2))
        cos2 = (gradient2 %*% (mx[,i])) / (norm_vec[i] * norm_grad2)
        p2 = p2 + (cos2-1)*gamma
        loglo = loglo + log(lo)
      }
      logp[r+1] = loglo + p2
      if((logp[r+1]-logp[r])/abs(logp[r]) <= 1e-5) break
    }
    classification = apply(z, 1, which.max)
    
    #check mode
    check_mode = GaussianMixtureMEM(x, pis, Mus, Sigmas)
    n_modes = check_mode$nmodes
    
    res = list(Mus = Mus, Sigmas = Sigmas, pis = pis, classification = classification, m = m, MuMclust = MuMclust, SigmaMclust = SigmaMclust, PiMclust = PiMclust, p2 = p2, n_modes = n_modes, x = x)
    
  }
}
