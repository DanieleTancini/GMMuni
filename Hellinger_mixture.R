#Function for mixture Hellinger algorithm

hellinger.dist.mod = function(mu1, mu2, sigma1, sigma2, w1, w2){
  hd = sqrt(w1+w2-(2*sqrt(w1*w2))*exp(-bhattacharyya.dist(mu1, mu2, sigma1, sigma2)))
  return(as.numeric(hd))
}

Hellingher_algorithm = function(mod, ep){
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
  
  #Algorithm
  IH = 1*(H > 0 & H <= ep)
  merg = matrix(0, G, G)
  V = 1:G
  E = as.vector(t(which(IH == 1, arr.ind = TRUE)))
  if(length(E)==0){
    merg = merg
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
  res = list(merg = merg, H = H)
}


