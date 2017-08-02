###################################################################################################################################
##---------------------------------------------------------- kinship ------------------------------------------------------------##
###################################################################################################################################
kinshipFn <- function(G, method="additive", use="all") {
  G[G==1]<-.5
  G[G==2]<-1
  n0 <- sum(G==0,na.rm=TRUE)
  nh <- sum(G==0.5,na.rm=TRUE)
  n1 <- sum(G==1,na.rm=TRUE)
  nNA <- sum(is.na(G))

  stopifnot(n0+nh+n1+nNA == length(G))

  if ( method == "dominant" ) {
    SNP <- matrix(as.double(rowMeans(G,na.rm=TRUE) > 0.5),nrow(G),ncol(G))
    G[!is.na(G) & (G == 0.5)] <- SNP[!is.na(G) & (G == 0.5)]
  }
  else if ( method == "recessive" ) {
    SNP <- matrix(as.double(rowMeans(G,na.rm=TRUE) < 0.5),nrow(G),ncol(G))
    G[!is.na(G) & (G == 0.5)] <- SNP[!is.na(G) & (G == 0.5)]
  }
  else if ( ( method == "additive" ) && ( nh > 0 ) ) {
    dG <- G
    rG <- G
    SNP <- matrix(as.double(rowMeans(G,na.rm=TRUE) > 0.5),nrow(G),ncol(G))
    dG[!is.na(G) & (G==0.5)] <- SNP[!is.na(G) & (G==0.5)]
    SNP <- matrix(as.double(rowMeans(G,na.rm=TRUE) < 0.5),nrow(G),ncol(G))
    rG[!is.na(G) & (G==0.5)] <- SNP[!is.na(G) & (G==0.5)]
    G <- rbind(dG,rG)
  }

  if ( use == "all" ) {
    rMean <- matrix(rowMeans(G,na.rm=TRUE),nrow(G),ncol(G))
    G[is.na(G)] <- rMean[is.na(G)]
  }
  else if ( use == "complete.obs" ) {
    G <- G[rowSums(is.na(G))==0,]
  }

  n <- ncol(G)
  K <- matrix(nrow=n,ncol=n)
  diag(K) <- 1

  for(i in 2:n) {
    for(j in 1:(i-1)) {
      x <- G[,i]*G[,j] + (1-G[,i])*(1-G[,j])
      K[i,j] <- sum(x,na.rm=TRUE)/sum(!is.na(x))
      K[j,i] <- K[i,j]
    }
  }
  return(K)
}
