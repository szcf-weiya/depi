epistasisFn<-function(G,Y,coveriate)  
        {
################################################################################################
## G             : The genotype data matrix                                                   ##
## Y             : The data vector of phenotype                                               ##
## coveriate     : A data matrix of covariate (First q pca of Kinship matrix)                 ##
## min.pair.geno : The minimum number of non-zero individuals in X12 = X1*X2                  ##
## max.pair.corr : A numeric value between 0 and 1 indicating the maximum Pearson correlation ##
##                 that two markers are allowed. If the correlation between a pair of markers ##
##                 exceeds this threshold, the pair is not tested.                            ##
################################################################################################
  library(combinat)

  ID <- colnames(G)
  r1<-combn(length(ID),2)[1,]
  r2<-combn(length(ID),2)[2,]

  PVal<-NULL
  for(k in 1:1)
         {
           x1<-G[,ID[r1[k]]]
           x2<-G[,ID[r2[k]]]
           x3<-x1*x2
           
           max.pair.corr = .98
           if(var(x1)!=0 & var(x2)!=0)# in some cases all x1 or x2 are same
           {
           if( abs(cor(x1,x2)) < max.pair.corr )
              {

                  min.pair.geno = 0
                  if(sum(x3>0)>=min.pair.geno)
                     {
                       Y<-as.vector(Y)
                       if(is.null(coveriate))
                           {
                            Z<- data.frame(cbind( x1 , x2 , x3))
                           } else {
                                    Z<-  data.frame(cbind(x1 , x2 , x3 , coveriate))
                                  }
                       model = Y ~ .
                       modelSummaryCof<-summary(lm(model, data = Z))$ coefficients
                       if(nrow(modelSummaryCof) == (ncol(Z)+1) )
                             {
                               P<-modelSummaryCof[2:4,4]
                               P<-signif(P, digits = 6)
                               PVal<-rbind(PVal,c(ID[r1[k]],ID[r2[k]], P))
                             }
                     } 
              }
           }      
         }
          PVAL<-data.frame(PVal)
         # colnames(PVAL)<-c("SNP1", "SNP2", "P1", "P2", "P")
          PVAL$P1<-as.numeric(as.vector(PVAL$P1))
          PVAL$P2<-as.numeric(as.vector(PVAL$P2))
          PVAL$P<-as.numeric(as.vector(PVAL$P))
return(PVAL)
         } 



     