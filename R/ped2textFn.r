ped2textFn<-function(pedFile
                     #,outFileName
                    )
  {
   Data<-pedFile[,-c(1:6)]
   n<-nrow(Data)
   p<-ncol(Data)

   G<-NULL
   k<-1
   for(i in 1:(p/2))
      {
        xx1<-paste0(Data[,k],Data[,k+1])
        xx1[xx1=="11"]<-0
        xx1[xx1=="12"]<-1
        xx1[xx1=="22"]<-2
        G<-cbind(G,as.numeric(xx1))
        k<-k+2
      }
   #write.table(G,"outFileName",row.names=F,col.names=F,quote=F,sep="\t")
   return(G)
  }