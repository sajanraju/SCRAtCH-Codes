# Author: Rafael Bargiela, PhD. Bangor University (UK). 2021

# Fold change calculation and t.test statistics
################################################
  # This function used the transformed, filtered, normalized and also imputed data returned from functions l2.filter.Norm and Imputation.
  # Here we get the fold change among two designated groups and calculate statistical significant differences among each or the proteins, based on 
  # t-test and moderate t-test based on limma package

    # M.Imputed: Previously filtered and normalized matrix (usually log-2 transformed). Usually, data was also previously imputed.
    # groups.vector: Vector assigning groups to each of the samples in matrix (columns). Only TWO different groups allowed.
    ## WARNING: Regard that groups.vector is not equal to the vector used on previous functions

FC.Stats<-function(M.Imputed,groups.vector){
  
  require(limma)
  
  G<-unique(groups.vector)      
  G1<-M[,grep(G[1],groups.vector,value=FALSE)]
  G2<-M[,grep(G[2],groups.vector,value=FALSE)]    
  G1av<-sapply(1:nrow(G1),function(x){mean(G1[x,],na.rm = TRUE)})
  G2av<-sapply(1:nrow(G2),function(x){mean(G2[x,],na.rm = TRUE)})
  names(G1av)<-rownames(G1)
  names(G2av)<-rownames(G2)
  FC<-G1av-G2av
  # T-test for each protein among two different groups ##
  # NOTE: assuming equal variances and different groups size
  n1<-sapply(1:nrow(G1),function(x){sum(!is.na(G1[x,]))}) # nr of cases for each protein on each group
  n2<-sapply(1:nrow(G2),function(x){sum(!is.na(G2[x,]))})  # If there aren't NA values n1 and n2 are equal for all proteins
  
  SSE1<-sapply(1:nrow(G1),function(x){sum((as.numeric(na.omit(G1[x,]))-G1av[x])^2)})
  SSE2<-sapply(1:nrow(G2),function(x){sum(as.numeric((na.omit(G2[x,]))-G2av[x])^2)})
  SSE<-SSE1+SSE2
  Sp2<-SSE/(n1+n2-2) # Pooled variance of the two groups (n1+n2-2 degrees of freedom)
  ESE<-sqrt((Sp2/n1)+(Sp2/n2))
  t.test<-FC/ESE
  p.values<-sapply(1:length(t.test),function(x){
    if(t.test[x]<0){
      2*pt(q=t.test[x],df=(n1[x]+n2[x]-2),lower.tail=TRUE)                
    }else{
      2*pt(q=t.test[x],df=(n1[x]+n2[x]-2),lower.tail=FALSE)
    }
  })
  # t.testP<-sapply(1:nrow(G1),function(x){t.test(na.omit(G1[x,]),na.omit(G2[x,]),var.equal=TRUE)$p.value})
  
  
  # Limma method
  Dat<-data.frame(M)
  Dat[is.na(Dat)]<-0
  groups<-factor(groups.vector)
  des<-model.matrix(~groups)
  fit<-lmFit(Dat,des)
  fit2<-eBayes(fit)
  res<-topTable(fit2,number=nrow(M))
  res<-res[rownames(M),]
  moderated.p.values<-res[,4]
  moderated.p.values.adj<-res[,5]
  FC.Stats<-matrix(c(G1av,G2av,FC,p.values,moderated.p.values,moderated.p.values.adj),nc=6,dimnames=list(rownames(M),c(paste(G[1],"log2 mean"),paste(G[2],"log2 mean"),"Fold change","t.test p-values","moderated-t.test p-values","moderated-t.test adj. p-values")))
  return(FC.Stats)
}
