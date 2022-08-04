# Author: Rafael Bargiela, PhD. Bangor University. 2021

# Function for log2 transformation, groups filtering and Median Normalization.
################################################################
# This function uses a raw intensities table and return filtered and normalized data on a data.frame
# Data is transformed to log2 values, then filtered to those proteins present in a minimum percentage of samples on each group.
# Finally, data is normalized over the median of each sample

      # P: Matrix with raw proteins intensities.
      # groups.vector: Vector of length equal to number of columns in P, with the group assignation for each column (sample)
      # min.count.per : Minimum count of samples on each group where each protein should be present (default 60%, set 0 to include all proteins)

l2.filter.Norm<-function(P,groups.vector,min.count.per=60){
      FM<-{}
      l2P<-log(P,base=2) # Transforming RAW intensities to log2
      G<-unique(groups.vector) # Getting the different groups
      
      Gnames<-vector("character")
      for(g in 1:length(G)){
        group.names<-grep(G[g],groups.vector,value=FALSE,perl=TRUE) # Searching group on Metadata matrix (samples are by row)
        M<-as.matrix(l2P[,group.names]) # Specific group log2 matrix
        M[is.infinite(M)]<-0
        
        # Minimum count of samples of the group where protein should be present (set as 60%)
        if(min.count.per==0){
          min.count<-1 
        }else{
          min.count<-ceiling(ncol(M)*min.count.per/100) 
        }
        for(r in 1:nrow(M)){
          if(sum(M[r,]!=0,na.rm=TRUE)>=min.count){
            Gnames<-append(Gnames,rownames(M)[r],length(Gnames))
          }
        }
      }
      filter.names<-unique(Gnames)
      FM$l2P.filtered<-l2P[filter.names,]  # Filtered matrix with proteins present in at least 60% samples of any group
      
      # Normalizing filtered proteins to the median of each sample
      FM$Median.Norm<-sapply(1:ncol(FM$l2P.filtered),function(x){
        log2<-as.matrix(FM$l2P.filtered[,x])
        log2[is.infinite(log2)]<-NA
        gMedian<-median(log2,na.rm=TRUE)
        log2-gMedian
      })
      rownames(FM$Median.Norm)<-rownames(FM$l2P.filtered)
      colnames(FM$Median.Norm)<-colnames(FM$l2P.filtered)
      return(FM)
  }
