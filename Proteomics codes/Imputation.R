# Author: Rafael Bargiela, PhD. Bangor University (UK). 2021

# Function for data imputation of proteomics data
#################################################
  # This function uses the filtered normalized proteomics data returned from Log2.filter.norm function. Based on the fact that there can't be 0 
  # or missing values on proteomics, it assigns then a value following the Probabilistic Minimum Imputation or QRILC methods

    # M.filt.Norm: Matrix with filtered and normalized data (usually returned from Log2.filter.norm function)
    # method: Inputation method. Values allowed are MinProb or QRILC
    # tune.sigma: constant to modify the SD based on the median of the sd of all proteins. If 1, SD for imputed values is equal to the sd of the non-zero values.
    # q: Probability of the minimum quantile to get the random distrubution based on SD for imputation
    
# WARNING: Imputation produce slightly different data for missing values each time you run it due to the random assignation, 
#           so data could look a bit different.

Imputation<-function(M.filt.Norm, method="MinProb", q = 0.01, tune.sigma = 1){
  
  require(tmvtnorm)
  
  if(grepl("MinProb|QRILC",method)==FALSE){stop("Incorrect method name, only 'MinProb' or 'QRILC' are accepted")}
  
  M<-M.filt.Norm # Previously filtered matrix with proteins with more than >60% of non-NA values
  
  if(method=="MinProb"){
    print("## DATA IMPUTATION METHOD: Probabilistic Minimum Imputation (MinProb) ##")    
    PsSD<-apply(M,1,sd,na.rm=TRUE) # Getting sd of each protein, excluding missing values
    SD.imp<-median(PsSD,na.rm=TRUE)*tune.sigma # SD to get the imputed values, based on median proteins SD 
    # and tuned by constant set on tune.sigma (if 1, just ProtSD median)
  }else{
    print("## DATA IMPUTATION METHOD: Quantile regression (QRILC) ##")
    require(tmvtnorm)
  }
  
  for(c in 1:ncol(M)){
    sample<-M[,c]
    sample[is.infinite(sample)]<-NA # if inf haven't been switch to NA yet
    if(method=="MinProb"){
      min.quantile<-quantile(sample,probs = q, na.rm=TRUE) # Getting minimum value for the q quantile
      dataset.to.impute.miss<-rnorm(length(sample),mean = min.quantile,sd=SD.imp) # Getting random distribution of values based on SD for imputation
      # and minimum quantile value as mean
    }else{
      per.NAs<-length(which(is.na(sample)))/length(sample)
      upper.q<-0.99
      # Estimating mean and SD from quantile regression based on original data 
      q.normal = qnorm(seq((per.NAs+0.001),(upper.q+0.001),(upper.q-per.NAs)/(upper.q*100)),mean = 0, sd = 1)
      q.sample = quantile(sample,probs = seq(0.001,(upper.q+0.001),0.01),na.rm = TRUE) 
      QR = lm(q.sample ~ q.normal) # quantile regression
      mean.qr = as.numeric(QR$coefficients[1])
      sd.qr = as.numeric(QR$coefficients[2])
      
      dataset.to.impute.miss<-rtmvnorm(n=length(sample), mean = mean.qr,sigma = sd.qr*tune.sigma,
                                       upper = qnorm((per.NAs+0.001),mean = mean.qr,sd = sd.qr),
                                       algorithm=c("gibbs"))
    }
    sample.NAs<-which(is.na(sample))  # Changing Missing values for values on the created dataset 
    sample[sample.NAs]<-dataset.to.impute.miss[sample.NAs]
    M[,c]<-sample
  }
  return(M)
}
