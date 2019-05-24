################################################################

# (Method name TBD)
#
# Diego Greta Lauren and Richard 
# 

#Diego Vitali


################################### ORTHOGONAL-DIFF FUNCTION

land77 <- function (x, align = "", ideal = "", trim = 5, criterion = "nparallel", rep = 500, rc = T, paral = F){
  library(OpenRepGrid) 
  library(nFactors)
  # I have implemented Cohen'RC correction as Default in case neither the "align-by-ideal"
  # or the align-by-loading" gets selected. 
  
  # align-by-loading
  # swap constructs on the basis of: 
  # first component loadings
  if ( align == "loading") {
    #message("Constructs were selectively swapped according their loading on the first component of SVD")
  x <- swapPoles(x, pos = alignByLoadings(x)$reversed[,1])
  rc <- F
  }
  # Align-by-ideal self
  else { 
    if (align == "ideal"){
      if (ideal >0) {
        message("Construct have been realigned according ideal self scorings")
    x<- alignByIdeal(x,ideal = ideal)
    rc <- F
      } else{stop(message("Please indicate position of ideal self: e.g. ideal=10"))}
    }
    else{
      if (rc == TRUE){
      #  message("Using Cohen's rc for element correlations")
      }
      else{
      # Cohen's RC is applied to Element correlations
      #message("Warning: Elements correlation is susceptible to construct reflection and no construct alignment method was defined. Using Cohen's rc instead which is invariant to construct reflection")
      rc<-TRUE
      }
    }
  }
  
  X<-getRatingLayer(x, trim = trim)
  
  X# prepare the dimension index report with the 3 measures: n.factor, PVAFF and Intensity
  differentiation<-matrix(0, nrow=1, ncol=3)
  integration<-matrix(0, nrow=1, ncol=3)
  
  # calculate ordination
  
  rows<-dim(X)[1]
  # store each contruct ordination score
  ord_scores<-c()
  for (i in seq(rows)){
    ord_scores[i]<-length(unique(X[i,]))*(max(X[i,])-min(X[1,]))
  }
  ordination<-mean(ord_scores)
  
########################## between-construct 
  
  r  <- cor(t(X)) # get construct correlation
  
  colnames(differentiation) <- c("N.Factors_BTW", "Intensity", "PVAFF_BTW")
  
  ev <- eigen(r, symmetric = T) # get eigenvalues by eigen value decomposition
  # only do parallel analysis if required
  if(paral == TRUE){
  ap <- parallel(subject=nrow(r),var=ncol(r), rep = rep,quantile = .05) # Parallel: distribution of the eigenvalues of correlation matrices of random uncorrelated standardized normal variables.
  nS <- nScree(x=ev$values, aparallel=ap$eigen$qevpea, cor = T)
  # at this point we should look at what criterion to extract factors? 
  # criterion = "noc"       # Optimal coordinates (this seems better to Rache')
  # criterion = "naf"       # Acceleration factor
  # criterion = "nparallel" # parallel analysis
  # criterion = "nkaiser"   # Kaiser criterion (def don't want this one)

  differentiation[1,1]<- as.numeric(nS$Components[paste(criterion)])
  } else{
    differentiation[1,1]<- NA
  }
  differentiation[1,2]<- round(indexIntensity(x)$c.int.mean, 2)
  differentiation[1,3]<- round(ev$values[1]^2/sum(ev$values^2),2)
  

######################### within-construct 
  if(rc==TRUE){
    r <- elementCor(x, rc = rc)
  }else{
    r  <- cor(X) # get element correlation
  }
  
  colnames(integration) <- c("N.Factors_WTH", "Ordination", "PVAFF_WTH")
  
  ev <- eigen(r, symmetric = T) # get eigenvalues
  # only do parallel analysis if required
  if(paral == TRUE){
  ap <- parallel(subject=nrow(r),var=ncol(r), rep= rep, quantile = .05)
  nS <- nScree(x=ev$values, aparallel=ap$eigen$qevpea, cor = T)
  
  integration[1,1]<- as.numeric(nS$Components[paste(criterion)])
  }else{
  integration[1,1]<- NA
  }
  # intensity here is calculated without RC as poles were swapped 
  # at the beginning of the procedure. 
  integration[1,2]<- ordination
  integration[1,3]<- round(ev$values[1]^2/sum(ev$values^2), 2)
  list<-list(criterion =  criterion, El.alignment = align, cohenRC=rc)
  output<-data.frame(list,differentiation,integration)
  return(output)
}
