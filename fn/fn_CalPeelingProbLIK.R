CalPeelingProbLIK <- function(allef, LIK, ped, counselee.id, nloci=1, mRate=0){
  ## Aim: calculate likelihood of the pedigree through peeling and also the 
  #  posterior probability before normalization
  ## Author: Gang Peng
  ##
  ## Input:
  ## allef: List. Allele frequency for each locus/gene. If there are two genes,
  #  two alleles (allele frequency is 0.1 and 0.9) for gene 1 
  #  and three alleles (allele frequncy is 0.2, 0.2 and 0.6) for gene 2.
  #  allef = list(c(0.1,0.9),c(0.2,0.2,0.6))
  ## LIK: likelihood for each kinds of genotype
  ## ped: pedigree structure, data frame, ID, Gender, FatherID, MotherID
  #  Gender: Male, 1; Female, 0
  ## counselee.id: the ID of sample who you want to estimate the posterior
  #  probability. If you want to estimate multiple samples at the same time,
  #  just set counselee.id as a vector of ID for all samples.
  ## nloci: number of loci
  ## mRate: muation rate. If all loci/genes have same mutation rate, you just
  #  need to set one value for mRate. If loci/geens have different mutation rate,
  #  set mRate as a vector of different mutation rate for each locus/gene.
  ##
  ## Output
  # a list containning pp and pedLIK
  # pp: posterior probability without normalization
  # pedLIK: pedigree likelihood
  
  #input check
  if(length(allef)!=nloci){
    print("ERROR: number of loci in allef and nloci don't match!")
    return (-1);
  }
  
  nAllele <- rep(0,nloci)
  for(i in 1:nloci){
    nAllele[i] <- length(allef[[i]])
  }
  
  genoProb <- NULL
  for(i in 1:nloci){
    tmp <- outer(allef[[i]],allef[[i]])
    tmp <- tmp*2
    for(j in 1:nAllele[i]){
      tmp[j,j] <- tmp[j,j]/2
    }
    genoProb <- c(genoProb,as.vector(tmp[!upper.tri(tmp)]) )
  }
  
  if(length(mRate)==1){
    mRate <- rep(mRate, nloci)
  }
  else{
    if(length(mRate)!=nloci){
      sprintf("Parameter mRate and nloci is not set correctly!")
      return
    }
  }
  
  likelihood <- as.vector(t(LIK))
  
  id <- ped$ID
  gender <- ped$Gender
  fid <- ped$FatherID
  mid <- ped$MotherID
  
  for(i in 1:length(gender)){
    if(is.na(gender[i])){
      gender[i] <- 0
    }
    else{
      if(gender[i]==0){
        gender[i] <- 2
      }
    }
  }
  
  nGenoAll <- 1
  for(i in 1:nloci){
    nGenoAll <- nGenoAll*(nAllele[i]+1)*nAllele[i]/2
  }
  
  posteriorProb <- 1:(length(counselee.id)*nGenoAll)
  pedLIK <- 0;
  rltPC <- .C("CalPeelingProbLIK", genoProb = as.double(genoProb), likelihood = as.double(likelihood),
              id = as.integer(id), gender = as.integer(gender), fid = as.integer(fid),
              mid = as.integer(mid), nSample = as.integer(length(id)), cid = as.integer(counselee.id), 
              ncid = as.integer(length(counselee.id)), nloci = as.integer(nloci), nAllele = as.integer(nAllele),
              mRate = as.double(mRate), posteriorProb = as.double(posteriorProb),
              pedLIK = as.double(pedLIK))
  
  return(log(rltPC$pedLIK))
}
