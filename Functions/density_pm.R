### Auxiliary functions ###

# Finds the peak
peaks <- function(x,partial=TRUE){
  if (partial){ #includes the first and the last element
    which(diff(c(TRUE,diff(x)>=0,FALSE))<0)
  }else {
    which(diff(diff(x)>=0)<0)+1
  }
}


#Function that finds the valleys
valleys <- function(x,partial=TRUE){
  if (partial){ #includes the first and the last element
    which(diff(c(FALSE,diff(x)>0,TRUE))>0)
  }else {
    which(diff(diff(x)>0)>0)+1
  }
}


#Function that calculates the probability masses
#ker: kernel density
#v: valleys
probability_mass <- function(ker,v){
  require(sfsmisc, quietly = TRUE)
  ker$y[which(ker$x<0)] = 0
  pm = c()
  for(j in 1:(length(v)-1)){
    pm[j] = integrate.xy(ker$x,ker$y,ker$x[v[j]],ker$x[v[j+1]], use.spline = FALSE)
  }
  pm = pm/sum(pm)
  return(pm)
}


#Function that tests if pm < value
#pm: probability masses
test_pm <- function(pm,value){ 
  p = c()
  num_pm = length(pm)
  for(j in 1:num_pm){
    if(pm[j]<value){
      p = c(p,j)
    }
  }
  return(p)
}



### Main functions ###

# Estimate the kernel density and calculate the probability masses
# obs : data set
# num.points : number of points for the estimation of the kernel density
density_pm <- function(obs, num.points, p.bw = p.bw, threshold.modes = threshold.modes){
  
  # fit model
  ker = density(obs, bw = p.bw, n = num.points)
  
  # find peaks
  p = peaks(ker$y)
  
  # find valleys
  v = valleys(ker$y)
  
  # probability masses
  pm = probability_mass(ker,v)
  num_pm = length(pm)
  
  # number of pm < threshold.modes
  p.del = test_pm(pm,threshold.modes)
    
  # delete modes with probability masses < threshold.modes
  for(j in 1:num_pm){
    if(j %in% p.del){
      p[j] = NA
      v[j+1] = NA
    }
  }  
  p = p[!is.na(p)]
  v = v[!is.na(v)]
  
  # probability masses (without the ones < threshold.modes)
  pm = probability_mass(ker,v)
  num_pm = length(pm)
  
  # number of pm<0.05
  p5 = test_pm(pm,0.05)
  
  # number of pm<0.10
  p10 = test_pm(pm,0.10)
  
  estimate = list(ker=ker,peaks=p,valleys=v,pm=pm,p5=p5,p10=p10)
  
  return(estimate)
}



# Estimate the kernel density and calculate the probability masses, and do bootstrap
# obs : data set
# num.points : number of points for the estimation of the kernel density
density_pm_boot <- function(obs, num.points, p.bw = p.bw, threshold.modes = threshold.modes){
  
  # fit model
  ker = density(obs,bw = p.bw, n = num.points)
  
  # find peaks
  p = peaks(ker$y)
  
  # find valleys
  v = valleys(ker$y)
  
  # probability masses
  pm = probability_mass(ker, v)
  num_pm = length(pm)
  
  # number of pm < threshold.modes
  p.del = test_pm(pm,threshold.modes)
  
  # delete modes with probability masses < threshold.modes
  for(j in 1:num_pm){
    if(j %in% p.del){
      p[j] = 0
      v[j+1] = 0
    }
  }  
  p = p[! p == 0]
  v = v[! v == 0]
  
  # probability masses (without the ones < threshold.modes)
  pm = probability_mass(ker,v)
  num_pm = length(pm)
  
  estimate = list(ker=ker,peaks=p,valleys=v,pm=pm)
  
  return(estimate)
}
