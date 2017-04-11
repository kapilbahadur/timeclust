# This function expands the MixSim package to provide additional functionality.
#
# MixSim is used to draw data from individual distributions, with user defined
# number of dimensions, number of observations, number of true clusters and 
# overlap of distributions (and thus clusters).
#
### Input variables ###
#
# dim - positive integer - number of dimensions of each observation (represented by columns)
# clust - positive integer - number of true clusters for output data
# size - positive integer - number of observation in the data (rows)
# overlap - positive double - overlap value to be passed to MixSim initialisation, smaller numbers
#           generate more separated clusters. Works reliably up to the value 0.75, after which 
#           the simulation times become long
# timestamp - boolean - indicate whether randomised timestamp data should be generated
# uniqe_timestamps - boolean - if TRUE, observations can have the same timestamp
# regular - boolean - if TRUE, timestamps will be spaced out equally, otherwise they will be random
# date_start - date of the format 'yyyy/mm/dd' (including apostrophes) - the earliest date point, from which
#             timestamps will be generated. Can remain NULL if timestamp is set to false.
# date_end - date of the format 'yyyy/mm/dd' (including apostrophes) - the latest date point, to which
#             timestamps will be generated.  Can remain NULL if timestamp is set to false.
# NA_val - positive integer - the number of missing values to generate in the data set. If 0 or NULL (default), 
#         no missing values will be generated
#
### The output ###
#
# The function returns the resulting data frame. It is thus necessary to save the output into a variable, e.g.:
# simulated_output <- create_data(...)
#
# The "dim", "clust", "size" and "overlap" values are passed to the MixSim function to generate the initial 
# complete dataset with the desired dimensions, observations, clustering and overlap. After the complete data
# is returned, NA values are generated and replace random entries, if the value of NA_val is specified to be
# greater than 0. A column of "true" cluster membership is also stored.
#
# If "timestamp" is set to TRUE, the next step will generate a vector of timestamp data, containing a timestamp 
# for each observation (row). These will be either completely randomly spaced, random but unique (i.e. one 
# observation per day) if "unique_timestamps" is set to TRUE, or regularly spaced if "regular" is set to TRUE. 
# The timestamps are then appended to the rest of the data as the first column. The entire data frame is then 
# sorted based on the timestamp data, in an ascending fashion.
# 


create_data <- function(dim=NULL, clust=NULL, size=NULL, overlap=0.01, NA_val=NULL, timestamp=FALSE, unique_timestamps=FALSE, regular=FALSE, date_start=NULL, date_end=NULL){

require(MixSim)
require(zoo)
  
if (is.null(dim))
    stop("dim: Select a valid integer number of dimensions (1 or more). ")
if (is.null(clust))
    stop("clust: Select a valid integer number of clusters (2 or more).")
if (is.null(size))
    stop("size: Select a valid integer number of observations (2 or more).")  

  
if (dim<1 | !all.equal(dim, as.integer(dim)))
  stop("dim: Select a valid integer number of dimensions (1 or more).")
  
if (clust<2 | !all.equal(clust, as.integer(clust)))
  stop("clust: Select a valid integer number of clusters (2 or more).")
  
if (overlap<0)
    stop("overlap: Please selet a non-negative value for overlap. Values over 0.75 can cause long simulation times.")

if (size<2 | !all.equal(size, as.integer(size)))
  stop("size: Select a valid integer number of observations (2 or more).")
  
if (timestamp==TRUE & is.null(date_start))
  stop("date_start: Select a valid stream start date in the form 'yyyy/mm/dd'.")
  
if (timestamp==TRUE & is.null(date_end))
  stop("date_end: Select a valid stream end date in the form 'yyyy/mm/dd'.")

if (!is.null(NA_val)) {
    if (!all.equal(NA_val, as.integer(NA_val)) | NA_val<0)
  stop("Select a valid integer number of missing values (argument NA_val).")  
  
  if (NA_val > dim*size)
  stop("NA_val: Selected number of missing values is greater than the number of data entries.")
}
  
repeat{
  Q <- MixSim(BarOmega = overlap, K = clust, p = dim)     # initialisation for MixSim  
  if (Q$fail == 0) break
}

A <- simdataset(n = size, Pi = Q$Pi, Mu = Q$Mu, S = Q$S, int = c(0, 1)) # output numerical values from MixSim

if(!is.null(NA_val)){
  na.which <- sample(1:(dim*size), NA_val)          # choose random positions for the given number of NA values
  na.convert <- rep(FALSE, length.out=dim*size)     
  na.convert[na.which] <- TRUE                      # create a grid indicating which positions are NA
  na.conv.matrix <- matrix(na.convert, nrow=size, ncol=dim)  # convert to matrix
  A$X[na.conv.matrix==TRUE] <- NA                   # replace with NAs
}

if(timestamp==FALSE) {
  
sim_data <- as.data.frame(cbind(A$X, A$id))         
colnames(sim_data)[ncol(sim_data)] <- "true_clust"  # combine output and label
} 

if(timestamp==TRUE) {

  if (regular==FALSE) {
    stamp <- sample(seq(as.Date(date_start), as.Date(date_end), by="day"), size, replace=!unique_timestamps) # pick random time points between start and end
  }
  
  if (regular==TRUE) {
    stamp <- seq(from=as.Date(date_start), to=as.Date(date_end),length.out = size) # or pick equally spaced time points between start and end
  }
  
  
sim_data <- as.data.frame(cbind(stamp, A$X, A$id))  # combine with data
colnames(sim_data)[1] <- "timestamp"              
colnames(sim_data)[ncol(sim_data)] <- "true_clust"  
sim_data$timestamp <- as.Date(sim_data$timestamp)   # fix format

sim_data <- sim_data[order(sim_data$timestamp),]
rownames(sim_data) <- 1:nrow(sim_data) 
}

return(sim_data)
}
