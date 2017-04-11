# The timeclust fucntion uses native R clustering fuctions to perform a 2 stage
# clustering of timestamp-labelled data. The two stages are used to cluster 
# numerical data and timestamp data respectively.
#
### Input variables ###
#
# clust_data - a data frame or a matrix containing real numerical values
# time_data - a vector of dates containing timestamps of the data
# clust_method - the method for clustering numerical data - either "hclust", "kmeans" or "kpod" (including quotation marks)
# k - positive integer - the desired number of clusters for the numerical data
# dist_m - the distance method for the numerical data, in case clust_method is set to "hclust" - any native R 
#           method is possible (e.g. "euclidean", "maximum", "manhattan", etc., including quotation marks)
# hclust_m - the linkage method for the numerical data, in case clust_method is set to "hclust" - any native R
#           method is possible (e.g. "single", "complete", "average", etc., including quotation marks)
# time_method - the method for clustering timestamp data - either "hclust", "kmeans" or "kpod" (including quotation marks)
# k_t - positive integer - the desired number of clusters for the timestamp data
# dist_m_t - the distance method for the timestamp data, in case time_method is set to "hclust" - any native R 
#           method is possible (e.g. "euclidean", "maximum", "manhattan", etc., including quotation marks)
# hclust_m_t - the linkage method for the timestamp data, in case time_method is set to "hclust" - any native R
#           method is possible (e.g. "single", "complete", "average", etc., including quotation marks)
# timefirst - boolean - if FALSE, numerical data is clustered first, then timestamp based clusters are determined within
#             the numerical clusters (and vice versa for TRUE)
#
### The output ###
#
# The function returns a data frame with two columns named "timeLabels" and "dataLabels" that need to be stored, e.g.:
# timeclust_output_clusters <- timeclust(...)
#
# The function uses built-in R functions hclust and kmeans to perform clustering. In case missing values are observed in the
# numerical data, the method used is k-POD clustering, from the kpodclustr package. The functions performs clustering in
# two stages. The method of clustering can be specified for either of these individually, as can be the desired number of 
# clusters, the distance method employed and the linkage used if hierarchical clustering is used. The order is depends on 
# whether the setting of the "timefirst" parameter. 
#
# If "timefirst" is set to FALSE, the numerical data saved in "clust_data" is clustered first using the clustering method 
# defined in "clust_method", with the number of clusters specified in "k". If hierarchical clustering is selected in 
# "clust_method", the parameter "dist_m" can be specified to choose a method for distance measurement (default: "euclidean"), 
# and "hclust_m" can be specified to select linkage method (default: "average"). Distance and linkage are passed to 
# the built-in hclust function, therefore any supported methods can be used. Subsequently, the "time_data" is clustered 
# using the clustering method defined in "time_method", with the number of clusters specified in "k_t". k_t clusters 
# are determined within each of the numerical clusters individually. Again, if hierarchical clustering is selected in 
# "time_method", the parameter "dist_m_t" can be specified to choose a method for distance measurement (default: "euclidean"), 
# and "hclust_m_t" can be specified to select linkage method (default: "average").
#
# Alternatively, if "timefirst" is set to TRUE, the timestamp data saved in "time_data" is clustered first using 
# the clustering method defined in "time_method", with the number of clusters specified in "k_t". If hierarchical clustering 
# is selected in "time_method", the parameter "dist_m_t" can be specified to choose a method for distance measurement 
# (default: "euclidean"), and "hclust_m_t" can be specified to select linkage method (default: "average"). As in the previous 
# case, distance and linkage are passed to the built-in hclust function, therefore any supported methods can be used.
# Afterwards, the numerical data stored in "clust_data" are clustered within each of the clusters from the previous clustering
# step. These are clustered into clusters stored in "k". If hierarchical clustering is chosen, distance and linkage can be 
# specified with "dist_m" and "hclust_m" as before.



timeclust <- function(clust_data=NULL, time_data=NULL, clust_method = NULL, k=NULL, dist_m="euclidean", hclust_m="average",
                       time_method=NULL, k_t=NULL, dist_m_t="euclidean", hclust_m_t="average", timefirst=FALSE){
  
  if (sum(is.na(clust_data))>0 && clust_method!="kpod") {
    require(kpodclustr)
    clust_method = "kpod"
    print("Missing values detected: setting clust_method to k-POD clustering.")
  }
  
  #  **** OPTION 1 ****
  if (timefirst==FALSE) {
    #cluster data first
    if (clust_method=="hclust"){
      dis <- dist(clust_data, method = dist_m)
      h_tree <- hclust(dis, method = hclust_m)
      cut_h <- cutree(h_tree, k)
      grp <- cut_h
    }
    
    if (clust_method=="kmeans") {
      k_grp <- kmeans(clust_data, k)
      grp <- k_grp$cluster
    }
    if (clust_method=="kpod") {
      k_grp <- kpod(as.matrix(clust_data), k)
      grp <- k_grp$cluster
    }
    
    #now cluster time
    
    out <- rep(0, length.out=(length(time_data))) # vector to store timestamp clustering results
    
    for (i in 1:k){
      k_t_sub <- k_t
      time_data_sub <- as.Date(time_data[grp==i])        # cycle through each cluster, then group timestamps within those clusters separately
      if (length(time_data_sub) < k_t) {
        k_t_sub <- length(time_data_sub)
      }
      if(k_t_sub == 1){
        out[grp==i] <- 1
      } else {
        if (time_method=="hclust"){
         dis_t <- dist(time_data_sub, method = dist_m_t)
          h_tree_t <- hclust(dis_t, method = hclust_m_t)
          cut_h_t <- cutree(h_tree_t, k_t_sub)
          grp_t <- cut_h_t
        }
        if (time_method=="kmeans") {
          k_grp_t <- kmeans(time_data_sub, k_t_sub)
          grp_t <- k_grp_t$cluster
        }
        out[grp==i] <- grp_t                      # save labels
      }
    }
    out_clusts <- as.data.frame(cbind(grp, out))
    colnames(out_clusts) <- c("dataLabels", "timeLabels")
  }
  
  #  **** OPTION 2 ****
  if (timefirst==TRUE){
    #cluster time first
    if (time_method=="hclust"){
      dis_t <- dist(as.Date(time_data), method = dist_m_t)
      h_tree_t <- hclust(dis_t, method = hclust_m_t)
      cut_h_t <- cutree(h_tree_t, k_t)
      grp_t <- cut_h_t
    }
    
    if (time_method=="kmeans") {
      k_grp_t <- kmeans(as.Date(time_data), k_t)
      grp_t <- k_grp_t$cluster
    }
    
    #now cluster data
    
    out <- rep(0, length.out=(nrow(clust_data)))
    
    for (i in 1:k_t){
      k_sub <- k
      clust_data_sub <- clust_data[grp_t==i,]
      if (nrow(clust_data_sub) < k) {
        k_sub <- nrow(clust_data_sub)
      }
      if (k_sub == 1) {
        out[grp_t==i] <- 1
      } else {
        if (clust_method=="hclust"){
          dis <- dist(clust_data_sub, method = dist_m)
          h_tree <- hclust(dis, method = hclust_m)
          cut_h <- cutree(h_tree, k_sub)
          grp <- cut_h
        }
        if (clust_method=="kmeans") {
          k_grp <- kmeans(clust_data_sub, k_sub)
          grp <- k_grp$cluster
        }
        if (clust_method=="kpod") {
          k_grp <- kpod(as.matrix(clust_data), k_sub)
          grp <- k_grp$cluster
        }
        out[grp_t==i] <- grp
      }
    }
    out_clusts <- as.data.frame(cbind(grp_t, out))
    colnames(out_clusts) <- c("timeLabels", "dataLabels")
  }
  
  # OUTPUT
  return(out_clusts)
}