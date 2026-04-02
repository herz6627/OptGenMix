#' Mean pairwise distance
#' 
#' Calculates the mean value for pairwise matrix, 
#' when supplied with weights that should be 
#' applied to the individuals 
#' 
#' @param sm  - pairwise matrix of values
#' @param w   - vector of weights
#' @return mean similarity 
#' @author Jason Bragg (jasongbragg@gmail.com)
#' @note
#' Sarah Herzog updated this function from the orignial code to run much faster for simple weighting by sub-setting the matrices.
#' @export

weighted_mean_of_pairwise_matrix <- function(sm, w=NULL) {
  
  if(!is.matrix(sm)) {stop("sm must be a matrix.")}
  
  wlen <- length(w)
  wsum <- sum(w)
  
  if (length(w) != nrow(sm)) {stop("Length of w must match number of rows in sm.")}
  
  # if weights are max 1, the calculation is simple
  if (all(w <= 1)) {
    
    sub_sm <- sm[w == 1, w == 1] # selected only weighted individuals
    return(mean(sub_sm[upper.tri(sub_sm, diag = F)])) # dont need to wory about the diagonal when an individual is only selected max once
    
  } else{
    
    # if weights are more complicated we need to adjust methods
    
    # subset matrix by removing unselected individuals (speeds up the function)
    idx <- which(w > 0)
    
    if (length(idx) < 2) return(NA_real_)
    
    sub_sm <- sm[idx, idx]
    sub_w  <- w[idx]
    
    W <- outer(sub_w, sub_w) # get matrix of weights
    weighted_mat <- sub_sm * W # multiply matrix by weights
    
    diag_term <- sum(diag(sub_sm) * (sub_w - 1) * sub_w / 2) # diagonal contribution (repeated sampling of an individual)
  
    off_diag_term <- sum(weighted_mat[upper.tri(weighted_mat)]) # off-diagonal (i < j)
    
    numerator <- off_diag_term + diag_term # sum of all trait values
    
    denominator <- sum(sub_w) * (sum(sub_w) - 1) / 2 # number of all pairwise combinations
    
    return(numerator / denominator) # get mean
    
  }
}
