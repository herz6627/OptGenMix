#' Single-objective simulated annealing
#' 
#' Performs single-objective simulated annealing to find mixture of individuals 
#' that maximize a function (e.g. a measure of genetic diversity).
#' 
#' @param dat Matrix. Trait data for first trait value (convert a vector or have a pairwise matrix).
#' @param measure Character. Objective function used to evaluate candidate solutions.
#' Determines both the required input data type and the optimization target.
#'
#' Supported options:
#' \describe{
#'   \item{"nei"}{(genotype_matrix). Maximizes expected heterozygosity. 
#'   Interpreted biologically as maximizing neutral genetic diversity.}
#'
#'   \item{"shannon"}{(genotype_matrix). Maximizes Shannon/Hill diversity (q = 0, 1, 2).
#'   Emphasizes allelic diversity, with sensitivity to rare alleles depending on \code{q}.}
#'
#'   \item{"allele_enrichment"}{(genotype_matrix). Maximizes allele capture relative to a baseline.
#'   Useful for enriching rare or target alleles.}
#'
#'   \item{"vector_weighted_mean"}{(numeric_vector). Maximizes the weighted mean.
#'   Typically used to increase average trait values or breeding values.}
#'
#'   \item{"negative_vector_weighted_mean"}{(numeric_vector). Minimizes the weighted mean.
#'   Equivalent to selecting for lower trait values.}
#'
#'   \item{"sum_squared_difference"}{(numeric_vector). Maximizes total squared deviation.
#'   Promotes increased variance or spread in trait values. Difference is between trait value and \code{disp}}
#'
#'   \item{"negative_sum_squared_difference"}{(numeric_vector). Minimizes total squared deviation.
#'   Promotes uniformity (reduced variance). Difference is between trait value and \code{disp}}
#'
#'   \item{"vector_diff_weighted_mean"}{(numeric_vector). Maximizes mean absolute deviation.
#'   Encourages dispersion from \code{disp}. Same calculation as "sum_squared_difference", but uses the absolute difference between \code{dat} and \code{disp}.}
#'
#'   \item{"negative_vector_diff_weighted_mean"}{(numeric_vector). Minimizes mean absolute deviation.
#'   Encourages clustering around \code{disp}.}
#'
#'   \item{"matrix_weighted_mean"}{(pairwise_matrix). Maximizes the weighted mean of pairwise values.
#'   When applied to a distance matrix, this increases overall divergence among selected individuals.}
#'
#'   \item{"negative_matrix_weighted_mean"}{(pairwise_matrix). Minimizes the weighted mean of pairwise values.
#'   When applied to a distance matrix, this promotes similarity among selected individuals.}
#' }
#' 
#' @inheritParams optimize_multi_objective
#' 
#' @details
#' The choice of \code{measure} must be compatible with the structure of the input data:
#' \itemize{
#'   \item \code{genotype_matrix}: loci × individuals or similar genetic representation
#'   \item \code{numeric_vector}: trait or score per individual
#'   \item \code{pairwise_matrix}: symmetric matrix of pairwise distances or similarities
#' }
#'
#' Some measures are simple sign inversions of others (prefixed with \code{"negative_"}),
#' allowing minimization objectives to be handled within a maximization framework.
#' Original code by Jason Bragg (jasongbragg@gmail.com) with modifications by Sarah Herzog (sherzog@mobot.org).
#' 
#' @return A list containing the results of the single-objective optimization chain.
#'
#' The returned object has the following components:
#'
#' \describe{
#'   \item{\code{weight}}{A matrix of dimension \code{max_steps x nrow(v1)}.
#'   Each row represents the weight vector at a given iteration of the chain.}
#'
#'   \item{\code{value}}{A numeric vector of length \code{max_steps}.
#'   The value of the first objective (\code{measure_1}) at each cooling step.}
#'
#'   \item{\code{accept}}{A logical vector of length \code{max_steps}.
#'   Indicates whether the proposed weights were accepted at each iteration.}
#' }
#'
#'
#' @return a list of output 
#' @author Jason Bragg (jasongbragg@gmail.com)
#' @export


optimize_single_objective <- function(dat,
                                      measure=NULL, 
                                      max_steps=10000, 
                                      N_t=NULL, 
                                      initial_weights=NULL, 
                                      weights_max=NULL, 
                                      weights_min=NULL, 
                                      max_t=1, 
                                      q=NULL, 
                                      p_depends_delta=FALSE, 
                                      disp=0) {
  
  ### parse arguments
  
  N_g = nrow(dat) # number of individuals
  v = dat
  
  # check for N_t or initial_weights
  if ( is.null(N_t) & is.null(initial_weights)) {stop("Provide either initial_weights or a number of target individuals (N_t)")}
  
  if ( !is.null(N_t) && !is.null(initial_weights) && N_t != sum(initial_weights)) {stop("Conflict in number of target individuals and initial weights provided")}

  # if N_t but no initial weights, generate initial weights
  if ( !is.null(N_t) & is.null(initial_weights)) {
    cat( "   Generating initial weights  \n" )
    initial_weights <- propose_initial_weights(N_g, N_t)
  }
  
  summary <- NULL

  # generate a value of objective measure for initial
  summary <- generate_measure_value(v, measure=measure, w=initial_weights, q=q, disp=disp)
  
  # if objective measure initial returns NULL, problem
  if (is.null(summary)) {
    proceed=FALSE
    cat( "   Initial weights returned null value for measure \n" )
  }
  
    ### run Simulated Annealing chain
    
    ### set up list to store chain
    weight  <- mat.or.vec(max_steps,N_g)
    value   <- mat.or.vec(max_steps,1)
    accept  <- mat.or.vec(max_steps,1)
    chain   <- list(weight=weight, value=value, accept=accept)
    chain$value[1]   <- summary
    chain$weight[1,] <- initial_weights
    
    # initialize some params to begin
    s <- 2
    weights <- initial_weights
    
    # MAIN CHAIN LOOP
    while ( s <= max_steps ) {
      
      proposed_weights <- propose_new_weights(weights, w_max=weights_max, w_min=weights_min)
      
      proposal_summary <- generate_measure_value(v, measure=measure, w=proposed_weights, q=q, disp=disp)
      
      temp             <- temp_scheduler(s, max_steps, max_t=max_t)
      
      accept_proposal  <- proposal_accept_reject(summary=summary, proposal_summary=proposal_summary, temp, p_depends_delta=p_depends_delta)
      
      if (accept_proposal) {
        weights <- proposed_weights
        summary <- proposal_summary
      } 
      
      chain$accept[s]  <- accept_proposal
      chain$value[s]   <- summary
      chain$weight[s,] <- weights
      
      s <- s + 1
      
      pfreq <- floor(max_steps / 20)
      if (s %% pfreq == 0) {cat("   up to step", s, "\n")}
      
    } # end MAIN CHAIN LOOP
    
    return(chain)

}

