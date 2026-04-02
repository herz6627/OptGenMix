#' Multi-objective simulated annealing
#' 
#' Performs multi-objective simulated annealing to find mixture of individuals 
#' that maximize a function (e.g. a measure of genetic diversity).
#' 
#' @param v1 Matrix. Trait data for first trait value (convert a vector or have a pairwise matrix).
#' @param v2 Matrix. Trait data for second trait value (convert a vector or have a pairwise matrix).
#' @param measure_1 Character. Objective function used to evaluate candidate solutions.
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
#' @param measure_2 Character. Secondary objective function used in
#' multi-objective optimization. Uses the same options as \code{measure_1};
#' see \code{measure_1} for full descriptions.
#' 
#' @param max_steps         Numeric. Cooling schedule, reduces temperature each iteration. Larger values give the optimizer time to test different options through slower cooling.
#' @param N_t               Integer. Number of individuals to select in the final "population."
#' @param initial_weights   Numeric. A vector of initial weights for each individual. If none provided, values are randomly assigned.
#' @param weights_max       Numeric. Maximum weights for each individual
#' @param weights_min       Numeric. Minimum weights for each individual
#' @param max_t             Numeric. Maximum temperature. At higher temperatures, the algorithm accepts solutions with worse fit to explore the search space and avoid local optima.
#' @param min_t             Numeric. Minimum temperature. Only relevant if nda = TRUE, as when nda = TRUE, the algorithm holds temp at min_t when temp < min_t. Allows for exploration of non-dominant solutions.
#' @param q                 Numeric. q value for Shannon diversity. Options are 0, 1, or 2.
#' @param p_depends_delta   Logical. If TRUE, accept depends on difference between new and old . If TRUE: if both scenarios are worse: sums the difference between the the proposals for each trait and the best proposal, then multiplies by cboth. If only one of the trait proposals is worse, multiplies by the corresponding c1/c2. If FALSe: c1 or c2 is used for the numerator. 
#' @param disp              Numeric. Displacement from zero. Used to modify the measure values, thus modifying accpetance probability.
#' @param c1                Numeric. Multiplier for numerator when proposal for value 1 worse, 2 better. Increasing (>1) this value adds resistance to worsening v1 whereas <1 lessens the resistance to accepting a worse v1.
#' @param c2                Numeric. Multiplier for numerator when proposal for value 2 worse, 1 better. Increasing (>1) this value adds resistance to worsening v2 whereas <1 lessens the resistance to accepting a worse v2.
#' @param cboth             Numeric. Multiplier for numerator when proposal for when measures for v1 and v2 are worse. Increasing (>1) this value adds resistance to worsening measures whereas <1 lessens the resistance to accepting a worse measures.
#' @param nda               Logical. if TRUE: creates a non-dominated archive. Allowing for non-dominated scenarios allows the algorithm to explore the Pareto Front more.
#' @param nd_samples        Numeric. Number of allowed non-dominated archive samples.
#' 
#' @details
#' The choice of \code{measure_1} must be compatible with the structure of the input data:
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
#' @return A list containing the results of the multi-objective optimization chain.
#'
#' The returned object has the following components:
#'
#' \describe{
#'   \item{\code{weight}}{A matrix of dimension \code{max_steps x nrow(v1)}.
#'   Each row represents the weight vector at a given iteration of the chain.}
#'
#'   \item{\code{value_1}}{A numeric vector of length \code{max_steps}.
#'   The value of the first objective (\code{measure_1}) at each cooling step.}
#'
#'   \item{\code{value_2}}{A numeric vector of length \code{max_steps}.
#'   The value of the second objective (\code{measure_2}) at each cooling step.}
#'
#'   \item{\code{accept}}{A logical vector of length \code{max_steps}.
#'   Indicates whether the proposed weights were accepted at each iteration.}
#'
#'   \item{\code{archive}}{(Optional; only if \code{nda = TRUE}).
#'   A list containing non-dominated solutions encountered during the run.
#'   It includes:
#'     \itemize{
#'       \item \code{archive_values}: A matrix with two columns (\code{value_1}, \code{value_2}) for resulting measure values for \code{v1} and \code{v2} respectively
#'       for each non-dominated solution.
#'       \item \code{archive_weights}: A matrix where each row is a weight vector corresponding
#'       to a non-dominated solution.
#'     }
#'   }
#' }
#'
#'
#' If \code{nda = TRUE}, the function also returns a Pareto archive of non-dominated
#' solutions discovered during the search. The optimization may terminate early if
#' \code{nd_samples} non-dominated solutions are found.
#' @author Jason Bragg (jasongbragg@gmail.com)
#' @export


optimize_multi_objective <- function(v1=NULL, 
                                     v2=NULL, 
                                     measure_1 = c("nei", "shannon", "allele_enrichment", "vector_weighted_mean",
                                                   "negative_vector_weighted_mean", "sum_squared_difference",
                                                   "negative_sum_squared_difference", "vector_diff_weighted_mean",
                                                   "negative_vector_diff_weighted_mean", "matrix_weighted_mean",
                                                   "negative_matrix_weighted_mean"),
                                     measure_2 = c("nei", "shannon", "allele_enrichment", "vector_weighted_mean",
                                                   "negative_vector_weighted_mean", "sum_squared_difference",
                                                   "negative_sum_squared_difference", "vector_diff_weighted_mean",
                                                   "negative_vector_diff_weighted_mean", "matrix_weighted_mean",
                                                   "negative_matrix_weighted_mean"), 
                                     max_steps=10000, 
                                     N_t=NULL, 
                                     initial_weights=NULL, 
                                     weights_max=NULL, 
                                     weights_min=NULL, 
                                     max_t=1, 
                                     min_t=0, 
                                     q=NULL, 
                                     p_depends_delta=FALSE, 
                                     disp=0, 
                                     c1=1, 
                                     c2=1, 
                                     cboth=1, 
                                     nda=FALSE, 
                                     nd_samples=100) {
  
  cat( "\n\n" )
  cat( "  Multi-objective optimization commencing. \n" )

  ### parse arguments
  measure_1 <- match.arg(measure_1)
  measure_2 <- match.arg(measure_2)
  
  # check for gt or sm
  if ( is.null(v1) ) {
    
    stop("Provide at least one matrix of values.")
  
    } else {
    N_g = nrow(v1) 
  }
  
  
  if ( is.null(v2) ) {
    cat( "   Both objective criteria will be applied to v1. \n" )
    v2=v1
  }
  
  # check for N_t or initial_weights
  if ( is.null(N_t) & is.null(initial_weights)) {stop("Provide either initial_weights or a number of target individuals (N_t)")}
  
  # if N_t but no initial weights, generate initial weights
  if ( !is.null(N_t) & is.null(initial_weights)) {
    cat( "   Generating initial weights.  \n" )
    initial_weights <- propose_initial_weights(N_g, N_t)
  }
  
  
  summary_1 <- NULL
  summary_2 <- NULL
  
  
  # generate a value of objective measure for initial
  summary_1 <- generate_measure_value(v1, measure=measure_1, w=initial_weights, q=q, disp=disp)
  summary_2 <- generate_measure_value(v2, measure=measure_2, w=initial_weights, q=q, disp=disp)
  
  if (is.null(summary_1) || is.null(summary_2)) {stop("Initial weights returned null value for a measure.")}   # if objective measure initial returns NULL, problem
  
  # if creating a non-dominated archive
  if (nda) {
    
    archive      <- list()
    nda_complete <- FALSE
    
  } else {
    
    nda_complete <- FALSE
    
  }
  
    ### run Simulated Annealing chain
    
    ### set up list to store chain
    weight  <- mat.or.vec(max_steps,N_g)
    value_1 <- mat.or.vec(max_steps,1)
    value_2 <- mat.or.vec(max_steps,1)
    accept  <- mat.or.vec(max_steps,1)
    chain   <- list(weight=weight, value_1=value_1, value_2=value_2, accept=accept)
    chain$value_1[1]   <- summary_1
    chain$value_2[1]   <- summary_2
    chain$weight[1,] <- initial_weights
    
    if (nda) {
      
      cat("   First sample added to non-dominated archive  \n")
      archive[[1]] <- list(value_1=summary_1, value_2=summary_2, weights=initial_weights)
      
    }
    
    
    # initialize some params to begin
    s <- 2
    weights <- initial_weights
    cat("   Commencing optimization \n")
    
    # MAIN CHAIN LOOP
    while ( s <= max_steps & !nda_complete ) {
      
      proposed_weights   <- propose_new_weights(weights, w_max=weights_max, w_min=weights_min)
      
      proposal_summary_1 <- generate_measure_value(v1, measure=measure_1, w=proposed_weights, q=q, disp=disp)
      proposal_summary_2 <- generate_measure_value(v2, measure=measure_2, w=proposed_weights, q=q, disp=disp)
      
      temp               <- temp_scheduler(s, max_steps, max_t=max_t)
      
      # hold at moderate temp to allow
      # exploration of nd solutions
      if (nda) {
        if (temp < min_t) {
          temp <- min_t
        }
      }
      
      accept_proposal    <- multi_accept_reject(summary_1=summary_1, summary_2=summary_2, 
                                                proposal_summary_1=proposal_summary_1, proposal_summary_2=proposal_summary_2, 
                                                temp, p_depends_delta=p_depends_delta, c1=c1, c2=c2, cboth=cboth)
      
      
      if (accept_proposal) {
        weights <- proposed_weights
        summary_1 <- proposal_summary_1
        summary_2 <- proposal_summary_2
      } 
      
      if (nda) { # moved this outside of accept_proposal so that we can explore un-accepted front space- Herzog
        archive <- reconcile_sample_nondominated_archive(summary_1, summary_2, weights, archive)
        if ( length(archive) >= nd_samples ) { nda_complete <- TRUE }
      }
      
      chain$accept[s]    <- accept_proposal
      chain$value_1[s]   <- summary_1
      chain$value_2[s]   <- summary_2
      chain$weight[s,]   <- weights
      
      s <- s + 1
      
      pfreq <- floor(max_steps / 20)
      if (s %% pfreq == 0) {cat("   Step:", s, "\n")}
      
    } # end MAIN CHAIN LOOP
    
    if (nda) {
      
      value_1 <- unlist(do.call("cbind", archive)[1,])
      value_2 <- unlist(do.call("cbind", archive)[2,])
      archive_values <- cbind(value_1, value_2)
      
      archive_weights          <- matrix( unlist(do.call("cbind", archive)[3,]),nrow=length(archive),byrow=TRUE)
      archive_matrices         <- list(archive_values=archive_values, archive_weights=archive_weights)
      chain[["archive"]]       <- archive_matrices
    }
    
    return(chain)
    
}

