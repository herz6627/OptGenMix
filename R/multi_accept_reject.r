#' Determines acceptance or rejection of a proposal
#' 
#' @param summary          Summary of current.
#' @param proposal_summary Summary of new proposal.
#' @param t                Current temperature.
#' @param p_depends_delta  If TRUE, accept depends on difference between new and old.
#' @param c1               Multiplier for numerator when proposal for value 1 worse, 2 better.
#' @param c2               Multiplier for numerator when proposal for value 2 worse, 1 better.
#' @param cboth            Multiplier for numerator when proposals for both values are worse.
#' @return Boolean, TRUE if proposal accepted.
#' @author Jason Bragg (jasongbragg@gmail.com)
#' @export

multi_accept_reject <- function(summary_1=summary_1, summary_2=summary_2, proposal_summary_1=proposal_summary_1, proposal_summary_2=proposal_summary_2, t, p_depends_delta=FALSE, c1=1, c2=1, cboth=1) {

   accept=FALSE
   
   # add catch for when t = 0
   if (t == 0) {t = 1e-10}  # avoid divide-by-zero

   # if better, accept
   if ( proposal_summary_1 > summary_1 & proposal_summary_2 > summary_2 ) {

      accept=TRUE

   } else {

      delta_1 <- summary_1 - proposal_summary_1 
      delta_2 <- summary_2 - proposal_summary_2 

      # both proposals worse
      if ( proposal_summary_1 <= summary_1 & proposal_summary_2 <= summary_2 ) {

         if (p_depends_delta) {
            numerator <- (delta_1 + delta_2) * cboth
         } else {
            numerator <- cboth
         }

      }

      # objective 1 better, objective 2 worse
      if ( proposal_summary_1 > summary_1 & proposal_summary_2 <= summary_2 ) {

         if (p_depends_delta) {
            numerator <- delta_2 * c2
         } else {
            numerator <- c2
         }

      }

      # objective 2 better, objective 1 worse
      if ( proposal_summary_1 <= summary_1 & proposal_summary_2 > summary_2 ) {

         if (p_depends_delta) {
            numerator <- delta_1 * c1
         } else {
            numerator <- c1
         }

      }

      p_accept_worse <- exp( ( -numerator ) / t  )

      #cat(p_accept_worse, "\n")
      if ( runif(1) < p_accept_worse ) {
         accept=TRUE
      }
   }

   return(accept)
}


