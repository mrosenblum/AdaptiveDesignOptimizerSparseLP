#' Return Default Solver
#'
#' @return A string of the solver type
#' @export
#'
#' @examples
#' get_default_solver()
get_default_solver = function() {
  if (requireNamespace("matlabr", quietly = TRUE) && matlabr::have_matlab() &&
        (requireNamespace("Rcplex", quietly = TRUE) ||
         requireNamespace("cplexAPI", quietly = TRUE) )){
    return("cplex")
  }
  if (requireNamespace("matlabr", quietly = TRUE)) {
    if (matlabr::have_matlab()) {
       return("matlab")
    }
  }
  if (requireNamespace("gurobi", quietly = TRUE)) {
    return("gurobi")
  }
  if (requireNamespace("Rglpk", quietly = TRUE)) {
    return("glpk")
  }
  stop("No solver found, please install a linear program solver. See README file for links to four solver options.")
}
