#' Return Default Solver
#'
#' @return A string of the solver type
#' @export
#'
#' @examples
#' get_default_solver()
get_default_solver = function() {
  if (requireNamespace("gurobi", quietly = TRUE)) {
    return("gurobi")
  }
  if (requireNamespace("Rglpk", quietly = TRUE)) {
    return("glpk")
  }
  if (requireNamespace("Rcplex", quietly = TRUE) ||
      requireNamespace("cplexAPI", quietly = TRUE) ) {
    return("cplex")
  }
  if (requireNamespace("matlabr", quietly = TRUE)) {
    if (matlabr::have_matlab()) {
      return("matlab")
    }
  }
  stop("No solver found, please install a solver!")
}
