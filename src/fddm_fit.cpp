// Log-likelihood and gradient function for the Ratcliff DDM
#include <Rcpp.h>
// [[Rcpp::plugins(cpp11)]]
#include <vector>

using std::vector;

//' @export fddmFit
class fddmFit {
  public:
    vector<double> rt {};
    vector<double> response {};
    double err_tol {0.000001};
    double switch_thresh {0.8};
    vector<double> likelihood {};
    int Nrt {};
    double rt0 {1e6};
    vector<double> parameters {0.0, 0.0, 0.0, 0.0, 0.0};
    fddmFit();
    fddmFit(const vector<double>& rt_vector,
             const SEXP& response_vector);
    fddmFit(const vector<double>& rt_vector,
             const SEXP& response_vector,
             const double& error_tolerance);
    fddmFit(const vector<double>& rt_vector,
             const SEXP& response_vector,
             const double& error_tolerance,
             const double& switching_threshold);
    double calc_loglik(const vector<double>& temp_params);
    vector<double> calc_gradient(const vector<double>& temp_params);
};

#include "fitting_helper_functions/class_methods.h"


RCPP_MODULE(fddmFit) {
  using namespace Rcpp;

  class_<fddmFit>( "fddmFit")
    .constructor("Empty constructor, requires response times and responses")
    .constructor<vector<double>, SEXP >("Constructor given response times and responses")
    .constructor<vector<double>, SEXP, double>("Constructor given response times, responses, and error tolerance")
    .constructor<vector<double>, SEXP, double, double>("Constructor given response times, responses, error tolerance, and switching threshold")
    .field("rt", &fddmFit::rt)
    .field("response", &fddmFit::response)
    .field("err_tol", &fddmFit::err_tol)
    .field("switch_thresh", &fddmFit::switch_thresh)
    .field("likelihood", &fddmFit::likelihood)
    .field("parameters", &fddmFit::parameters)
    .method("calculate_loglik", &fddmFit::calc_loglik)
    .method("calculate_gradient", &fddmFit::calc_gradient)
  ;
}
