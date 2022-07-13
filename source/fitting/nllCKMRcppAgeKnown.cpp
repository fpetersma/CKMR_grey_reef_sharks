#include <Rcpp.h>
using namespace Rcpp;
// #include <cmath>
// #include <math.h>

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// =============================================================================
// PAIR PROBABILITY RCPP
// -----------------------------------------------------------------------------
// [[Rcpp::export]]
double nllPOPCKMRcppAgeKnown(List dat, List par) {
  // ===========================================================================
  // 1. EXTRACT DATA OBJECTS
  // ---------------------------------------------------------------------------
  const CharacterVector s1 = dat["s1"];         // sex of individual 1
  const CharacterVector s2 = dat["s2"];         // sex of individual 2
  const IntegerVector a1 = dat["a1"];           // length of individual 1
  const IntegerVector a2 = dat["a2"];           // length of individual 2
  const IntegerVector c1 = dat["c1"];           // capture year of individual 1
  const IntegerVector c2 = dat["c2"];           // capture year of individual 2
  const LogicalVector pair_found = 
    dat["pair_found"];
  const IntegerVector cov_combo_freq = 
    dat["cov_combo_freq"];
  
  const int alpha_m = dat["alpha_m"];           // male age of maturity
  const int alpha_f = dat["alpha_f"];           // female age of maturity
  const int max_age = dat["max_age"];           // maximum age
  const int t0 = dat["t0"];                     // reference  year for abundance
  const int n = dat["n"];                       // number of observations
  // const double vbgf_l_inf = dat["vbgf_l_inf"];
  // const double vbgf_k = dat["vbgf_k"];
  // const double vbgf_t0 = dat["vbgf_t0"];
  
  // Add parameters to be kept as constants here
  const double r = exp(double(dat["r"]));             // growth parameter
  // const double sigma_vbgf = exp(double(dat["sigma_vbgf"]));
  // const double phi = exp(double(dat["phi"])) /        // survival parameter
  //   (1.0 + exp(double(dat["phi"])));                  
  
  
  // ===========================================================================
  
  // ===========================================================================
  // 2. EXTRACT INDIVIDUAL PARAMETERS
  // ---------------------------------------------------------------------------
  const double phi = exp(double(par["phi"])) /        // survival parameter
    (1.0 + exp(double(par["phi"])));
  const double N_t0_m = exp(double(par["N_t0_m"]));   // male abundance
  const double N_t0_f = exp(double(par["N_t0_f"]));   // female abundance
  // const double sigma_vbgf = exp(double(par["sigma_vbgf"]));
  // const double r = exp(double(par["r"]));             // growth parameter
  
  // std::cout << "r: " << r << std::endl;
  // ===========================================================================
  
  // ===========================================================================
  // 3. DERIVE THE NEGATIVE LOG LIKELIHOOD
  // ---------------------------------------------------------------------------
  double nll = 0;
  
  for (int index = 0; index < n; index++) { 
    
    // Extract birth year of parent and offspring
    int y1 = c1[index] - a1[index];
    int y2 = c2[index] - a2[index];
    
    double prob = 0.0; 
    
    // Probability stays zero if offspring was born maturity of parent + 1
    if (y2 >= y1 + alpha_f + 1) {
      
      if (s1[index] == "F") {
        // What are the birthing  years?
        NumericVector birthing_years = NumericVector::create(13, 15, 17, 19);

        // Was mother at birthing age when offspring was born?
        if (is_true(any(birthing_years == y2 - y1))) {
          
          // Derive the ERRO in the year before the birth year of the offspring
          prob = 1.0 / (N_t0_f * pow(r, y2 - 1 - t0));
          // std::cout << "Female abundance in year " << y2 << " is: " << N_t0_f * pow(r, y2 - 1 - t0) << std::endl;
          
          // Account for survival of parent i if j was born after c1 
          if (c1[index] < y2) {
            std::cout << "parent i if j was born after c1 + 1! " << std::endl;
            prob *=  pow(phi, c1[index] - y2);
          }
        }
        
      }
      if (s1[index] == "M") {
        // Derive the ERRO in the year before the birth year of the offspring
        prob = 1.0 / (N_t0_m * pow(r, y2 - 1 - t0));

        // Account for survival of parent i if j was born after c1 + 1
        if (c1[index] + 1 < y2) {
          std::cout << "parent i if j was born after c1 + 1! " << std::endl;
          prob *=  pow(phi, c1[index] - y2 - 1);
        }
      }
      
      // std::cout << "out: " << out << std::endl;
    }
    
    if (pair_found[index] == false) {
      prob = 1 - prob; 
    }
    
    // std::cout << "prob: " << prob << std::endl;
    
    if (prob == 0) {
      prob = pow(10, -60);
    }
    
    // Update the negative log-likelihood
    nll -= log(prob) * cov_combo_freq[index];
  }
  // Return the negative log-likelihood
  return nll;
} 
// =============================================================================
