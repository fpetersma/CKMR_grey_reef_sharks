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
  const CharacterVector kinship =               // So far, can be either U, S, 
    dat["kinship"];                             // or PO/OP
  const IntegerVector cov_combo_freq = 
    dat["cov_combo_freq"];
  
  const int alpha_m = dat["alpha_m"];           // male age of maturity
  const int alpha_f = dat["alpha_f"];           // female age of maturity
  const int max_age = dat["max_age"];           // maximum age
  const int t0 = dat["t0"];                     // reference  year for abundance
  const int n = dat["n"];                       // number of observations
  
  // Add parameters to be kept as constants here
  // const double r = exp(double(dat["r"]));             // growth parameter
  // const double sigma_l = exp(double(dat["sigma_l"]));
  const double phi = exp(double(dat["phi"])) /        // survival parameter
    (1.0 + exp(double(dat["phi"])));

  // ===========================================================================
  
  // ===========================================================================
  // 2. EXTRACT INDIVIDUAL PARAMETERS
  // // ---------------------------------------------------------------------------
  // const double phi = exp(double(par["phi"])) /        // survival parameter
  // (1.0 + exp(double(par["phi"])));
  const double N_t0_m = exp(double(par["N_t0_m"]));   // male abundance
  const double N_t0_f = exp(double(par["N_t0_f"]));   // female abundance
  // const double sigma_l = exp(double(par["sigma_l"]));
  const double r = exp(double(par["r"]));             // growth parameter
  
  // std::cout << "r: " << r << std::endl;
  // ===========================================================================
  
  // ===========================================================================
  // 3. DERIVE THE NEGATIVE LOG LIKELIHOOD
  // ---------------------------------------------------------------------------
  double nll = 0;
  
  for (int index = 0; index < n; index++) { 
    
    // Extract birth year of individual 1 and individual 2
    int y1 = c1[index] - a1[index];
    int y2 = c2[index] - a2[index];
    
    double prob = 0.0; 
    
    // Probability stays zero if offspring was born before maturity of parent
    if ((s1[index] == "F") & (y2 >= y1 + alpha_f) & 
        (a1[index] + y2 - c1[index] <= max_age)) {
      prob = 1.0 / (N_t0_f * pow(r, y2 - t0));
      // std::cout << "Female abundance in year " << y2 << " is: " << N_t0_f * pow(r, y2 - 1 - t0) << std::endl;
        
      // Account for survival of parent i if j was born after c1 
      if (c1[index] < y2) {
        // std::cout << "parent i if j was born after c1 + 1! " << std::endl;
        prob *=  pow(phi, y2 - c1[index]);
      }
    }
    if ((s1[index] == "M") & (y2 >= y1 + alpha_m) & 
        (a1[index] + y2 - c1[index] <= max_age)) {
      // Derive the ERRO in the year before the birth year of the offspring
      prob = 1.0 / (N_t0_m * pow(r, y2 - t0));
      
      // Account for survival of parent i if j was born after c1 + 1
      if (c1[index] < y2) {
        // std::cout << "parent i if j was born after c1 + 1! " << std::endl;
        prob *=  pow(phi, y2 - c1[index]);
      }
    }
    
    double prob12 = prob;
    // :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    // Now the other way around, i.e., 1 is the offspring and 2 is the parent
    // :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    prob = 0.0;

    // Probability stays zero if offspring was born maturity of parent
    if ((s2[index] == "F") & (y1 >= y2 + alpha_f) &
        (a2[index] + y1 - c2[index] <= max_age)) {
      // Derive the ERRO in the year before the birth year of the offspring
      prob = 1.0 / (N_t0_f * pow(r, y1 - t0));
      // std::cout << "Female abundance in year " << y2 << " is: " << N_t0_f * pow(r, y2 - 1 - t0) << std::endl;

      // Account for survival of parent j if i was born after c2
      if (c2[index] < y1) {
        // std::cout << "parent j if i was born after c2 + 1! " << std::endl;
        prob *=  pow(phi, y1 - c2[index]);
      }
    }
    if ((s2[index] == "M") & (y1 >= y2 + alpha_m) &
        (a2[index] + y1 - c2[index] <= max_age)) {
      // Derive the ERRO in the year before the birth year of the offspring
      prob = 1.0 / (N_t0_m * pow(r, y1 - t0));

      // Account for survival of parent i if j was born after c1 + 1
      if (c2[index] < y1) {
        // std::cout << "parent i if j was born after c1 + 1! " << std::endl;
        prob *=  pow(phi, y1 - c2[index]);
      }
    }
    prob += prob12;
    
    if (kinship[index] != "PO/OP") {
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
