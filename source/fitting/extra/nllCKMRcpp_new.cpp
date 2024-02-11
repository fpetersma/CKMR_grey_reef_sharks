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
double nllPOPCKMRcpp(List dat, List par) {
  // ===========================================================================
  // 1. EXTRACT DATA OBJECTS
  // ---------------------------------------------------------------------------
  const CharacterVector s1 = dat["s1"];         // sex of individual 1
  const CharacterVector s2 = dat["s2"];         // sex of individual 2
  const IntegerVector l1 = dat["l1"];           // length of individual 1
  const IntegerVector l2 = dat["l2"];           // length of individual 2
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
  const double vbgf_l_inf = dat["vbgf_l_inf"];
  const double vbgf_k = dat["vbgf_k"];
  const double vbgf_a0 = dat["vbgf_a0"];
  
  // Add parameters to be kept as constants here
  const double r = exp(double(dat["r"]));             // growth parameter
  const double sigma_vbgf = exp(double(dat["sigma_vbgf"]));
  const double phi = exp(double(dat["phi"])) /        // survival parameter
    (1.0 + exp(double(dat["phi"])));
  
  
  // ===========================================================================
  
  // ===========================================================================
  // 2. EXTRACT INDIVIDUAL PARAMETERS
  // ---------------------------------------------------------------------------
  // const double phi = exp(double(par["phi"])) /        // survival parameter
  //   (1.0 + exp(double(par["phi"])));
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
  
  IntegerVector ages = seq(0, max_age);
  int n_ages = ages.size();
  
  for (int index = 0; index < n; index++) { 
    
    double prob1 = 0.0; 
    
    // Loop over the ages for individual 1
    // .........................................................................
    for (int age_index_1 = 0; age_index_1 < n_ages; age_index_1++) {
      int a1 = ages[age_index_1];
      int y1 = c1[index] - a1;          // Extract birth year based on a1 and c1
      
      double prob2 = 0.0;
      
      // Loop over the ages for individual 2
      // .......................................................................
      for (int age_index_2 = 0; age_index_2 < n_ages; age_index_2++) {
        int a2 = ages[age_index_2];
        int y2 = c2[index] - a2;        // Extract birth year based on a2 and c2
        
        // Set probability for comparison to 0 to start
        double prob = 0.0;
        
        // :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
        //               1 is the parent, 2 is the offspring
        // :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
        
        // Derivation if parent is female
        // .....................................................................
        if (s1[index] == "F") {
          if (y2 - 1 <= y1 + alpha_f) {
            // Derive probability based on ERRO in y2 - 1
            double prob12 = 1.0 / (N_t0_f * pow(r, y2 - 1 - t0));
            
            // Correct for survival of mother, if required
            if (c1[index] < y2) {
              prob12 *= pow(phi, c1[index] - y2);    
            } 
            
            // Add prob12 to prob
            prob += prob12;
          }
        }
        // Derivation is parent if male  
        // .....................................................................  
        else if (s1[index] == "M") {
          if (y2 - 1 <= y1 + alpha_m) {
            // Derive probability based on ERRO in y2 - 1
            double prob12 = 1.0 / (N_t0_m * pow(r, y2 - 1 - t0));
            
            // Correct for survival of father, if required
            if (c1[index] < y2 - 1) {
              prob12 *=  pow(phi, c1[index] - y2 - 1);
            } 
            
            // Add prob12 to prob
            prob += prob12;
          }
        } 
        // What if the sex is neither?
        // .....................................................................
        else {
          stop("Sex has to be either 'F' (female) or 'M' (male)!");
        }
        
        // :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
        //               2 is the parent, 1 is the offspring
        // :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
        
        // // Derivation if the parent is female
        // // .....................................................................
        // if (s2[index] == "F") {
        //   if (y1 - 1 <= y2 + alpha_f) {
        //     // Derive probability based on ERRO in y1 - 1
        //     double prob21 = 1.0 / (N_t0_f * pow(r, y1 - 1 - t0));
        //     
        //     // Correct for survival of mother, if required
        //     if (c2[index] < y1) {
        //       prob21 *= pow(phi, c2[index] - y1);    
        //     } 
        //     
        //     // Add prob21 to prob
        //     prob += prob21;
        //   }
        // }
        // // Derivation is the parent if male  
        // // .....................................................................  
        // else if (s2[index] == "M") {
        //   if (y1 - 1 <= y2 + alpha_m) {
        //     // Derive probability based on ERRO in y1 - 1
        //     double prob21 = 1.0 / (N_t0_m * pow(r, y1 - 1 - t0));
        //     
        //     // Correct for survival of father, if required
        //     if (c2[index] < y1 - 1) {
        //       prob21 *=  pow(phi, c2[index] - y1 - 1);
        //     } 
        //     
        //     // Add prob21 to prob
        //     prob += prob21;
        //   } 
        // }
        // // What if the sex is neither?
        // // .....................................................................
        // else {
        //   stop("Sex has to be either 'F' (female) or 'M' (male)!");
        // }
        
        // :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
        //        What is the pdf of age given length for indiv 2?
        // :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
        
        // Find expected age based on observed length and VBGF parameters
        double a2_exp = vbgf_a0 - log(1 - l2[index] / vbgf_l_inf) / vbgf_k;
        // std::cout << "a2_exp: " << a2_exp << std::endl;
        
        // Multiply prob by probability density of individual 2's age
        prob *= R::dnorm(a2, a2_exp, sigma_vbgf, false);
        
        // std::cout << "prob: " << prob << std::endl;
        
        // Add prob to prob2
        prob2 += prob;
      } // End of a2 loop
      
      // :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      //         What is the pdf of age given length for indiv 1?
      // :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      
      // Find expected age based on observed length and VBGF parameters
      double a1_exp = vbgf_a0 - log(1 - l1[index] / vbgf_l_inf) / vbgf_k;
      
      // Multiply prob2 by probability density of individual 1's age
      prob2 *= R::dnorm(a1, a1_exp, sigma_vbgf, false);
      
      // Add prob2 to prob1
      prob1 += prob2;
    } // End of a1 loop
    
    
    // std::cout << "prob1: " << prob1 << std::endl;
    
    // Check if the kinship for this pair was observed, or not
    if (pair_found[index] == false) {
      prob1 = 1 - prob1; 
    }

    
    // Set prob2 to a very small non-zero value to avoid underflow
    if (prob1 == 0) {
      prob1 = pow(10, -60);
    }
    
    // Update the negative log-likelihood
    nll -= log(prob1) * cov_combo_freq[index];
    
  }
  // Return the negative log-likelihood
  return nll;
} 
// =============================================================================