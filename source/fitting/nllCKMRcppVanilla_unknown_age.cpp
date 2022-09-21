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
// Non-vectorised probability mass functions used in the likelihood
// -----------------------------------------------------------------------------
// Discrete for of the normal distribution 
double pDiscreteNorm(int x, 
                     double mu, 
                     double sigma) {
  double out = R::pnorm(x + 0.5, mu, sigma, true, false) - 
    R::pnorm(x - 0.5, mu, sigma, true, false);
  
  return out;
}

// Probability mass function for length given age
double fLengthGivenAge(int l, 
                       int a, 
                       double sigma_l, 
                       int l_inf, 
                       double k,
                       double a_0) {
  double expected_length = l_inf * (1 - std::exp(-k * (a - a_0)));
  
  double out = pDiscreteNorm(l, expected_length, sigma_l);
  
  return out;
}

// Probability mass function for sampled age
double fSampledAge(int a, 
                   double p_geom, 
                   int max_age) {
  double out = R::dgeom(a, p_geom, false) / 
    R::pgeom(max_age + 1, p_geom, true, false);
  
  return out;
}

// Probability mass function for sampled length
double fSampledLength(int l, 
                      double sigma_l,
                      int max_age,
                      double p_geom,
                      int l_inf, 
                      double k,
                      double a_0) {
  
  double prob_mass = 0.0;
  
  for (int a = 0; a <= max_age; a++) {
    prob_mass += fLengthGivenAge(l, a, sigma_l, l_inf, k, a_0) *
      fSampledAge(a, p_geom, max_age);
  }
  
  return prob_mass;
}

// Probability mass function for age given length
double fAgeGivenLength(int a, 
                       int l, 
                       double sigma_l, 
                       double p_geom,
                       int max_age,
                       int l_inf, 
                       double k, 
                       double a_0) {
  
  double prob_length_given_age = fLengthGivenAge(l, a, sigma_l, l_inf, k, a_0);
  
  double prob_sampled_age = fSampledAge(a, p_geom, max_age);
  
  double prob_sampled_length = fSampledLength(l, sigma_l, max_age, p_geom, 
                                              l_inf, k, a_0);
  
  double out = prob_length_given_age * prob_sampled_age / prob_sampled_length;
  
  return out;
}

// -----------------------------------------------------------------------------

// =============================================================================
// PAIR PROBABILITY RCPP
// -----------------------------------------------------------------------------
// [[Rcpp::export]]
double nllPOPCKMRcppAgeUnknown(List dat, List par) {
  // ===========================================================================
  // 1. EXTRACT DATA OBJECTS
  // ---------------------------------------------------------------------------
  const CharacterVector s1 = dat["s1"];         // sex of individual 1
  const CharacterVector s2 = dat["s2"];         // sex of individual 2
  const IntegerVector l1 = dat["l1"];           // length of individual 1
  const IntegerVector l2 = dat["l2"];           // length of individual 2
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
  const double vbgf_l_inf = dat["vbgf_l_inf"];
  const double vbgf_k = dat["vbgf_k"];
  const double vbgf_a0 = dat["vbgf_a0"];
  
  // Add parameters to be kept as constants here
  // const double r = std::exp(double(dat["r"]));             // growth parameter
  const double sigma_l = std::exp(double(dat["sigma_l"]));
  const double phi = std::exp(double(dat["phi"])) /        // survival parameter
    (1.0 + std::exp(double(dat["phi"])));
  
  // ===========================================================================
  
  // ===========================================================================
  // 2. EXTRACT INDIVIDUAL PARAMETERS
  // // ---------------------------------------------------------------------------
  // const double phi = std::exp(double(par["phi"])) /        // survival parameter
  // (1.0 + std::exp(double(par["phi"])));
  const double N_t0_m = std::exp(double(par["N_t0_m"]));   // male abundance
  const double N_t0_f = std::exp(double(par["N_t0_f"]));   // female abundance
  // const double sigma_l = std::exp(double(par["sigma_l"]));
  const double r = std::exp(double(par["r"]));             // growth parameter
  
  // std::cout << "r: " << r << std::endl;
  // ===========================================================================
  
  // ===========================================================================
  // 3. DERIVE THE NEGATIVE LOG LIKELIHOOD
  // ---------------------------------------------------------------------------
  double nll = 0;
  
  for (int index = 0; index < n; index++) { 
    
    double prob = 0.0;
    
    // Loop over the ages
    for (int a1 = 0; a1 <= max_age; a1++) {
      // Extract birth year of individual 2
      int y1 = c1[index] - a1;
      
      double prob1 = 0.0;
      
      for (int a2 = 0; a2 <= max_age; a2++) {
        // Extract birth year of individual 2
        int y2 = c2[index] - a2;
        
        // :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
        // First the initial direction, i.e., 1 is the parent and 2 is the offspring
        // :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
        double prob12 = 0.0; 
        
        // Probability stays zero if offspring was born before maturity of parent
        if ((s1[index] == "F") & (y2 >= y1 + alpha_f) & 
            (a1 + y2 - c1[index] <= max_age)) {
          prob12 = 1.0 / (N_t0_f * std::pow(r, y2 - t0));
          // std::cout << "Female abundance in year " << y2 << " is: " << N_t0_f * std::pow(r, y2 - 1 - t0) << std::endl;
          
          // Account for survival of parent i if j was born after c1 
          if (c1[index] < y2) {
            // std::cout << "parent i if j was born after c1 + 1! " << std::endl;
            prob12 *=  std::pow(phi, y2 - c1[index]);
          }
        }
        if ((s1[index] == "M") & (y2 >= y1 + alpha_m) & 
            (a1 + y2 - c1[index] <= max_age)) {
          // Derive the ERRO in the year before the birth year of the offspring
          prob12 = 1.0 / (N_t0_m * std::pow(r, y2 - t0));
          
          // Account for survival of parent i if j was born after c1 + 1
          if (c1[index] < y2) {
            // std::cout << "parent i if j was born after c1 + 1! " << std::endl;
            prob12 *=  std::pow(phi, y2 - c1[index]);
          }
        }

        // :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
        // Now the other way around, i.e., 1 is the offspring and 2 is the parent
        // :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
        double prob21 = 0.0;
        
        // Probability stays zero if offspring was born maturity of parent
        if ((s2[index] == "F") & (y1 >= y2 + alpha_f) &
            (a2 + y1 - c2[index] <= max_age)) {
          // Derive the ERRO in the year before the birth year of the offspring
          prob21 = 1.0 / (N_t0_f * std::pow(r, y1 - t0));
          // std::cout << "Female abundance in year " << y2 << " is: " << N_t0_f * std::pow(r, y2 - 1 - t0) << std::endl;
          
          // Account for survival of parent j if i was born after c2
          if (c2[index] < y1) {
            // std::cout << "parent j if i was born after c2 + 1! " << std::endl;
            prob21 *=  std::pow(phi, y1 - c2[index]);
          }
        }
        if ((s2[index] == "M") & (y1 >= y2 + alpha_m) &
            (a2 + y1 - c2[index] <= max_age)) {
          // Derive the ERRO in the year before the birth year of the offspring
          prob21 = 1.0 / (N_t0_m * std::pow(r, y1 - t0));
          
          // Account for survival of parent i if j was born after c1 + 1
          if (c2[index] < y1) {
            // std::cout << "parent i if j was born after c1 + 1! " << std::endl;
            prob21 *=  std::pow(phi, y1 - c2[index]);
          }
        }
        double prob2 = prob12 + prob21;
        
        prob2 *= fAgeGivenLength(a2, l2[index], sigma_l, 1 - phi, max_age, 
                                 vbgf_l_inf, vbgf_k, vbgf_a0);
        
        // // Find expected length based on age and VBGF parameters
        // double l2_exp = vbgf_l_inf * (1 - std::exp(-vbgf_k * (a2 - vbgf_a0)));
        // // Find expected length based on age and VBGF parameters
        // double l1_exp = vbgf_l_inf * (1 - std::exp(-vbgf_k * (a1 - vbgf_a0)));

        // double a2_exp = vbgf_a0 - log(1 - l2[index] / vbgf_l_inf) / vbgf_k;
        // double a1_exp = vbgf_a0 - log(1 - l1[index] / vbgf_l_inf) / vbgf_k;
        
        // // Correct prob2 for f(a2)
        // prob2 *= R::dnorm(a2, a2_exp, sigma_l, false);
        // // Correct prob2 for f(a2)
        // prob2 *= R::dnorm(a1, a1_exp, sigma_l, false);
        
        // Add prob2 to prob_1
        prob1 += prob2;
      }
      
      prob1 *= fAgeGivenLength(a1, l1[index], sigma_l, 1 - phi, max_age, 
                               vbgf_l_inf, vbgf_k, vbgf_a0);
      // // Find expected length based on age and VBGF parameters
      // double l1_exp = vbgf_l_inf * (1 - std::exp(-vbgf_k * (a1 - vbgf_a0)));
      // 
      // // Correct prob2 for f(a2)
      // prob1 *= R::dnorm(l1_exp, l1[index], sigma_l, false);
      // 
      prob += prob1;
    }
    
    if (kinship[index] != "PO/OP") {
      prob = 1 - prob; 
    }
    
    // std::cout << "prob: " << prob << std::endl;
    
    if (prob == 0) {
      prob = std::pow(10, -60);
    }
    
    // Update the negative log-likelihood
    nll -= std::log(prob) * cov_combo_freq[index];
  }
  // Return the negative log-likelihood
  return nll;
} 
// =============================================================================
