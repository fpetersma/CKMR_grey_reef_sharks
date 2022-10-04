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
double fSampledLengthFast(int l, 
                          double sigma_l,
                          int max_age,
                          NumericVector probs_sampled_ages,
                          int l_inf, 
                          double k,
                          double a_0) {
  
  double prob_mass = 0.0;
  
  for (int a = 0; a <= max_age; a++) {
    prob_mass += fLengthGivenAge(l, a, sigma_l, l_inf, k, a_0) *
      probs_sampled_ages[a];
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
double fAgeGivenLengthFast(int a, 
                           int l, 
                           double sigma_l, 
                           NumericVector probs_sampled_ages,
                           int max_age,
                           int l_inf, 
                           double k, 
                           double a_0) {
  
  double prob_length_given_age = fLengthGivenAge(l, a, sigma_l, l_inf, k, a_0);
  
  double prob_sampled_length = fSampledLengthFast(l, sigma_l, max_age, 
                                                  probs_sampled_ages, 
                                                  l_inf, k, a_0);
  
  double out = prob_length_given_age * probs_sampled_ages[a] / 
    prob_sampled_length;
  
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
  // DO SOME THINGS TO OPTIMISE THE CODE
  // ---------------------------------------------------------------------------
  // Derive fSampledAge() for all possible ages
  NumericVector probs_sampled_ages (max_age + 1);
  for (int a = 0; a <= max_age; a++) {
    probs_sampled_ages[a] = fSampledAge(a, 1 - phi, max_age);
    // std::cout << "prob of " << a << ": " << probs_sampled_ages[a] << std::endl;
  }

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
      
      double prob12_1 = 0.0;
      
      // max_age_2 depends on c2, y1, and the alpha for the gender of 1
      int max_age_2 = c2[index] - y1;
      // corrected for alpha below
      if (s1[index] == "F") {
        max_age_2 -= alpha_f;
      } else {
        max_age_2 -= alpha_m;
      }

      for (int a2 = 0; a2 <= max_age_2; a2++) {
        // Extract birth year of individual 2
        int y2 = c2[index] - a2;
        
        // :::::::s::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
        // First the initial direction, i.e., 1 is the parent and 2 is the offspring
        // :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
        double prob12_2 = 0.0; 
        
        // Probability stays zero if offspring was born before maturity of parent
        if ((s1[index] == "F") & (a1 + y2 - c1[index] <= max_age)) {
          prob12_2 = 1.0 / (N_t0_f * std::pow(r, y2 - t0));

          // Account for survival of parent i if j was born after c1 
          if (c1[index] < y2) {
            prob12_2 *=  std::pow(phi, y2 - c1[index]);
          }
        }
        if ((s1[index] == "M") & (a1 + y2 - c1[index] <= max_age)) {
          // Derive the ERRO in the year before the birth year of the offspring
          prob12_2 = 1.0 / (N_t0_m * std::pow(r, y2 - t0));
          
          // Account for survival of parent i if j was born after c1 + 1
          if (c1[index] < y2) {
            prob12_2 *=  std::pow(phi, y2 - c1[index]);
          }
        }
        // Fast version
        prob12_1 += prob12_2 * fAgeGivenLengthFast(a2, l2[index], sigma_l,
                                                   probs_sampled_ages, max_age,
                                                   vbgf_l_inf, vbgf_k, vbgf_a0);
        
        // Slow version
        // prob12_1 += prob12_2 * fAgeGivenLength(a2, l2[index], sigma_l,
        //                                        1 - phi, max_age, vbgf_l_inf,
        //                                        vbgf_k, vbgf_a0);
        
      }
      // Fast version
      prob += prob12_1 * fAgeGivenLengthFast(a1, l1[index], sigma_l,
                                             probs_sampled_ages, max_age,
                                             vbgf_l_inf, vbgf_k, vbgf_a0);
      // Slow version
      // prob += prob12_1 * fAgeGivenLength(a1, l1[index], sigma_l,
      //                                    1 - phi, max_age, vbgf_l_inf,
      //                                    vbgf_k, vbgf_a0);
    }
    
    // Loop over the ages
    for (int a2 = 0; a2 <= max_age; a2++) {
      // Extract birth year of individual 2
      int y2 = c2[index] - a2;
      
      double prob21_1 = 0.0;
      
      // max_age_1 depends on c1, y2, and the alpha for the gender of 2
      int max_age_1 = c1[index] - y2;
      // corrected for alpha below
      if (s1[index] == "F") {
        max_age_1 -= alpha_f;
      } else {
        max_age_1 -= alpha_m;
      }
      
      for (int a1 = 0; a1 <= max_age_1; a1++) {
        
        // Extract birth year of individual 1
        int y1 = c1[index] - a1;
        
        // :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
        // First the initial direction, i.e., 1 is the parent and 2 is the offspring
        // :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
        double prob21_2 = 0.0; 
        
        // Probability stays zero if offspring was born before maturity of parent
        if ((s2[index] == "F") & (a2 + y1 - c2[index] <= max_age)) {
          prob21_2 = 1.0 / (N_t0_f * std::pow(r, y1 - t0));
          
          // Account for survival of parent i if j was born after c1 
          if (c2[index] < y1) {
            prob21_2 *=  std::pow(phi, y1 - c2[index]);
          }
        }
        if ((s2[index] == "M") & (a2 + y1 - c2[index] <= max_age)) {
          // Derive the ERRO in the year before the birth year of the offspring
          prob21_2 = 1.0 / (N_t0_m * std::pow(r, y1 - t0));
          
          // Account for survival of parent i if j was born after c1 + 1
          if (c2[index] < y1) {
            prob21_2 *=  std::pow(phi, y1 - c2[index]);
          }
        }
        // Fast version
        prob21_1 += prob21_2 * fAgeGivenLengthFast(a1, l1[index], sigma_l,
                                                   probs_sampled_ages, max_age,
                                                   vbgf_l_inf, vbgf_k, vbgf_a0);
        // Slow version
        // prob21_1 += prob21_2 * fAgeGivenLength(a1, l1[index], sigma_l,
        //                                        1 - phi, max_age, vbgf_l_inf,
        //                                        vbgf_k, vbgf_a0);
      }
      // Fast version
      prob += prob21_1 * fAgeGivenLengthFast(a2, l2[index], sigma_l,
                                             probs_sampled_ages, max_age,
                                             vbgf_l_inf, vbgf_k, vbgf_a0);
      // Slow version
      // prob += prob21_1 * fAgeGivenLength(a2, l2[index], sigma_l,
      //                                    1 - phi, max_age, vbgf_l_inf,
      //                                    vbgf_k, vbgf_a0);

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
