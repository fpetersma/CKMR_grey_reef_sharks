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
double nllPOPCKMRcppAgeUnknownGestation(List dat, List par) {
  // ===========================================================================
  // 1. EXTRACT DATA OBJECTS
  // ---------------------------------------------------------------------------
  const IntegerVector s1 = dat["s1"];         // sex of individual 1 (1=F, 0=M)
  const IntegerVector s2 = dat["s2"];         // sex of individual 2
  const IntegerVector l1 = dat["l1"];           // length of individual 1
  const IntegerVector l2 = dat["l2"];           // length of individual 2
  const IntegerVector c1 = dat["c1"];           // capture year of individual 1
  const IntegerVector c2 = dat["c2"];           // capture year of individual 2
  const IntegerVector kinship =               // So far, can be either U=0, PO/OP=1, 
    dat["kinship"];                             // or U=2.
  const IntegerVector cov_combo_freq = 
    dat["cov_combo_freq"];
  
  const int alpha_m = dat["alpha_m"];           // male age of maturity
  const int alpha_f = dat["alpha_f"];           // female age of maturity
  const int max_age = dat["max_age"];           // maximum age
  const int max_length = dat["max_length"];     // maximum possible length
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
  }
  
  // Derive fAgeGivenLength for all age-length combos
  NumericMatrix age_length_prob_matrix ((max_age + 1), (max_length + 1));
  for (int age = 0; age <= max_age; age++) {
    for (int length = 0; length <= max_length; length++) {
      age_length_prob_matrix(age, length) = fAgeGivenLengthFast(age, 
                             length, 
                             sigma_l,
                             probs_sampled_ages, 
                             max_age,
                             vbgf_l_inf, 
                             vbgf_k, 
                             vbgf_a0);
    }
  }
  // std::cout << "Age-Length probability matrix (age = 13, length = 171): " << age_length_prob_matrix(13, 171) << std::endl;
  // std::cout << "Age-Length probability matrix (age = 10, length = 150): " << age_length_prob_matrix(10, 150) << std::endl;
  // std::cout << "Age-Length probability matrix (age = 5, length = 136): " << age_length_prob_matrix(3, 136) << std::endl;
  // 
  
  // // Below, to avoid getting negative indices, add 40 to the indices
  // // Derive values N_t0_f back in time for 0:19 years of population growth r
  // NumericVector N_t0_f_vector (40);
  // for (int i = -40; i < 0; i++) {
  //   N_t0_f_vector[i + 40] = N_t0_f * std::pow(r, i);
  // }
  // // Derive values N_t0_m back in time for 0:19 years of population growth r
  // NumericVector N_t0_m_vector (40);
  // for (int i = -40; i < 0; i++) {
  //   N_t0_f_vector[i + 40] = N_t0_m * std::pow(r, i);
  // }
  // 
  // // Derive values for phi for 0:10 years of survival
  // NumericVector phi_vector (10);
  // for (int i = 0; i < 10; i++) {
  //   phi_vector[i] = std::pow(phi, i);
  // }
  
  // ===========================================================================
  // 3. DERIVE THE NEGATIVE LOG LIKELIHOOD
  // ---------------------------------------------------------------------------
  double nll = 0;
  
  // Loop over all comparisonss
  for (int index = 0; index < n; index++) { 
    
    // Start with probability zero, and add to this
    double prob = 0.0;
    
    // :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    // Indiv 1 is the potential parent and indiv 2 the offspring
    // :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    
    // :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    // Start of faster but more complicated version
    // :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    
    // // For the max age for individual 1, derive the probabilities for all age for indiv 2
    // int lowest_y1 = c1[index] - max_age;
    // 
    // // max_age_2 depends on c2, y1, and the alpha for the gender of 1
    // int max_age_2 = c2[index] - lowest_y1;
    // // corrected for alpha below
    // if (s1[index] == 1) {
    //   max_age_2 -= alpha_f;
    // } else {
    //   max_age_2 -= alpha_m;
    // }
    // 
    // // Create vector to store probabilities
    // NumericVector probs_a2_given_a1 (max_age_2 + 1);
    // 
    // // Derive the probabilities for every a2, and store in vector
    // for (int a2 = 0; a2 <= max_age_2; a2++) {
    //   // Extract birth year of individual 2
    //   int y2 = c2[index] - a2;
    //   
    //   double prob_a2 = 0.0;
    //   
    //   // If sex is female
    //   if (s1[index] == 1) {
    //     prob_a2 = 1.0 / (N_t0_f * std::pow(r, y2 - t0));
    //     // prob12_2 = 1.0 / (N_t0_f_vector[y2 - t0] + 40);
    //     
    //     // Account for survival of parent i if j was born after c1
    //     if (c1[index] < y2) {
    //       prob_a2 *= std::pow(phi, y2 - c1[index]);
    //       // prob12_2 *= phi_vector[y2 - c1[index]];
    //     }
    //   }
    //   // If sex is male
    //   if (s1[index] == 0) {
    //     // Derive the ERRO in the year before the birth year of the offspring
    //     prob_a2 = 1.0 / (N_t0_m * std::pow(r, y2 - t0));
    //     // prob12_2 = 1.0 / (N_t0_m_vector[y2 - t0] + 40);
    //     
    //     // Account for survival of parent i if j was born after c1 + 1
    //     if (c1[index] < y2) {
    //       prob_a2 *= std::pow(phi, y2 - c1[index]);
    //       // prob12_2 *= phi_vector[y2 - c1[index]];
    //     }
    //   }
    //   // 'Look up table'-version
    //   probs_a2_given_a1[a2] = prob_a2 * age_length_prob_matrix(a2, l2[index]);
    //   // // Conventional version
    //   // probs_a2_given_a1[a2] = prob_a2 * fAgeGivenLengthFast(a2, l2[index], sigma_l,
    //   //                                            probs_sampled_ages, max_age,
    //   //                                            vbgf_l_inf, vbgf_k, vbgf_a0);
    // }
    // 
    // // std::cout << "probs_a2_given_a1: " << probs_a2_given_a1 << std::endl;
    // 
    // int size_vector1 = probs_a2_given_a1.size();
    // 
    // NumericVector temp_1 (size_vector1);
    // for (int a2 = 0; a2 < size_vector1; a2++) {
    //   int a1 = a2 + max_age - max_age_2;
    //   int y1 = c1[index] - a1;
    //   
    //   // Set to zero if parent was not mature in birth year offspring
    //   // Sum using a loop to be able to ignore probabilities where birth year
    //   // of parent is before the earliest birth year according to age of offspring.
    //   double special_sum = 0.0;
    //   for (int i = 0; i <= a2; i++) {
    //     // Only add probability if the combination of a2 and a1 does not imply
    //     // that the parent was older than max_age in the birth year of offspring
    //     if (y1 >= c2[index] - i - max_age) {
    //       special_sum += probs_a2_given_a1[i];
    //     }
    //   }
    //   temp_1[a2] = special_sum;
    //   
    //   // correct for f(a1) = f(a2 + max_age - max_age_2)
    //   // 'Look up table'-version
    //   temp_1[a2] *= age_length_prob_matrix(a1, l1[index]);
    //   // // Conventional version
    //   // temp_1[a2] *= fAgeGivenLengthFast(a1,
    //   //                                   l1[index], sigma_l,
    //   //                                   probs_sampled_ages, max_age,
    //   //                                   vbgf_l_inf, vbgf_k, vbgf_a0);
    // }
    // // std::cout << "temp_1: " << temp_1 << std::endl;
    // // std::cout << "prob new: " << sum(temp_1) << std::endl;
    // 
    // 
    // // Add to the main probability before comparing the other direction
    // prob += sum(temp_1);
    
    // :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    // Start of slower but less complicated version
    // :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    
    // Loop over the ages
    for (int a1 = 0; a1 <= max_age; a1++) {
      // Extract birth year of individual 2
      int y1 = c1[index] - a1;

      double prob12_1 = 0.0;

      // max_age_2 depends on c2, y1, and the alpha for the gender of 1
      int max_age_2 = c2[index] - y1 - 1; // -1 for gestation
      // corrected for alpha below
      if (s1[index] == 1) {
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
        if ((s1[index] == 1) & (a1 + y2 - c1[index] <= max_age)) {
          // Derive the ERRO in the year before the birth year of the offspring
          prob12_2 = 1.0 / (N_t0_f * std::pow(r, y2 - t0 - 1)); // subtract 1 for one year gestation

          // Account for survival of parent i if j was born after c1
          if (c1[index] < y2) {
            prob12_2 *= std::pow(phi, y2 - c1[index]); // don't subtract here
          } else {
            // Account for survival of mom from mating to birth (1 year)
            prob12_2 *= phi;
            }
        }
        if ((s1[index] == 0) & (a1 + y2 - c1[index] - 1 <= max_age)) {
          // Derive the ERRO in the year before the birth year of the offspring
          prob12_2 = 1.0 / (N_t0_m * std::pow(r, y2 - t0 - 1));

          // Account for survival of parent i if j was born after c1 + 1
          if (c1[index] + 1 < y2) {
            prob12_2 *= std::pow(phi, y2 - c1[index] - 1);
          }
        }
        // std::cout << "prob12_2 for a1=" << a1 << " and a2=" << a2 << " is: " << prob12_2 << std::endl;
        // Look up table-version
        prob12_1 += prob12_2 * age_length_prob_matrix(a2, l2[index]);
        // // Fast version
        // prob12_1 += prob12_2 * fAgeGivenLengthFast(a2, l2[index], sigma_l,
        //                                            probs_sampled_ages, max_age,
        //                                            vbgf_l_inf, vbgf_k, vbgf_a0);
        // std::cout << "Corrected prob12_2 for a1=" << a1 << " and a2=" << a2 << " is: " << prob12_2 * fAgeGivenLengthFast(a2, l2[index], sigma_l,
        //                                                                                probs_sampled_ages, max_age,
        //                                                                                vbgf_l_inf, vbgf_k, vbgf_a0) << std::endl;

      }
      // std::cout << "prob12_1 for a1=" << a1 << " is: " << prob12_1 << std::endl;
      // Look up table-version
      prob += prob12_1 * age_length_prob_matrix(a1, l1[index]);
      // // Fast version
      // prob += prob12_1 * fAgeGivenLengthFast(a1, l1[index], sigma_l,
      //                                        probs_sampled_ages, max_age,
      //                                        vbgf_l_inf, vbgf_k, vbgf_a0);
      // std::cout << "Corrected prob12_1 for a1=" << a1 << " is: " << prob12_1 * fAgeGivenLengthFast(a1, l1[index], sigma_l,
      //                                                                    probs_sampled_ages, max_age,
      //                                                                    vbgf_l_inf, vbgf_k, vbgf_a0) << std::endl;
    }
    
    
    
    
    // :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    // Indiv 2 is the potential parent and indiv 1 the offspring
    // :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    
    // :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    // Start of faster but more complicated version
    // :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    
    // // For the max age for individual 2, derive the probabilities for all age for indiv 1
    // int lowest_y2 = c2[index] - max_age;
    // 
    // // max_age_1 depends on c2, y1, and the alpha for the gender of 2
    // int max_age_1 = c1[index] - lowest_y2;
    // // corrected for alpha below
    // if (s2[index] == 1) {
    //   max_age_1 -= alpha_f;
    // } else {
    //   max_age_1 -= alpha_m;
    // }
    // 
    // // Create vector to store probabilities
    // NumericVector probs_a1_given_a2 (max_age_1 + 1);
    // 
    // // Derive the probabilities for every a1, and store in vector
    // for (int a1 = 0; a1 <= max_age_1; a1++) {
    //   // Extract birth year of individual 1
    //   int y1 = c1[index] - a1;
    // 
    //   double prob_a1 = 0.0;
    // 
    //   // If parent sex is female
    //   if (s2[index] == 1) {
    //     prob_a1 = 1.0 / (N_t0_f * std::pow(r, y1 - t0));
    //     // prob12_2 = 1.0 / (N_t0_f_vector[y2 - t0] + 40);
    // 
    //     // Account for survival of parent i if j was born after c1
    //     if (c2[index] < y1) {
    //       prob_a1 *= std::pow(phi, y1 - c2[index]);
    //       // prob12_2 *= phi_vector[y2 - c1[index]];
    //     }
    //   }
    //   // If parent sex is male
    //   if (s2[index] == 0) {
    //     // Derive the ERRO in the year before the birth year of the offspring
    //     prob_a1 = 1.0 / (N_t0_m * std::pow(r, y1 - t0));
    //     // prob12_2 = 1.0 / (N_t0_m_vector[y2 - t0] + 40);
    // 
    //     // Account for survival of parent i if j was born after c1 + 1
    //     if (c2[index] < y1) {
    //       prob_a1 *= std::pow(phi, y1 - c2[index]);
    //       // prob12_2 *= phi_vector[y2 - c1[index]];
    //     }
    //   }
    //   // 'Look up table'-version
    //   probs_a1_given_a2[a1] = prob_a1 * age_length_prob_matrix(a1, l1[index]);
    //   // // Conventional version
    //   // probs_a1_given_a2[a1] = prob_a1 * fAgeGivenLengthFast(a1, l1[index], sigma_l,
    //   //                                                       probs_sampled_ages, max_age,
    //   //                                                       vbgf_l_inf, vbgf_k, vbgf_a0);
    // }
    // 
    // // std::cout << "probs_a1_given_a2: " << probs_a1_given_a2 << std::endl;
    // 
    // int size_vector2 = probs_a1_given_a2.size();
    // 
    // NumericVector temp_2 (size_vector2);
    // for (int a1 = 0; a1 < size_vector2; a1++) {
    //   int a2 = a1 + max_age - max_age_1;
    //   int y2 = c2[index] - a2;
    // 
    //   // Set to zero if parent was not mature in birth year offspring
    //   // Sum using a loop to be able to ignore probabilities where birth year
    //   // of parent is before the earliest birth year according to age of offspring.
    //   double special_sum = 0.0;
    //   for (int i = 0; i <= a1; i++) {
    //     // Only add probability if the combination of a2 and a1 does not imply
    //     // that the parent was older than max_age in the birth year of offspring
    //     if (y2 >= c1[index] - i - max_age) {
    //       special_sum += probs_a1_given_a2[i];
    //     }
    //   }
    //   temp_2[a1] = special_sum;
    // 
    //   // correct for f(a1) = f(a2 + max_age - max_age_2)
    //   // 'Look up table'-version
    //   temp_2[a1] *= age_length_prob_matrix(a2, l2[index]);
    //   // // Conventional version
    //   // temp_2[a1] *= fAgeGivenLengthFast(a2,
    //   //                                   l2[index], sigma_l,
    //   //                                   probs_sampled_ages, max_age,
    //   //                                   vbgf_l_inf, vbgf_k, vbgf_a0);
    // }
    // // std::cout << "temp_2: " << temp_2 << std::endl;
    // // std::cout << "prob new: " << sum(temp_2) << std::endl;
    // 
    // prob += sum(temp_2);
    
    // :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    // Start of slower but less complicated version
    // :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    
    // Loop over the ages, the other way around now
    for (int a2 = 0; a2 <= max_age; a2++) {
      // Extract birth year of individual 2
      int y2 = c2[index] - a2;

      double prob21_1 = 0.0;

      // max_age_1 depends on c1, y2, and the alpha for the gender of 2
      int max_age_1 = c1[index] - y2 - 1; // -1 for gestation
      // corrected for alpha below
      if (s1[index] == 1) {
        max_age_1 -= alpha_f;
      } else {
        max_age_1 -= alpha_m;
      }

      for (int a1 = 0; a1 <= max_age_1; a1++) {

        // Extract birth year of individual 1
        int y1 = c1[index] - a1;

        // :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
        // Now the other direction, i.e., 2 is the parent and 1 is the offspring
        // :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
        double prob21_2 = 0.0;

        // Probability stays zero if offspring was born before maturity of parent
        if ((s2[index] == 1) & (a2 + y1 - c2[index] <= max_age)) {
          prob21_2 = 1.0 / (N_t0_f * std::pow(r, y1 - t0 - 1)); // subtract 1 for one year gestation

          // Account for survival of parent i if j was born after c1
          if (c2[index] < y1) {
            prob21_2 *=  std::pow(phi, y1 - c2[index]);
          } else {
            // Account for survival of mom from mating to birth (1 year)
            prob21_2 *= phi;
          }
        }
        if ((s2[index] == 0) & (a2 + y1 - c2[index] - 1<= max_age)) {
          // Derive the ERRO in the year before the birth year of the offspring
          prob21_2 = 1.0 / (N_t0_m * std::pow(r, y1 - t0 - 1));

          // Account for survival of parent i if j was born after c1 + 1
          if (c2[index] + 1 < y1) {
            prob21_2 *=  std::pow(phi, y1 - c2[index] - 1);
          }
        }
        // Look up table version
        prob21_1 += prob21_2 * age_length_prob_matrix(a1, l1[index]);
        // // Fast version
        // prob21_1 += prob21_2 * fAgeGivenLengthFast(a1, l1[index], sigma_l,
        //                                            probs_sampled_ages, max_age,
        //                                            vbgf_l_inf, vbgf_k, vbgf_a0);
        // std::cout << "Corrected prob21_2 for a2=" << a2 << " and a1=" << a1 << " is: " << prob21_2 * fAgeGivenLengthFast(a1, l1[index], sigma_l,
        //                                                                                probs_sampled_ages, max_age,
        //                                                                                vbgf_l_inf, vbgf_k, vbgf_a0) << std::endl;
      }
      // Look up table version
      prob += prob21_1 * age_length_prob_matrix(a2, l2[index]);
      // // Fast version
      // prob += prob21_1 * fAgeGivenLengthFast(a2, l2[index], sigma_l,
      //                                        probs_sampled_ages, max_age,
      //                                        vbgf_l_inf, vbgf_k, vbgf_a0);
      // std::cout << "corrected prob21_1 for a2=" << a2 << ": " << prob21_1 * fAgeGivenLengthFast(a2, l2[index], sigma_l,
      //                                                   probs_sampled_ages, max_age,
      //                                                   vbgf_l_inf, vbgf_k, vbgf_a0)<< std::endl;
    }

    // std::cout << "prob final : " << prob << std::endl;
    
    if (kinship[index] != 1) {
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
