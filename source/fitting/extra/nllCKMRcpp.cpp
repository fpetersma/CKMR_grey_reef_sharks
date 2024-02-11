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
  const double vbgf_t0 = dat["vbgf_t0"];

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
  
  // const double g0 = exp(double(par["logit_g0"])) / 
  //   (exp(double(par["logit_g0"])) + 1.0);
  // const double beta_r = exp(double(par["log_beta_r"]));
  // const double sd_r = exp(double(par["log_sd_r"]));
  // const double mu_s = exp(double(par["log_mu_s"]));
  // const double sd_s = exp(double(par["log_sd_s"]));
  // ===========================================================================
  
  // ===========================================================================
  // 3. DERIVE THE NEGATIVE LOG LIKELIHOOD
  // ---------------------------------------------------------------------------
  double nll = 0;
  
  for (int index = 0; index < n; index++) { 
    // Derive potential ages for parent
    IntegerVector ages_parent = seq(0, max_age);
    int n_ages_parent = ages_parent.size();
    
    double output_parent_level = 0;
    
    // Loop over ages for parent
    for (int index_age_parent = 0; index_age_parent < n_ages_parent; index_age_parent++) {
      int a1 = ages_parent[index_age_parent];
      int y1 = c1[index] - a1;
      
      int max_age_offspring = 0; // place holder to allow compiling
      // Set maximum age for the offspring
      if (s1[index] == "F") {
        max_age_offspring = c2[index] - c1[index] - alpha_f + a1 - 1;
      } 
      if (s1[index] == "M") {
        max_age_offspring = c2[index] - c1[index] - alpha_m + a1 - 1;
      }
      
      // Probability equals zero if max_age_offspring is negative
      if (max_age_offspring < 0) {
        output_parent_level += 0;
      } else {
        // Derive potential ages for offspring
        IntegerVector ages_offspring = seq(0, max_age_offspring);
        int n_ages_offspring = ages_offspring.size();
        
        // Create output for prob for every potential age of the offspring
        double output_offspring_level = 0;
        
        // std::cout << "max_age_offspring: " << max_age_offspring << std::endl;
        
        if (s1[index] == "F") {
          // Loop over ages for offspring
          for (int index_age_offspring = 0; index_age_offspring < n_ages_offspring / 2; index_age_offspring++) {
            
            // testing!
            int correct_offspring_age_index = index_age_offspring * 2;
            
            int a2 = ages_offspring[correct_offspring_age_index];
            int y2 =  c2[index] - a2;
            
            double out = 0.0; // place holder to allow compiling
            
            if (s1[index] == "F") {
              // Derive the ERRO in the year before the birth year of the offspring
              out = 1.0 / (N_t0_f * pow(r, y2 - 1 - t0));
              // std::cout << "Female abundance in year " << y2 << " is: " << N_t0_f * pow(r, y2 - 1 - t0) << std::endl;
              // out = 1.0 / (N_t0_f * pow(1.0 + r, y2 - 1 - t0)); # old, no need to add the 1 to r I think
              
              // Account for survival of parent i if j was born after c1 
              if (c1[index] < y2) {
                out *=  pow(phi, c1[index] - y2);
              }
            } 
            // if (s1[index] == "M") {
            //   // Derive the ERRO in the year before the birth year of the offspring
            //   out = 1.0 / (N_t0_m * pow(r, y2 - 1 - t0));
            //   // out = 1.0 / (N_t0_m * pow(1.0 + r, y2 - 1 - t0)); # old, no need to add the 1 to r I think
            //   
            //   // Account for survival of parent i if j was born after c1 + 1
            //   if (c1[index] < y2 - 1) {
            //     out *=  pow(phi, c1[index] - y2 - 1);
            //   }
            // } 
            
            // Find expected length baharased on age and VBGF parameters
            double l2_exp = vbgf_l_inf * (1 - exp(-vbgf_k * (a2 - vbgf_t0)));
            
            // Multiply above by probability density of offspring age
            out *= R::dnorm(l2[index], l2_exp, sigma_vbgf, false);
            
            // Add to output_offspring_level
            output_offspring_level += out;
            
            // std::cout << "out: " << out << std::endl;
          }
        } else if (s1[index] == "M") {
          // Loop over ages for offspring
          for (int index_age_offspring = 0; index_age_offspring < n_ages_offspring; index_age_offspring++) {
            
            int a2 = ages_offspring[index_age_offspring];
            int y2 =  c2[index] - a2;
            
            double out = 0.0; // place holder to allow compiling
            
            // if (s1[index] == "F") {
            //   // Derive the ERRO in the year before the birth year of the offspring
            //   out = 1.0 / (N_t0_f * pow(r, y2 - 1 - t0));
            //   std::cout << "Female abundance in year " << y2 << " is: " << N_t0_f * pow(r, y2 - 1 - t0) << std::endl;
            //   // out = 1.0 / (N_t0_f * pow(1.0 + r, y2 - 1 - t0)); # old, no need to add the 1 to r I think
            //   
            //   // Account for survival of parent i if j was born after c1 
            //   if (c1[index] < y2) {
            //     out *=  pow(phi, c1[index] - y2);
            //   }
            // } 
            if (s1[index] == "M") {
              // Derive the ERRO in the year before the birth year of the offspring
              out = 1.0 / (N_t0_m * pow(r, y2 - 1 - t0));
              // out = 1.0 / (N_t0_m * pow(1.0 + r, y2 - 1 - t0)); # old, no need to add the 1 to r I think
              
              // Account for survival of parent i if j was born after c1 + 1
              if (c1[index] < y2 - 1) {
                out *=  pow(phi, c1[index] - y2 - 1);
              }
            } 
            
            // Find expected length baharased on age and VBGF parameters
            double l2_exp = vbgf_l_inf * (1 - exp(-vbgf_k * (a2 - vbgf_t0)));
            
            // Multiply above by probability density of offspring age
            out *= R::dnorm(l2[index], l2_exp, sigma_vbgf, false);
            
            // Add to output_offspring_level
            output_offspring_level += out;
            
            // std::cout << "out: " << out << std::endl;
          }
        }
      
        
        // // Sum the probabilities for every potential offspring age
        // double output_offspring_level = output_offspring_level.sum(); 
        
        // Find expected length based on age and VBGF parameters
        double l1_exp = vbgf_l_inf * (1 - exp(-vbgf_k * (a1 - vbgf_t0)));
        
        // Multiply above by probability density of parent age
        output_offspring_level *= R::dnorm(l1[index], l1_exp, sigma_vbgf, false);
        
        // Add to output_parent_level
        output_parent_level += output_offspring_level;
        
        // std::cout << "output_offspring_level: " << output_offspring_level << std::endl;
      }
    }
    
    // Sum the probabilities for every potential offspring age
    double prob = output_parent_level;
    
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
