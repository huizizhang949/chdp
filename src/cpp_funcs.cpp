#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix compute_probability_cpp(NumericMatrix Y_d, NumericVector Beta,
                                      NumericMatrix mu, NumericMatrix phi,
                                      NumericVector Q_J, NumericVector t_d,
                                      NumericVector t_star_J, NumericVector sigma_star_2_J) {

  int G = Y_d.nrow(); // Number of rows in Y_d (G)
  int C = Y_d.ncol(); // Number of columns in Y_d (C)
  int J = mu.nrow(); // Number of rows in mu and phi (J)

  // Create a matrix to store the result (C x J)
  NumericMatrix probabilities(C, J);

  // Loop over each column (1:C)
  for (int c = 0; c < C; ++c) {

    // Initialize a vector to store log probabilities (LP) for the current column
    NumericVector LP(J);

    // Loop over each j (1:J)
    for (int j = 0; j < J; ++j) {

      // Compute the first part: sum of dnbinom values
      double sum_dnbinom = 0;
      for (int g = 0; g < G; ++g) {
        double mu_val = mu(j, g) * Beta[c];
        double phi_val = phi(j, g);
        double y_val = Y_d(g, c);

        // Compute the negative binomial log probability for each g
        sum_dnbinom += R::dnbinom_mu(y_val, phi_val, mu_val, true);
      }

      // Compute the second part: Gaussian term
      double logQ = log(Q_J[j]);
      double t_diff = t_d[c] - t_star_J[j];
      double gaussian_term = -pow(t_diff, 2) / (2 * sigma_star_2_J[j]);

      // Compute the total log probability for this j
      LP[j] = sum_dnbinom + logQ + gaussian_term;
    }

    // Compute the normalizing constant (nc) as the negative max(LP)
    double nc = -max(LP);

    // Compute the probabilities (P) and normalize
    NumericVector P(J);
    double sum_exp_LP = 0;

    for (int j = 0; j < J; ++j) {
      P[j] = exp(LP[j] + nc);
      sum_exp_LP += P[j];
    }

    // Normalize the probabilities and store them in the result matrix
    for (int j = 0; j < J; ++j) {
      probabilities(c, j) = P[j] / sum_exp_LP;
    }
  }

  return probabilities;
}
