#include <Rcpp.h>
using namespace Rcpp;


//============================================================================//
//' R wrapper for the \code{CalcMatchFactor} C++ function.
//'
//' @description
//'   Call the \code{CalcMatchFactor} C++ function.
//'
//' @param ms_u,ms_l
//'   A numeric matrix. Mass values, weighted and original intensities (see the
//'   \code{\link{LibrarySearch}} function).
//' @param is_identity
//'   A logical value. If \code{TRUE}, the Identity algorithm is used to
//'   calculate the match factor. If \code{FALSE}, the Similarity algorithm is
//'   used.
//' @param min_mz,max_mz
//'   An integer value. Boundaries of the m/z range.
//' @param is_reverse_search
//'   A logical value. If \code{TRUE}, reverse search is performed and RMF is
//'   returned.
//'
//' @details
//'   All m/z values should be unique and sorted in the accessing order.
//'   Original intensities of mass spectral peaks should be normalized; the base
//'   peak should have intensity of 999. Any checks of the input data are not
//'   performed (i.e., if m/z values are not sorted in the ascending order, the
//'   calculated match factor is going to be incorrect).
//'
//' @return
//'   A numeric value. The match factor calculated between a pair of mass
//'   spectra.
//'
//' @noRd
//'
// [[Rcpp::export(".CalcMatchFactorC")]]

double CalcMatchFactor(NumericMatrix ms_u,
                       NumericMatrix ms_l,
                       bool is_identity,
                       int min_mz = -1,
                       int max_mz = -1,
                       bool is_reverse_search = false,
                       bool is_debugging = false) {

  // tolerance
  const double tol = sqrt(DBL_EPSILON);
  // Rprintf("tol = %e\n", tol);

  // counters
  int i_u = 0;
  int i_l = 0;
  const int n_peaks_u = ms_u.nrow();
  const int n_peaks_l = ms_l.nrow();

  // column indexes
  const int mz = 0;
  const int intst = 1;
  const int w1 = 2; // weighted intensity for the cosine similarity
  const int w2 = 3; // weighted intensity for the the second term

  // finding starting m/z value
  if (min_mz > 0) {
    while (tol < min_mz - ms_u(i_u, mz) && i_u < n_peaks_u) i_u++;
    while (tol < min_mz - ms_l(i_l, mz) && i_l < n_peaks_l) i_l++;
  } else if (ms_u(i_u, mz) - ms_l(i_l, mz) > tol) {
    while (ms_u(i_u, mz) - ms_l(i_l, mz) > tol && i_l < n_peaks_l) i_l++;
  } else {
    while (tol < ms_l(i_l, mz) - ms_u(i_u, mz) && i_u < n_peaks_u) i_u++;
  }

  // if (is_debugging) {
  //   Rprintf("min_mz = %d; i_u = %d (out of %d); i_l = %d (out of %d)\n",
  //           min_mz, i_u + 1, n_peaks_u, i_l + 1, n_peaks_l);
  // }

  // at least one mass spectrum does not contain peaks in the set m/z range
  if (i_u >= n_peaks_u || i_l >= n_peaks_l) {
    return 0;
  }

  // Identity algorithm:
  //   mf = (term1 * n1 + term2 * n2) / (n1 + n2)
  //   term1 = cos_a^2 = sum_ul^2 / (sum_uu * sum_ll)
  //   term2 = r = sum_rm / sum_m

  // Similarity algorithm:
  //   mf = cos_a^2 = sum_ul^2 / (sum_uu * sum_ll)

  double sum_ul = 0;
  double sum_uu = 0;
  double sum_ll = 0;
  double sum_rm = 0;
  int sum_m = 0;

  int n1 = 0; // the number of peaks in both library and unknown spectrum
  int n2 = 0; // the number of ratios of intensities

  double prev_w2_u = 0;
  double prev_w2_l = 0;

  bool should_calc_term2 = false;

  bool should_continue = true;
  while (should_continue) {

    enum {kBoth, kUnknown, kLibrary}; // Which spectra do contain this peak?
    int calc_type;
    if (i_u >= n_peaks_u) {
      calc_type = kLibrary;
    } else if (i_l >= n_peaks_l) {
      calc_type = kUnknown;
    } else if (ms_u(i_u, mz) - ms_l(i_l, mz) > tol) {
      calc_type = kLibrary;
    } else if (tol < ms_l(i_l, mz) - ms_u(i_u, mz)) {
      calc_type = kUnknown;
    } else {
      if (is_reverse_search && ms_l(i_l, intst) < tol) {
        calc_type = kUnknown;
      } else {
        calc_type = kBoth;
      }
    }

    if (calc_type == kBoth) {
      if (ms_u(i_u, intst) > 1 + tol || ms_l(i_l, intst) > 1 + tol) {
        // Rprintf("mz = %d\n", ms_u(i_u, mz));
        sum_ul += ms_u(i_u, w1) * ms_l(i_l, w1);
        sum_uu += ms_u(i_u, w1) * ms_u(i_u, w1);
        sum_ll += ms_l(i_l, w1) * ms_l(i_l, w1);
        if (ms_u(i_u, intst) > tol && ms_l(i_l, intst) > tol) {
          n1++;
        }

        if (is_identity) {
          if (should_calc_term2) {
            double r_num = ms_u(i_u, w2) * prev_w2_l;
            double r_denom = prev_w2_u * ms_l(i_l, w2);
            if (r_num > tol && r_denom > tol) {
              if (r_denom > r_num) sum_rm += ms_u(i_u, mz) * r_num / r_denom;
              else sum_rm += ms_u(i_u, mz) * r_denom / r_num;
              sum_m += ms_u(i_u, mz);
              n2++;
            }
          }
          should_calc_term2 = true;
          prev_w2_u = ms_u(i_u, w2);
          prev_w2_l = ms_l(i_l, w2);
        }
      }
      i_u++;
      i_l++;

    } else if (calc_type == kUnknown) {
      if (ms_u(i_u, intst) > 1) {
        should_calc_term2 = false;
        if (!is_reverse_search) {
          sum_uu += ms_u(i_u, w1) * ms_u(i_u, w1);
        }
      }
      i_u++;

    } else { // i.e., 'calc_type == kLibrary'
      if (ms_l(i_l, intst) > 1) {
        should_calc_term2 = false;
        sum_ll += ms_l(i_l, w1) * ms_l(i_l, w1);
      }
      i_l++;
    }

    if (i_u >= n_peaks_u && i_l >= n_peaks_l) {
      should_continue = false;
    }

    if (max_mz > 0) {
      if ((i_u >= n_peaks_u || ms_u(i_u, mz) - max_mz > tol) &&
          (i_l >= n_peaks_l || ms_l(i_l, mz) - max_mz > tol)){
        should_continue = false;
      }
    }
  }

  if (n1 == 0) {
    return 0;
  }

  // if (is_debugging) {
  //   double term1 = sum_ul * sum_ul / (sum_uu * sum_ll);
  //   Rprintf("term1 = %.1lf; n1 = %d", 1000.0 * term1 - 0.5, n1);
  //   Rprintf(" (sum_ul = %.1lf; sum_uu = %.1lf; sum_ll = %.1lf)",
  //           sum_ul, sum_uu, sum_ll);
  //   if (sum_m > 0) {
  //     double term2 = sum_rm / sum_m;
  //     Rprintf("\n");
  //     Rprintf("term2 = %.1lf; n2 = %d", 1000.0 * term2 - 0.5, n2);
  //     Rprintf("; mf = %.1lf\n",
  //             1000.0 * (term1 * n1 + term2 * n2) / (n1 + n2) - 0.5);
  //   } else {
  //     Rprintf("; mf = %.1lf\n", 1000.0 * term1 - 0.5);
  //   }
  // }

  double term1 = sum_ul * sum_ul / (sum_uu * sum_ll); // i.e., cos_a^2;
  if (sum_m > 0) { // if 'sum_m = 0', then 'n2 = 0' and 'term2 = 0'
    return 1000.0 * (term1 * n1 + (sum_rm / sum_m) * n2) / (n1 + n2) - 0.5;
  } else {
    return 1000.0 * term1 - 0.5;
  }
}




