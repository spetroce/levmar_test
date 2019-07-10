#include <vector>
#include <array>
#include <cmath>
#include "levmar/levmar.h"
#include "mio/altro/rand.h"
#include "mio/altro/error.h"

#ifndef LM_DBL_PREC
#error Demo program assumes that levmar has been compiled with double precision, see LM_DBL_PREC!
#endif

#define NUM_EXP_PARAM 3

/*
y = f(x) = a*exp(-bx) + c = p[0]*exp(-p[1]*x) + p[2]
traditionally:
  m is the number of training examples
  n is the number of features (x0 x1 x2... xn)
  y is the column vector of all training labels with dimensionality m (ie. f(x))
  x is the column vector of all the feature inputs with dimensionality n+1 (ie. the known variables of f(x))
    so x and y are vectors known as the training data
  θ is the parameter column vector with dimensionality n + 1 (ie. THE THINGS WE ARE TRYING TO ESTIMATE)

levmar uses slightly different notation:
  m is the parameter vector dimension. For the problem here, m = 3.
  n is the measurement vector dimension. (english translation: the number of training examples)
  x is the "measurement vector". (english translation: this is y, aka f(x))
  p is the parameter vector with dimension m.
*/


void ExpFuncGenerateData(const size_t n, std::vector<double> &y_data_vec,
    std::vector<double> &x_data_vec, const std::array<double, NUM_EXP_PARAM> p,
    const bool with_rand_noise) {
  y_data_vec.resize(n);
  x_data_vec.resize(n);
  const double mean = 0.0, std_dev = 0.01;
  mio::CRandNum<double> rand_num(mean, std_dev);
  for (size_t i = 0; i < n; ++i) {
    double x = i + 1;
    x_data_vec[i] = x;
    y_data_vec[i] = p[0]*std::exp(-p[1]*x) + p[2];
    if (with_rand_noise) {
      y_data_vec[i] += rand_num.GetRandNum();
    }
  }
}


void ExpFunc(double *p, double *x, int m, int n, void *user_data) {
  std::vector<double> &x_data_vec = *static_cast< std::vector<double>* >(user_data);
  STD_INVALID_ARG_E(x_data_vec.size() == n)
  STD_INVALID_ARG_E(m == NUM_EXP_PARAM)
  for (size_t i = 0; i < n; ++i)
    x[i] = p[0]*std::exp(-p[1]*x_data_vec[i]) + p[2]; // f(x) = a*e^(-bx)+c
}


void JacExpFunc(double *p, double *jac, int m, int n, void *user_data) {
  std::vector<double> &x_data_vec = *static_cast< std::vector<double>* >(user_data);
  STD_INVALID_ARG_E(x_data_vec.size() == n)
  STD_INVALID_ARG_E(m == NUM_EXP_PARAM)
  double temp;
  // solve jacobian for each training sample
  for (size_t i = 0, j = 0; i < n; ++i) {
    temp = std::exp(-p[1]*x_data_vec[i]);
    jac[j++] = temp; // df_da = e^{-bx}
    jac[j++] = -p[0]*x_data_vec[i]*temp; // df_db = -a*x*e^{-bx}
    jac[j++] = 1.0; // df_dc = 1
  }
}


void ExpLevMarTest() {
  const size_t n = 1000, // number of training samples
               m = NUM_EXP_PARAM; // number of parameters to estimate
  double opts[LM_OPTS_SZ], // opts[4] = minim. options [\mu, \epsilon1, \epsilon2, \epsilon3]
         info[LM_INFO_SZ];
  const std::array<double, NUM_EXP_PARAM> p_actual = {5.0, 0.1, 1.0};
  std::array<double, NUM_EXP_PARAM> p = {0.0, 0.0, 0.0};
  std::vector<double> y_data_vec, x_data_vec;
  ExpFuncGenerateData(n, y_data_vec, x_data_vec, p_actual, true);
  const bool default_opt = true;
  opts[0] = LM_INIT_MU; // scale factor for initial \mu
  opts[1] = 1E-15; // stopping threshold for ||J^T e||_inf
  opts[2] = 1E-15; // stopping threshold for ||Dp||_2
  opts[3] = 1E-20; // stopping threshold for ||e||_2
  void *adata = static_cast<void*>(&x_data_vec);
  int ret = dlevmar_der(ExpFunc, JacExpFunc, p.data(), y_data_vec.data(), m, n, 10000,
                        default_opt ? NULL : opts, info, NULL, NULL, adata);
  printf("Levenberg-Marquardt returned in %g iter, reason %g, sumsq %g [%g]\n",
    info[5], info[6], info[1], info[0]);
  printf("Best fit parameters: %.7g %.7g %.7g\n", p[0], p[1], p[2]);
  printf("   Values should be: %.7g %.7g %.7g\n", p_actual[0], p_actual[1], p_actual[2]);
}


int main(int argc, char *argv[]) {
  ExpLevMarTest();
  return 0;
}
