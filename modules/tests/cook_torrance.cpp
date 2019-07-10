#include <stdio.h>
#include <iostream>
#include <array>
#include "mio/altro/types.hpp"
#include "mio/math/math.hpp"
#include "mio/altro/rand.hpp"
#include "mio/altro/algorithm.hpp"
#include "cookTorrance.h"
#include "levmar.h"


/*** Beckmann Distribution ***/

void LMBeckmann(double *p, double *x, int m, int n, void *user_data){
  std::vector<double> &phi_vec = *static_cast<std::vector<double>*>(user_data);
  STD_INVALID_ARG_E(phi_vec.size() == n)
  printf("LMBeckmann() : p[0] = %f\n", p[0]);

  for(size_t i = 0; i < n; ++i)
    x[i] = BeckmannDistribution(p[0], phi_vec[i]);
}


void LMBeckmannJac(double *p, double *jac, int m, int n, void *user_data){
  std::vector<double> &phi_vec = *static_cast<std::vector<double>*>(user_data);
  STD_INVALID_ARG_E(phi_vec.size() == n)

  // Solve jacobian for each training sample
  for(size_t i = 0, j = 0; i < n; ++i)
    jac[j++] = BeckmannDistribution_d_dmu(p[0], phi_vec[i]);
}


void LevMarTestBeckmann(){
  // Generate some data
  const size_t num_train_sample = 100;
  std::vector<double> y_data_vec(num_train_sample), phi_vec(num_train_sample);
  mio::FillRange(phi_vec, sm::DegToRad(1.0), sm::DegToRad(89.0));
  mio::CRandNum<double> rand_num(0.0, 0.01); //mean, stddev
  const double mu = 0.5; // Try to estimate this parameter

  for(size_t i = 0; i < num_train_sample; ++i){
    y_data_vec[i] = BeckmannDistribution(mu, phi_vec[i]);
    y_data_vec[i] += rand_num.GetRandNum(); // Add some noise
  }

  // Estimate parameters with levmar
  const size_t num_param = 1,
               max_iter = 1000;
  std::array<double, num_param> p = {1.0}; // Initial value for mu parameter
  double opts[LM_OPTS_SZ], info[LM_INFO_SZ];
  const bool default_opt = true;
  opts[0] = 1E-40;
  opts[1] = 1E-40;
  opts[2] = 1E-40;
  opts[3] = 1E-40;

  int ret = dlevmar_der( LMBeckmann, LMBeckmannJac, p.data(), y_data_vec.data(), num_param, num_train_sample, max_iter,
                         default_opt ? NULL : opts, info, NULL, NULL, static_cast<void*>(&phi_vec) );
  printf("Levenberg-Marquardt returned in %g iter, reason %g, sumsq %g [%g]\n", info[5], info[6], info[1], info[0]);
  printf("Best fit parameter: %.7g\n", p[0]);
}


/*** Fresnel ***/

struct FresnelData{
  std::vector<double> theta_vec;
  double err_a, err_gamma;
};


void LMFresnel(double *p, double *x, int m, int n, void *user_data){
  FresnelData &fresnel_data = *static_cast<FresnelData*>(user_data);
  STD_INVALID_ARG_E(fresnel_data.theta_vec.size() == n)
  printf("LMFresnel() : p[0] = %f, p[1] = %f\n", p[0], p[1]);

  for(size_t i = 0; i < n; ++i)
    x[i] = FresnelSimple(p[0], p[1], fresnel_data.theta_vec[i]);
    //x[i] = FresnelMetalApprox(p[0], p[1], fresnel_data.theta_vec[i], fresnel_data.err_a, fresnel_data.err_gamma);
}


void LMFresnelJac(double *p, double *jac, int m, int n, void *user_data){
  FresnelData &fresnel_data = *static_cast<FresnelData*>(user_data);
  STD_INVALID_ARG_E(fresnel_data.theta_vec.size() == n)

  // Solve jacobian for each training sample
  double dF_dn, dF_dk;
  for(size_t i = 0, j = 0; i < n; ++i){
    FresnelSimple_deriv(p[0], p[1], fresnel_data.theta_vec[i], dF_dn, dF_dk);
    //FresnelMetalApprox_deriv(p[0], p[1], fresnel_data.theta_vec[i], dF_dn, dF_dk);
    jac[j++] = dF_dn;
    jac[j++] = dF_dk;
  }
}


// Try to estimate the n and k parameters of the Fresnel equation for reflectivity
void LevMarTestFresnel(){
  // Generate some data
  const size_t num_train_sample = 25;
  std::vector<double> y_data_vec(num_train_sample);
  FresnelData fresnel_data; // This struct contains static variable values (values that will not be estimated)
  fresnel_data.theta_vec.resize(num_train_sample);
  mio::FillRange(fresnel_data.theta_vec, sm::DegToRad(1.0), sm::DegToRad(89.0));
  //fresnel_data.err_a = 2.0*p_actual[0]; // err_a and err_gamma are parameters of FresnelMetalApprox()
  //fresnel_data.err_gamma = 7;
  mio::CRandNum<double> rand_num(0.0, 0.01); //mean, stddev
  // These are the parameters we will try to estimate
  const double n = 1.5,
               k = 5.0;

  for(size_t i = 0; i < num_train_sample; ++i){
    y_data_vec[i] = FresnelSimple(n, k, fresnel_data.theta_vec[i]);
    //y_data_vec[i] = FresnelMetalApprox(n, k, fresnel_data.theta_vec[i], fresnel_data.err_a, fresnel_data.err_gamma);
    y_data_vec[i] += rand_num.GetRandNum(); // Add some noise
  }

  // Estimate parameters with levmar
  const size_t num_param = 2,
               max_iter = 1000;
  std::array<double, num_param> p = {1.0, 1.0}; // Initial values for n and k
  double opts[LM_OPTS_SZ], info[LM_INFO_SZ];
  const bool default_opt = true;
  opts[0] = 1E-40;
  opts[1] = 1E-40;
  opts[2] = 1E-40;
  opts[3] = 1E-40;

  int ret = dlevmar_der( LMFresnel, LMFresnelJac, p.data(), y_data_vec.data(), num_param, num_train_sample, max_iter,
                         default_opt ? NULL : opts, info, NULL, NULL, static_cast<void*>(&fresnel_data) );
  printf("Levenberg-Marquardt returned in %g iter, reason %g, sumsq %g [%g]\n", info[5], info[6], info[1], info[0]);
  printf("Best fit parameters: %.7g %.7g\n", p[0], p[1]);
}


/*** Cook Torrance ***/

struct CookTorranceData{
  vec3d_t H, L, N, V;
};


void LMCookTorrance(double *p, double *x, int m, int n, void *user_data){
  std::vector<CookTorranceData> &ctd_vec = *static_cast< std::vector<CookTorranceData>* >(user_data);
  STD_INVALID_ARG_E(ctd_vec.size() == n)
  printf("LMCookTorrance() : p[0] = %f, p[1] = %f, p[2] = %f\n", p[0], p[1], p[2]);

  for(size_t i = 0; i < n; ++i)
    x[i] = CookTorrance(ctd_vec[i].H, ctd_vec[i].L, ctd_vec[i].N, ctd_vec[i].V, p[0], p[1], p[2]);
}


void LMCookTorranceJac(double *p, double *jac, int m, int n, void *user_data){
  std::vector<CookTorranceData> &ctd_vec = *static_cast< std::vector<CookTorranceData>* >(user_data);
  STD_INVALID_ARG_E(ctd_vec.size() == n)

  // Solve jacobian for each training sample
  double dR_dn, dR_dk, dR_dmu;
  for(size_t i = 0, j = 0; i < n; ++i){
    CookTorrance_deriv(ctd_vec[i].H, ctd_vec[i].L, ctd_vec[i].N, ctd_vec[i].V, p[0], p[1], p[2], dR_dn, dR_dk, dR_dmu);
    jac[j++] = dR_dn;
    jac[j++] = dR_dk;
    jac[j++] = dR_dmu;
  }
}


void LevMarTestCookTorrance(){
  // Generate some data
  const size_t num_train_sample = 1000;
  std::vector<CookTorranceData> ctd_vec(num_train_sample);
  std::vector<double> y_data_vec(num_train_sample), theta_vec(num_train_sample);
  const vec3d_t V(0.0, 0.0, 1.0),
                L( 0.0, std::sin( sm::DegToRad(45) ), std::cos( sm::DegToRad(45) ) );
  const vec3d_t H = sm::HalfwayVectorNorm3(L, V);
  mio::FillRange(theta_vec, sm::DegToRad(1.0), sm::DegToRad(89.0));
  // These are the parameters we will try to estimate
  const double n = 1.5,
               k = 5.0,
               mu = 0.5;

  for(size_t i = 0; i < num_train_sample; ++i){
    ctd_vec[i].H = H;
    ctd_vec[i].L = L;
    //simulate flat plate that is rotation, camera and illuminator are stationary
    ctd_vec[i].N = vec3d_t( 0.0, std::sin(theta_vec[i]), std::cos(theta_vec[i]) );
    ctd_vec[i].V = V;
  }
  for(size_t i = 0; i < num_train_sample; ++i)
    y_data_vec[i] = CookTorrance(ctd_vec[i].H, ctd_vec[i].L, ctd_vec[i].N, ctd_vec[i].V, n, k, mu);

  // Estimate parameters with levmar
  const size_t num_param = 3,
               max_iter = 1000;
  std::array<double, num_param> p = {1.45, 4.95, 0.6}; //initial values for n, k, mu
  double opts[LM_OPTS_SZ], info[LM_INFO_SZ];
  const bool default_opt = true;
  opts[0] = 1E-10;
  opts[1] = 1E-15;
  opts[2] = 1E-15;
  opts[3] = 1E-20;

  int ret = dlevmar_der( LMCookTorrance, LMCookTorranceJac, p.data(), y_data_vec.data(), num_param, num_train_sample, max_iter,
                         default_opt ? NULL : opts, info, NULL, NULL, static_cast<void*>(&ctd_vec) );
  printf("Levenberg-Marquardt returned in %g iter, reason %g, sumsq %g [%g]\n", info[5], info[6], info[1], info[0]);
  printf("Best fit parameters: %.7g %.7g %.7g\n", p[0], p[1], p[2]);
}


int main(int argc, char *argv[]){
  LevMarTestBeckmann();
  //LevMarTestFresnel();
  //LevMarTestCookTorrance();

  return 0;
}

