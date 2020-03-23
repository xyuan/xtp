/*
 *            Copyright 2009-2019 The VOTCA Development Team
 *                       (http://www.votca.org)
 *
 *      Licensed under the Apache License, Version 2.0 (the "License")
 *
 * You may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *              http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

#pragma once
#ifndef _VOTCA_XTP_GW_H
#define _VOTCA_XTP_GW_H

#include "votca/xtp/logger.h"
#include <fstream>
#include <votca/xtp/orbitals.h>
#include <votca/xtp/rpa.h>
#include <votca/xtp/sigma_base.h>
#include <votca/xtp/threecenter.h>
namespace votca {
namespace xtp {

class GW {
 public:
  GW(Logger& log, TCMatrix_gwbse& Mmn, const Eigen::MatrixXd& vxc,
     const Eigen::VectorXd& dft_energies)
      : _log(log),
        _Mmn(Mmn),
        _vxc(vxc),
        _dft_energies(dft_energies),
        _rpa(log, Mmn){};

  struct options {
    Index homo;
    Index qpmin;
    Index qpmax;
    Index rpamin;
    Index rpamax;
    double eta = 1e-3;
    double g_sc_limit = 1e-5;
    Index order = 12;  // only needed for complex integration sigma default:12
    double alpha =
        0.1;  // only needed for complex integration of sigma: default 0.1
    Index g_sc_max_iterations = 50;
    double gw_sc_limit = 1e-5;
    Index gw_sc_max_iterations = 50;
    double shift = 0;
    double ScaHFX = 0.0;
    std::string sigma_integration = "ppm";
    Index reset_3c = 5;  // how often the 3c integrals in iterate should be
                         // rebuild
    std::string qp_solver = "grid";
    double qp_solver_alpha = 0.75;
    Index qp_grid_steps = 601;       // Number of grid points
    double qp_grid_spacing = 0.005;  // Spacing of grid points in Ha
    std::string quadrature_scheme;  // Kind of Gaussian-quadrature scheme to use
    Index qp_training_points =
        5;  // Number of starting training point to use for Kernel Regression method
    double qp_spread =
        1.0;  // Spread of laplacian kernel for Kernel Regression method
    double qp_mae_tol = 5e-4; //To prove that we learn the curve we use the MAE 
    double qp_fixedpoint_tol = 1e-3; // Tolerance to be reach when usign Atkin method for fixed point solver
    double qp_grid_hartree = 0.5;  // How many hartree on the left and right of
                                   // the pre-shooted center of the grid
    std::string qp_test_points = "both"; //Other options are even and odd
  };

  void configure(const options& opt);

  Eigen::VectorXd getGWAResults() const;
  // Calculates the diagonal elements up to self consistency
  void CalculateGWPerturbation();

  // Calculated offdiagonal elements as well
  void CalculateHQP();

  Eigen::MatrixXd getHQP() const;

  // Diagonalize QP particle Hamiltonian
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> DiagonalizeQPHamiltonian()
      const;

  void PlotSigma(std::string filename, Index steps, double spacing,
                 std::string states) const;

 private:
  Index _qptotal;

  Eigen::MatrixXd _Sigma_x;
  Eigen::MatrixXd _Sigma_c;
  Eigen::MatrixXd _Sigma_c_i;

  options _opt;

  std::unique_ptr<Sigma_base> _sigma = nullptr;
  Logger& _log;
  TCMatrix_gwbse& _Mmn;
  const Eigen::MatrixXd& _vxc;
  const Eigen::VectorXd& _dft_energies;

  RPA _rpa;
  // small class which calculates f(w) with and df/dw(w)
  // f=Sigma_c(w)+offset-w
  // offset= e_dft+Sigma_x-Vxc
  class QPFunc {
   public:
    QPFunc(Index gw_level, const Sigma_base& sigma, double offset)
        : _gw_level(gw_level), _offset(offset), _sigma_c_func(sigma){};
    std::pair<double, double> operator()(double frequency) const {
      std::pair<double, double> value;

      value.first =
          _sigma_c_func.CalcCorrelationDiagElement(_gw_level, frequency);
      value.second = _sigma_c_func.CalcCorrelationDiagElementDerivative(
          _gw_level, frequency);
      value.first += (_offset - frequency);
      value.second -= 1.0;
      return value;
    }
    double value(double frequency) const {
      return _sigma_c_func.CalcCorrelationDiagElement(_gw_level, frequency) +
             _offset - frequency;
    }
    double deriv(double frequency) const {
      return _sigma_c_func.CalcCorrelationDiagElementDerivative(_gw_level,
                                                                frequency) +
             _offset - frequency;
    }

   private:
    Index _gw_level;
    double _offset;
    const Sigma_base& _sigma_c_func;
  };

  double SolveQP_Bisection(double lowerbound, double f_lowerbound,
                           double upperbound, double f_upperbound,
                           const QPFunc& f) const;

  Eigen::VectorXd Laplacian_Kernel(double x1, Eigen::VectorXd x2,
                                   double sigma) const;

  double CalcHomoLumoShift(Eigen::VectorXd frequencies) const;
  Eigen::VectorXd ScissorShift_DFTlevel(
      const Eigen::VectorXd& dft_energies) const;
  void PrintQP_Energies(const Eigen::VectorXd& qp_diag_energies) const;
  void PrintGWA_Energies() const;

  Eigen::VectorXd SolveQP(const Eigen::VectorXd& frequencies) const;
  boost::optional<double> SolveQP_Grid(double intercept0, double frequency0,
                                       Index gw_level) const;

  boost::optional<double> SolveQP_Grid_reduced_interval(double intercept0,
                                                        double frequency0,
                                                        Index gw_level) const;
  boost::optional<double> SolveQP_Regression(double intercept0,
                                             double frequency0,
                                             Index gw_level) const;
  boost::optional<double> SolveQP_SelfConsistent(double intercept0,
                                                 double frequency0,
                                                 Index gw_level) const;
  boost::optional<double> SolveQP_FixedPoint(double intercept0,
                                             double frequency0,
                                             Index gw_level) const;
  boost::optional<double> SolveQP_Linearisation(double intercept0,
                                                double frequency0,
                                                Index gw_level) const;
  bool Converged(const Eigen::VectorXd& e1, const Eigen::VectorXd& e2,
                 double epsilon) const;
  void ExportCorrelationDiags(const Eigen::VectorXd& frequencies) const;
};
}  // namespace xtp
}  // namespace votca

#endif /* _VOTCA_XTP_BSE_H */
