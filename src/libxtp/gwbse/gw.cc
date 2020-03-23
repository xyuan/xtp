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

#include "votca/xtp/rpa.h"
#include "votca/xtp/sigma_ci.h"
#include "votca/xtp/sigma_exact.h"
#include "votca/xtp/sigma_ppm.h"
#include <fstream>
#include <iostream>
#include <votca/xtp/IndexParser.h>
#include <votca/xtp/gw.h>
#include <votca/xtp/newton_rapson.h>
#include <votca/xtp/selfconsistentsolver.h>
namespace votca {
namespace xtp {

void GW::configure(const options& opt) {
  _opt = opt;
  _qptotal = _opt.qpmax - _opt.qpmin + 1;
  _rpa.configure(_opt.homo, _opt.rpamin, _opt.rpamax);
  if (_opt.sigma_integration == "exact") {
    _sigma = std::make_unique<Sigma_Exact>(Sigma_Exact(_Mmn, _rpa));
  } else if (_opt.sigma_integration == "ci") {
    _sigma = std::make_unique<Sigma_CI>(Sigma_CI(_Mmn, _rpa));
  } else if (_opt.sigma_integration == "ppm") {
    _sigma = std::make_unique<Sigma_PPM>(Sigma_PPM(_Mmn, _rpa));
  }

  XTP_LOG(Log::error, _log)
      << TimeStamp() << " Using " << _opt.sigma_integration
      << " for Correlation part of self-energy" << std::flush;
  Sigma_base::options sigma_opt;
  sigma_opt.order = _opt.order;
  sigma_opt.alpha = _opt.alpha;
  sigma_opt.homo = _opt.homo;
  sigma_opt.qpmax = _opt.qpmax;
  sigma_opt.qpmin = _opt.qpmin;
  sigma_opt.rpamin = _opt.rpamin;
  sigma_opt.rpamax = _opt.rpamax;
  sigma_opt.eta = _opt.eta;
  sigma_opt.quadrature_scheme = _opt.quadrature_scheme;
  _sigma->configure(sigma_opt);
  _Sigma_x = Eigen::MatrixXd::Zero(_qptotal, _qptotal);
  _Sigma_c = Eigen::MatrixXd::Zero(_qptotal, _qptotal);
}

double GW::CalcHomoLumoShift(Eigen::VectorXd frequencies) const {
  double DFTgap = _dft_energies(_opt.homo + 1) - _dft_energies(_opt.homo);
  double QPgap = frequencies(_opt.homo + 1 - _opt.qpmin) -
                 frequencies(_opt.homo - _opt.qpmin);
  return QPgap - DFTgap;
}

Eigen::MatrixXd GW::getHQP() const {
  return _Sigma_x + _Sigma_c - _vxc +
         Eigen::MatrixXd(
             _dft_energies.segment(_opt.qpmin, _qptotal).asDiagonal());
}

Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> GW::DiagonalizeQPHamiltonian()
    const {
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> qpdiag(getHQP());
  PrintQP_Energies(qpdiag.eigenvalues());
  return qpdiag;
}

void GW::PrintGWA_Energies() const {
  Eigen::VectorXd gwa_energies = getGWAResults();
  double shift = CalcHomoLumoShift(gwa_energies);

  XTP_LOG(Log::error, _log)
      << (boost::format(
              "  ====== Perturbative quasiparticle energies (Hartree) ====== "))
             .str()
      << std::flush;
  XTP_LOG(Log::error, _log)
      << (boost::format("   DeltaHLGap = %1$+1.6f Hartree") % shift).str()
      << std::flush;

  for (Index i = 0; i < _qptotal; i++) {
    std::string level = "  Level";
    if ((i + _opt.qpmin) == _opt.homo) {
      level = "  HOMO ";
    } else if ((i + _opt.qpmin) == _opt.homo + 1) {
      level = "  LUMO ";
    }

    XTP_LOG(Log::error, _log)
        << level
        << (boost::format(" = %1$4d DFT = %2$+1.4f VXC = %3$+1.4f S-X = "
                          "%4$+1.4f S-C = %5$+1.4f GWA = %6$+1.4f") %
            (i + _opt.qpmin) % _dft_energies(i + _opt.qpmin) % _vxc(i, i) %
            _Sigma_x(i, i) % _Sigma_c(i, i) % gwa_energies(i))
               .str()
        << std::flush;
  }
  return;
}

void GW::PrintQP_Energies(const Eigen::VectorXd& qp_diag_energies) const {
  Eigen::VectorXd gwa_energies = getGWAResults();
  XTP_LOG(Log::error, _log)
      << TimeStamp() << " Full quasiparticle Hamiltonian  " << std::flush;
  XTP_LOG(Log::error, _log)
      << (boost::format(
              "  ====== Diagonalized quasiparticle energies (Hartree) "
              "====== "))
             .str()
      << std::flush;
  for (Index i = 0; i < _qptotal; i++) {
    std::string level = "  Level";
    if ((i + _opt.qpmin) == _opt.homo) {
      level = "  HOMO ";
    } else if ((i + _opt.qpmin) == _opt.homo + 1) {
      level = "  LUMO ";
    }
    XTP_LOG(Log::error, _log)
        << level
        << (boost::format(" = %1$4d PQP = %2$+1.6f DQP = %3$+1.6f ") %
            (i + _opt.qpmin) % gwa_energies(i) % qp_diag_energies(i))
               .str()
        << std::flush;
  }
  return;
}

Eigen::VectorXd GW::ScissorShift_DFTlevel(
    const Eigen::VectorXd& dft_energies) const {
  Eigen::VectorXd shifted_energies = dft_energies;
  shifted_energies.segment(_opt.homo + 1, dft_energies.size() - _opt.homo - 1)
      .array() += _opt.shift;
  return shifted_energies;
}

void GW::CalculateGWPerturbation() {

  _Sigma_x = (1 - _opt.ScaHFX) * _sigma->CalcExchangeMatrix();
  XTP_LOG(Log::error, _log)
      << TimeStamp() << " Calculated Hartree exchange contribution"
      << std::flush;
  // dftenergies has size aobasissize
  // rpaenergies/Mmn have size rpatotal
  // gwaenergies/frequencies have size qptotal
  // homo index is relative to dft_energies
  XTP_LOG(Log::error, _log)
      << TimeStamp() << " Scissor shifting DFT energies by: " << _opt.shift
      << " Hrt" << std::flush;
  Eigen::VectorXd dft_shifted_energies = ScissorShift_DFTlevel(_dft_energies);
  _rpa.setRPAInputEnergies(
      dft_shifted_energies.segment(_opt.rpamin, _opt.rpamax - _opt.rpamin + 1));
  Eigen::VectorXd frequencies =
      dft_shifted_energies.segment(_opt.qpmin, _qptotal);
  for (Index i_gw = 0; i_gw < _opt.gw_sc_max_iterations; ++i_gw) {

    if (i_gw % _opt.reset_3c == 0 && i_gw != 0) {
      _Mmn.Rebuild();
      XTP_LOG(Log::info, _log)
          << TimeStamp() << " Rebuilding 3c integrals" << std::flush;
    }
    _sigma->PrepareScreening();
    XTP_LOG(Log::info, _log)
        << TimeStamp() << " Calculated screening via RPA" << std::flush;
    XTP_LOG(Log::info, _log) << TimeStamp() << " Solving QP equations \n"
                             << std::flush;
    frequencies = SolveQP(frequencies);

    if (_opt.gw_sc_max_iterations > 1) {

      Eigen::VectorXd rpa_energies_old = _rpa.getRPAInputEnergies();
      _rpa.UpdateRPAInputEnergies(_dft_energies, frequencies, _opt.qpmin);
      XTP_LOG(Log::error, _log)
          << TimeStamp() << " GW_Iteration:" << i_gw
          << " Shift[Hrt]:" << CalcHomoLumoShift(frequencies) << std::flush;
      if (Converged(_rpa.getRPAInputEnergies(), rpa_energies_old,
                    _opt.gw_sc_limit)) {
        XTP_LOG(Log::error, _log)
            << TimeStamp() << " Converged after " << i_gw + 1
            << " GW iterations." << std::flush;
        break;
      } else if (i_gw == _opt.gw_sc_max_iterations - 1) {
        XTP_LOG(Log::error, _log)
            << TimeStamp()
            << " WARNING! GW-self-consistency cycle not converged after "
            << _opt.gw_sc_max_iterations << " iterations." << std::flush;
        XTP_LOG(Log::error, _log)
            << TimeStamp() << "      Run continues. Inspect results carefully!"
            << std::flush;
        break;
      }
    }
  }
  _Sigma_c.diagonal() = _sigma->CalcCorrelationDiag(frequencies);
  PrintGWA_Energies();
}

Eigen::VectorXd GW::getGWAResults() const {
  return _Sigma_x.diagonal() + _Sigma_c.diagonal() - _vxc.diagonal() +
         _dft_energies.segment(_opt.qpmin, _qptotal);
}

Eigen::VectorXd GW::SolveQP(const Eigen::VectorXd& frequencies) const {
  const Eigen::VectorXd intercepts =
      _dft_energies.segment(_opt.qpmin, _qptotal) + _Sigma_x.diagonal() -
      _vxc.diagonal();
  Eigen::VectorXd frequencies_new = frequencies;
  Eigen::Array<bool, Eigen::Dynamic, 1> converged =
      Eigen::Array<bool, Eigen::Dynamic, 1>::Zero(_qptotal);

#pragma omp parallel for schedule(dynamic)
  for (Index gw_level = 0; gw_level < _qptotal; ++gw_level) {
    double initial_f = frequencies[gw_level];
    double intercept = intercepts[gw_level];
    boost::optional<double> newf;
    if (_opt.qp_solver == "fixedpoint") {
      newf = SolveQP_FixedPoint(intercept, initial_f, gw_level);
    }
    if (_opt.qp_solver == "scf") {
      newf = SolveQP_SelfConsistent(intercept, initial_f, gw_level);
    }
    if (_opt.qp_solver == "regression") {
      newf = SolveQP_Regression(intercept, initial_f, gw_level);
    }
    if (_opt.qp_solver == "grid") {
      newf = SolveQP_Grid(intercept, initial_f, gw_level);
    }
    if (_opt.qp_solver == "grid_ps") {
      newf = SolveQP_Grid_reduced_interval(intercept, initial_f, gw_level);
    }
    if (newf) {
      frequencies_new[gw_level] = newf.value();
      converged[gw_level] = true;
    } else {
      newf = SolveQP_Grid(intercept, initial_f, gw_level);
      if (newf) {
        frequencies_new[gw_level] = newf.value();
        converged[gw_level] = true;
      } else {
        newf = SolveQP_Linearisation(intercept, initial_f, gw_level);
        if (newf) {
          frequencies_new[gw_level] = newf.value();
        }
      }
    }
  }

  if ((converged == false).any()) {
    std::vector<Index> states;
    for (Index s = 0; s < converged.size(); s++) {
      if (!converged[s]) {
        states.push_back(s);
      }
    }
    IndexParser rp;
    XTP_LOG(Log::error, _log) << TimeStamp() << " Not converged PQP states are:"
                              << rp.CreateIndexString(states) << std::flush;
    XTP_LOG(Log::error, _log)
        << TimeStamp() << " Increase the grid search interval" << std::flush;
  }
  return frequencies_new;
}

boost::optional<double> GW::SolveQP_Linearisation(double intercept0,
                                                  double frequency0,
                                                  Index gw_level) const {
  boost::optional<double> newf = boost::none;

  double sigma = _sigma->CalcCorrelationDiagElement(gw_level, frequency0);
  double dsigma_domega =
      _sigma->CalcCorrelationDiagElementDerivative(gw_level, frequency0);
  double Z = 1.0 - dsigma_domega;
  if (std::abs(Z) > 1e-9) {
    newf = frequency0 + (intercept0 - frequency0 + sigma) / Z;
  }
  return newf;
}

boost::optional<double> GW::SolveQP_Grid(double intercept0, double frequency0,
                                         Index gw_level) const {
  const double range =
      _opt.qp_grid_spacing * double(_opt.qp_grid_steps - 1) / 2.0;
  boost::optional<double> newf = boost::none;
  QPFunc fqp(gw_level, *_sigma.get(), intercept0);
  double freq_prev = frequency0 - range;
  double targ_prev = fqp.value(freq_prev);
  Index numbersofcalls = 1;  // This is for pre-shooting
  Index numbersofcalls_bisection = 0;
  double qp_energy = 0.0;
  double gradient_max = std::numeric_limits<double>::max();
  bool pole_found = false;
  for (Index i_node = 1; i_node < _opt.qp_grid_steps; ++i_node) {
    double freq = freq_prev + _opt.qp_grid_spacing;
    double targ = fqp.value(freq);
    numbersofcalls += 1;
    if (targ_prev * targ < 0.0) {  // Sign change
    std::pair<double,Index> fandcount = SolveQP_Bisection_c(freq_prev, targ_prev, freq, targ, fqp);
      double f = fandcount.first;
      numbersofcalls_bisection += fandcount.second;
      //double f = SolveQP_Bisection(freq_prev, targ_prev, freq, targ, fqp);
      double gradient = std::abs(fqp.deriv(f));
      numbersofcalls++;
      if (gradient < gradient_max) {
        qp_energy = f;
        gradient_max = gradient;
        pole_found = true;
      }
    }
    freq_prev = freq;
    targ_prev = targ;
  }

  if (pole_found) {
    newf = qp_energy;
  }
  XTP_LOG(Log::error, _log)
      << "Level " << gw_level << " Sigma evaluations " << numbersofcalls+numbersofcalls_bisection << "\n"
      << std::flush;
  numbersofcalls = 0;
  return newf;
}

boost::optional<double> GW::SolveQP_Grid_reduced_interval(
    double intercept0, double frequency0, Index gw_level) const {

  boost::optional<double> newf = boost::none;
  QPFunc fqp(gw_level, *_sigma.get(), intercept0);
  // To reduce the frequency window where to find the fixed point we compute
  // the interesection between the linearization of fqp =
  // sigma_c+intercept0-frequency0 around the frequency0 and the line fqp = 0.
  // In other words this is the solution of the linearization problem. We use
  // this as starting point. We shift +- 0.5 Hartree from that as default
  double initial_targ_prev = fqp.value(frequency0);
  double initial_targ_prev_div = fqp.deriv(frequency0);
  Index numbersofcalls = 2;  // This is for pre-shooting
  Index numbersofcalls_bisection = 0;
  double initial_f = frequency0;
  frequency0 = initial_f + (initial_targ_prev) / (1.0 - initial_targ_prev_div);

  const double range = _opt.qp_grid_hartree;
  double freq_prev = frequency0 - range;
  double targ_prev = fqp.value(freq_prev);
  double qp_energy = 0.0;
  double gradient_max = std::numeric_limits<double>::max();
  bool pole_found = false;

  for (Index i_node = 1; i_node < _opt.qp_grid_steps; ++i_node) {

    double freq =
        freq_prev + 2.0 * range / _opt.qp_grid_steps;  // _opt.qp_grid_spacing;
    double targ = fqp.value(freq);
    numbersofcalls++;

    if (targ_prev * targ < 0.0) {  // Sign change
      std::pair<double,Index> fandcount = SolveQP_Bisection_c(freq_prev, targ_prev, freq, targ, fqp);
      double f = fandcount.first;
      numbersofcalls_bisection += fandcount.second;
      //double f = SolveQP_Bisection(freq_prev, targ_prev, freq, targ, fqp);
      double gradient = std::abs(fqp.deriv(f));
      numbersofcalls++;
      if (gradient < gradient_max) {
        qp_energy = f;
        gradient_max = gradient;
        pole_found = true;
      }
    }
    freq_prev = freq;
    targ_prev = targ;
  }

  if (pole_found) {
    newf = qp_energy;
  }
  XTP_LOG(Log::error, _log)
      << "Level " << gw_level << " Sigma evaluations " << numbersofcalls+numbersofcalls_bisection << "\n"
      << std::flush;
  return newf;
  numbersofcalls = 0;
}

boost::optional<double> GW::SolveQP_Regression(double intercept0,
                                               double frequency0,
                                               Index gw_level) const {

  boost::optional<double> newf = boost::none;
  QPFunc fqp(gw_level, *_sigma.get(), intercept0);

  // To reduce the frequency window where to find the fixed point we compute
  // the interesection between the linearization of fqp =
  // sigma_c+intercept0-frequency0 around the frequency0 and the line fqp = 0.
  // In other words this is the solution of the linearization problem. We use
  // this as starting point. We shift +- 0.5 Hartree from that as default

  double initial_targ_prev = fqp.value(frequency0);
  double initial_targ_prev_div = fqp.deriv(frequency0);

  const double range = _opt.qp_grid_hartree;
  double initial_f = frequency0;
  frequency0 = initial_f + (initial_targ_prev) / (1.0 - initial_targ_prev_div);

  double freq_prev = frequency0 - range;
  double targ_prev = fqp.value(freq_prev);

  Index numbersofcalls = 3;

  double qp_energy = 0.0;

  bool pole_found = false;
  bool mae_test_pass = false;
  double mae_final = 0;
  Eigen::VectorXd alphas;
  Index num_max = _opt.qp_training_points;
  Eigen::VectorXd frequencies;

  double delta = (2.0 * range) / ((1.0 * num_max - 1.0));

  Eigen::VectorXd freq_training_i(num_max);
  Eigen::VectorXd sigma_training_i(num_max);

  std::vector<Eigen::VectorXd> sigma;
  std::vector<Eigen::VectorXd> freq;

  for (Index j = 0; j < freq_training_i.size(); ++j) {
    freq_training_i(j) = freq_prev + j * delta;
    sigma_training_i(j) =
        fqp.value(freq_training_i(j)) + freq_training_i(j) - intercept0;
    numbersofcalls++;
  }

  freq.push_back(freq_training_i);
  sigma.push_back(sigma_training_i);

  Index step = 1;
  double error_fp = 0;
  double nc_en = 0;
  Index test_size = freq_training_i.size() - 1;
  while (step < 20) {

    double mae = 0;

    Eigen::VectorXd freq_training = freq[step - 1];

    Eigen::VectorXd sigma_training = sigma[step - 1];

    Eigen::MatrixXd kernel(freq_training.size(), sigma_training.size());

    if (_opt.qp_test_points == "both") {
      test_size = freq_training.size() - 1;
    } else {
      test_size = (freq_training.size() - 1) / 2;
    }

    Eigen::VectorXd freq_test((freq_training.size() - 1) / 2);
    Eigen::VectorXd sigma_test((sigma_training.size() - 1) / 2);

    delta *= 0.5;

    for (Index j = 0; j < freq_test.size(); ++j) {
      if (_opt.qp_test_points == "odd") {
        freq_test(j) = freq_training(2 * j + 1) - delta;
      } else if (_opt.qp_test_points == "even") {
        freq_test(j) = freq_training(2 * j) + delta;
      } else {
        freq_test(j) = freq_training(j) + delta;
      }
      sigma_test(j) = fqp.value(freq_test(j)) + freq_test(j) - intercept0;
      numbersofcalls++;
    }

    for (Index i = 0; i < freq_training.size(); ++i) {
      kernel.row(i) =
          Laplacian_Kernel(freq_training(i), freq_training, _opt.qp_spread);
    }
    kernel.diagonal().array() += 1e-8;
    alphas = kernel.colPivHouseholderQr().solve(sigma_training);

    Index num_test = freq_test.size();
    for (Index t = 0; t < freq_test.size() - 1; t++) {
      mae +=
          std::abs(Laplacian_Kernel(freq_test(t), freq_training, _opt.qp_spread)
                       .dot(alphas) -
                   sigma_test(t));
    }

    mae /= (double)num_test;
    mae_final = mae;
    if (mae < _opt.qp_mae_tol) {
      frequencies = freq_training;
      mae_test_pass = true;
      break;
    } else {

      Eigen::VectorXd temp_f(freq_training.size() + freq_test.size());

      Eigen::VectorXd temp_s(sigma_training.size() + sigma_test.size());

      temp_f << freq_training, freq_test;

      temp_s << sigma_training, sigma_test;

      freq.push_back(temp_f);
      sigma.push_back(temp_s);

      step++;
    }
  }
  if (mae_test_pass == false) {
    XTP_LOG(Log::error, _log) << " MAE test failed for Level: " << gw_level
                              << " Number of points increased" << std::flush;
  } else {
    XTP_LOG(Log::debug, _log) << " MAE test passed for Level: " << gw_level
                              << " Step needed: " << step << std::flush;

    // // Stupid Fixed point solver
    // double p0 = frequency0;
    // Index fps = 1;

    // while (fps <= 50000) {
    //   double p = Laplacian_Kernel(p0, frequencies,
    //   _opt.qp_spread).dot(alphas) +
    //              intercept0;
    //   if (std::abs(p - p0) < _opt.qp_fixedpoint_tol) {
    //     qp_energy = p;
    //     pole_found = true;
    //     break;
    //   }
    //   error_fp = std::abs(p - p0);
    //   nc_en = p;
    //   fps++;
    //   p0 = p;
    // }

    // Aitkin Method
    double p0 = frequency0;
    for (Index i = 0; i < 100; i++) {
      double p1 =
          Laplacian_Kernel(p0, frequencies, _opt.qp_spread).dot(alphas) +
          intercept0;
      double p2 =
          Laplacian_Kernel(p1, frequencies, _opt.qp_spread).dot(alphas) +
          intercept0;
      double d = (p2 - p1) - (p1 - p0);
      if (std::abs(d) < 1e-16) {
        XTP_LOG(Log::error, _log)
            << "This denominator should not be so small" << std::flush;
        break;
      }
      double ap = p2 - ((p2 - p1) * (p2 - p1)) / d;
      if (std::abs(ap - p2) < _opt.qp_fixedpoint_tol) {
        qp_energy = ap;
        pole_found = true;
        break;
      }
      error_fp = std::abs(ap - p0);
      nc_en = ap;
      p0 = ap;
    }
  }
  XTP_LOG(Log::error, _log)
      << "Level " << gw_level << " Sigma evaluations " << numbersofcalls << "\n"
      << std::flush;
  if (pole_found) {
    newf = qp_energy;
  } else {

    XTP_LOG(Log::error, _log) << " Fixed point not found for " << gw_level
                              << " going to grid evaluation. "
                              << "Error: " << error_fp
                              << "Not converged energy " << nc_en << std::flush;
  }

  numbersofcalls = 0;
  return newf;
}  // namespace xtp

boost::optional<double> GW::SolveQP_SelfConsistent(double intercept0,
                                                   double frequency0,
                                                   Index gw_level) const {
  boost::optional<double> newf = boost::none;
  QPFunc f(gw_level, *_sigma.get(), intercept0);
  Selfconsistentsolver<QPFunc> scf = Selfconsistentsolver<QPFunc>(
      _opt.g_sc_max_iterations, _opt.g_sc_limit, _opt.qp_solver_alpha);
  double freq_new = scf.FindRoot(f, frequency0);
  if (scf.getInfo() == Selfconsistentsolver<QPFunc>::success) {
    newf = freq_new;
  }
  return newf;
}

boost::optional<double> GW::SolveQP_FixedPoint(double intercept0,
                                               double frequency0,
                                               Index gw_level) const {
  boost::optional<double> newf = boost::none;
  QPFunc f(gw_level, *_sigma.get(), intercept0);
  NewtonRapson<QPFunc> newton = NewtonRapson<QPFunc>(
      _opt.g_sc_max_iterations, _opt.g_sc_limit, _opt.qp_solver_alpha);
  double freq_new = newton.FindRoot(f, frequency0);
  if (newton.getInfo() == NewtonRapson<QPFunc>::success) {
    newf = freq_new;
  }
  return newf;
}

// https://en.wikipedia.org/wiki/Bisection_method
std::pair<double,Index> GW::SolveQP_Bisection_c(double lowerbound, double f_lowerbound,
                             double upperbound, double f_upperbound,
                             const QPFunc& f) const {

  Index callf = 0;
  if (f_lowerbound * f_upperbound > 0) {
    throw std::runtime_error(
        "Bisection needs a postive and negative function value");
  }
  double zero = 0.0;
  while (true) {
    double c = 0.5 * (lowerbound + upperbound);
    if (std::abs(upperbound - lowerbound) < _opt.g_sc_limit) {
      zero = c;
      break;
    }
    double y_c = f.value(c);
    callf++;
    if (std::abs(y_c) < _opt.g_sc_limit) {
      zero = c;
      break;
    }
    if (y_c * f_lowerbound > 0) {
      lowerbound = c;
      f_lowerbound = y_c;
    } else {
      upperbound = c;
      f_upperbound = y_c;
    }
  }
  return std::make_pair(zero,callf);
}

double GW::SolveQP_Bisection(double lowerbound, double f_lowerbound,
                             double upperbound, double f_upperbound,
                             const QPFunc& f) const {

  Index callf = 0;
  if (f_lowerbound * f_upperbound > 0) {
    throw std::runtime_error(
        "Bisection needs a postive and negative function value");
  }
  double zero = 0.0;
  while (true) {
    double c = 0.5 * (lowerbound + upperbound);
    if (std::abs(upperbound - lowerbound) < _opt.g_sc_limit) {
      zero = c;
      break;
    }
    double y_c = f.value(c);
    callf++;
    if (std::abs(y_c) < _opt.g_sc_limit) {
      zero = c;
      break;
    }
    if (y_c * f_lowerbound > 0) {
      lowerbound = c;
      f_lowerbound = y_c;
    } else {
      upperbound = c;
      f_upperbound = y_c;
    }
  }
  XTP_LOG(Log::error, _log) << "Bisection evaluations " << callf << std::flush;
  return zero;
}


Eigen::VectorXd GW::Laplacian_Kernel(double x1, Eigen::VectorXd x2,
                                     double sigma) const {
  Eigen::VectorXd kernel(x2.size());
  for (Index j = 0; j < x2.size(); ++j) {
    kernel(j) = std::exp(std::abs(x1 - x2(j)) / sigma);
  }
  return kernel;
}

bool GW::Converged(const Eigen::VectorXd& e1, const Eigen::VectorXd& e2,
                   double epsilon) const {
  Index state = 0;
  bool energies_converged = true;
  double diff_max = (e1 - e2).cwiseAbs().maxCoeff(&state);
  if (diff_max > epsilon) {
    energies_converged = false;
  }
  XTP_LOG(Log::info, _log) << TimeStamp() << " E_diff max=" << diff_max
                           << " StateNo:" << state << std::flush;
  return energies_converged;
}

void GW::CalculateHQP() {
  Eigen::VectorXd diag_backup = _Sigma_c.diagonal();
  _Sigma_c = _sigma->CalcCorrelationOffDiag(getGWAResults());
  _Sigma_c.diagonal() = diag_backup;
}

void GW::PlotSigma(std::string filename, Index steps, double spacing,
                   std::string states) const {

  Eigen::VectorXd frequencies =
      _rpa.getRPAInputEnergies().segment(_opt.qpmin - _opt.rpamin, _qptotal);

  std::vector<Index> state_inds;
  IndexParser rp;
  std::vector<Index> parsed_states = rp.CreateIndexVector(states);
  for (Index gw_level : parsed_states) {
    if (gw_level >= _opt.qpmin && gw_level <= _opt.qpmax) {
      state_inds.push_back(gw_level);
    }
  }
  XTP_LOG(Log::error, _log)
      << TimeStamp() << " PQP(omega) written to '" << filename
      << "' for states " << rp.CreateIndexString(state_inds) << std::flush;

  const Index num_states = state_inds.size();

  const Eigen::VectorXd intercept =
      _dft_energies.segment(_opt.qpmin, _qptotal) + _Sigma_x.diagonal() -
      _vxc.diagonal();
  Eigen::MatrixXd mat = Eigen::MatrixXd::Zero(steps, 2 * num_states);
#pragma omp parallel for schedule(dynamic)
  for (Index grid_point = 0; grid_point < steps; grid_point++) {
    const double offset =
        ((double)grid_point - ((double)(steps - 1) / 2.0)) * spacing;
    for (Index i = 0; i < num_states; i++) {
      const Index gw_level = state_inds[i];
      const double omega = frequencies(gw_level) + offset;
      double sigma = _sigma->CalcCorrelationDiagElement(gw_level, omega);
      mat(grid_point, 2 * i) = omega;
      mat(grid_point, 2 * i + 1) = sigma + intercept[gw_level];
    }
  }

  std::ofstream out;
  out.open(filename);
  for (Index i = 0; i < num_states; i++) {
    const Index gw_level = state_inds[i];
    out << boost::format("#%1$somega_%2$d\tE_QP(omega)_%2$d") %
               (i == 0 ? "" : "\t") % gw_level;
  }
  out << std::endl;
  boost::format numFormat("%+1.6f");
  Eigen::IOFormat matFormat(Eigen::StreamPrecision, 0, "\t", "\n");
  out << numFormat % mat.format(matFormat) << std::endl;
  out.close();
}

}  // namespace xtp
}  // namespace votca
