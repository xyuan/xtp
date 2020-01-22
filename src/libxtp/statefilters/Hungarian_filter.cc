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

#include "Hungarian_filter.h"
#include <votca/xtp/aomatrix.h>
#include <fstream>
#include <iostream>
namespace votca {
namespace xtp {

void Hungarian_filter::Initialize(const tools::Property& options) {
  _threshold = options.ifExistsReturnElseThrowRuntimeError<double>(".");
}

void Hungarian_filter::Info(Logger& log) const {
  if (_threshold == 0.0) {
    XTP_LOG(Log::error, log)
        << "Using assignment filter with no threshold " << std::flush;
  } else {
    XTP_LOG(Log::error, log) << "Using assignment filter with threshold "
                             << _threshold << std::flush;
  }
}

Eigen::VectorXd Hungarian_filter::CalculateOverlap(const Orbitals& orb,
                                                   QMStateType type) const {
  AOOverlap S_ao;
  S_ao.Fill(orb.SetupDftBasis());
  Eigen::MatrixXd coeffs = CalcAOCoeffs(orb, type);
  Index basis = orb.getBasisSetSize();
  Index nn = orb.BSESinglets().eigenvalues().size();
  Eigen::MatrixXd overlap = Eigen::MatrixXd::Zero(nn,nn);
  Eigen::VectorXd optimalassignment = Eigen::VectorXd::Zero(nn);
#pragma omp parallel for schedule(dynamic)
  for (Index i = 0; i < coeffs.cols(); i++) {
    for (Index j = 0; j < _laststatecoeffs.cols(); j++) {
      Eigen::VectorXd _laststatecoeff = _laststatecoeffs.col(j);
      {
        Eigen::VectorXd state = coeffs.col(i).head(basis * basis);
        Eigen::Map<const Eigen::MatrixXd> mat(state.data(), basis, basis);
        Eigen::VectorXd laststate = _laststatecoeff.head(basis * basis);
        Eigen::Map<const Eigen::MatrixXd> lastmat(laststate.data(), basis,
                                                  basis);

        overlap(j,i) = (mat * S_ao.Matrix() * lastmat.transpose())
                         .cwiseProduct(S_ao.Matrix())
                         .sum();
      }
      if (!orb.getTDAApprox()) {
        Eigen::VectorXd state = coeffs.col(i).tail(basis * basis);
        Eigen::Map<const Eigen::MatrixXd> mat(state.data(), basis, basis);
        Eigen::VectorXd laststate = _laststatecoeff.tail(basis * basis);
        Eigen::Map<const Eigen::MatrixXd> lastmat(laststate.data(), basis,
                                                  basis);

        overlap(j,i) -= (mat * S_ao.Matrix() * lastmat.transpose())
                             .cwiseProduct(S_ao.Matrix())
                             .sum();
      }
    }
  }
  std::string file_in = "overlapmatrix.txt";
  std::cout << "Running Optimal Assignment Python script" << std::endl;
  std::ofstream overlapmatrix(file_in);
  overlapmatrix << overlap.cwiseAbs2();
  overlapmatrix.close();
  std::string file_out = "assignment_out.txt";
   
  std::string command =
      "python optimalassigment.py " + file_out + " " + _laststate; //std::to_string(ss);
  std::system(command.c_str());
  std::ifstream ifs(file_out);
  
  std::vector<double> values;
  double val;
  int i = 0;
  while (ifs >> val) {
    optimalassignment(i) = val;
    i += 1;
  }
  return optimalassignment;
}

Eigen::MatrixXd Hungarian_filter::CalcExcitonAORepresentation(
    const Orbitals& orb, QMStateType type) const {
  Eigen::MatrixXd coeffs;
  Index nostates = orb.NumberofStates(type);
  Index bse_cmax = orb.getBSEcmax();
  Index bse_cmin = orb.getBSEcmin();
  Index bse_vmax = orb.getBSEvmax();
  Index bse_vmin = orb.getBSEvmin();
  Index bse_vtotal = bse_vmax - bse_vmin + 1;
  Index bse_ctotal = bse_cmax - bse_cmin + 1;
  Index basis = orb.getBasisSetSize();
  Index bse_size_ao = basis * basis;
  auto occlevels = orb.MOs().eigenvectors().block(
      0, bse_vmin, orb.MOs().eigenvectors().rows(), bse_vtotal);
  auto virtlevels = orb.MOs().eigenvectors().block(
      0, bse_cmin, orb.MOs().eigenvectors().rows(), bse_ctotal);

  if (orb.getTDAApprox()) {
    coeffs.resize(bse_size_ao, nostates);
  } else {
    coeffs.resize(2 * bse_size_ao, nostates);
  }
#pragma omp parallel for schedule(dynamic)
  for (Index i = 0; i < nostates; i++) {
    {
      Eigen::VectorXd exciton;
      if (type == QMStateType::Singlet) {
        exciton = orb.BSESinglets().eigenvectors().col(i);
      } else {
        exciton = orb.BSETriplets().eigenvectors().col(i);
      }
      Eigen::Map<const Eigen::MatrixXd> mat(exciton.data(), bse_ctotal,
                                            bse_vtotal);
      const Eigen::MatrixXd aomatrix =
          occlevels * mat.transpose() * virtlevels.transpose();
      coeffs.col(i).head(bse_size_ao) =
          Eigen::Map<const Eigen::VectorXd>(aomatrix.data(), bse_size_ao);
    }
    if (!orb.getTDAApprox()) {
      Eigen::VectorXd exciton;
      if (type == QMStateType::Singlet) {
        exciton = orb.BSESinglets().eigenvectors2().col(i);
      } else {
        exciton = orb.BSETriplets().eigenvectors2().col(i);
      }
      Eigen::Map<const Eigen::MatrixXd> mat(exciton.data(), bse_ctotal,
                                            bse_vtotal);
      const Eigen::MatrixXd aomatrix =
          occlevels * mat.transpose() * virtlevels.transpose();
      coeffs.col(i).tail(bse_size_ao) =
          Eigen::Map<const Eigen::VectorXd>(aomatrix.data(), bse_size_ao);
    }
  }
  return coeffs;
}

Eigen::MatrixXd Hungarian_filter::CalcAOCoeffs(const Orbitals& orb,
                                               QMStateType type) const {
  Eigen::MatrixXd coeffs;
  if (type.isSingleParticleState()) {
    if (type == QMStateType::DQPstate) {
      coeffs = orb.CalculateQParticleAORepresentation();
    } else {
      coeffs = orb.MOs().eigenvectors();
    }
  } else {
    coeffs = CalcExcitonAORepresentation(orb, type);
  }
  return coeffs;
}

void Hungarian_filter::UpdateHist(const Orbitals& orb, QMState state) {
  Eigen::MatrixXd aocoeffs = CalcAOCoeffs(orb, state.Type());
  Index offset = 0;
  if (state.Type() == QMStateType::DQPstate) {
    offset = orb.getGWAmin();
  }
  _laststatecoeffs = aocoeffs;
  _laststate = state.ToString();
  //_laststatecoeffs = aocoeffs.col(state.StateIdx() - offset);
}

std::vector<Index> Hungarian_filter::CalcIndeces(const Orbitals& orb,
                                                 QMStateType type) const {
  Index offset = 0;
  if (type.isGWState()) {
    offset = orb.getGWAmin();
  }
  Eigen::MatrixXd Overlap = CalculateOverlap(orb, type);
  return ReduceAndSortIndecesUp(Overlap, offset, _threshold);
}

void Hungarian_filter::WriteToCpt(CheckpointWriter& w) {
  w(_laststatecoeffs, "laststatecoeffs");
  w(_threshold, "threshold");
}

void Hungarian_filter::ReadFromCpt(CheckpointReader& r) {
  r(_laststatecoeffs, "laststatecoeffs");
  r(_threshold, "threshold");
}

}  // namespace xtp
}  // namespace votca