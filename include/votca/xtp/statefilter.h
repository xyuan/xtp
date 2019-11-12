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
#ifndef VOTCA_XTP_STATEFILTER_H
#define VOTCA_XTP_STATEFILTER_H

#include <votca/xtp/logger.h>
#include <votca/xtp/orbitals.h>
#include <votca/xtp/populationanalysis.h>
#include <votca/xtp/qmfragment.h>
#include <votca/xtp/qmstate.h>

namespace votca {
namespace xtp {
/**
 *  \brief  Filters from a spectrum of states the state, which fullfills certain
 * criteria
 *
 *
 */

class Statefilter {

 public:
  void Initialize(tools::Property& options);
  void setLogger(Logger* log) { _log = log; }
  void setInitialState(const QMState& state) { _statehist.push_back(state); }
  void PrintInfo() const;
  
  QMState InitialState() const { return _statehist[0]; }
  QMState CalcStateAndUpdate(const Orbitals& orbitals);
  QMState CalcState(const Orbitals& orbitals) const;

  void WriteToCpt(CheckpointWriter& w) const;

  void ReadFromCpt(CheckpointReader& r);

 private:
  std::vector<int> OscFilter(const Orbitals& orbitals) const;
  std::vector<int> LocFilter(const Orbitals& orbitals) const;
  std::vector<int> DeltaQFilter(const Orbitals& orbitals) const;
  std::vector<int> OverlapFilter(const Orbitals& orbitals) const;
  std::vector<int> OverlapFilterBSE(const Orbitals& orbitals) const;
  std::vector<int> WassersteinFilter(const Orbitals& orbitals) const;
 //std::vector<int> VariationalFilter(const Orbitals& orbitals) const;
 std::vector<int> DensityFilter(const Orbitals& orbitals) const;
std::vector<int> MAE(const Orbitals& orbitals) const;
std::vector<int> MSE(const Orbitals& orbitals) const;
  Eigen::VectorXd CalculateOverlap(const Orbitals& orbitals) const;
  Eigen::VectorXd CalculateOverlapBSE(const Orbitals& orbitals) const;
  Eigen::VectorXd CalculateDNorm(const Orbitals& orbitals) const;
  Eigen::VectorXd CalculateWassersteinNorm(const Orbitals& orbitals) const;
  void UpdateLastCoeff(const Orbitals& orbitals);
  void UpdateLastCoeff_matrix(const Orbitals& orbitals);
  void UpdateLastDmat(const Orbitals& orbitals);
  void UpdateLastBSE_R(const Orbitals& orbitals);
  void UpdateLastBSE_AR(const Orbitals& orbitals);
  void UpdateLastBSE_energy(const Orbitals& orbitals);
  void calculateCube(const Orbitals& orbitals,Eigen::MatrixXd mat, std::string fileout) const;
  Eigen::MatrixXd CalcOrthoCoeffs(const Orbitals& orbitals) const;

  std::vector<int> CollapseResults(
      std::vector<std::vector<int> >& results) const;
  std::vector<int> ComparePairofVectors(std::vector<int>& vec1,
                                        std::vector<int>& vec2) const;

  
  Eigen::VectorXd EvaluateBasisAtPosition(const AOBasis& dftbasis,
                                                 const Eigen::Vector3d& pos) const;
  Logger* _log;

  std::vector<QMState> _statehist;

  bool _use_oscfilter = false;
  double _oscthreshold = 0.0;

  bool _use_overlapfilter = false;
  bool _use_overlapfilter_bse = false;
  bool _use_variationalfilter = false;
  bool _use_densityfilter =false;
  bool _use_maefilter = false;
  bool _use_msefilter = false;
  bool _use_wasserstein = false;
  
  double _overlapthreshold = 0.0;
  double _variationalthreshold = 0.0;
  double _dmatthreshold = 0.0;
  
  Eigen::VectorXd _laststatecoeff; 
  Eigen::MatrixXd _laststatecoeff_mat; 
  Eigen::VectorXd _lastbse_R;
  Eigen::VectorXd _lastbse_AR;
  double _lastbseenergy;
  Eigen::MatrixXd _lastdmat;
  
  
  bool _use_localisationfilter = false;
  std::vector<QMFragment<BSE_Population> > _fragment_loc;
  double _loc_threshold = 0.0;

  bool _use_dQfilter = false;
  std::vector<QMFragment<BSE_Population> > _fragment_dQ;
  double _dQ_threshold = 0.0;
  double _padding = 3.5;
  int _xsteps = 40;
  int _ysteps = 40;
  int _zsteps = 40;
};
}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_STATEFILTER_H
