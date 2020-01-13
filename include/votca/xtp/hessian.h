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
#ifndef VOTCA_XTP_HESSIAN_H
#define VOTCA_XTP_HESSIAN_H

#include <stdio.h>
#include <votca/xtp/gwbseengine.h>
#include <votca/xtp/logger.h>
#include <votca/xtp/qmatom.h>
#include <votca/xtp/segment.h>

namespace votca {
namespace xtp {

class StateTracker;

class Hessian {
 public:
  Forces(GWBSEEngine& gwbse_engine, const StateTracker& tracker)
      : _gwbse_engine(gwbse_engine), _tracker(tracker){};

  void Initialize(tools::Property& options);
  void Calculate(const Orbitals& orbitals);

  void setLog(Logger* pLog) { _pLog = pLog; }

  const Eigen::MatrixX3d& GetHessian() const { return _hessian; };
  void Report() const;

 private:
  Eigen::Matrix3d NumHessianCentral(Orbitals orbitals, Index atom_index_a, Index atom_index_b);
  
  double _displacementt;
  
  GWBSEEngine& _gwbse_engine;
  const StateTracker& _tracker;
  
  Eigen::MatrixXd _hessian;
  Logger* _pLog;
};

}  // namespace xtp
}  // namespace votca
#endif  // VOTCA_XTP_FORCES_H
