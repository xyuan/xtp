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
#ifndef _VOTCA_XTP_TRANSIENTABSORPTION_H
#define _VOTCA_XTP_TRANSIENTABSORPTION_H

#include <stdio.h>

#include <votca/xtp/logger.h>
#include <votca/xtp/orbitals.h>
#include <votca/xtp/qmstate.h>
#include <votca/xtp/qmtool.h>
#include <votca/xtp/eigen.h>



namespace votca {
namespace xtp {

class TransientAbsorption : public QMTool {
 public:
  TransientAbsorption() = default;

  ~TransientAbsorption() override = default;

  std::string Identify() final { return "transientabsorption"; }

  void Initialize(tools::Property& options) override;

  bool Evaluate() override;

 private:
  std::string _orbfile;
  std::string _output_file;
  Logger _log;
  Orbitals _orbitals;
  std::vector<Eigen::Vector3d> _transition_dipoles;

  void CalcSingletTransitionDipole();
  void CalcInterlevelDipoles();
  
};

}  // namespace xtp
}  // namespace votca

#endif