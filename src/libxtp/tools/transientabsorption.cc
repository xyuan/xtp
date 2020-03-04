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

#include "transientabsorption.h"
#include <iostream>
#include <numeric>
#include <stdio.h>
#include <vector>
#include <votca/tools/constants.h>
#include <votca/tools/elements.h>
#include "votca/xtp/aomatrix.h"
#include "votca/xtp/aomatrix3d.h"

namespace votca {
namespace xtp {

void TransientAbsorption::Initialize(tools::Property& options) {
  std::string key = "options." + Identify();

  _orbfile =
      options.ifExistsReturnElseThrowRuntimeError<std::string>(key + ".input");
  _output_file =
      options.ifExistsReturnElseThrowRuntimeError<std::string>(key + ".output");

  _orbitals.ReadFromCpt(_orbfile);
  XTP_LOG(Log::info, _log) << "Reading serialized QM data from " << _orbfile
                           << std::flush;
}

bool TransientAbsorption::Evaluate() {

  CalcSingletTransitionDipole();

  return true;
}



void TransientAbsorption::CalcInterlevelDipoles(){
  const Eigen::MatrixXd& dft_orbitals = _orbitals.MOs().eigenvectors;
  AOBasis basis = SetupDftBasis();

  AODipole dft_dipole;
  dft_dipole.Fill(basis);

}

/*
 * Computes transition dipoles between excited states (singlets only)
 */
void TransientAbsorption::CalcSingletTransitionDipole() {

  Index numofstates = _orbitals.BSESinglets.eigenvalues().size();
  Index ctotal = _orbitals.getBSEctotal();
  Index vtotal = _orbitals.getBSEvtotal();

  // create space to store the results
  _transition_dipoles.resize(0);
  _transition_dipoles.reserve(numofstates * (numofstates + 1) / 2);

  for (Index i_level = 0; i_level < numofstates - 1; i_level++)
    for (Index i_exc = i_level; i_exc < numofstates; i_exc++) {

      Eigen::VectorXd coeffs_level =
          _orbitals.BSESinglets.eigenvectors().col(i_level);
      Eigen::VectorXd coeffs_exc =
          _orbitals.BSESinglets.eigenvectors().col(i_exc);
      if (!_orbitals.BSESinglets.eigenvectors2().empty()) {  // check if TDA is
                                                             // used or not
        coeffs_level += _orbitals.BSESinglets.eigenvectors2().col(i_level);
        coeffs_exc += _orbitals.BSESinglets.eigenvectors2().col(i_exc);
      }
      // map the long row (array) of coeffs to a matrix form
      Eigen::Map<Eigen::MatrixXd> mat_exc(coeffs_exc.data(), ctotal, vtotal);
      Eigen::Map<Eigen::MatrixXd> mat_level(coeffs_level.data(), ctotal,
                                            vtotal);
      // creat a place to store the result
      Eigen::Vector3d tdipole = Eigen::Vector3d::Zero();

      for (Index i = 0; i < 3; i++) {
        // Compute the transition dipole moment seperately for the x, y and z
        // coord
      }
    }
}

}  // namespace xtp
}  // namespace votca