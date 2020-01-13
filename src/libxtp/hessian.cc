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

#include <boost/format.hpp>

#include <votca/tools/elements.h>

#include "votca/xtp/statetracker.h"
#include <votca/xtp/atom.h>
#include <votca/xtp/hessian.h>

namespace votca {
namespace xtp {

using std::flush;
void Hessian::Initialize(tools::Property& options) {

  std::vector<std::string> choices = {"central"};
  _hessian_method =
      options.ifExistsAndinListReturnElseThrowRuntimeError<std::string>(
          ".method", choices);
  _displacementt = options.ifExistsReturnElseReturnDefault<double>(
      ".displacement", 0.001);  // Angstrom
  _displacementt *= tools::conv::ang2bohr;

  return;
}

void Hessian::Calculate(const Orbitals& orbitals) {

  Index natoms = orbitals.QMAtoms().size();
  _hessian = Eigen::MatrixXd::Zero(3 * natoms, 3 * natoms);
  for (Index atom_index_a = 0; atom_index_a < natoms; atom_index++) {
    for (Index atom_index_b = 0; atom_index_b < natoms; atom_index++) {
      XTP_LOG(Log::debug, *_pLog)
          << "HESSIAN--DEBUG working on atoms " << atom_index_a << '\t'
          << atom_index_b << flush;
      Eigen::MatrixX3d ab_hessian = Eigen::Matrix3d::Zero();
      // Calculate Hessian on this atom
      if (_hessian_method == "central") {
        ab_hessian = NumHessianCentral(orbitals, atom_index_a, atom_index_b);
      }
      _hessian.block<3, 3>(3 * atom_index_a, 3 * atom_index_b) = ab_hessian;
    }
  }
  return;
}

void Hessian::Report() const {

  Index natoms = orbitals.QMAtoms().size();
  XTP_LOG(Log::error, *_pLog)
      << (boost::format(" ---- Hessian (Hartree/Bohr^2)   ")).str() << flush;
  XTP_LOG(Log::error, *_pLog)
      << (boost::format("      %1$s differences   ") % _hessian_method).str()
      << flush;
  XTP_LOG(Log::error, *_pLog)
      << (boost::format("      displacement %1$1.4f Angstrom   ") %
          (_displacementt * tools::conv::bohr2ang))
             .str()
      << flush;
  XTP_LOG(Log::error, *_pLog) << (boost::format(" H(a,b) ")).str() << flush;

  for (Index i = 0; i < natoms; i++) {
    for (Index j = 0; i < natoms; i++) {
      XTP_LOG(Log::error, *_pLog)
          << _hessian.block<3, 3>(3 * i, 3 * j) << flush;
    }
  }
  return;
}

Eigen::Matrix3d Hessian::NumHessianCentral(Orbitals orbitals,
                                           Index atom_index_a,
                                           Index atom_index_b) {
  Eigen::Matrix3d hessian = Eigen::Matrix3d::Zero();
  // Calculate the case when no displacement has been done
  _gwbse_engine.ExcitationEnergies(orbitals);
  double e = orbitals.getTotalStateEnergy(_tracker.CalcState(orbitals));

  const Eigen::Vector3d current_pos_a =
      orbitals.QMAtoms()[atom_index_a].getPos();
  const Eigen::Vector3d current_pos_b =
      orbitals.QMAtoms()[atom_index_b].getPos();
  for (Index i_cart = 0; i_cart < 3; i_cart++) {
    for (Index j_cart = 0; j_cart < 3; j_cart++) {
      XTP_LOG(Log::debug, *_pLog)
          << "HESSIAN--DEBUG \t Cartesian component " << i_cart
          << '\t Cartesian component' << j_cart << flush;
      Eigen::Vector3d displacement_vec = Eigen::Vector3d::Zero();
      displacement_vec_a[i_cart] = _displacementt;
      displacement_vec_b[j_cart] = _displacementt;

      // For each (i,j) I need 8 numbers. Some combination of them will give the
      // hessian for atom a and b

      // E(a_i+h,b_j)
      Eigen::Vector3d pos_displaced_a = current_pos_a + displacement_vec_a;
      orbitals.QMAtoms()[atom_index_a].setPos(pos_displaced_a);
      _gwbse_engine.ExcitationEnergies(orbitals);
      double h = orbitals.getTotalStateEnergy(_tracker.CalcState(orbitals));
      // E(a_i+h,b_j+h)
      Eigen::Vector3d pos_displaced_b = current_pos_b + displacement_vec_b;
      orbitals.QMAtoms()[atom_index_b].setPos(pos_displaced_b);
      _gwbse_engine.ExcitationEnergies(orbitals);
      double i = orbitals.getTotalStateEnergy(_tracker.CalcState(orbitals));
      // E(a_i+h,b_j-h)
      Eigen::Vector3d pos_displaced_b = current_pos_b - displacement_vec_b;
      orbitals.QMAtoms()[atom_index_b].setPos(pos_displaced_b);
      _gwbse_engine.ExcitationEnergies(orbitals);
      double g = orbitals.getTotalStateEnergy(_tracker.CalcState(orbitals));

      // E(a_i-h,b_j)
      Eigen::Vector3d pos_displaced_a = current_pos_a - displacement_vec_a;
      orbitals.QMAtoms()[atom_index_a].setPos(pos_displaced_a);
      orbitals.QMAtoms()[atom_index_b].setPos(current_pos_b);
      _gwbse_engine.ExcitationEnergies(orbitals);
      double b = orbitals.getTotalStateEnergy(_tracker.CalcState(orbitals));
      // E(a_i-h,b_j+h)
      Eigen::Vector3d pos_displaced_b = current_pos_b + displacement_vec_b;
      orbitals.QMAtoms()[atom_index_b].setPos(pos_displaced_b);
      _gwbse_engine.ExcitationEnergies(orbitals);
      double c = orbitals.getTotalStateEnergy(_tracker.CalcState(orbitals));
      // E(a_i-h,b_j-h)
      Eigen::Vector3d pos_displaced_b = current_pos_b - displacement_vec_b;
      orbitals.QMAtoms()[atom_index_b].setPos(pos_displaced_b);
      _gwbse_engine.ExcitationEnergies(orbitals);
      double a = orbitals.getTotalStateEnergy(_tracker.CalcState(orbitals));

      // E(a_i,b_j+h)
      orbitals.QMAtoms()[atom_index_a].setPos(current_pos_a);
      Eigen::Vector3d pos_displaced_b = current_pos_b + displacement_vec_b;
      orbitals.QMAtoms()[atom_index_b].setPos(pos_displaced_b);
      _gwbse_engine.ExcitationEnergies(orbitals);
      double f = orbitals.getTotalStateEnergy(_tracker.CalcState(orbitals));
      // E(a_i,b_j-h)
      Eigen::Vector3d pos_displaced_b = current_pos_b - displacement_vec_b;
      orbitals.QMAtoms()[atom_index_b].setPos(pos_displaced_b);
      _gwbse_engine.ExcitationEnergies(orbitals);
      double d = orbitals.getTotalStateEnergy(_tracker.CalcState(orbitals));

      // Calculate the Hessian
      hessian(i_cart, i_cart) =
          0.25 * (b - 2 * e + h) / (_displacementt * _displacementt);
      hessian(j_cart, j_cart) =
          0.25 * (d - 2 * e + f) / (_displacementt * _displacementt);
      hessian(i_cart, j_cart) =
          0.25 * (a - g - c - i) / (_displacementt * _displacementt);
      // restore original coordinate into orbital
      orbitals.QMAtoms()[atom_index_a].setPos(current_pos_a);
      orbitals.QMAtoms()[atom_index_b].setPos(current_pos_b);
    }
  }
  return hessian;
}

}  // namespace xtp
}  // namespace votca
