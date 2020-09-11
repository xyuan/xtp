
#pragma once
#ifndef VOTCA_XTP_LIBINTSANDBOX_PRIVATE_H
#define VOTCA_XTP_LIBINTSANDBOX_PRIVATE_H

// Point Libint to the data folder with the basisset info
#define DATADIR std::string(getenv("VOTCASHARE")) + "/xtp/data"

// VOTCA includes
#include <votca/tools/constants.h>

// Local VOTCA includes
#include "votca/xtp/ERIs.h"
#include "votca/xtp/convergenceacc.h"
#include "votca/xtp/ecpaobasis.h"
#include "votca/xtp/logger.h"
#include "votca/xtp/qmmolecule.h"
#include "votca/xtp/qmtool.h"
#include "votca/xtp/staticsite.h"
#include "votca/xtp/vxc_grid.h"
#include "votca/xtp/vxc_potential.h"

// Libint include
#include <libint2.hpp>

// VOTCA includes
#include <votca/tools/property.h>

namespace votca {
namespace xtp {

class LibintSandbox final : public QMTool {
 public:
  LibintSandbox() = default;

  ~LibintSandbox() final = default;

  std::string Identify() final { return "libintsandbox"; }

  void Initialize(const tools::Property& user_options) final;
  bool Evaluate() final;

 private:
  using real_t = libint2::scalar_type;
  using Matrix =
      Eigen::Matrix<real_t, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
  Logger _log;
  Index _numofelectrons = 0;

  // Libint
  Matrix compute_1body_ints(
      const std::vector<libint2::Shell>& shells, libint2::Operator obtype,
      const std::vector<libint2::Atom>& atoms = std::vector<libint2::Atom>());
  size_t nbasis(const std::vector<libint2::Shell>& shells);
  size_t max_nprim(const std::vector<libint2::Shell>& shells);
  int max_l(const std::vector<libint2::Shell>& shells);
  std::vector<size_t> map_shell_to_basis_function(
      const std::vector<libint2::Shell>& shells);

  // VOTCA
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_LIBINTSANDBOX_PRIVATE_H