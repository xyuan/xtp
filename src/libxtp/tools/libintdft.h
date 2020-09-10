
#pragma once
#ifndef VOTCA_XTP_LIBINTDFT_PRIVATE_H
#define VOTCA_XTP_LIBINTDFT_PRIVATE_H

// VOTCA includes
#include <votca/tools/constants.h>

// Local VOTCA includes
#include "votca/xtp/ERIs.h"
#include "votca/xtp/convergenceacc.h"
#include "votca/xtp/ecpaobasis.h"
#include "votca/xtp/logger.h"
#include "votca/xtp/qmtool.h"
#include "votca/xtp/staticsite.h"
#include "votca/xtp/vxc_grid.h"
#include "votca/xtp/vxc_potential.h"

// VOTCA includes
#include <votca/tools/property.h>vo

namespace votca {
namespace xtp {

class LibintDFT final : public QMTool {
 public:
  LibintDFT() = default;

  ~LibintDFT() final = default;

  std::string Identify() final { return "libintdft"; }

  void Initialize(const tools::Property& user_options) final;
  bool Evaluate() final;

 private:
  Logger _log;
  Index _numofelectrons = 0;

  // AO Matrices
  AOOverlap _dftAOoverlap;
  void SetupInvariantMatrices();
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_LIBINTDFT_PRIVATE_H