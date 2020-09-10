#include "libintdft.h"

// Third party includes
#include <boost/algorithm/string.hpp>
#include <fstream>
#include <libint2.hpp>

// Local VOTCA includes
#include "votca/xtp/molden.h"
#include "votca/xtp/orbitals.h"
#include <votca/tools/constants.h>

namespace votca {
namespace xtp {

void LibintDFT::Initialize(const tools::Property& user_options) {
  tools::Property options =
      LoadDefaultsAndUpdateWithUserOptions("xtp", user_options);

  _job_name = options.ifExistsReturnElseReturnDefault<std::string>("job_name",
                                                                   _job_name);
}

bool LibintDFT::Evaluate() {
  _log.setReportLevel(Log::current_level);
  _log.setMultithreading(true);
  _log.setCommonPreface("\n... ...");

  XTP_LOG(Log::error, _log) << "Starting Libint DFT Calculation " << std::flush;
  libint2::initialize();

  // Load atoms
  std::ifstream input_file(_job_name + ".xyz");
  std::vector<libint2::Atom> atoms = libint2::read_dotxyz(input_file);

  // Load basisset
  libint2::BasisSet obs("def2-tzvp", atoms);

  // Get Number of electrons
  for (Index i = 0; i < atoms.size(); ++i) {
    _numofelectrons += atoms[i].atomic_number;
  }

  SetupInvariantMatrices();

  XTP_LOG(Log::error, _log)
      << TimeStamp() << " Total number of electrons: " << _numofelectrons
      << std::flush;

  /*
  Prepare(orb.QMAtoms());
  Mat_p_Energy H0 = SetupH0(orb.QMAtoms());
  tools::EigenSystem MOs;
  MOs.eigenvalues() = Eigen::VectorXd::Zero(H0.cols());
  MOs.eigenvectors() = Eigen::MatrixXd::Zero(H0.rows(), H0.cols());
  Vxc_Potential<Vxc_Grid> vxcpotential = SetupVxc(orb.QMAtoms());
  ConfigOrbfile(orb);

  if (_with_guess) {
    XTP_LOG(Log::error, *_pLog)
        << TimeStamp() << " Reading guess from orbitals object/file" << flush;
    MOs = orb.MOs();
    MOs.eigenvectors() = OrthogonalizeGuess(MOs.eigenvectors());
  } else {
    XTP_LOG(Log::error, *_pLog)
        << TimeStamp() << " Setup Initial Guess using: " << _initial_guess
        << flush;
    if (_initial_guess == "independent") {
      MOs = IndependentElectronGuess(H0);
    } else if (_initial_guess == "atom") {
      MOs = ModelPotentialGuess(H0, orb.QMAtoms(), vxcpotential);
    } else {
      throw std::runtime_error("Initial guess method not known/implemented");
    }
  }

  Eigen::MatrixXd Dmat = _conv_accelerator.DensityMatrix(MOs);
  XTP_LOG(Log::info, *_pLog)
      << TimeStamp() << " Guess Matrix gives N=" << std::setprecision(9)
      << Dmat.cwiseProduct(_dftAOoverlap.Matrix()).sum() << " electrons."
      << flush;

  XTP_LOG(Log::error, *_pLog) << TimeStamp() << " STARTING SCF cycle" << flush;
  XTP_LOG(Log::error, *_pLog)
      << " ----------------------------------------------"
         "----------------------------"
      << flush;

  for (Index this_iter = 0; this_iter < _max_iter; this_iter++) {
    XTP_LOG(Log::error, *_pLog) << flush;
    XTP_LOG(Log::error, *_pLog) << TimeStamp() << " Iteration " << this_iter + 1
                                << " of " << _max_iter << flush;

    Mat_p_Energy e_vxc = vxcpotential.IntegrateVXC(Dmat);
    XTP_LOG(Log::info, *_pLog)
        << TimeStamp() << " Filled DFT Vxc matrix " << flush;

    Mat_p_Energy ERIs = CalculateERIs(Dmat);
    Eigen::MatrixXd H = H0.matrix() + ERIs.matrix() + e_vxc.matrix();
    double Eone = Dmat.cwiseProduct(H0.matrix()).sum();
    double Etwo = 0.5 * ERIs.energy() + e_vxc.energy();
    double exx = 0.0;
    if (_ScaHFX > 0) {
      Mat_p_Energy EXXs = CalcEXXs(MOs.eigenvectors(), Dmat);
      XTP_LOG(Log::info, *_pLog)
          << TimeStamp() << " Filled DFT Electron exchange matrix" << flush;
      H -= 0.5 * _ScaHFX * EXXs.matrix();
      exx = -_ScaHFX / 4 * EXXs.energy();
    }
    Etwo += exx;
    double totenergy = Eone + H0.energy() + Etwo;
    XTP_LOG(Log::info, *_pLog) << TimeStamp() << " Single particle energy "
                               << std::setprecision(12) << Eone << flush;
    XTP_LOG(Log::info, *_pLog) << TimeStamp() << " Two particle energy "
                               << std::setprecision(12) << Etwo << flush;
    XTP_LOG(Log::info, *_pLog)
        << TimeStamp() << std::setprecision(12) << " Local Exc contribution "
        << e_vxc.energy() << flush;
    if (_ScaHFX > 0) {
      XTP_LOG(Log::info, *_pLog)
          << TimeStamp() << std::setprecision(12)
          << " Non local Ex contribution " << exx << flush;
    }
    XTP_LOG(Log::error, *_pLog) << TimeStamp() << " Total Energy "
                                << std::setprecision(12) << totenergy << flush;

    Dmat = _conv_accelerator.Iterate(Dmat, H, MOs, totenergy);

    PrintMOs(MOs.eigenvalues(), Log::info);

    XTP_LOG(Log::info, *_pLog) << "\t\tGAP "
                               << MOs.eigenvalues()(_numofelectrons / 2) -
                                      MOs.eigenvalues()(_numofelectrons / 2 - 1)
                               << flush;

    if (_conv_accelerator.isConverged()) {
      XTP_LOG(Log::error, *_pLog)
          << TimeStamp() << " Total Energy has converged to "
          << std::setprecision(9) << _conv_accelerator.getDeltaE()
          << "[Ha] after " << this_iter + 1
          << " iterations. DIIS error is converged up to "
          << _conv_accelerator.getDIIsError() << flush;
      XTP_LOG(Log::error, *_pLog)
          << TimeStamp() << " Final Single Point Energy "
          << std::setprecision(12) << totenergy << " Ha" << flush;
      XTP_LOG(Log::error, *_pLog) << TimeStamp() << std::setprecision(12)
                                  << " Final Local Exc contribution "
                                  << e_vxc.energy() << " Ha" << flush;
      if (_ScaHFX > 0) {
        XTP_LOG(Log::error, *_pLog)
            << TimeStamp() << std::setprecision(12)
            << " Final Non Local Ex contribution " << exx << " Ha" << flush;
      }

      Mat_p_Energy EXXs = CalcEXXs(MOs.eigenvectors(), Dmat);
      exx = -1.0 / 4.0 * EXXs.energy();
      XTP_LOG(Log::error, *_pLog) << TimeStamp() << std::setprecision(12)
                                  << " EXX energy " << exx << " Ha" << flush;

      XTP_LOG(Log::error, *_pLog)
          << TimeStamp() << std::setprecision(12) << " Final EXX Total energy "
          << totenergy - e_vxc.energy() + (1.0 - _ScaHFX) * exx << " Ha"
          << flush;

      PrintMOs(MOs.eigenvalues(), Log::error);
      orb.setQMEnergy(totenergy);
      orb.MOs() = MOs;
      CalcElDipole(orb);
      break;
    } else if (this_iter == _max_iter - 1) {
      XTP_LOG(Log::error, *_pLog)
          << TimeStamp() << " DFT calculation has not converged after "
          << _max_iter
          << " iterations. Use more iterations or another convergence "
             "acceleration scheme."
          << std::flush;
      return false;
    }
  }
  return true;
}
  */

  libint2::finalize();

  XTP_LOG(Log::error, _log) << std::endl;
  return true;
}

void LibintDFT::SetupInvariantMatrices() {

  _dftAOoverlap.Fill(_dftbasis);

  XTP_LOG(Log::info, *_pLog)
      << TimeStamp() << " Filled DFT Overlap matrix." << flush;

  _conv_opt.numberofelectrons = _numofelectrons;
  _conv_accelerator.Configure(_conv_opt);
  _conv_accelerator.setLogger(_pLog);
  _conv_accelerator.setOverlap(_dftAOoverlap, 1e-8);
  _conv_accelerator.PrintConfigOptions();

  if (_four_center_method == "RI") {
    // prepare invariant part of electron repulsion integrals
    _ERIs.Initialize(_dftbasis, _auxbasis);
    XTP_LOG(Log::info, *_pLog)
        << TimeStamp() << " Inverted AUX Coulomb matrix, removed "
        << _ERIs.Removedfunctions() << " functions from aux basis" << flush;
    XTP_LOG(Log::error, *_pLog)
        << TimeStamp()
        << " Setup invariant parts of Electron Repulsion integrals " << flush;
  } else {

    if (_four_center_method == "cache") {

      XTP_LOG(Log::info, *_pLog)
          << TimeStamp() << " Calculating 4c integrals. " << flush;
      _ERIs.Initialize_4c_small_molecule(_dftbasis);
      XTP_LOG(Log::error, *_pLog)
          << TimeStamp() << " Calculated 4c integrals. " << flush;
    }

    if (_with_screening && _four_center_method == "direct") {
      XTP_LOG(Log::info, *_pLog)
          << TimeStamp() << " Calculating 4c diagonals. " << flush;
      _ERIs.Initialize_4c_screening(_dftbasis, _screening_eps);
      XTP_LOG(Log::info, *_pLog)
          << TimeStamp() << " Calculated 4c diagonals. " << flush;
    }
  }

  return;
}

}  // namespace xtp
}  // namespace votca