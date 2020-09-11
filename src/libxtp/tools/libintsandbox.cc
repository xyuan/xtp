#include "libintsandbox.h"

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

void LibintSandbox::Initialize(const tools::Property& user_options) {
  tools::Property options =
      LoadDefaultsAndUpdateWithUserOptions("xtp", user_options);

  _job_name = options.ifExistsReturnElseReturnDefault<std::string>("job_name",
                                                                   _job_name);
}

bool LibintSandbox::Evaluate() {
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

  XTP_LOG(Log::error, _log)
      << TimeStamp() << " Total number of electrons: " << _numofelectrons
      << std::flush;

  libint2::finalize();

  XTP_LOG(Log::error, _log) << std::endl;
  return true;
}

}  // namespace xtp
}  // namespace votca