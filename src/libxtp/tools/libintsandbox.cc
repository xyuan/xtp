#include "libintsandbox.h"

// Third party includes
#include <boost/algorithm/string.hpp>
#include <fstream>

// Local VOTCA includes
#include "votca/xtp/molden.h"
#include "votca/xtp/orbitals.h"
#include "votca/xtp/orbreorder.h"
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

  // Setup VOTCA basis
  QMMolecule atoms("", 0);
  atoms.LoadFromFile(_job_name + ".xyz");
  BasisSet bs;
  bs.Load("def2-tzvp");
  AOBasis basis;
  basis.Fill(bs, atoms);

  // Convert to Libint Basis
  std::vector<libint2::Shell> libintBasis = make_libint_basis(basis);

  // Computer Overlap with LIBINT
  libint2::initialize();

  auto S = compute_1body_ints(libintBasis, libint2::Operator::overlap);

  XTP_LOG(Log::error, _log)
      << "\n\tOverlap Integrals Libint BEFORE reordering:\n";
  XTP_LOG(Log::error, _log) << S << std::endl;

  Eigen::MatrixXd S_Xd = static_cast<Eigen::MatrixXd>(S);

  OrbReorder reorder(_libint_reorder, _libint_multipliers);

  reorder.reorderRowsAndCols(S_Xd, basis);

  XTP_LOG(Log::error, _log)
      << "\n\tOverlap Integrals Libint AFTER reordering:\n";
  XTP_LOG(Log::error, _log) << S_Xd << std::endl;

  // Compute Overlap with VOTCA
  AOOverlap dftAOoverlap;
  dftAOoverlap.Fill(basis);

  XTP_LOG(Log::error, _log) << "\n\tOverlap Integrals VOTCA:\n";
  XTP_LOG(Log::error, _log) << dftAOoverlap.Matrix() << std::endl;

  XTP_LOG(Log::error, _log) << std::endl;

  // Compare
  XTP_LOG(Log::error, _log)
      << "Matrices are approximately equal: "
      << S_Xd.isApprox(dftAOoverlap.Matrix(), 1e-5) << std::endl;
  libint2::finalize();
  return true;
}

size_t LibintSandbox::nbasis(const std::vector<libint2::Shell>& shells) {
  size_t n = 0;
  for (const auto& shell : shells) n += shell.size();
  return n;
}

size_t LibintSandbox::max_nprim(const std::vector<libint2::Shell>& shells) {
  size_t n = 0;
  for (auto shell : shells) n = std::max(shell.nprim(), n);
  return n;
}

int LibintSandbox::max_l(const std::vector<libint2::Shell>& shells) {
  int l = 0;
  for (auto shell : shells)
    for (auto c : shell.contr) l = std::max(c.l, l);
  return l;
}

LibintSandbox::Matrix LibintSandbox::compute_1body_ints(
    const std::vector<libint2::Shell>& shells, libint2::Operator obtype,
    const std::vector<libint2::Atom>& atoms) {
  using libint2::Engine;
  using libint2::Operator;
  using libint2::Shell;

  const auto n = nbasis(shells);
  LibintSandbox::Matrix result(n, n);

  // construct the overlap integrals engine
  libint2::Engine engine(obtype, max_nprim(shells), max_l(shells), 0);
  // nuclear attraction ints engine needs to know where the charges sit ...
  // the nuclei are charges in this case; in QM/MM there will also be classical
  // charges
  if (obtype == libint2::Operator::nuclear) {
    std::vector<std::pair<real_t, std::array<real_t, 3>>> q;
    for (const auto& atom : atoms) {
      q.push_back({static_cast<real_t>(atom.atomic_number),
                   {{atom.x, atom.y, atom.z}}});
    }
    engine.set_params(q);
  }

  auto shell2bf = map_shell_to_basis_function(shells);

  // buf[0] points to the target shell set after every call  to engine.compute()
  const auto& buf = engine.results();

  // loop over unique shell pairs, {s1,s2} such that s1 >= s2
  // this is due to the permutational symmetry of the real integrals over
  // Hermitian operators: (1|2) = (2|1)
  for (auto s1 = 0; s1 != shells.size(); ++s1) {

    auto bf1 = shell2bf[s1];  // first basis function in this shell
    auto n1 = shells[s1].size();

    for (auto s2 = 0; s2 <= s1; ++s2) {

      auto bf2 = shell2bf[s2];
      auto n2 = shells[s2].size();

      // compute shell pair
      engine.compute(shells[s1], shells[s2]);

      // "map" buffer to a const Eigen Matrix, and copy it to the corresponding
      // blocks of the result
      Eigen::Map<const Matrix> buf_mat(buf[0], n1, n2);
      result.block(bf1, bf2, n1, n2) = buf_mat;
      if (s1 != s2)  // if s1 >= s2, copy {s1,s2} to the corresponding {s2,s1}
                     // block, note the transpose!
        result.block(bf2, bf1, n2, n1) = buf_mat.transpose();
    }
  }

  return result;
}

std::vector<size_t> LibintSandbox::map_shell_to_basis_function(
    const std::vector<libint2::Shell>& shells) {
  std::vector<size_t> result;
  result.reserve(shells.size());

  size_t n = 0;
  for (auto shell : shells) {
    result.push_back(n);
    n += shell.size();
  }

  return result;
}

std::vector<libint2::Shell> LibintSandbox::make_libint_basis(
    const AOBasis& aobasis) {
  std::vector<libint2::Shell> shells;

  for (const auto& shell : aobasis) {
    libint2::svector<real_t> decays;
    libint2::svector<libint2::Shell::real_t> contractions;
    const Eigen::Vector3d& pos = shell.getPos();
    for (const auto& primitive : shell) {
      decays.push_back(primitive.getDecay());
      contractions.push_back(primitive.getContraction());
    }

    shells.push_back({decays,
                      {// compute integrals in sphericals
                       {static_cast<int>(shell.getL()), true, contractions}},
                      // Atomic Coordinates
                      {{pos[0], pos[1], pos[2]}}});
  }
  return shells;
}  // namespace xtp

}  // namespace xtp
}  // namespace votca