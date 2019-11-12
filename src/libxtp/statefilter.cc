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

#include "votca/xtp/aomatrix.h"
#include <numeric>
#include <votca/xtp/statefilter.h>
#include <eigen3/Eigen/src/Core/Matrix.h>
#include <c++/5/cmath>

namespace votca {
namespace xtp {
using std::flush;

void Statefilter::Initialize(tools::Property& options) {
  if (options.exists("oscillator_strength")) {
    _use_oscfilter = true;
    _oscthreshold = options.ifExistsReturnElseThrowRuntimeError<double>(
        "oscillator_strength");
  }
  if (options.exists("overlap")) {
    _use_overlapfilter = true;
    _overlapthreshold =
        options.ifExistsReturnElseReturnDefault<double>("overlap", 0.0);
  }
  if (options.exists("overlapbse")) {
    _use_overlapfilter_bse = true;
    _overlapthreshold =
        options.ifExistsReturnElseReturnDefault<double>("overlapbse", 0.0);
  }
  if (options.exists("mae")) {
    _use_maefilter = true;
    _overlapthreshold =
        options.ifExistsReturnElseReturnDefault<double>("mae", 0.0);
  }
  if (options.exists("mse")) {
    _use_msefilter = true;
    _overlapthreshold =
        options.ifExistsReturnElseReturnDefault<double>("mse", 0.0);
  }
  
  if (options.exists("wasserstein")) {
    _use_wasserstein = true;
    _overlapthreshold =
        options.ifExistsReturnElseReturnDefault<double>("wasserstein", 0.0);
  }
  
  if (options.exists("density")) {
    _use_densityfilter = true;
    _dmatthreshold =
        options.ifExistsReturnElseReturnDefault<double>("density", 0.0);
  }
  if (options.exists("localisation")) {
    _use_localisationfilter = true;
    std::string indices =
        options.ifExistsReturnElseThrowRuntimeError<std::string>(
            "localisation.fragment");
    QMFragment<BSE_Population> reg =
        QMFragment<BSE_Population>("Fragment", 0, indices);
    _loc_threshold = options.ifExistsReturnElseThrowRuntimeError<double>(
        "localisation.threshold");
    _fragment_loc.push_back(reg);
  }

  if (options.exists("charge_transfer")) {
    _use_dQfilter = true;
    std::string indices =
        options.ifExistsReturnElseThrowRuntimeError<std::string>(
            "charge_transfer.fragment");
    QMFragment<BSE_Population> reg =
        QMFragment<BSE_Population>("Fragment", 0, indices);
    _dQ_threshold = options.ifExistsReturnElseThrowRuntimeError<double>(
        "charge_transfer.threshold");
    _fragment_dQ.push_back(reg);
  }
  
  if (_use_dQfilter && _use_localisationfilter) {
    throw std::runtime_error(
        "Cannot use localisation and charge_transfer filter at the same time.");
  } 
}

void Statefilter::PrintInfo() const {
  XTP_LOG(logDEBUG, *_log) << "Initial state: " << _statehist[0].ToString()
                           << flush;
  if (_statehist.size() > 1) {
    XTP_LOG(logDEBUG, *_log)
        << "Last state: " << _statehist.back().ToString() << flush;
  }
  if (_use_oscfilter) {
    XTP_LOG(logDEBUG, *_log) << "Using oscillator strength filter with threshold "
                             << _oscthreshold << flush;
  }
  if (_use_overlapfilter) {
    if (_overlapthreshold == 0.0) {
      XTP_LOG(logDEBUG, *_log)
          << "Using overlap filter with no threshold " << flush;
    } else {
      XTP_LOG(logDEBUG, *_log)
          << "Using overlap filter with threshold " << _overlapthreshold << flush;
    }
  }
  if (_use_overlapfilter_bse) {
    if (_overlapthreshold == 0.0) {
      XTP_LOG(logDEBUG, *_log)
          << "Using overlap filter BSE with no threshold " << flush;
    } else {
      XTP_LOG(logDEBUG, *_log)
          << "Using overlap filter BSE with threshold " << _overlapthreshold << flush;
    }
  }
    
  if (_use_densityfilter) {
    if (_dmatthreshold == 0.0) {
      XTP_LOG(logDEBUG, *_log)
          << "Using density filter with no threshold " << flush;
    } else {
      XTP_LOG(logDEBUG, *_log)
          << "Using density filter with threshold " << _dmatthreshold << flush;
    }
  }
  if (_use_localisationfilter) {
    XTP_LOG(logDEBUG, *_log)
        << "Using localisation filter for fragment" << _fragment_loc[0].name()
        << " with cutoff " << _loc_threshold << flush;
  }
  if (_use_dQfilter) {
    XTP_LOG(logDEBUG, *_log)
        << "Using Delta Q filter for fragment" << _fragment_dQ[0].name()
        << "with cutoff  " << _dQ_threshold << flush;
  }
  if (_use_oscfilter && _use_dQfilter) {
    XTP_LOG(logDEBUG, *_log) << "WARNING: filtering for optically active CT "
                                "transition - might not make sense... "
                             << flush;
  }
  if (_use_dQfilter + _use_oscfilter + _use_localisationfilter +
          _use_overlapfilter + _use_densityfilter + _use_overlapfilter_bse + _use_maefilter + _use_msefilter + _use_wasserstein<
      1) {
    XTP_LOG(logDEBUG, *_log) << "WARNING: No filter is used " << flush;
  }
  if (_use_maefilter) {
    XTP_LOG(logDEBUG, *_log) << "WARNING: I am combining BSE-Overlap and Density. It might not work as you expect"
                             << flush;
    
    
  }
   if (_use_msefilter) {
    XTP_LOG(logDEBUG, *_log) << "WARNING: I am combining BSE-Overlap and Density. It might not work as you expect"
                             << flush;
    
    
  }
  if (_use_wasserstein) {
    XTP_LOG(logDEBUG, *_log) << "WARNING: I am using Wasserstein" << flush;
  }
}

std::vector<int> Statefilter::ComparePairofVectors(
    std::vector<int>& vec1, std::vector<int>& vec2) const {
  std::vector<int> result(std::min(vec1, vec2));
  std::vector<int>::iterator it;
  std::sort(vec1.begin(), vec1.end());
  std::sort(vec2.begin(), vec2.end());
  it = std::set_intersection(vec1.begin(), vec1.end(), vec2.begin(), vec2.end(),
                             result.begin());
  result.resize(it - result.begin());
  return result;
}

std::vector<int> Statefilter::CollapseResults(
    std::vector<std::vector<int> >& results) const {
  if (results.size() == 1) {
    return results[0];
  } else {
    std::vector<int> result = results[0];
    for (unsigned i = 1; i < results.size(); i++) {
      result = ComparePairofVectors(result, results[i]);
    }
    return result;
  }
}


Eigen::VectorXd Statefilter::EvaluateBasisAtPosition(const AOBasis& dftbasis,
                                                 const Eigen::Vector3d& pos) const {

  // get value of orbitals at each gridpoint
  Eigen::VectorXd tmat = Eigen::VectorXd::Zero(dftbasis.AOBasisSize());
  for (const AOShell& shell : dftbasis) {
    const double decay = shell.getMinDecay();
    const Eigen::Vector3d& shellpos = shell.getPos();
    Eigen::Vector3d dist = shellpos - pos;
    double distsq = dist.squaredNorm();
    // if contribution is smaller than -ln(1e-10), calc density
    if ((decay * distsq) < 20.7) {
      Eigen::VectorBlock<Eigen::VectorXd> tmat_block =
          tmat.segment(shell.getStartIndex(), shell.getNumFunc());
      shell.EvalAOspace(tmat_block, pos);
    }
  }
  return tmat;
}



QMState Statefilter::CalcState(const Orbitals& orbitals) const {

  if (_use_dQfilter + _use_oscfilter + _use_localisationfilter +
          _use_overlapfilter + _use_densityfilter +_use_overlapfilter_bse + _use_maefilter + _use_msefilter + _use_wasserstein
           <
      1) {
    return _statehist[0];
  }

  std::vector<std::vector<int> > results;
  if (_use_oscfilter) {
    results.push_back(OscFilter(orbitals));
  }
  if (_use_localisationfilter) {
    results.push_back(LocFilter(orbitals));
  }
  if (_use_overlapfilter) {
    results.push_back(OverlapFilter(orbitals));
  }
  if (_use_overlapfilter_bse) {
    results.push_back(OverlapFilterBSE(orbitals));
  }
  
  if (_use_densityfilter) {
    results.push_back(DensityFilter(orbitals));
  }
  if (_use_dQfilter) {
    results.push_back(DeltaQFilter(orbitals));
  }
  if (_use_maefilter ) {
    results.push_back(MAE(orbitals));
  }
  if (_use_msefilter ) {
    results.push_back(MSE(orbitals));
  }
  if (_use_wasserstein ) {
    results.push_back(WassersteinFilter(orbitals));
  }
  std::vector<int> result = CollapseResults(results);
  QMState state;
  if (result.size() < 1) {
    state = _statehist.back();
    XTP_LOG(logDEBUG, *_log)
        << "No State found by filter using last state: " << state.ToString()
        << flush;
  } else {
    state = QMState(_statehist.back().Type(), result[0], false);
    XTP_LOG(logDEBUG, *_log) << "Next State is: " << state.ToString() << flush;
  }
  return state;
}

QMState Statefilter::CalcStateAndUpdate(const Orbitals& orbitals) {
  QMState result = CalcState(orbitals);
  _statehist.push_back(result);
  if (_use_overlapfilter) {
    UpdateLastCoeff(orbitals);
  }
  if (_use_overlapfilter_bse) {
    UpdateLastBSE_R(orbitals);
    UpdateLastCoeff_matrix(orbitals);
  }
  if(_use_densityfilter){
    UpdateLastDmat(orbitals);
  }
  if (_use_maefilter ) {
    UpdateLastBSE_R(orbitals);
    UpdateLastCoeff_matrix(orbitals);
    UpdateLastDmat(orbitals);
    UpdateLastBSE_energy(orbitals);
  }
  if (_use_msefilter ) {
    UpdateLastBSE_R(orbitals);
    UpdateLastCoeff_matrix(orbitals);
    UpdateLastDmat(orbitals);
    UpdateLastBSE_energy(orbitals);
  }
  if (_use_wasserstein ) {
      UpdateLastDmat(orbitals);
      UpdateLastBSE_energy(orbitals);
  }
return result;
}

std::vector<int> Statefilter::OscFilter(const Orbitals& orbitals) const {
  Eigen::VectorXd oscs = orbitals.Oscillatorstrengths();
  std::vector<int> indexes;
  for (int i = 0; i < oscs.size(); i++) {
    if (oscs[i] > _oscthreshold) indexes.push_back(i);
  }
  return indexes;
}

std::vector<int> Statefilter::LocFilter(const Orbitals& orbitals) const {
  std::vector<int> indexes;
  Lowdin low;
  std::vector<QMFragment<BSE_Population> > loc = _fragment_loc;
  low.CalcChargeperFragment(loc, orbitals, _statehist[0].Type());
  const Eigen::VectorXd& popE = loc[0].value().E;
  const Eigen::VectorXd& popH = loc[0].value().H;
  for (int i = 0; i < popE.size(); i++) {
    if (popE[i] > _loc_threshold && popH[i] > _loc_threshold) {
      indexes.push_back(i);
    }
  }
  return indexes;
}

std::vector<int> Statefilter::DeltaQFilter(const Orbitals& orbitals) const {
  std::vector<int> indexes;
  Lowdin low;
  std::vector<QMFragment<BSE_Population> > loc = _fragment_dQ;
  low.CalcChargeperFragment(loc, orbitals, _statehist[0].Type());
  Eigen::VectorXd dq = (loc[0].value().H - loc[0].value().E).cwiseAbs();

  for (int i = 0; i < dq.size(); i++) {
    if (dq[i] > _dQ_threshold) {
      indexes.push_back(i);
    }
  }
  return indexes;
}



// This overlap is wrong: it should be something like coeff_new * overlap * coeff_old
Eigen::VectorXd Statefilter::CalculateOverlap(const Orbitals& orbitals) const {
  Eigen::MatrixXd coeffs = CalcOrthoCoeffs(orbitals);
  Eigen::VectorXd overlap = (coeffs * _laststatecoeff).cwiseAbs2();
  return overlap;
}


Eigen::VectorXd Statefilter::CalculateOverlapBSE(const Orbitals& orbitals) const {
    //Define some useful integers as the numeber of states, the number of occ and virt level
    int nostates=orbitals.NumberofStates(_statehist[0].Type());
    int v =  orbitals.getBSEvmax()+1;
    int c = orbitals.getBSEcmax()+1;
    //Degine the AO overlap which is not step dependent unless you do geo. optimization  
    AOOverlap S_ao ;
    S_ao.Fill(orbitals.SetupDftBasis()); 
    // Get the whole MO matrix at this step (n)
    Eigen::MatrixXd MO_mat = orbitals.MOs().eigenvectors();
    //Compute C as MO_mat(step n-1) * S_ao * MO_mat(step n)
    Eigen::MatrixXd C = _laststatecoeff_mat.transpose() * S_ao.Matrix() * MO_mat;
    
    //Get the diagonal blocks (to be updated with off diagonal block and beyond TDA overlap
    Eigen::MatrixXd M_vv = C.block(0,0,v,v);
    Eigen::MatrixXd M_cc = C.block(v,v,c-v,c-v);
    // Let's define the Reference vector containing the singlet amplitudes at step n-1
    Eigen::VectorXd Ref = _lastbse_R;
    // I have to reshape this vector in a matrix form
    Eigen::MatrixXd A_vc_old = Eigen::Map<Eigen::MatrixXd>(Ref.data(),v,c-v);
    //Compute A_vc * M_cc'
    Eigen::MatrixXd two = A_vc_old *  M_cc;
    
    // Inizialise the BSE overlap vector     
    Eigen::VectorXd overlap_bse=Eigen::VectorXd::Zero(nostates);
    for(int i=0;i<nostates;i++){
        QMState state(_statehist[0].Type(),i,false);
        Eigen::VectorXd a = orbitals.BSESinglets().eigenvectors().col(state.Index());
        Eigen::MatrixXd A_vc = Eigen::Map<Eigen::MatrixXd>(a.data(),v,c-v); 
        
        Eigen::MatrixXd one =  M_vv * A_vc;
        
        Eigen::MatrixXd O = two.transpose() * one;
        double ov = (O).trace();
        overlap_bse(i) = ov; 
         
    }
    std::cout << " \n Testing overlap Singlet \n " << overlap_bse.cwiseAbs() << std::endl;
    
    return overlap_bse.cwiseAbs();
    
}

Eigen::VectorXd Statefilter::CalculateDNorm(const Orbitals& orbitals) const{
    
    int nostates=orbitals.NumberofStates(_statehist[0].Type());
    
    Eigen::VectorXd norm=Eigen::VectorXd::Zero(nostates);
    for(int i=0;i<nostates;i++){
        QMState state(_statehist[0].Type(),i,false);
        Eigen::MatrixXd diff = (orbitals.DensityMatrixFull(state)-_lastdmat);
        norm(i) = diff.norm()/(_lastdmat.norm());
    }
    return norm;    
}

    
//Eigen::VectorXd Statefilter::CalculateWassersteinNorm(const Orbitals& orbitals) const{
//    int nostates = orbitals.NumberofStates(_statehist[0].Type());
//    Eigen::VectorXd wasserstein_norm = Eigen::VectorXd::Zero(nostates);
//    std::ofstream txtout;
//    std::cout << "\n Density Matrix dimension " << _lastdmat.rows() << "\t" << _lastdmat.cols() << std::endl;
//    std::string cmd_mkdir = "mkdir Density";
//    std::system(cmd_mkdir.c_str());
//    std::string path = "./Density/" ;
//    std::cout << "Writing density matrix files in "<< path <<  std::endl;
//    //This generates the file for the old density matrix
//    std::string file_name = path+"lastdensitymatrix.txt" ;
//    txtout.open(file_name);
//    txtout << _lastdmat ;
//    txtout.close();
//    //This loop generates the files for all the new states (the ones to be compared with the old density)
//    for(int i=0;i<nostates;i++){
//        QMState state(_statehist[0].Type(),i,false);
//        file_name = path+"density_s"+std::to_string(i+1)+".txt";
//        txtout.open(file_name);
//        txtout <<  orbitals.DensityMatrixFull(state);
//        txtout.close();
//    }
//    
//    std::string file_out = "wasserstein_out.txt";
//    std::cout << "Running Wasserstein Python script" << std::endl;
//    
//    std::string command = "python wasserstein.py " + std::to_string(nostates) + " " + file_out ;
//    std::system(command.c_str());
//    
//    std::cout << "Printing Wasserstein measure in " << file_out << std::endl   ;
//    std::ifstream ifs( file_out );
//
//    std::vector< double > values;
//    double val;
//    int i = 0;
//    while( ifs >> val ){
//        wasserstein_norm(i) = val;
//        i += 1;
//    }
//    
//    return wasserstein_norm;                
//}


void Statefilter::calculateCube(const Orbitals& orbitals, Eigen::MatrixXd mat, std::string fileout) const {
  // generate cube grid
  const QMMolecule& atoms = orbitals.QMAtoms();
  std::pair<Eigen::Vector3d, Eigen::Vector3d> minmax = atoms.CalcSpatialMinMax();  
  double xstart = minmax.first.x() - _padding;
  double xstop = minmax.second.x() + _padding;
  double ystart = minmax.first.y() - _padding;
  double ystop = minmax.second.y() + _padding;
  double zstart = minmax.first.z() - _padding;
  double zstop = minmax.second.z() + _padding;

  double xincr = (xstop - xstart) / double(_xsteps);
  double yincr = (ystop - ystart) / double(_ysteps);
  double zincr = (zstop - zstart) / double(_zsteps);
  
  // load DFT basis set (element-wise information) from xml file
  BasisSet dftbs;
  dftbs.Load(orbitals.getDFTbasisName());
  // fill DFT AO basis by going through all atoms
  AOBasis dftbasis;
  dftbasis.Fill(dftbs, orbitals.QMAtoms());
  
  
  std::ofstream out;
  out.open(fileout);
  // eval density at cube grid points
  for (int ix = 0; ix <= _xsteps; ix++) {
    double x = xstart + double(ix) * xincr;
    for (int iy = 0; iy <= _ysteps; iy++) {
      double y = ystart + double(iy) * yincr;
      int Nrecord = 0;
      for (int iz = 0; iz <= _zsteps; iz++) {
        double z = zstart + double(iz) * zincr;
        Nrecord++;
        Eigen::Vector3d pos(x, y, z);
        Eigen::VectorXd tmat = EvaluateBasisAtPosition(dftbasis, pos);
        double value = 0.0;
        value = (tmat.transpose() * mat * tmat).value();
        
//        if (Nrecord == 6 || iz == _zsteps) {
//          out << boost::format("%1$E \n") % value;
//          Nrecord = 0;
//        } else {
          out << boost::format("%1$E ") % value;
//        }
      }  // z-component
    }    // y-component
  }  // x-component

  out.close();
  
  return;
}

Eigen::VectorXd Statefilter::CalculateWassersteinNorm(const Orbitals& orbitals) const{
    int nostates = orbitals.NumberofStates(_statehist[0].Type());
    Eigen::VectorXd wasserstein_norm = Eigen::VectorXd::Zero(nostates);
    
    std::string cmd_mkdir = "mkdir Density";
    std::system(cmd_mkdir.c_str());
    std::string path = "./Density/" ;
    std::cout << "Writing density files in "<< path <<  std::endl;
    //This generates the file for the old density 
    std::string file_name = path+"lastdensitymatrix.txt" ;
    //Calculate old density
    calculateCube(orbitals,_lastdmat,file_name);
    //This loop generates the files for all the new states (the ones to be compared with the old density)
    for(int i=0;i<nostates;i++){
        QMState state(_statehist[0].Type(),i,false);
        file_name = path+"density_s"+std::to_string(i+1)+".txt";
        calculateCube(orbitals,orbitals.DensityMatrixFull(state),file_name);
    }
    
    std::string file_out = "wasserstein_out.txt";
    std::cout << "Running Wasserstein Python script" << std::endl;
    
    std::string command = "python wasserstein.py " + std::to_string(nostates) + " " + file_out ;
    std::system(command.c_str());
    
    std::cout << "Printing Wasserstein measure in " << file_out << std::endl   ;
    std::ifstream ifs( file_out );

    std::vector< double > values;
    double val;
    int i = 0;
    while( ifs >> val ){
        wasserstein_norm(i) = val;
        i += 1;
    }
    
    return wasserstein_norm;                
}


Eigen::MatrixXd Statefilter::CalcOrthoCoeffs(const Orbitals& orbitals) const {
  QMStateType type = _statehist[0].Type();
  Eigen::MatrixXd coeffs;
  if (type.isSingleParticleState()) {
    if (type == QMStateType::DQPstate) {
      coeffs = orbitals.CalculateQParticleAORepresentation();
    } else {
      coeffs = orbitals.MOs().eigenvectors();
    }
  } else {
    //throw std::runtime_error("Overlap for excitons not implemented yet");
      std::cout << "Please be careful with this overlap" << std::endl;
  }
  return coeffs;
}

void Statefilter::UpdateLastCoeff(const Orbitals& orbitals) {
  Eigen::MatrixXd ortho_coeffs = CalcOrthoCoeffs(orbitals);
  int offset = 0;
  if (_statehist[0].Type() == QMStateType::DQPstate) {
    offset = orbitals.getGWAmin();
  }
  _laststatecoeff = ortho_coeffs.col(_statehist.back().Index() - offset);
}


void Statefilter::UpdateLastCoeff_matrix(const Orbitals& orbitals) {
     _laststatecoeff_mat = orbitals.MOs().eigenvectors();
  }


void Statefilter::UpdateLastDmat(const Orbitals& orbitals){
    _lastdmat=orbitals.DensityMatrixFull(_statehist.back());
}

void Statefilter::UpdateLastBSE_R(const Orbitals& orbitals){   
    _lastbse_R = 
           orbitals.BSESinglets().eigenvectors().col(_statehist.back().Index());
}

void Statefilter::UpdateLastBSE_AR(const Orbitals& orbitals){
    _lastbse_AR = 
          orbitals.BSESinglets().eigenvectors2().col(_statehist.back().Index());
}

void Statefilter::UpdateLastBSE_energy(const Orbitals& orbitals){
    _lastbseenergy = orbitals.BSESinglets().eigenvalues()(_statehist.back().Index());
}

std::vector<int> Statefilter::OverlapFilter(const Orbitals& orbitals) const {
  std::vector<int> indexes;
  if (_statehist.size() <= 1) {
    indexes = std::vector<int>{_statehist[0].Index()};
    return indexes;
  }

  Eigen::VectorXd Overlap = CalculateOverlap(orbitals);
  int validelements = Overlap.size();
  for (int i = 0; i < Overlap.size(); i++) {
    if (Overlap(i) < _overlapthreshold) {
      validelements--;
    }
  }

  std::vector<int> index = std::vector<int>(Overlap.size());
  std::iota(index.begin(), index.end(), 0);
  std::stable_sort(index.begin(), index.end(), [&Overlap](int i1, int i2) {
    return Overlap[i1] > Overlap[i2];
  });

  int offset = 0;
  if (_statehist[0].Type() == QMStateType::DQPstate) {
    offset = orbitals.getGWAmin();
  }

  for (int i : index) {
    if (int(indexes.size()) == validelements) {
      break;
    }
    indexes.push_back(i + offset);
  }
  return indexes;
}

std::vector<int> Statefilter::OverlapFilterBSE(const Orbitals& orbitals) const {
  std::vector<int> indexes;
  if (_statehist.size() <= 1) {
    indexes = std::vector<int>{_statehist[0].Index()};
    return indexes;
  }

  Eigen::VectorXd Overlap = CalculateOverlapBSE(orbitals);
  int validelements = Overlap.size();
  for (int i = 0; i < Overlap.size(); i++) {
    if (Overlap(i) < _overlapthreshold) {
      validelements--;
    }
  }

  std::vector<int> index = std::vector<int>(Overlap.size());
  std::iota(index.begin(), index.end(), 0);
  std::stable_sort(index.begin(), index.end(), [&Overlap](int i1, int i2) {
    return Overlap[i1] > Overlap[i2];
  });

  int offset = 0;
  if (_statehist[0].Type() == QMStateType::DQPstate) {
    offset = orbitals.getGWAmin();
  }

  for (int i : index) {
    if (int(indexes.size()) == validelements) {
      break;
    }
    indexes.push_back(i + offset);
  }
  return indexes;
}

std::vector<int> Statefilter::DensityFilter(const Orbitals& orbitals) const { 
  std::vector<int> indexes;
  if (_statehist.size() <= 1) {
    indexes = std::vector<int>{_statehist[0].Index()};
    return indexes;
  }
  
  Eigen::VectorXd dnorm = CalculateDNorm(orbitals);
  
  int validelements = dnorm.size();
  for (int i = 0; i < dnorm.size(); i++) {
    if (dnorm(i) > _dmatthreshold) {
      validelements--;
    }
  }
  std::vector<int> index = std::vector<int>(dnorm.size());
  std::iota(index.begin(), index.end(), 0);
  std::stable_sort(index.begin(), index.end(), [&dnorm](int i1, int i2) {
    return dnorm[i1] < dnorm[i2];
  });
  int offset = 0;
  if (_statehist[0].Type() == QMStateType::DQPstate) {
    offset = orbitals.getGWAmin();
  }

  for (int i : index) {
    if (int(indexes.size()) == validelements) {
      break;
    }
    indexes.push_back(i + offset);
  }
  return indexes;
}



std::vector<int> Statefilter::WassersteinFilter(const Orbitals& orbitals) const { 
  std::vector<int> indexes;
  if (_statehist.size() <= 1) {
    indexes = std::vector<int>{_statehist[0].Index()};
    return indexes;
  }
  
  Eigen::VectorXd ddnorm = CalculateWassersteinNorm(orbitals);
  Eigen::VectorXd ones = Eigen::VectorXd::Ones(ddnorm.size()); 
  Eigen::VectorXd en_m = orbitals.BSESinglets().eigenvalues() - (_lastbseenergy * ones) ;

  Eigen::VectorXd dnorm = ddnorm.array();// + penalty_num;
   
  int validelements = dnorm.size();
//  for (int i = 0; i < dnorm.size(); i++) {
//    if (en_m(i) > 0.1*en_m.maxCoeff()) {
//      validelements--;
//    }
//  }
    
  std::cout << "\n Wasserstein distance on Excited State Density matrix \n" << std::endl;
  std::cout << ddnorm <<  std::endl;
     
  std::vector<int> index = std::vector<int>(dnorm.size());
  std::iota(index.begin(), index.end(), 0);
  std::stable_sort(index.begin(), index.end(), [&dnorm](int i1, int i2) {
    return dnorm[i1] < dnorm[i2];
  });
  int offset = 0;
  if (_statehist[0].Type() == QMStateType::DQPstate) {
    offset = orbitals.getGWAmin();
  }

  for (int i : index) {
    if (int(indexes.size()) == validelements) {
      break;
    }
    indexes.push_back(i + offset);
  }
  return indexes;
}



std::vector<int> Statefilter::MAE(const Orbitals& orbitals) const {
    std::vector<int> indexes;
  if (_statehist.size() <= 1) {
    indexes = std::vector<int>{_statehist[0].Index()};
    return indexes;
  }
    
  Eigen::VectorXd dnorm = CalculateDNorm(orbitals);
  
  Eigen::VectorXd overlap = CalculateOverlapBSE(orbitals);
    Eigen::VectorXd ones = Eigen::VectorXd::Ones(dnorm.size());
  
  Eigen::VectorXd en_m = orbitals.BSESinglets().eigenvalues() - (_lastbseenergy * ones) ;
  Eigen::VectorXd en_p = orbitals.BSESinglets().eigenvalues() + (_lastbseenergy * ones) ;    
  
  Eigen::VectorXd en = en_m.array()/en_p.array();
  
  Eigen::VectorXd d = (dnorm + (-overlap + ones) + en.cwiseAbs())/3.;
    
  Eigen::MatrixXd result = Eigen::MatrixXd::Zero(dnorm.size(),3);
  result.col(0) = en;
  result.col(1) = dnorm;
  result.col(2) = overlap;
  
  std::cout << "\n !s3 " << orbitals.BSESinglets().eigenvalues()(2) << std::endl;
  std::cout << "\n !s4 " << orbitals.BSESinglets().eigenvalues()(3) << std::endl;
  
  std::cout << "\n !f3 " << orbitals.Oscillatorstrengths()(2) << std::endl;
  std::cout << "\n !f4 " << orbitals.Oscillatorstrengths()(3) << std::endl;
  
  std::cout << "\n Energy Density Overlap \n" << std::endl;
  std::cout << result << std::endl;
  
  std::cout << "\n !den3 " << en(2) << std::endl;
  std::cout << "\n !den4 " << en(3) << std::endl;
  
  std::cout << "\n !ddens3 " << dnorm(2) << std::endl;
  std::cout << "\n !ddens4 " << dnorm(3) << std::endl;
  
  std::cout << "\n !overlap3 " << overlap(2) << std::endl;
  std::cout << "\n !overlap4 " << overlap(3) << std::endl;
  
  int validelements = d.size();
  
  std::vector<int> index = std::vector<int>(d.size());
  std::iota(index.begin(), index.end(), 0);
  std::stable_sort(index.begin(), index.end(), [&d](int i1, int i2) {
    return d[i1] < d[i2];
  });
  int offset = 0;
  if (_statehist[0].Type() == QMStateType::DQPstate) {
    offset = orbitals.getGWAmin();
  }

  for (int i : index) {
    if (int(indexes.size()) == validelements) {
      break;
    }
    indexes.push_back(i + offset);
  }
  return indexes;
} 


std::vector<int> Statefilter::MSE(const Orbitals& orbitals) const {
    std::vector<int> indexes;
  if (_statehist.size() <= 1) {
    indexes = std::vector<int>{_statehist[0].Index()};
    return indexes;
  }
    
  Eigen::VectorXd dnorm = CalculateDNorm(orbitals);
  Eigen::VectorXd overlap = CalculateOverlapBSE(orbitals);
  Eigen::VectorXd ones = Eigen::VectorXd::Ones(dnorm.size());
  Eigen::VectorXd en = orbitals.BSESinglets().eigenvalues() - _lastbseenergy * ones ;

  Eigen::VectorXd d = (dnorm.cwiseAbs2() + (-overlap + ones).cwiseAbs2() + en.cwiseAbs2())/3.;
  
  
  Eigen::MatrixXd result = Eigen::MatrixXd::Zero(dnorm.size(),3);
  result.col(0) = en;
  result.col(1) = dnorm;
  result.col(2) = overlap;
  
  
  std::cout << "\n !s3 " << orbitals.BSESinglets().eigenvalues()(2) << std::endl;
  std::cout << "\n !s4 " << orbitals.BSESinglets().eigenvalues()(3) << std::endl;
  
  std::cout << "\n !f3 " << orbitals.Oscillatorstrengths()(2) << std::endl;
  std::cout << "\n !f4 " << orbitals.Oscillatorstrengths()(3) << std::endl;
  
  std::cout << "\n Energy Density Overlap \n" << std::endl;
  std::cout << result << std::endl;
  
  std::cout << "\n !den3 " << en(2) << std::endl;
  std::cout << "\n !den4 " << en(3) << std::endl;
  
  std::cout << "\n !ddens3 " << dnorm(2) << std::endl;
  std::cout << "\n !ddens4 " << dnorm(3) << std::endl;
  
  std::cout << "\n !overlap3 " << overlap(2) << std::endl;
  std::cout << "\n !overlap4 " << overlap(3) << std::endl;
  
  
  int validelements = d.size();
 
  std::vector<int> index = std::vector<int>(d.size());
  std::iota(index.begin(), index.end(), 0);
  std::stable_sort(index.begin(), index.end(), [&d](int i1, int i2) {
    return d[i1] < d[i2];
  });
  int offset = 0;
  if (_statehist[0].Type() == QMStateType::DQPstate) {
    offset = orbitals.getGWAmin();
  }

  for (int i : index) {
    if (int(indexes.size()) == validelements) {
      break;
    }
    indexes.push_back(i + offset);
  }
  return indexes;
} 


void Statefilter::WriteToCpt(CheckpointWriter& w) const {
  std::vector<std::string> statehiststring;
  statehiststring.reserve(_statehist.size());
  for (const QMState& s : _statehist) {
    statehiststring.push_back(s.ToString());
  }
  w(statehiststring, "statehist");
  w(_use_oscfilter, "oscfilter");
  w(_oscthreshold, "oscthreshold");

  w(_use_overlapfilter, "overlapfilter");
  w(_use_overlapfilter_bse, "overlapfilterbse");
  w(_use_maefilter, "maefilter");
  w(_use_msefilter, "msefilter");
  w(_use_wasserstein, "wasserstein");
  w(_overlapthreshold, "overlapthreshold");
  w(_lastbse_R, "lastbse_R");
  w(_lastbse_AR, "lastbse_AR");
  w(_lastbseenergy, "lastbseenergies");
  
  
  w(_use_densityfilter, "densityfilter");
  w(_dmatthreshold, "dmatthreshold");
  w(_lastdmat, "lastdmat");

  w(_use_localisationfilter, "localisationfilter");
  w(_loc_threshold, "locthreshold");
  w(int(_fragment_loc.size()), "loc_fragments");
  for (unsigned i = 0; i < _fragment_loc.size(); i++) {
    CheckpointWriter ww = w.openChild("fragment_loc_" + std::to_string(i));
    _fragment_loc[i].WriteToCpt(ww);
  }

  w(_use_dQfilter, "dQfilter");
  w(_dQ_threshold, "dQthreshold");
  w(int(_fragment_dQ.size()), "dQ_fragments");
  for (unsigned i = 0; i < _fragment_dQ.size(); i++) {
    CheckpointWriter ww = w.openChild("fragment_dQ_" + std::to_string(i));
    _fragment_dQ[i].WriteToCpt(ww);
  }
}

void Statefilter::ReadFromCpt(CheckpointReader& r) {
  std::vector<std::string> statehiststring;
  r(statehiststring, "statehist");
  _statehist.clear();
  _statehist.reserve(statehiststring.size());
  for (const std::string& s : statehiststring) {
    _statehist.push_back(QMState(s));
  }
  r(_use_oscfilter, "oscfilter");
  r(_oscthreshold, "oscthreshold");
 
//  r(_use_variationalfilter, "variationalfilter");
//  r(_variationalthreshold, "variationalthreshold");
  r(_use_overlapfilter, "overlapfilter");
  r(_overlapthreshold, "overlapthreshold");
  r(_laststatecoeff, "laststatecoeff");
  r(_use_overlapfilter_bse, "overlapfilterbse");
  r(_use_maefilter, "maefilter");
  r(_use_msefilter, "msefilter");
  r(_use_wasserstein, "wasserstein");
  r(_lastbse_R, "lastbse_R");
  r(_lastbse_AR, "lastbse_AR");
  
  r(_use_densityfilter, "densityfilter");
  r(_dmatthreshold, "dmatthreshold");
  r(_lastdmat, "lastdmat");
  r(_lastbseenergy, "lastbseenergies");
  
  r(_use_localisationfilter, "localisationfilter");
  r(_loc_threshold, "locthreshold");
  _fragment_loc.clear();
  int loc_size = 0;
  r(loc_size, "loc_fragments");
  _fragment_loc.resize(loc_size);
  for (unsigned i = 0; i < _fragment_loc.size(); i++) {
    CheckpointReader rr = r.openChild("fragment_loc_" + std::to_string(i));
    _fragment_loc[i].ReadFromCpt(rr);
  }

  r(_use_dQfilter, "dQfilter");
  r(_dQ_threshold, "dQthreshold");
  _fragment_dQ.clear();
  int dQ_size = 0;
  r(dQ_size, "dQ_fragments");
  _fragment_dQ.resize(dQ_size);
  for (unsigned i = 0; i < _fragment_dQ.size(); i++) {
    CheckpointReader rr = r.openChild("fragment_dQ_" + std::to_string(i));
    _fragment_dQ[i].ReadFromCpt(rr);
  }
}

}  // namespace xtp
}  // namespace votca
