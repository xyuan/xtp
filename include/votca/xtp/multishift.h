///*
// *            Copyright 2009-2019 The VOTCA Development Team
// *                       (http://www.votca.org)
// *
// *      Licensed under the Apache License, Version 2.0 (the "License")
// *
// * You may not use this file except in compliance with the License.
// * You may obtain a copy of the License at
// *
// *              http://www.apache.org/licenses/LICENSE-2.0
// *
// * Unless required by applicable law or agreed to in writing, software
// * distributed under the License is distributed on an "AS IS" BASIS,
// * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// * See the License for the specific language governing permissions and
// * limitations under the License.
// */

#pragma once
#ifndef MULTISHIFT_H
#define MULTISHIFT_H

#include <votca/tools/property.h>
#include <fstream>
#include <votca/xtp/eigen.h>
#include <votca/xtp/logger.h>

namespace votca{
namespace xtp{
class Multishift{
public:
    Multishift(){};
    
    void setBasisSize(double basis_size);
    
    Eigen::VectorXcd CBiCG(Eigen::MatrixXcd A, Eigen::VectorXcd b);
    
    Eigen::VectorXcd DoMultishift(Eigen::MatrixXcd A, Eigen::VectorXcd b, std::complex<double> w);

    void testMultishift();
    
private:

    int _basis_size;
    
    std::vector<std::complex<double>> _a;
    
    std::vector<std::complex<double>> _b;
    
    std::vector<Eigen::VectorXcd> _r;
    
};
}
}
#endif