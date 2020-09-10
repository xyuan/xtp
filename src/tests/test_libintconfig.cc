/*
 * Copyright 2009-2020 The VOTCA Development Team (http://www.votca.org)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */
#define BOOST_TEST_MAIN

#define BOOST_TEST_MODULE libint_config

// override for default location of basissets for Libint
#define DATADIR std::string(XTP_TEST_DATA_FOLDER) + "/libintconfig"

#include <iostream>
#include <fstream>

// Third party includes
#include <boost/test/unit_test.hpp>
#include <libint2.hpp>


// VOTCA includes
#include <votca/tools/eigenio_matrixmarket.h>
#include <votca/tools/filesystem.h>


using namespace votca;

BOOST_AUTO_TEST_SUITE(libint_config)

BOOST_AUTO_TEST_CASE(overlap_integrals) {

libint2::initialize();

std::string xyzfilename = std::string(XTP_TEST_DATA_FOLDER) + "/libintconfig/hydrogen.xyz";
std::ifstream input_file(xyzfilename);
std::vector<libint2::Atom> atoms = libint2::read_dotxyz(input_file);

libint2::BasisSet obs("def2-tzvp", atoms);


std::cout << std::string(XTP_TEST_DATA_FOLDER) << std::endl;
std::copy(std::begin(obs), std::end(obs),
          std::ostream_iterator<libint2::Shell>(std::cout, "\n"));


libint2::finalize();

}


BOOST_AUTO_TEST_SUITE_END()
