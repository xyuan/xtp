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
#include <stdio.h>
#include <votca/tools/constants.h>
#include <votca/tools/elements.h>
#include <votca/xtp/aobasis.h>
#include <votca/xtp/cubefile_writer.h>
#include <votca/xtp/orbitals.h>
#include <boost/format.hpp>

namespace votca {
namespace xtp {

void TransientAbsorption::Initialize(tools::Property& options) {
    std::string key = "options." + Identify();

    std::cout << options.ifExistsReturnElseThrowRuntimeError<std::string>(key + ".message") << std::endl;

}

bool TransientAbsorption::Evaluate(){
    return true;
}

} // namespace xtp
} // namespace votca