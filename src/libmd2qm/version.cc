/* 
 * Copyright 2009 The VOTCA Development Team (http://www.votca.org)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
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

#include <votca/tools/version.h>
#include <iostream>
#include "version.h"
#include "votca_config.h"

namespace votca { namespace md2qm {

#ifdef HGVERSION
  static const std::string version_str = VERSION " " HGVERSION " (compiled " __DATE__ ", " __TIME__ ")";
#else
  static const std::string version_str = VERSION "(compiled " __DATE__ ", " __TIME__ ")";
#endif

const std::string &MD2QMVersionStr()
{
    return version_str;
}

void HelpTextHeader(const std::string &tool_name)
{
    std::cout
         << "==================================================\n"
         << "========   VOTCA (http://www.votca.org)   ========\n"
         << "==================================================\n\n"
	 << "please submit bugs to " PACKAGE_BUGREPORT "\n\n" 
	 << tool_name << ", version " << votca::md2qm::MD2QMVersionStr() 
         << "\nvotca_tools, version " << votca::tools::ToolsVersionStr() 
         << "\n\n";
}

}}

