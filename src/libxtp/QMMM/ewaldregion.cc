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

#include <votca/xtp/ewaldregion.h>
#include <votca/xtp/polarregion.h>
#include <votca/xtp/qmregion.h>
#include <votca/xtp/staticregion.h>

namespace votca {
namespace xtp {

void EwaldRegion::Initialize(const tools::Property&) {}

void EwaldRegion::Evaluate(std::vector<std::unique_ptr<Region> >&) {}

double EwaldRegion::charge() const { return 0; }

void EwaldRegion::AppendResult(tools::Property&) const {}

void EwaldRegion::Reset() { return; }
double EwaldRegion::InteractwithQMRegion(const QMRegion&) { return 0.0; }
double EwaldRegion::InteractwithPolarRegion(const PolarRegion&) { return 0.0; }
double EwaldRegion::InteractwithStaticRegion(const StaticRegion&) {
  throw std::runtime_error("Ewald Region and static region cannot coexist.");
  return 0.0;
}

void EwaldRegion::WritePDB(csg::PDBWriter&) const {}

void EwaldRegion::WriteToCpt(CheckpointWriter& w) const {
  w(_id, "id");
  w(identify(), "type");
  w(_size, "size");
}

void EwaldRegion::ReadFromCpt(CheckpointReader& r) {
  r(_id, "id");
  r(_size, "size");
}

}  // namespace xtp
}  // namespace votca
