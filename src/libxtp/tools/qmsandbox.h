/*
 *            Copyright 2009-2020 The VOTCA Development Team
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

#pragma once
#ifndef _VOTCA_XTP_QMSANDBOX_H
#define _VOTCA_XTP_QMSANDBOX_H

#include <votca/xtp/logger.h>
#include <votca/xtp/qmtool.h>

#include <votca/xtp/aobasis.h>
#include <votca/xtp/aomatrix.h>
#include <votca/xtp/qmpackagefactory.h>

namespace votca {
namespace xtp {

class QMSandbox : public QMTool {
 public:
  QMSandbox() = default;
  ~QMSandbox() override = default;

  std::string Identify() override { return "qmsandbox"; }

  void Initialize(tools::Property& user_options) override;
  bool Evaluate() override;

 private:
  std::string _orbfile;
  std::string _espfile;
  std::string _mpsfiled;
  std::string _mpsfileds;
  std::string _mpsfileq;
  std::string _mpsfileqs;
};

void QMSandbox::Initialize(tools::Property& user_options) {

  // get pre-defined default options from VOTCASHARE/xtp/xml/qmsandbox.xml
  LoadDefaults("xtp");
  // update options with user specified input
  UpdateWithUserOptions(user_options);

  _orbfile =
      _options.ifExistsReturnElseThrowRuntimeError<std::string>(".orbfile");
  _espfile =
      _options.ifExistsReturnElseThrowRuntimeError<std::string>(".espfile");
  _mpsfiled =
      _options.ifExistsReturnElseThrowRuntimeError<std::string>(".dipole");

  _mpsfileds = _options.ifExistsReturnElseThrowRuntimeError<std::string>(
      ".dipole_split");
  _mpsfileq =
      _options.ifExistsReturnElseThrowRuntimeError<std::string>(".quadrupole");
  _mpsfileqs = _options.ifExistsReturnElseThrowRuntimeError<std::string>(
      ".quadrupole_split");
}

bool QMSandbox::Evaluate() { return true; }

}  // namespace xtp
}  // namespace votca

#endif
