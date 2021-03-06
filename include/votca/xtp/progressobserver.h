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
/// For an earlier history see ctp repo commit
/// 77795ea591b29e664153f9404c8655ba28dc14e9

#pragma once
#ifndef VOTCA_XTP_PROGRESSOBSERVER_H
#define VOTCA_XTP_PROGRESSOBSERVER_H

// Standard includes
#include <vector>

// Third party includes
#include <boost/interprocess/sync/file_lock.hpp>
#include <boost/program_options.hpp>

// VOTCA includes
#include <votca/tools/mutex.h>
#include <votca/tools/property.h>

namespace votca {
namespace xtp {

class QMThread;

template <typename JobContainer>
class ProgObserver {

  using Job = typename std::iterator_traits<
      typename JobContainer::iterator>::value_type;

  using Result = typename Job::JobResult;

 public:
  void InitCmdLineOpts(const boost::program_options::variables_map &optsMap);
  void InitFromProgFile(std::string progFile, QMThread &thread);
  ProgObserver::Job *RequestNextJob(QMThread &thread);
  void ReportJobDone(Job &job, Result &res, QMThread &thread);

  void SyncWithProgFile(QMThread &thread);
  void LockProgFile(QMThread &thread);
  void ReleaseProgFile(QMThread &thread);

  std::string GenerateHost();
  std::string GenerateTime();

 private:
  std::string _lockFile = "";
  std::string _progFile = "";
  Index _cacheSize = -1;
  JobContainer _jobs;

  std::vector<Job *> _jobsToProc;
  std::vector<Job *> _jobsToSync;

  using iterator = typename JobContainer::iterator;
  iterator _metajit;
  using iterator_vec = typename std::vector<Job *>::iterator;
  iterator_vec _nextjit;
  tools::Mutex _lockThread;
  std::unique_ptr<boost::interprocess::file_lock> _flock;

  std::map<std::string, bool> _restart_hosts;
  std::map<std::string, bool> _restart_stats;
  bool _restartMode = false;
  Index _jobsReported = 0;

  bool _moreJobsAvailable = false;
  Index _startJobsCount = 0;
  Index _maxJobs = 0;
};

}  // namespace xtp
}  // namespace votca
#endif  // VOTCA_XTP_PROGRESSOBSERVER_H
