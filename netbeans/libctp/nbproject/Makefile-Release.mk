#
# Generated Makefile - do not edit!
#
# Edit the Makefile in the project folder instead (../Makefile). Each target
# has a -pre and a -post target defined where you can add customized code.
#
# This makefile implements configuration specific macros and targets.


# Environment
MKDIR=mkdir
CP=cp
GREP=grep
NM=nm
CCADMIN=CCadmin
RANLIB=ranlib
CC=gcc
CCC=g++
CXX=g++
FC=gfortran
AS=as

# Macros
CND_PLATFORM=GNU-Linux-x86
CND_DLIB_EXT=so
CND_CONF=Release
CND_DISTDIR=dist
CND_BUILDDIR=build

# Include project Makefile
include Makefile

# Object Directory
OBJECTDIR=${CND_BUILDDIR}/${CND_CONF}/${CND_PLATFORM}

# Object Files
OBJECTFILES= \
	${OBJECTDIR}/_ext/484457893/apolarsite.o \
	${OBJECTDIR}/_ext/484457893/calculatorfactory.o \
	${OBJECTDIR}/_ext/238600121/eoutersphere.o \
	${OBJECTDIR}/_ext/238600121/idft.o \
	${OBJECTDIR}/_ext/484457893/fragment.o \
	${OBJECTDIR}/_ext/484457893/job.o \
	${OBJECTDIR}/_ext/484457893/molecule.o \
	${OBJECTDIR}/_ext/484457893/orbitals.o \
	${OBJECTDIR}/_ext/484457893/parallelpaircalc.o \
	${OBJECTDIR}/_ext/484457893/parallelxjobcalc.o \
	${OBJECTDIR}/_ext/484457893/polarseg.o \
	${OBJECTDIR}/_ext/484457893/polarsite.o \
	${OBJECTDIR}/_ext/484457893/polartop.o \
	${OBJECTDIR}/_ext/484457893/progressobserver.o \
	${OBJECTDIR}/_ext/484457893/qmapplication.o \
	${OBJECTDIR}/_ext/484457893/qmcalculator.o \
	${OBJECTDIR}/_ext/484457893/qmdatabase.o \
	${OBJECTDIR}/_ext/484457893/qmmachine.o \
	${OBJECTDIR}/_ext/484457893/qmnblist.o \
	${OBJECTDIR}/_ext/484457893/qmpackagefactory.o \
	${OBJECTDIR}/_ext/648834637/gaussian.o \
	${OBJECTDIR}/_ext/648834637/nwchem.o \
	${OBJECTDIR}/_ext/648834637/turbomole.o \
	${OBJECTDIR}/_ext/484457893/qmpair.o \
	${OBJECTDIR}/_ext/484457893/qmtool.o \
	${OBJECTDIR}/_ext/484457893/segment.o \
	${OBJECTDIR}/_ext/484457893/segmenttype.o \
	${OBJECTDIR}/_ext/484457893/statesaversqlite.o \
	${OBJECTDIR}/_ext/484457893/toolfactory.o \
	${OBJECTDIR}/_ext/484457893/topology.o \
	${OBJECTDIR}/_ext/484457893/version.o \
	${OBJECTDIR}/_ext/484457893/version_nb.o \
	${OBJECTDIR}/_ext/484457893/xinductor.o \
	${OBJECTDIR}/_ext/484457893/xinteractor.o \
	${OBJECTDIR}/_ext/484457893/xjob.o \
	${OBJECTDIR}/_ext/484457893/xmapper.o


# C Compiler Flags
CFLAGS=

# CC Compiler Flags
CCFLAGS=
CXXFLAGS=

# Fortran Compiler Flags
FFLAGS=

# Assembler Flags
ASFLAGS=

# Link Libraries and Options
LDLIBSOPTIONS=

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS}
	"${MAKE}"  -f nbproject/Makefile-${CND_CONF}.mk ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/liblibctp.a

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/liblibctp.a: ${OBJECTFILES}
	${MKDIR} -p ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}
	${RM} ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/liblibctp.a
	${AR} -rv ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/liblibctp.a ${OBJECTFILES} 
	$(RANLIB) ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/liblibctp.a

${OBJECTDIR}/_ext/484457893/apolarsite.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libctp/apolarsite.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/484457893
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/484457893/apolarsite.o ../../src/libctp/apolarsite.cc

${OBJECTDIR}/_ext/484457893/calculatorfactory.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libctp/calculatorfactory.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/484457893
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/484457893/calculatorfactory.o ../../src/libctp/calculatorfactory.cc

${OBJECTDIR}/_ext/238600121/eoutersphere.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libctp/calculators/eoutersphere.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/238600121
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/238600121/eoutersphere.o ../../src/libctp/calculators/eoutersphere.cc

${OBJECTDIR}/_ext/238600121/idft.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libctp/calculators/idft.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/238600121
	${RM} $@.d
	$(COMPILE.c) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/238600121/idft.o ../../src/libctp/calculators/idft.cc

${OBJECTDIR}/_ext/484457893/fragment.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libctp/fragment.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/484457893
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/484457893/fragment.o ../../src/libctp/fragment.cc

${OBJECTDIR}/_ext/484457893/job.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libctp/job.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/484457893
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/484457893/job.o ../../src/libctp/job.cc

${OBJECTDIR}/_ext/484457893/molecule.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libctp/molecule.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/484457893
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/484457893/molecule.o ../../src/libctp/molecule.cc

${OBJECTDIR}/_ext/484457893/orbitals.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libctp/orbitals.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/484457893
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/484457893/orbitals.o ../../src/libctp/orbitals.cc

${OBJECTDIR}/_ext/484457893/parallelpaircalc.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libctp/parallelpaircalc.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/484457893
	${RM} $@.d
	$(COMPILE.c) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/484457893/parallelpaircalc.o ../../src/libctp/parallelpaircalc.cc

${OBJECTDIR}/_ext/484457893/parallelxjobcalc.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libctp/parallelxjobcalc.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/484457893
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/484457893/parallelxjobcalc.o ../../src/libctp/parallelxjobcalc.cc

${OBJECTDIR}/_ext/484457893/polarseg.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libctp/polarseg.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/484457893
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/484457893/polarseg.o ../../src/libctp/polarseg.cc

${OBJECTDIR}/_ext/484457893/polarsite.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libctp/polarsite.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/484457893
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/484457893/polarsite.o ../../src/libctp/polarsite.cc

${OBJECTDIR}/_ext/484457893/polartop.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libctp/polartop.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/484457893
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/484457893/polartop.o ../../src/libctp/polartop.cc

${OBJECTDIR}/_ext/484457893/progressobserver.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libctp/progressobserver.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/484457893
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/484457893/progressobserver.o ../../src/libctp/progressobserver.cc

${OBJECTDIR}/_ext/484457893/qmapplication.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libctp/qmapplication.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/484457893
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/484457893/qmapplication.o ../../src/libctp/qmapplication.cc

${OBJECTDIR}/_ext/484457893/qmcalculator.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libctp/qmcalculator.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/484457893
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/484457893/qmcalculator.o ../../src/libctp/qmcalculator.cc

${OBJECTDIR}/_ext/484457893/qmdatabase.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libctp/qmdatabase.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/484457893
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/484457893/qmdatabase.o ../../src/libctp/qmdatabase.cc

${OBJECTDIR}/_ext/484457893/qmmachine.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libctp/qmmachine.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/484457893
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/484457893/qmmachine.o ../../src/libctp/qmmachine.cc

${OBJECTDIR}/_ext/484457893/qmnblist.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libctp/qmnblist.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/484457893
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/484457893/qmnblist.o ../../src/libctp/qmnblist.cc

${OBJECTDIR}/_ext/484457893/qmpackagefactory.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libctp/qmpackagefactory.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/484457893
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/484457893/qmpackagefactory.o ../../src/libctp/qmpackagefactory.cc

${OBJECTDIR}/_ext/648834637/gaussian.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libctp/qmpackages/gaussian.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/648834637
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/648834637/gaussian.o ../../src/libctp/qmpackages/gaussian.cc

${OBJECTDIR}/_ext/648834637/nwchem.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libctp/qmpackages/nwchem.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/648834637
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/648834637/nwchem.o ../../src/libctp/qmpackages/nwchem.cc

${OBJECTDIR}/_ext/648834637/turbomole.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libctp/qmpackages/turbomole.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/648834637
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/648834637/turbomole.o ../../src/libctp/qmpackages/turbomole.cc

${OBJECTDIR}/_ext/484457893/qmpair.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libctp/qmpair.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/484457893
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/484457893/qmpair.o ../../src/libctp/qmpair.cc

${OBJECTDIR}/_ext/484457893/qmtool.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libctp/qmtool.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/484457893
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/484457893/qmtool.o ../../src/libctp/qmtool.cc

${OBJECTDIR}/_ext/484457893/segment.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libctp/segment.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/484457893
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/484457893/segment.o ../../src/libctp/segment.cc

${OBJECTDIR}/_ext/484457893/segmenttype.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libctp/segmenttype.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/484457893
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/484457893/segmenttype.o ../../src/libctp/segmenttype.cc

${OBJECTDIR}/_ext/484457893/statesaversqlite.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libctp/statesaversqlite.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/484457893
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/484457893/statesaversqlite.o ../../src/libctp/statesaversqlite.cc

${OBJECTDIR}/_ext/484457893/toolfactory.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libctp/toolfactory.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/484457893
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/484457893/toolfactory.o ../../src/libctp/toolfactory.cc

${OBJECTDIR}/_ext/484457893/topology.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libctp/topology.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/484457893
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/484457893/topology.o ../../src/libctp/topology.cc

${OBJECTDIR}/_ext/484457893/version.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libctp/version.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/484457893
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/484457893/version.o ../../src/libctp/version.cc

${OBJECTDIR}/_ext/484457893/version_nb.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libctp/version_nb.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/484457893
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/484457893/version_nb.o ../../src/libctp/version_nb.cc

${OBJECTDIR}/_ext/484457893/xinductor.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libctp/xinductor.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/484457893
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/484457893/xinductor.o ../../src/libctp/xinductor.cc

${OBJECTDIR}/_ext/484457893/xinteractor.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libctp/xinteractor.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/484457893
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/484457893/xinteractor.o ../../src/libctp/xinteractor.cc

${OBJECTDIR}/_ext/484457893/xjob.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libctp/xjob.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/484457893
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/484457893/xjob.o ../../src/libctp/xjob.cc

${OBJECTDIR}/_ext/484457893/xmapper.o: nbproject/Makefile-${CND_CONF}.mk ../../src/libctp/xmapper.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/484457893
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/484457893/xmapper.o ../../src/libctp/xmapper.cc

# Subprojects
.build-subprojects:

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r ${CND_BUILDDIR}/${CND_CONF}
	${RM} ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/liblibctp.a

# Subprojects
.clean-subprojects:

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc