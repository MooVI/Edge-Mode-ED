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
	${OBJECTDIR}/XYZ_dists.o \
	${OBJECTDIR}/XYZ_mean.o \
	${OBJECTDIR}/XYZ_stagger.o \
	${OBJECTDIR}/decay_XYZ.o \
	${OBJECTDIR}/decay_ising.o \
	${OBJECTDIR}/exact_diag.o \
	${OBJECTDIR}/exact_diag_basis.o \
	${OBJECTDIR}/exact_diag_basis_J2.o \
	${OBJECTDIR}/exact_diag_basis_f.o \
	${OBJECTDIR}/exact_diag_basis_xyz.o \
	${OBJECTDIR}/exact_diag_mean.o \
	${OBJECTDIR}/ising_dists.o \
	${OBJECTDIR}/ising_dists_bulk.o \
	${OBJECTDIR}/ising_extra_terms_dists.o \
	${OBJECTDIR}/ising_kron_dists.o \
	${OBJECTDIR}/ising_kron_general.o \
	${OBJECTDIR}/kit_hubbard.o \
	${OBJECTDIR}/kronecker_product.o \
	${OBJECTDIR}/length_XYZ.o \
	${OBJECTDIR}/length_ising.o \
	${OBJECTDIR}/parafermion.o \
	${OBJECTDIR}/parafermion_basis.o \
	${OBJECTDIR}/parafermion_dists.o \
	${OBJECTDIR}/parafermion_mean.o \
	${OBJECTDIR}/parafermion_overlap.o \
	${OBJECTDIR}/timedecay_ising.o


# C Compiler Flags
CFLAGS=

# CC Compiler Flags
CCFLAGS=-m64 -lmpfr -lgmp -march=native -Ofast -pthread -no-pie
CXXFLAGS=-m64 -lmpfr -lgmp -march=native -Ofast -pthread -no-pie

# Fortran Compiler Flags
FFLAGS=

# Assembler Flags
ASFLAGS=

# Link Libraries and Options
LDLIBSOPTIONS=-L$$HOME/lib ../../Plotting/Plotter/dist/Release/GNU-Linux-x86/libplotter.a ../../NumericalMethods/NumericalMethods/dist/Release/GNU-Linux-x86/libnumericalmethods.a

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS}
	"${MAKE}"  -f nbproject/Makefile-${CND_CONF}.mk ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/tp11_1

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/tp11_1: ../../Plotting/Plotter/dist/Release/GNU-Linux-x86/libplotter.a

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/tp11_1: ../../NumericalMethods/NumericalMethods/dist/Release/GNU-Linux-x86/libnumericalmethods.a

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/tp11_1: ${OBJECTFILES}
	${MKDIR} -p ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}
	${LINK.cc} -o ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/tp11_1 ${OBJECTFILES} ${LDLIBSOPTIONS} -lgmp -lmpfr -lyaml-cpp

${OBJECTDIR}/XYZ_dists.o: XYZ_dists.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -I../.. -I.. -I$$HOME/include -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/XYZ_dists.o XYZ_dists.cpp

${OBJECTDIR}/XYZ_mean.o: XYZ_mean.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -I../.. -I.. -I$$HOME/include -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/XYZ_mean.o XYZ_mean.cpp

${OBJECTDIR}/XYZ_stagger.o: XYZ_stagger.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -I../.. -I.. -I$$HOME/include -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/XYZ_stagger.o XYZ_stagger.cpp

${OBJECTDIR}/decay_XYZ.o: decay_XYZ.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -I../.. -I.. -I$$HOME/include -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/decay_XYZ.o decay_XYZ.cpp

${OBJECTDIR}/decay_ising.o: decay_ising.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -I../.. -I.. -I$$HOME/include -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/decay_ising.o decay_ising.cpp

${OBJECTDIR}/exact_diag.o: exact_diag.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -I../.. -I.. -I$$HOME/include -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/exact_diag.o exact_diag.cpp

${OBJECTDIR}/exact_diag_basis.o: exact_diag_basis.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -I../.. -I.. -I$$HOME/include -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/exact_diag_basis.o exact_diag_basis.cpp

${OBJECTDIR}/exact_diag_basis_J2.o: exact_diag_basis_J2.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -I../.. -I.. -I$$HOME/include -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/exact_diag_basis_J2.o exact_diag_basis_J2.cpp

${OBJECTDIR}/exact_diag_basis_f.o: exact_diag_basis_f.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -I../.. -I.. -I$$HOME/include -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/exact_diag_basis_f.o exact_diag_basis_f.cpp

${OBJECTDIR}/exact_diag_basis_xyz.o: exact_diag_basis_xyz.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -I../.. -I.. -I$$HOME/include -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/exact_diag_basis_xyz.o exact_diag_basis_xyz.cpp

${OBJECTDIR}/exact_diag_mean.o: exact_diag_mean.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -I../.. -I.. -I$$HOME/include -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/exact_diag_mean.o exact_diag_mean.cpp

${OBJECTDIR}/ising_dists.o: ising_dists.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -I../.. -I.. -I$$HOME/include -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/ising_dists.o ising_dists.cpp

${OBJECTDIR}/ising_dists_bulk.o: ising_dists_bulk.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -I../.. -I.. -I$$HOME/include -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/ising_dists_bulk.o ising_dists_bulk.cpp

${OBJECTDIR}/ising_extra_terms_dists.o: ising_extra_terms_dists.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -I../.. -I.. -I$$HOME/include -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/ising_extra_terms_dists.o ising_extra_terms_dists.cpp

${OBJECTDIR}/ising_kron_dists.o: ising_kron_dists.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -I../.. -I.. -I$$HOME/include -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/ising_kron_dists.o ising_kron_dists.cpp

${OBJECTDIR}/ising_kron_general.o: ising_kron_general.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -I../.. -I.. -I$$HOME/include -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/ising_kron_general.o ising_kron_general.cpp

${OBJECTDIR}/kit_hubbard.o: kit_hubbard.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -I../.. -I.. -I$$HOME/include -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/kit_hubbard.o kit_hubbard.cpp

${OBJECTDIR}/kronecker_product.o: kronecker_product.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -I../.. -I.. -I$$HOME/include -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/kronecker_product.o kronecker_product.cpp

${OBJECTDIR}/length_XYZ.o: length_XYZ.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -I../.. -I.. -I$$HOME/include -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/length_XYZ.o length_XYZ.cpp

${OBJECTDIR}/length_ising.o: length_ising.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -I../.. -I.. -I$$HOME/include -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/length_ising.o length_ising.cpp

${OBJECTDIR}/parafermion.o: parafermion.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -I../.. -I.. -I$$HOME/include -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/parafermion.o parafermion.cpp

${OBJECTDIR}/parafermion_basis.o: parafermion_basis.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -I../.. -I.. -I$$HOME/include -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/parafermion_basis.o parafermion_basis.cpp

${OBJECTDIR}/parafermion_dists.o: parafermion_dists.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -I../.. -I.. -I$$HOME/include -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/parafermion_dists.o parafermion_dists.cpp

${OBJECTDIR}/parafermion_mean.o: parafermion_mean.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -I../.. -I.. -I$$HOME/include -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/parafermion_mean.o parafermion_mean.cpp

${OBJECTDIR}/parafermion_overlap.o: parafermion_overlap.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -I../.. -I.. -I$$HOME/include -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/parafermion_overlap.o parafermion_overlap.cpp

${OBJECTDIR}/timedecay_ising.o: timedecay_ising.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -I../.. -I.. -I$$HOME/include -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/timedecay_ising.o timedecay_ising.cpp

# Subprojects
.build-subprojects:

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r ${CND_BUILDDIR}/${CND_CONF}
	${RM} ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/tp11_1

# Subprojects
.clean-subprojects:

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
