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
CND_CONF=Debug
CND_DISTDIR=dist
CND_BUILDDIR=build

# Include project Makefile
include Makefile

# Object Directory
OBJECTDIR=${CND_BUILDDIR}/${CND_CONF}/${CND_PLATFORM}

# Object Files
OBJECTFILES= \
	${OBJECTDIR}/exact_diag.o \
	${OBJECTDIR}/exact_diag_basis.o \
	${OBJECTDIR}/exact_diag_basis_J2.o \
	${OBJECTDIR}/exact_diag_basis_f.o \
	${OBJECTDIR}/exact_diag_basis_xyz.o \
	${OBJECTDIR}/parafermion.o \
	${OBJECTDIR}/parafermion_basis.o \
	${OBJECTDIR}/parafermion_overlap.o


# C Compiler Flags
CFLAGS=

# CC Compiler Flags
CCFLAGS=-lgmp -lmpfr -msse4
CXXFLAGS=-lgmp -lmpfr -msse4

# Fortran Compiler Flags
FFLAGS=

# Assembler Flags
ASFLAGS=

# Link Libraries and Options
LDLIBSOPTIONS=../../Plotting/Plotter/dist/Release/GNU-Linux-x86/libplotter.a ../../NumericalMethods/NumericalMethods/dist/Release/GNU-Linux-x86/libnumericalmethods.a

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS}
	"${MAKE}"  -f nbproject/Makefile-${CND_CONF}.mk ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/tp11_1

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/tp11_1: ../../Plotting/Plotter/dist/Release/GNU-Linux-x86/libplotter.a

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/tp11_1: ../../NumericalMethods/NumericalMethods/dist/Release/GNU-Linux-x86/libnumericalmethods.a

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/tp11_1: ${OBJECTFILES}
	${MKDIR} -p ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}
	${LINK.cc} -o ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/tp11_1 ${OBJECTFILES} ${LDLIBSOPTIONS}

${OBJECTDIR}/exact_diag.o: exact_diag.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../.. -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/exact_diag.o exact_diag.cpp

${OBJECTDIR}/exact_diag_basis.o: exact_diag_basis.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../.. -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/exact_diag_basis.o exact_diag_basis.cpp

${OBJECTDIR}/exact_diag_basis_J2.o: exact_diag_basis_J2.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../.. -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/exact_diag_basis_J2.o exact_diag_basis_J2.cpp

${OBJECTDIR}/exact_diag_basis_f.o: exact_diag_basis_f.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../.. -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/exact_diag_basis_f.o exact_diag_basis_f.cpp

${OBJECTDIR}/exact_diag_basis_xyz.o: exact_diag_basis_xyz.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../.. -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/exact_diag_basis_xyz.o exact_diag_basis_xyz.cpp

${OBJECTDIR}/parafermion.o: parafermion.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../.. -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/parafermion.o parafermion.cpp

${OBJECTDIR}/parafermion_basis.o: parafermion_basis.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../.. -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/parafermion_basis.o parafermion_basis.cpp

${OBJECTDIR}/parafermion_overlap.o: parafermion_overlap.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../.. -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/parafermion_overlap.o parafermion_overlap.cpp

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
