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
CND_PLATFORM=MinGW-Windows
CND_DLIB_EXT=dll
CND_CONF=Debug
CND_DISTDIR=dist
CND_BUILDDIR=build

# Include project Makefile
include Makefile

# Object Directory
OBJECTDIR=${CND_BUILDDIR}/${CND_CONF}/${CND_PLATFORM}

# Object Files
OBJECTFILES= \
	${OBJECTDIR}/ElementModule.o \
	${OBJECTDIR}/MethodModule.o \
	${OBJECTDIR}/main.o


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
	"${MAKE}"  -f nbproject/Makefile-${CND_CONF}.mk ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/fem_fortran.exe

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/fem_fortran.exe: ${OBJECTFILES}
	${MKDIR} -p ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}
	${LINK.f} -o ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/fem_fortran ${OBJECTFILES} ${LDLIBSOPTIONS}

${OBJECTDIR}/ElementModule.o: ElementModule.f90 
	${MKDIR} -p ${OBJECTDIR}
	$(COMPILE.f) -g -o ${OBJECTDIR}/ElementModule.o ElementModule.f90

${OBJECTDIR}/MethodModule.o: MethodModule.f90 
	${MKDIR} -p ${OBJECTDIR}
	$(COMPILE.f) -g -o ${OBJECTDIR}/MethodModule.o MethodModule.f90

${OBJECTDIR}/main.o: main.f90 
	${MKDIR} -p ${OBJECTDIR}
	$(COMPILE.f) -g -o ${OBJECTDIR}/main.o main.f90

# Subprojects
.build-subprojects:

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r ${CND_BUILDDIR}/${CND_CONF}
	${RM} ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/fem_fortran.exe
	${RM} *.mod

# Subprojects
.clean-subprojects:

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
