#@section LICENSE
#Copyright (c) 2013-2017, Regents of the University of California
#All rights reserved.
#Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
#1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
#2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
#THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


AC_PREREQ([2.60])

AC_INIT([seismio], [0.90], [YFCUI@SDSC.EDU])
AM_INIT_AUTOMAKE

AC_CONFIG_MACRO_DIR([m4])

# Checks for programs.
AC_PROG_FC
AC_PROG_FC_C_O
AC_FC_WRAPPERS
F77=$FC
dnl AC_PROG_F77
dnl AC_PROG_F77_C_O
#AC_LANG(Fortran)
#ACX_MPI([], [AC_MSG_ERROR([could not find mpi library for Fortran])])
#AC_PROG_RANLIB

# Default configure
#LT_INIT
# Turn off shared libraries during beta-testing, since they make the build process take too long.
LT_INIT([disable-shared])



AC_ARG_WITH(phdf5, AS_HELP_STRING([--with-phdf5], [prefix for installed HDF5 with --enable-parallel]),S_HDF5_HOME=$withval,S_HDF5_HOME=)


AC_ARG_WITH(pnetcdf, AS_HELP_STRING([--with-pnetcdf], [prefix for installed PNetCDF]),S_PNETCDF_HOME=$withval,S_PNETCDF_HOME=)


AC_ARG_WITH(adios, AS_HELP_STRING([--with-adios], [prefix for installed ADIOS]),S_ADIOS_HOME=$withval,S_ADIOS_HOME=)

AC_ARG_WITH(mxml, AS_HELP_STRING([--with-mxml], [prefix for installed MXML below 2.7 Adios 1.9.0 required]),S_MXML_HOME=$withval,S_MXML_HOME=)

#HDF5_CPPFLAGS="${PHDF5_CPPFLAGS}"
#HDF5_CFLAGS="${PHDF5_CFLAGS}"
#HDF5_LDFLAGS="${PHDF5_LDFLAGS}"
#HDF5_LIBS="${PHDF5_LIBS}"
#HDF5_LIB=/home/xyz/Tools/phdf5/lib/libhdf5hl_fortran.la
#S_HDF5_HOME=/home/xyz/Tools/phdf5
#S_PNETCDF_HOME=/home/xyz/Tools/PnetCDF
#S_ADIOS_HOME=/home/xyz/Tools/adios
#S_ADIOS_HOME=$(shell adios_config -c -f)

AM_CONDITIONAL([HAVE_HDF5], [test -n "$S_HDF5_HOME"])
AM_CONDITIONAL([HAVE_NETCDF], [test -n "$S_PNETCDF_HOME"])
AM_CONDITIONAL([HAVE_ADIOS], [test -n "$S_ADIOS_HOME"])
#if test -n "$S_HDF5_HOME"; then
#    AC_SUBST(S_HDF5_HOME)
#fi
AC_SUBST(S_HDF5_HOME)
AC_SUBST(S_PNETCDF_HOME)
AC_SUBST(S_ADIOS_HOME)
AC_SUBST(S_MXML_HOME)

# Checks for libraries.

# Checks for header files.

# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.

#AC_CONFIG_HEADERS([src/interfaces/include/seism_io.h
#                   src/interfaces/include/seism_io.fh
#                   src/interfaces/include/seism_io_consts.h
#                   src/interfaces/include/seism_io_consts.fh])

AC_CONFIG_FILES([src/interfaces/Makefile
                 Makefile])
AC_OUTPUT
