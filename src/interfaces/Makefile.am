#@section LICENSE
#Copyright (c) 2013-2017, Regents of the University of California
#All rights reserved.
#Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
#1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
#2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
#THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

AUTOMAKE_OPTIONS = subdir-objects
#CC=mpicc
#FC=mpif90

#RANLIB
#noinst_LIBRARIES=libseismio.a
#libtool
#noinst_LTLIBRARIES=libseismio.la

#for install lib file
lib_LTLIBRARIES=libseismio.la
#libseismio_a_SOURCES = seism_io.c seism_io_fortran.c \
#	seism_io_mpiio.c seism_io_phdf5.c seism_io_read_binary.c \
#	seism_io_pnetcdf.c seism_io_adios.c
libseismio_la_SOURCES = seism_io.c seism_io_fortran.c \
	seism_io_mpiio.c seism_io_phdf5.c seism_io_read_binary.c \
	seism_io_pnetcdf.c seism_io_adios.c \
../f-io-libs/f-mpiio-lib/initials.f90 \
../f-io-libs/f-mpiio-lib/set_filetype_style.f90 \
../f-io-libs/f-mpiio-lib/openfile.f90 \
../f-io-libs/f-mpiio-lib/readfile.f90 \
../f-io-libs/f-mpiio-lib/writefile.f90 \
../f-io-libs/f-mpiio-lib/readmesh.f90 \
../c-io-libs/c-mpiio-lib/c_mpiio.c

AM_CFLAGS = -Dsun -lmpich
AM_FCFLAGS =

if HAVE_HDF5
HDF5_SRC=../f-io-libs/f-phdf5-lib/setup.f90 \
../f-io-libs/f-phdf5-lib/find_start_count.f90 \
../f-io-libs/f-phdf5-lib/h_initials.f90 \
../f-io-libs/f-phdf5-lib/h_openfile.f90 \
../f-io-libs/f-phdf5-lib/h_readfile.f90 \
../f-io-libs/f-phdf5-lib/h_writefile.f90 \
../c-io-libs/c-phdf5-lib/c_phdf5.c
HDF5_CFLAGS = -I@S_HDF5_HOME@/include -L@S_HDF5_HOME@/lib -lhdf5 -lhdf5_fortran
HDF5_DEFINEMACRO = -DS_PHDF5
endif

if HAVE_NETCDF
NETCDF_SRC=../f-io-libs/f-pnetcdf-lib/p_initials.f90 \
../f-io-libs/f-pnetcdf-lib/p_openfile.f90 \
../f-io-libs/f-pnetcdf-lib/p_readfile.f90 \
../f-io-libs/f-pnetcdf-lib/p_writefile.f90 \
../c-io-libs/c-pnetcdf-lib/c_pnetcdf.c
NETCDF_CFLAGS = -I@S_PNETCDF_HOME@/include -L@S_PNETCDF_HOME@/lib -lpnetcdf
NETCDF_DEFINEMACRO = -DS_PNETCDF
endif

if HAVE_ADIOS
ADIOS_SRC=../f-io-libs/f-adios-lib/a_initials.f90 \
../f-io-libs/f-adios-lib/a_openfile.f90 \
../f-io-libs/f-adios-lib/a_readfile.f90 \
../f-io-libs/f-adios-lib/a_writefile.f90 \
../f-io-libs/f-adios-lib/a_finalize.f90 \
../c-io-libs/c-adios-lib/c_adios.c
ADIOS_CFLAGS = -I@S_ADIOS_HOME@/include -L@S_ADIOS_HOME@/lib -ladios -ladiosf -L@S_MXML_HOME@/lib -lmxml
ADIOS_DEFINEMACRO = -DS_ADIOS
endif

libseismio_la_SOURCES += $(HDF5_SRC) $(NETCDF_SRC) $(ADIOS_SRC)
AM_CFLAGS +=  $(HDF5_DEFINEMACRO) $(NETCDF_DEFINEMACRO) $(ADIOS_DEFINEMACRO)  $(HDF5_CFLAGS) $(NETCDF_CFLAGS) $(ADIOS_CFLAGS)
AM_FCFLAGS +=  $(HDF5_DEFINEMACRO) $(NETCDF_DEFINEMACRO) $(ADIOS_DEFINEMACRO)  $(HDF5_CFLAGS) $(NETCDF_CFLAGS) $(ADIOS_CFLAGS)

#includedir =

#../f-io-libs/f-pnetcdf-lib/setup.f90 \
#../f-io-libs/f-pnetcdf-lib/find_start_count.f90 \
#../f-io-libs/f-adios-lib/setup.f90 \
#../f-io-libs/f-adios-lib/find_start_count.f90


#libseismio_la_LIBADD = ../f-io-libs/f-mpiio-lib/libf_mpiio.la ../c-io-libs/c-mpiio-lib/libc_mpiio.la ../f-io-libs/f-phdf5-lib/libf_hdf5.la ../c-io-libs/c-phdf5-lib/libc_hdf5.la ../f-io-libs/f-pnetcdf-lib/libf_netcdf.la ../c-io-libs/c-pnetcdf-lib/libc_netcdf.la ../f-io-libs/f-adios-lib/libf_adios.la ../c-io-libs/c-adios-lib/libc_adios.la

#libseismio_la_LDADD =

#libseismio_la_dir = @prefix@/lib
#libseismio_ladir = ./
#libseismio_la_HEADERS = seism_io.h seism_io_consts.h

#for install include files
include_HEADERS =  seism_io.h seism_io_consts.h seism_io.fh seism_io_consts.fh

# $(top_builddir)/lib/
# -I$(prefix)/include
# -L$(prefix)/lib -lseismio


#-I@S_ADIOS_HOME@/include -I@S_HDF5_HOME@/include -I@S_PNETCDF_HOME@/include
#AM_FCFLAGS = -I@S_ADIOS_HOME@/include -I@S_PNETCDF_HOME@/include -I@S_HDF5_HOME@/include
AUTOMAKE_OPTIONS=foreign
