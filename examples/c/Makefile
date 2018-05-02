#@section LICENSE
#Copyright (c) 2013-2017, Regents of the University of California
#All rights reserved.
#Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
#1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
#2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
#THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

CC = mpicc
FC = mpif90

CFLAGS = -std=c99 -g -Dsun -DS_PHDF5 -DS_PNETCDF -DS_ADIOS

HOME_SEISMIO=/path/to/seismio
HOME_PHDF5=/path/to/phdf5-1.8.19
HOME_PNETCDF=/path/to/pnetcdf-1.7.0
HOME_ADIOS=/path/to/adios-1.9.0
HOME_MXML=/path/to/mxml-2.11

INCDIR  = -I$(HOME_ADIOS)/include -I$(HOME_PHDF5)/include -I$(HOME_PNETCDF)/include -I$(HOME_MXML)/include -I$(HOME_SEISMIO)/include
OBJECTS	= testInterface.o
LIB	= -L$(HOME_ADIOS)/lib -L$(HOME_PHDF5)/lib -L$(HOME_PNETCDF)/lib -L$(HOME_MXML)/lib -L$(HOME_SEISMIO)/lib $(HOME_SEISMIO)/lib/libseismio.a -lm -ldl -lhdf5 -lhdf5_fortran -lpnetcdf -ladios -ladiosf -lmxml -lpthread

testInterface:	$(OBJECTS)
	$(FC) $(CFLAGS) $(INCDIR) -o	testInterface	$(OBJECTS)	$(LIB)

testInterface.o:	testInterface.c
	$(CC) $(CFLAGS) $(INCDIR) -c -o testInterface.o	testInterface.c

clean:
	rm -rf *.o testInterface
