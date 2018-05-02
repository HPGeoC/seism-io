/**
@section LICENSE
Copyright (c) 2013-2017, Regents of the University of California
All rights reserved.
Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/
#ifndef SEISM_CONSTS_H
#define SEISM_CONSTS_H

#define MB1 (1048576)
#define GB1 (1073741824)
#define GB10 (10737418240)
#define GB100 (107374182400)
#define TB1 (1099511627776)
// SEISM METHODS - C
#define SEISM_ADIOS 0   // USE ADIOS
#define SEISM_PHDF5 1   // USE PHDF5
#define SEISM_PNETC 2   // USE PNETCDF
#define SEISM_MPIIO 3   // USE RAW MPI-IO
#define SEISM_DEFAULT 3
// SEISM USER APP LANGUAGES
#define SEISM_APP_NOT_SET (-1)
#define SEISM_APP_C 0
#define SEISM_APP_F 1

// constants from seism_io.h
// set DEBUG to be 0 to disable printf statements
extern const int  DEBUG;
// set end of string EOS
extern const char EOS;
// OPTION LIST and SIZE
//extern const char const *opt[];
extern const int n_opt;
// end OPTION LIST
extern int SEISM_APP_LANG;

extern const char seismParamsFname[];

#endif
