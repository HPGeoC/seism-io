/**
@section LICENSE
Copyright (c) 2013-2017, Regents of the University of California
All rights reserved.
Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/
#ifndef SEISM_IO_C_H
#define SEISM_IO_C_H

void seism_createRegularGrid(int *nbgx, int *nedx, int *nskpx,
                             int *nbgy, int *nedy, int *nskpy, int *nbgz, int *nedz, int *nskpz,
                             int *gridID, int *err);
void seism_init(MPI_Comm *tcomm, int *rank, int *tcoords, int *tmaxdim,
                int *tnx, int *tny, int *tnz,
                int *tnxt, int *tnyt, int *tnzt,
                int *tghostx, int *tghosty, int *tghostz,
                int *tpx, int *tpy, int *tpz,
                char *t_seism_method, int *err);
void seism_file_open(char *fname, char *fmode, int *write_step, char *fdata,
                     int *psID, int *fpID, int *err);
void seism_readMesh(char *fname, int *nvar, void *mediaData, int *order,
                    int *err);
void seism_write(int *seism_f, void *var, int *err);
void seism_read(int *seism_f, void *var, int *err);
void seism_file_close(int *seism_f, int *err);
void seism_finalize(int *err);

#endif /*SEISM_IO_C_H*/
