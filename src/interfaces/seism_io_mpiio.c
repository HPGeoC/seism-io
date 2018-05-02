/**
@section LICENSE
Copyright (c) 2013-2017, Regents of the University of California
All rights reserved.
Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "seism_io_globals.h"
#include "seism_io_consts.h"
#include "c_io_libs.h"
#include <stdint.h>
#include <inttypes.h>
/* initialization function for MPI-IO
 * recproc, nprec, and start2dim are set by the
   initializer.
   ps  : PointSet for which we compute the recproc, nprec etc.
 */
void seism_initials_mpiio(PointSet *ps)
{

    int readwrite_style = (ps->isGrid ? 1 : 2);

    c_initials(&readwrite_style, &mpi_rank,
           &(ps->nz), ps->locz,
           &(ps->ny), ps->locy,
           &(ps->nx), ps->locx,
           zbounds,
           ybounds,
           xbounds,
           &(ps->recproc), &(ps->nprec), &(ps->start2dim));
    //printf("nprec=%d\n",ps->nprec);
    ps->prec = (int*)malloc(sizeof(int)*(ps->nprec)*maxdim);

    return;
}

/* Opens file for MPI-IO
   fp      : FileProfile
   err     : error code pointer. 0 is no error.
 */
void seism_file_open_mpiio(FileProfile *fp, int *err)
{

    PointSet *ps = fp->pointSet;
    if(!ps)
    {
        printf("ERROR! SEISM_IO: no PS in FP!\n");
        *err=-1;
        return;
    }

    int readwrite_style = (ps->isGrid ? 1 : 2);
    SpecificVars_mpiio *svars
        = (SpecificVars_mpiio*)malloc(sizeof(SpecificVars_mpiio));

    int64_t nrec;
    int  icount;

    MPI_File fh;
    MPI_Datatype filetype;
    MPI_Offset disp;

    //Uncommeted and adjusted fcomm and fdatatype by Hui 2016.07
    c_openfile(&mpi_rank, &(ps->recproc), comm, fp->filename, &readwrite_style, fp->filemode,
           &ghostx, &ghosty, &ghostz,
           &(ps->nz), ps->locz,
           &(ps->ny), ps->locy,
           &(ps->nx), ps->locx,
           &maxdim,
           zbounds,
           ybounds,
           xbounds,
           &(fp->datatype), &(fp->write_step),
           &(ps->nprec), ps->prec,
           &nrec, &disp, &filetype, &icount, &fh);

    svars->nrec = nrec;
    svars->disp = disp;
    svars->filetype = filetype;
    svars->icount = icount;
    svars->fh = fh;
    fp->specificVars = svars;
    fp->buffer = (void*)malloc(sizeof(fp->datatype)*(ps->nprec)*(fp->write_step));

    return;
}

/* Writes out the variable to the file identified with fp
   fp    : FileProfile pointer
   var   : void pointer to a multi-dimensional variable
   err   : error code
 */
void seism_write_mpiio(FileProfile *fp, void *tvar, int *err)
{

    int i,j,k;
    float *var = (float*)tvar;

    i = nxt+2*ghostx;
    j = nyt+2*ghosty;
    k = nzt+2*ghostz;
    PointSet *ps = fp->pointSet;
    SpecificVars_mpiio *sv = (SpecificVars_mpiio*)(fp->specificVars);

    c_writefile(&mpi_rank, &(ps->recproc), &(fp->datatype), &(fp->write_step),
            &(ps->nprec), ps->prec,
            &(sv->nrec), &(sv->disp), &(sv->filetype), &(sv->icount), &(sv->fh),
            fp->buffer, tvar, &i, &j, &k);

    return;
}

/* Copy from above seism_write_mpiio() . add by Hui 2016.08
   Read the file identified with fp to the variable
   fp    : FileProfile pointer
   var   : void pointer to a multi-dimensional variable
   err   : error code
 */
void seism_read_mpiio(FileProfile *fp, void *tvar, int *err)
{

    int i,j,k;
    float *var = (float*)tvar;

    i = nxt+2*ghostx;
    j = nyt+2*ghosty;
    k = nzt+2*ghostz;

    PointSet *ps = fp->pointSet;
    SpecificVars_mpiio *sv = (SpecificVars_mpiio*)(fp->specificVars);

    c_readfile(&mpi_rank, &(ps->recproc), &(fp->datatype), &(fp->write_step),
           &(ps->nprec), ps->prec,
           &(sv->nrec), &(sv->disp), &(sv->filetype), &(sv->icount), &(sv->fh),
           //fp->buffer, tvar, &(nx+2*ghostx), &(ny+2*ghosty), &(nz+2*ghostz));
           fp->buffer, tvar, &i, &j, &k);

    return;
}


/* Close MPIIO file
   fp    : FileProfile pointer
   err   : error code
 */
void seism_file_close_mpiio(FileProfile *fp, int *err)
{
    SpecificVars_mpiio *sv = fp->specificVars;
    c_closefile(&(sv->fh));
    return;
}
