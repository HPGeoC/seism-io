/**
@section LICENSE
Copyright (c) 2013-2017, Regents of the University of California
All rights reserved.
Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/
#ifdef S_PNETCDF

#include <stdio.h>
#include <stdlib.h>
#include "seism_io_globals.h"
#include "seism_io_consts.h"
#include "c_io_libs.h"

/* seism_init_pnetcdf function reads in the parameters text file and
 * //?creates groups and variables for PHDF5 library.
 *
 * No  is necessary.
 */
void seism_initials_pnetcdf(PointSet *ps)
{

    int readwrite_style = (ps->isGrid ? 1 : 2);
    c_p_initials(&readwrite_style, &mpi_rank,
             &(ps->nz), ps->locz,
             &(ps->ny), ps->locy,
             &(ps->nx), ps->locx,
             zbounds,
             ybounds,
             xbounds,
             &(ps->recproc), &(ps->nprec), &(ps->start2dim));

    if((ps->start2dim) >= 0)
    {
        ps->start2 = (MPI_Offset*)malloc(sizeof(MPI_Offset) *4 *(ps->start2dim));
    }
    ps->prec = (int*)malloc(sizeof(int)*(ps->nprec)*maxdim);

    return;
}

/* seism_file_open_pnetcdf function opens the PHDF5 file collectively,
 * prepares the dataspace, dataset and file to write
 *
 * Inputs:
   vname:  variable name to be the dataset
   fname:  PHDF5 filename
   comm:   MPI_Comm
 *
 * Outputs:
   file_id:  file handle for the opened file
   dset_id:  dataset id for the data to be written
   datatype: native datatype of the dataset
 */

void seism_file_open_pnetcdf(FileProfile *fp, int *err)
{

    PointSet *ps = fp->pointSet;
    if(!ps)
    {
        printf("ERROR! SEISM_IO: no PS in FP!\n");
        *err=-1;
        return;
    }

    int readwrite_style = (ps->isGrid ? 1 : 2);
    SpecificVars_pnetcdf *sv
        = (SpecificVars_pnetcdf*)malloc(sizeof(SpecificVars_pnetcdf));

    //if(DEBUG>0 && mpi_rank%400==0) printf("call c_openfile with: rank=%d, recproc=%d, fname=%s, rstyle=%d, fmode=%s\nnxyz=%dx%dx%d, locxyz[010]=%d,%d,%d, maxdim=%d\n",
    //    mpi_rank, ps->recproc, fp->filename, readwrite_style, fp->filemode,
    //    ps->nx, ps->ny, ps->nz, ps->locx[0], ps->locy[1], ps->locz[0], maxdim);
    //if(DEBUG>0 && mpi_rank%400==0) printf("%d) xyzbounds=[%d,%d]x[%d,%d]x[%d,%d]\tnprec=%d\n",
    //  mpi_rank,xbounds[0],xbounds[1],ybounds[0],ybounds[1],zbounds[0],zbounds[1],
    //  ps->nprec);

    //set default attribute
    //add by Hui  2016.08
    char attr[16]="type";
    char attr_name[64]="typename";
    char dimXname[16]="X";
    char dimYname[16]="Y";
    char dimZname[16]="Z";
    char dimTname[16]="T";
    char varname[16]="SeismIOVar";

    c_p_openfile(&mpi_rank, &(ps->recproc), comm, fp->filename, &readwrite_style, fp->filemode,
             &ghostx, &ghosty, &ghostz,
             &(ps->nz), ps->locz,
             &(ps->ny), ps->locy,
             &(ps->nx), ps->locx,
             &maxdim, zbounds, ybounds, xbounds,
             &(fp->datatype), &(fp->write_step),
             &(ps->nprec), ps->prec,
             attr, attr_name,
             dimXname, dimYname, dimZname, dimTname, varname,
             sv->start, sv->count, ps->start2, &(ps->start2dim), &sv->ntime, &sv->xvarid, &sv->icount, &sv->ncid);

    fp->specificVars = sv;
    fp->buffer = (void*)malloc(sizeof(fp->datatype)*(ps->nprec)*(fp->write_step));

    return;
}

/* seism_write_pnetcdf writes the data in the buffer to the file
 * specified by the dataset dset_id
 *
 * Inputs:
   vname:     variable name to be written (a file must have opened beforehand)
 */
void seism_write_pnetcdf(FileProfile *fp, void *tvar, int *err)
{

    int i,j,k;
    i = nxt+2*ghostx;
    j = nyt+2*ghosty;
    k = nzt+2*ghostz;
    PointSet *ps = fp->pointSet;
    SpecificVars_pnetcdf *sv = (SpecificVars_pnetcdf*)(fp->specificVars);

    c_p_writefile(&mpi_rank, &(ps->recproc), &(fp->datatype), &(fp->write_step),
              &(ps->nprec), ps->prec,
              sv->start, sv->count, ps->start2, &ps->start2dim, &sv->ntime, &sv->xvarid, &sv->icount, &sv->ncid,
              fp->buffer, tvar, &i, &j, &k
              );

    //if(DEBUG>0 && mpi_rank%400==0) printf("%d) c_p_writefile returns: icount/ws=%d/%d\n",mpi_rank,sv->icount,fp->write_step);

    return;
}

/* Copy from above seism_write_pnetcdf() . add by Hui 2016.08
   Read the file identified with fp to the variable
   fp    : FileProfile pointer
   var   : void pointer to a multi-dimensional variable
   err   : error code
 */
void seism_read_pnetcdf(FileProfile *fp, void *tvar, int *err)
{

    int i,j,k;
    i = nxt+2*ghostx;
    j = nyt+2*ghosty;
    k = nzt+2*ghostz;
    PointSet *ps = fp->pointSet;
    SpecificVars_pnetcdf *sv = (SpecificVars_pnetcdf*)(fp->specificVars);

    c_p_readfile(&mpi_rank, &(ps->recproc), &(fp->datatype), &(fp->write_step),
             &(ps->nprec), ps->prec,
             sv->start, sv->count, ps->start2, &ps->start2dim, &sv->ntime, &sv->xvarid, &sv->icount, &sv->ncid,
             fp->buffer, tvar, &i, &j, &k
             );

    //if(DEBUG>0 && mpi_rank%400==0) printf("%d) c_p_readfile returns: icount/ws=%d/%d\n",mpi_rank,sv->icount,fp->write_step);

    return;
}

//add by Hui  2016.08
void seism_file_close_pnetcdf(FileProfile *fp, int *err)
{
    SpecificVars_pnetcdf *sv = (SpecificVars_pnetcdf*)fp->specificVars;
    c_p_closefile(&(sv->ncid));
    return;
}
#endif
