/**
@section LICENSE
Copyright (c) 2013-2017, Regents of the University of California
All rights reserved.
Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/
#ifdef S_PHDF5

#include <stdio.h>
#include <stdlib.h>
#include "hdf5.h"
#include "seism_io_globals.h"
#include "seism_io_consts.h"
#include "c_io_libs.h"

const int DATATYPELEN = 6;
const char const *datatypes[] = {"integer", "char", "real", "float", "double", "real*8"};
hid_t *H5types;

/* seism_init_phdf5 function reads in the parameters text file and
 * //?creates groups and variables for PHDF5 library.
 *
 * No  is necessary.
 */
void seism_initials_phdf5(PointSet *ps)
{
    /*
       void seism_init_phdf5(MPI_Comm comm){
     */
    H5types = (hid_t*)malloc(sizeof(hid_t)*DATATYPELEN);
    H5types[0] = H5T_NATIVE_INT;
    H5types[1] = H5T_NATIVE_CHAR;
    H5types[2] = H5T_NATIVE_FLOAT;
    H5types[3] = H5T_NATIVE_FLOAT;
    H5types[4] = H5T_NATIVE_DOUBLE;
    H5types[5] = H5T_NATIVE_LDOUBLE;

    int readwrite_style = (ps->isGrid ? 1 : 2);
    c_h_initials(&readwrite_style, &mpi_rank,
             &(ps->nz), ps->locz,
             &(ps->ny), ps->locy,
             &(ps->nx), ps->locx,
             zbounds,
             ybounds,
             xbounds,
             &(ps->recproc), &(ps->nprec), &(ps->start2dim));
    //printf("nprec=%d\n",ps->nprec);
    if((ps->start2dim) >= 0)
    {
        ps->start2 = (MPI_Offset*)malloc(sizeof(MPI_Offset) *4 *(ps->start2dim));
    }
    ps->prec = (int*)malloc(sizeof(int)*(ps->nprec)*maxdim);

    return;
}

/* seism_file_open_phdf5 function opens the PHDF5 file collectively,
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
//(FileProfile *fp, int *err)
void seism_file_open_phdf5(FileProfile *fp, int *err)
{

    PointSet *ps = fp->pointSet;
    if(!ps)
    {
        printf("ERROR! SEISM_IO: no PS in FP!\n");
        *err=-1;
        return;
    }
    //printf("PS found!\n");
    //printf("PS[0],PS[5]=%d,%d\n",ps->prec[0],ps->prec[5]);
    int readwrite_style = (ps->isGrid ? 1 : 2);
    SpecificVars_phdf5 *sv
        = (SpecificVars_phdf5*)malloc(sizeof(SpecificVars_phdf5));

    //set default dataset name
    char dsetname[64] = "seismiodataset";
    char attr[64]="type";
    char attr_name[64]="typename";

    c_h_openfile(&mpi_rank, &(ps->recproc), comm, fp->filename, &readwrite_style, fp->filemode,
             &ghostx, &ghosty, &ghostz,
             &(ps->nz), ps->locz,
             &(ps->ny), ps->locy,
             &(ps->nx), ps->locx,
             &maxdim, zbounds, ybounds, xbounds,
             &(fp->datatype), &(fp->write_step),
             &(ps->nprec), ps->prec,
             dsetname, attr, attr_name,
             sv->start, sv->count, ps->start2, &(ps->start2dim), &sv->ntime, &sv->icount,
             sv->dimsf, &sv->iopen, &sv->type_id, &sv->plist_id, &sv->dataspace, &sv->dset_id, &sv->memspace, &sv->file_id); //,

    fp->specificVars = sv;
    fp->buffer = (void*)malloc(sizeof(fp->datatype)*(ps->nprec)*(fp->write_step));

    return;
}

/* seism_write_phdf5 writes the data in the buffer to the file
 * specified by the dataset dset_id
 *
 * Inputs:
   vname:     variable name to be written (a file must have opened beforehand)
 */
void seism_write_phdf5(FileProfile *fp, void *tvar, int *err)
{

    int i,j,k;
    i = nxt+2*ghostx;
    j = nyt+2*ghosty;
    k = nzt+2*ghostz;
    PointSet *ps = fp->pointSet;
    SpecificVars_phdf5 *sv = (SpecificVars_phdf5*)(fp->specificVars);

    c_h_writefile(&mpi_rank, &(ps->recproc), &(fp->datatype), &(fp->write_step),
              &(ps->nprec), ps->prec,
              sv->start, sv->count, ps->start2, &ps->start2dim, &sv->ntime, &sv->icount,
              sv->dimsf, &sv->iopen, &sv->type_id, &sv->plist_id, &sv->dataspace, &sv->dset_id, &sv->memspace, &sv->file_id,
              fp->buffer, tvar, &i, &j, &k
              );

    return;
}

/* Copy from above seism_write_phdf5() . add by Hui 2016.08
   Read the file identified with fp to the variable
   fp    : FileProfile pointer
   var   : void pointer to a multi-dimensional variable
   err   : error code
 */
void seism_read_phdf5(FileProfile *fp, void *tvar, int *err)
{

    int i,j,k;
    i = nxt+2*ghostx;
    j = nyt+2*ghosty;
    k = nzt+2*ghostz;
    PointSet *ps = fp->pointSet;
    SpecificVars_phdf5 *sv = (SpecificVars_phdf5*)(fp->specificVars);

    c_h_readfile(&mpi_rank, &(ps->recproc), &(fp->datatype), &(fp->write_step),
             &(ps->nprec), ps->prec,
             sv->start, sv->count, ps->start2, &ps->start2dim, &sv->ntime, &sv->icount,
             sv->dimsf, &sv->iopen, &sv->type_id, &sv->plist_id, &sv->dataspace, &sv->dset_id, &sv->memspace, &sv->file_id,
             fp->buffer, tvar, &i, &j, &k
             );

    return;
}

//add by Hui 2016.08
void seism_file_close_phdf5(FileProfile *fp, int *err)
{

    PointSet *ps = fp->pointSet;
    if(mpi_rank!= ps->recproc) return;
    SpecificVars_phdf5 *sv = (SpecificVars_phdf5*)fp->specificVars;
    c_h_closefile(&(sv->plist_id), &(sv->dataspace), &(sv->dset_id), &(sv->memspace), &(sv->file_id));
    return;
}
#endif
