/**
@section LICENSE
Copyright (c) 2013-2017, Regents of the University of California
All rights reserved.
Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/
#ifdef S_ADIOS
#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include "adios.h"
#include "adios_read.h"
#include "seism_io_globals.h"
#include "seism_io_consts.h"
#include "c_io_libs.h"


void seism_initials_adios(FileProfile *fp)
{

    PointSet *ps = fp->pointSet;
    int readwrite_style = (ps->isGrid ? 1 : 2);
    c_a_initials(&readwrite_style, &mpi_rank,
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

// calls file open for adios
void seism_file_open_adios(FileProfile *fp, int *err)
{

    PointSet *ps = fp->pointSet;
    if(!ps) { printf("ERROR! SEISM_IO: no PS in FP!\n"); *err=-1; return; }

    SpecificVars_adios *sv
        = (SpecificVars_adios*)malloc(sizeof(SpecificVars_adios));

    sv->icount= 0;
    sv->iopen = 0;

    //set the size of ADIOS buffer in MB
    sv->adios_buffer_size = (sizeof(fp->datatype)*(ps->nprec)/1024/1024)+1;
    if( sv->adios_buffer_size < 10) sv->adios_buffer_size=10;

    int readwrite_style = (ps->isGrid ? 1 : 2);

    //set default attribute
    char attr[16]="SeismIOtype";
    char attr_name[64]="typename";
    char grpname[64]="SeismIOgroup";
    char mdname[64]="MPI";
    char dimXname[16]="X";
    char dimYname[16]="Y";
    char dimZname[16]="Z";
    char dimTname[16]="T";
    char varname[16]="SeismIOVar";

    ///TODO: adios open moved to write() when write function valid
    if(strcmp(fp->filemode,"r") == 0)
    {

        c_a_openfile(&mpi_rank, &(ps->recproc), comm, fp->filename, &readwrite_style, fp->filemode,
                 &ghostx, &ghosty, &ghostz,
                 &(ps->nz), ps->locz,
                 &(ps->ny), ps->locy,
                 &(ps->nx), ps->locx,
                 &maxdim, zbounds, ybounds, xbounds,
                 &(fp->datatype), &(fp->write_step),
                 &(ps->nprec), ps->prec,
                 attr, attr_name, grpname, mdname, &sv->adios_buffer_size,
                 dimXname, dimYname, dimZname, dimTname,
                 sv->start, sv->count, ps->start2, &(ps->start2dim), &sv->ntime, varname, &sv->icount,
                 &sv->iopen,&sv->m_adios_group,&sv->adios_handle);
    }

    fp->specificVars = sv;
    fp->buffer = (void*)malloc(sizeof(fp->datatype)*(ps->nprec)*(fp->write_step));

    return;
}

// calls adios_write function for the given variable name and file pointer
// modified by Hui Zhou
void seism_write_adios(FileProfile *fp, void *tvar, int *err)
{

    int i,j,k;
    i = nxt+2*ghostx;
    j = nyt+2*ghosty;
    k = nzt+2*ghostz;
    PointSet *ps = fp->pointSet;
    SpecificVars_adios *sv = (SpecificVars_adios*)(fp->specificVars);
    int readwrite_style = (ps->isGrid ? 1 : 2);

    //set default attribute
    char attr[16]="SeismIOtype";
    char attr_name[64]="typename";
    char grpname[64]="SeismIOgroup";
    char mdname[64]="MPI";
    char dimXname[16]="X";
    char dimYname[16]="Y";
    char dimZname[16]="Z";
    char dimTname[16]="T";
    char varname[16]="SeismIOVar";

    c_a_openfile(&mpi_rank, &(ps->recproc), comm, fp->filename, &readwrite_style, fp->filemode,
             &ghostx, &ghosty, &ghostz,
             &(ps->nz), ps->locz,
             &(ps->ny), ps->locy,
             &(ps->nx), ps->locx,
             &maxdim, zbounds, ybounds, xbounds,
             &(fp->datatype), &(fp->write_step),
             &(ps->nprec), ps->prec,
             attr, attr_name, grpname, mdname, &sv->adios_buffer_size,
             dimXname, dimYname, dimZname, dimTname,
             sv->start, sv->count, ps->start2, &(ps->start2dim), &sv->ntime, varname, &sv->icount,
             &sv->iopen,&sv->m_adios_group,&sv->adios_handle);

    c_a_writefile(&mpi_rank, &(ps->recproc), &(fp->datatype), &(fp->write_step),
              &(ps->nprec), ps->prec,
              sv->start, sv->count, ps->start2, &ps->start2dim, &sv->ntime, varname, &sv->icount,
              &sv->iopen, &sv->adios_handle,
              fp->buffer, tvar, &i, &j, &k
              );

    //prevent some rank close adios_handle in advance.
    MPI_Barrier(*comm);

    //if(mpi_rank == ps->recproc) c_a_closefile(&(sv->adios_handle));
    if(mpi_rank == ps->recproc)
    {
        if(sv->icount == fp->write_step)
        {
            c_a_closefile(&(sv->adios_handle));
            sv->icount = 0;
        }
    }
    return;
}

// calls adios_read function for the given variable name and file pointer
void seism_read_adios(FileProfile *fp, void *tvar, int *err)
{

    int i,j,k;
    i = nxt+2*ghostx;
    j = nyt+2*ghosty;
    k = nzt+2*ghostz;
    PointSet *ps = fp->pointSet;
    SpecificVars_adios *sv = (SpecificVars_adios*)(fp->specificVars);
    int readwrite_style = (ps->isGrid ? 1 : 2);

    //set default attribute
    char attr[16]="SeismIOtype";
    char attr_name[64]="typename";
    char grpname[64]="SeismIOgroup";
    char mdname[64]="MPI";
    char dimXname[16]="X";
    char dimYname[16]="Y";
    char dimZname[16]="Z";
    char dimTname[16]="T";
    char varname[16]="SeismIOVar";

    c_a_readfile(&mpi_rank, &(ps->recproc), &(fp->datatype), &(fp->write_step),
             &(ps->nprec), ps->prec,
             sv->start, sv->count, ps->start2, &ps->start2dim, &sv->ntime, varname, &sv->icount,
             &sv->iopen, &sv->adios_handle,
             fp->buffer, tvar, &i, &j, &k
             );

    return;
}

// add by Hui 2016.08
void seism_file_close_adios(FileProfile *fp, int *err)
{

    SpecificVars_adios *sv = (SpecificVars_adios*)(fp->specificVars);

    ///TODO: adios close moved to write() when write function valid
    if(strcmp(fp->filemode, "r") == 0) {
        c_a_read_closefile(&(sv->adios_handle));
        sv->icount = 0;
    }

    return;
}
#endif
