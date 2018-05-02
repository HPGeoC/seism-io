/**
@section LICENSE
Copyright (c) 2013-2017, Regents of the University of California
All rights reserved.
Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/
#include <mpi.h>
#include <string.h>
#include <stdio.h>
//#include "c_io_libs.h"
#include "adios.h"

extern int DEBUG;
extern int SEISM_APP_LANG; //Not Set:-1 ; C:0 ; Fortran:1

/* all pointers */
/* Notes:
    1. reference fortran in lower case
    2. when linking place c library first
 */

void c_a_initials(int* write_style,int* rank,int* nk,int* kout,int* nj,int* jout,int* ni,int*iout,int* kbnd,int* jbnd,int* ibnd,
      int* recproc,int* nprec,int* start2dim)
{

    a_initials(write_style,rank,nk,kout,nj,jout,ni,iout,kbnd,jbnd,ibnd,recproc,nprec,start2dim);

}

void c_a_openfile(int* rank,int* recproc,MPI_Comm * comm,char* filename, int *write_style,char *rw,
      int* ghostx, int* ghosty, int* ghostz,
      int* nk,int* kout,int* nj,int*jout,int* ni,int* iout,int* maxdim,
      int* kbnd,int* jbnd,int* ibnd,MPI_Datatype* datatype,int* write_step,
      int *nprec, int* prec,
      char* attr,char*attr_name, char*grpname,char* mdname,int *adios_buffer_size,
      char*dimXname,char*dimYname,char*dimZname,char*dimTname,
      MPI_Offset*start,MPI_Offset*count,MPI_Offset*start2,int* start2dim,int*ntime,char*varname,int*icount,
      int*iopen,long*m_adios_group,long*adios_handle)
{

    MPI_Fint fcomm= MPI_Comm_c2f(*comm);
    MPI_Fint fdatatype = MPI_Type_c2f(*datatype);

    int filenameLen  = strlen(filename);
    int attrLen      = strlen(attr);
    int attr_nameLen = strlen(attr_name);
    int grpnameLen   = strlen(grpname);
    int mdnameLen    = strlen(mdname);
    int dimXnameLen  = strlen(dimXname);
    int dimYnameLen  = strlen(dimYname);
    int dimZnameLen  = strlen(dimZname);
    int dimTnameLen  = strlen(dimTname);
    int varnameLen   = strlen(varname);

    /* pass pointers to fortran codes */
    a_openfile(rank,recproc,&fcomm,filename, write_style,&(rw[0]),
          ghostx, ghosty, ghostz,
          nk,kout,nj,jout,ni,iout,maxdim,
          kbnd,jbnd,ibnd,&fdatatype,write_step,
          nprec, prec,
          attr,attr_name, grpname,mdname,adios_buffer_size,
          dimXname,dimYname,dimZname,dimTname,
          start,count,start2,start2dim,ntime,varname,icount,iopen,m_adios_group,adios_handle,
          &filenameLen,&attrLen,&attr_nameLen,&grpnameLen,&mdnameLen,
          &dimXnameLen,&dimYnameLen,&dimZnameLen,&dimTnameLen,&varnameLen);

    // switch the order between C and Fortran
    // when C call this seismIO Lib
    // modified by Hui Zhou 2016.08
    if(SEISM_APP_LANG == 0)
    {
    /* switch the order between C and Fortran */
    int i,j;
    for (i=0; i<*nprec; i++) {
        j=*(prec+2*(*nprec)+i);
        *(prec+2*(*nprec)+i)= *(prec+i);
        *(prec+i)=j;
    }
    }else{ // FORTRAN APP
        ;     //do nothing
    }

}

void c_a_writefile(int* rank,int* recproc,
       MPI_Datatype* datatype,int* write_step,int* nprec,int* prec,
       MPI_Offset* start,MPI_Offset* count,MPI_Offset* start2,int* start2dim,int *ntime,char*varname,int*icount,
       int*iopen,long*adios_handle,
       void*buffer,void*var,int*nx,int*ny,int*nz)
{
    MPI_Fint fdatatype = MPI_Type_c2f(*datatype);

    // reverse dimensions of var
    // reverse dimensions of var when C call this seismIO Lib
    // add by Hui Zhou 2016.08
    if(SEISM_APP_LANG == 0)
    {
        int swaptmp= *nx;
        *nx = *nz;
        *nz = swaptmp;
    }else{ // FORTRAN APP
        ; //do nothing
    }
    // reverse dimensions of var

    int varnameLen = strlen(varname);

    if (*datatype==MPI_INT) {
    a_writefile_integer(rank,recproc,
               &fdatatype,write_step,nprec,prec,
               start,count,start2,start2dim,ntime,varname,icount,iopen,adios_handle,
               (int*)buffer,(int*)var,nx,ny,nz,&varnameLen);
    }

    if (*datatype==MPI_FLOAT) {
    a_writefile_real(rank,recproc,
            &fdatatype,write_step,nprec,prec,
            start,count,start2,start2dim,ntime,varname,icount,iopen,adios_handle,
            (float*)buffer,(float*)var,nx,ny,nz,&varnameLen);
    }
    if (*datatype==MPI_DOUBLE) {
    a_writefile_double(rank,recproc,
              &fdatatype,write_step,nprec,prec,
              start,count,start2,start2dim,ntime,varname,icount,iopen,adios_handle,
              (double*)buffer,(double*)var,nx,ny,nz,&varnameLen);
    }
}



void c_a_readfile(int* rank,int* recproc,
      MPI_Datatype* datatype,int* write_step,int* nprec,int* prec,
      MPI_Offset* start,MPI_Offset* count,MPI_Offset* start2,int* start2dim,int *ntime,char*varname,int*icount,
      int*iopen,long*adios_handle,
      void*buffer,void*var,int*nx,int*ny,int*nz)
{
    MPI_Fint fdatatype = MPI_Type_c2f(*datatype);

    // reverse dimensions of var
    // reverse dimensions of var when C call this seismIO Lib
    // add by Hui Zhou 2016.08
    if(SEISM_APP_LANG == 0)
    {
        int swaptmp= *nx;
        *nx = *nz;
        *nz = swaptmp;
    }else{ // FORTRAN APP
        ; //do nothing
    }
    // reverse dimensions of var

    int varnameLen = strlen(varname);

    if (*datatype==MPI_INT) {
    a_readfile_integer(rank,recproc,
              &fdatatype,write_step,nprec,prec,
              start,count,start2,start2dim,ntime,varname,icount,iopen,adios_handle,
              (int*)buffer,(int*)var,nx,ny,nz,&varnameLen);
    }

    if (*datatype==MPI_FLOAT) {
    a_readfile_real(rank,recproc,
               &fdatatype,write_step,nprec,prec,
               start,count,start2,start2dim,ntime,varname,icount,iopen,adios_handle,
               (float*)buffer,(float*)var,nx,ny,nz,&varnameLen);
    }
    if (*datatype==MPI_DOUBLE) {
    a_readfile_double(rank,recproc,
             &fdatatype,write_step,nprec,prec,
             start,count,start2,start2dim,ntime,varname,icount,iopen,adios_handle,
             (double*)buffer,(double*)var,nx,ny,nz,&varnameLen);
    }
}

//add by Hui 2016.08
void c_a_closefile(int64_t* adios_handle)
{

    a_closefile(adios_handle);

    return;
}

//add by Hui 2016.08
void c_a_read_closefile(int64_t* adios_handle)
{

    a_read_closefile(adios_handle);

    return;
}
