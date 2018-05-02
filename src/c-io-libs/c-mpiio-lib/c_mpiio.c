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
#include <stdint.h>
#include <inttypes.h>
/*
   MPI_Fint MPI_Type_c2f(MPI_Datatype datatype)
   MPI_Fint MPI_Comm_c2f(MPI_Comm comm)
 */

extern int DEBUG;
extern int SEISM_APP_LANG; //Not Set:-1 ; C:0 ; Fortran:1

/* all pointers */
void c_initials(int* write_style,int* rank,int* nk,int* kout,int* nj,int* jout,int* ni,int*iout,int* kbnd,int* jbnd,int* ibnd,
        int* recproc,int* nprec,int* start2dim){

    initials(write_style,rank,nk,kout,nj,jout,ni,iout,kbnd,jbnd,ibnd,recproc,nprec,start2dim);

}

void c_openfile(int* rank,int* recproc,MPI_Comm * comm,char* filename, int *write_style,char *rw,
        int* ghostx, int* ghosty, int* ghostz,
        int* nk,int* kout,int* nj,int*jout,int* ni,int* iout,int* maxdim,
        int* kbnd,int* jbnd,int* ibnd,MPI_Datatype* datatype,int* write_step,
        int *nprec, int* prec,
        int64_t* nrec,MPI_Offset* disp,int* filetype,int*icount,int*fh){

    MPI_Fint fcomm= MPI_Comm_c2f(*comm);
    MPI_Fint fdatatype = MPI_Type_c2f(*datatype);

    int filenameLen = strlen(filename);

    /* pass pointers to fortran codes */
    openfile(rank,recproc,&fcomm,filename, write_style,&(rw[0]),
                ghostx, ghosty, ghostz,
                nk,kout,nj,jout,ni,iout,maxdim,
                kbnd,jbnd,ibnd,&fdatatype,write_step,
                nprec, prec,
                nrec,disp,filetype,icount,fh,&filenameLen);

    // switch the order between C and Fortran
    // when C call this seismIO Lib
    // modified by Hui Zhou  2016.08
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
        ;             //do nothing
    }

    return;
}

void c_writefile(int* rank,int* recproc,
         MPI_Datatype* datatype,int* write_step,int* nprec,int* prec,
         int64_t* nrec,MPI_Offset* disp,int* filetype,int*icount,int*fh,
         void*buffer,void* var,int *nx,int* ny,int* nz){

    MPI_Fint fdatatype = MPI_Type_c2f(*datatype);
    int lang = SEISM_APP_LANG;
    // reverse dimensions of var
    // reverse dimensions of var when C call this seismIO Lib
    // add by Hui Zhou  2016.08
    if(SEISM_APP_LANG == 0)
    {
        int swaptmp= *nx;
        *nx = *nz;
        *nz = swaptmp;
    }
    else{ // FORTRAN APP
        ; //do nothing
    }

    if (*datatype==MPI_INT) {
        writefile_integer(rank,recproc,
                         &fdatatype,write_step,nprec,prec,
                         nrec, disp,filetype,icount,fh,
                         (int*)buffer,(int*)var,nx,ny,nz);
    }
    if (*datatype==MPI_FLOAT) {
        writefile_real(rank,recproc,
                          &fdatatype,write_step,nprec,prec,
                          nrec, disp,filetype,icount,fh,
                          (float*)buffer,(float*)var,nx,ny,nz);
    }
    if (*datatype==MPI_DOUBLE) {
        writefile_double(rank,recproc,
                        &fdatatype,write_step,nprec,prec,
                        nrec, disp,filetype,icount,fh,
                        (double*)buffer,(double*)var,nx,ny,nz);
    }
    return;
}


void c_readfile(int* rank,int* recproc,
        MPI_Datatype* datatype,int* write_step,int* nprec,int* prec,
        int64_t* nrec,MPI_Offset* disp,int* filetype,int*icount,int*fh,
        void*buffer,void* var,int *nx,int* ny,int* nz){

    MPI_Fint fdatatype = MPI_Type_c2f(*datatype);
    int lang = SEISM_APP_LANG;
    // reverse dimensions of var
    // reverse dimensions of var when C call this seismIO Lib
    // add by Hui Zhou 2016.08
    if(SEISM_APP_LANG == 0)
    {
        int swaptmp= *nx;
        *nx = *nz;
        *nz = swaptmp;
    }
    else{ // FORTRAN APP
        ; //do nothing
    }

    if (*datatype==MPI_INT) {
        readfile_integer(rank,recproc,
                        &fdatatype,write_step,nprec,prec,
                        nrec, disp,filetype,icount,fh,
                        (int*)buffer,(int*)var,nx,ny,nz);
    }
    if (*datatype==MPI_FLOAT) {
        readfile_real(rank,recproc,
                         &fdatatype,write_step,nprec,prec,
                         nrec, disp,filetype,icount,fh,
                         (float*)buffer,(float*)var,nx,ny,nz);
    }
    if (*datatype==MPI_DOUBLE) {
        readfile_double(rank,recproc,
                           &fdatatype,write_step,nprec,prec,
                           nrec, disp,filetype,icount,fh,
                           (double*)buffer,(double*)var,nx,ny,nz);
    }
    return;
}

void c_closefile(int*fh)
{

    closefile(fh);

    return;
}
