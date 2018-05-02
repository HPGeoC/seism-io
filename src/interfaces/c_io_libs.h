/**
@section LICENSE
Copyright (c) 2013-2017, Regents of the University of California
All rights reserved.
Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/
#ifndef CIOLIBS_H
#define CIOLIBS_H
#include <stdint.h>
#include <inttypes.h>
//MPIIO
extern void c_initials(int* write_style,int* rank,int* nk,int* kout,int* nj,int* jout,int* ni,int*iout,int* kbnd,int* jbnd,int* ibnd,
        int* recproc,int* nprec,int* start2dim);
extern void c_openfile(int* rank,int* recproc,MPI_Comm * comm,char* filename, int *write_style,char *rw,
        int* ghostx, int* ghosty, int* ghostz,
        int* nk,int* kout,int* nj,int*jout,int* ni,int* iout,int* maxdim,
        int* kbnd,int* jbnd,int* ibnd,MPI_Datatype* datatype,int* write_step,
        int *nprec, int* prec,
        int64_t* nrec, MPI_Offset* disp,int* filetype,int*icount,int*fh);
extern void c_writefile(int* rank,int* recproc,
                MPI_Datatype* datatype,int* write_step,int* nprec,int* prec,
                int64_t* nrec,MPI_Offset* disp,int* filetype,int*icount,int*fh,
        void*buffer,void* var,int *nx,int* ny,int* nz);
extern void c_readfile(int* rank,int* recproc,
                MPI_Datatype* datatype,int* write_step,int* nprec,int* prec,
                int64_t* nrec,MPI_Offset* disp,int* filetype,int*icount,int*fh,
        void*buffer,void* var,int *nx,int* ny,int* nz);
extern void c_closefile(int*fh);


//PHDF5
#ifdef S_PHDF5
extern void c_h_initials(int* write_style,int* rank,int* nk,int* kout,int* nj,int* jout,int* ni,
        int*iout,int* kbnd,int* jbnd,int* ibnd,
	int* recproc,int* nprec,int* start2dim);
extern void c_h_openfile(int* rank,int* recproc,MPI_Comm * comm,char* filename, int *write_style,char *rw,
        int* ghostx, int* ghosty, int* ghostz,
        int* nk,int* kout,int* nj,int*jout,int* ni,int* iout,int* maxdim,
        int* kbnd,int* jbnd,int* ibnd,MPI_Datatype* datatype,int* write_step,
        int *nprec, int* prec,
        char* dsetname,char* attr,char* attr_name,
        MPI_Offset*start,MPI_Offset*count,MPI_Offset*start2,int* start2dim, int*ntime,int*icount,hsize_t*   dimsf,int*iopen,hid_t* type_id,
        hid_t* plist_id,hid_t* dataspace,hid_t* dset_id,hid_t* memspace,hid_t* file_id);
extern void c_h_writefile(int* rank,int* recproc,
        MPI_Datatype* datatype,int* write_step,int* nprec,int* prec,
        MPI_Offset* start,MPI_Offset* count,MPI_Offset* start2,int*start2dim,int *ntime,int* icount,hsize_t* dimsf,int* iopen,hid_t*type_id,
        hid_t*plist_id,hid_t*dataspace,hid_t*dset_id,hid_t*memspace,hid_t*file_id,
        void*buffer,void*var,int*nx,int*ny,int*nz);
extern void c_h_readfile(int* rank,int* recproc,
        MPI_Datatype* datatype,int* write_step,int* nprec,int* prec,
        MPI_Offset* start,MPI_Offset* count,MPI_Offset* start2,int*start2dim,int *ntime,int* icount,hsize_t* dimsf,int* iopen,hid_t*type_id,
        hid_t*plist_id,hid_t*dataspace,hid_t*dset_id,hid_t*memspace,hid_t*file_id,
        void*buffer,void*var,int*nx,int*ny,int*nz);
extern void c_h_closefile(hid_t* plist_id,hid_t*dataspace,hid_t*dset_id,hid_t*memspace,hid_t*file_id);
#endif /*PHDF5*/

//PNETCDF
#ifdef S_PNETCDF
extern void c_p_initials(int* write_style,int* rank,int* nk,int* kout,int* nj,int* jout,int* ni,int*iout,int* kbnd,int* jbnd,int* ibnd,
    int* recproc,int* nprec, int* start2dim);
extern void c_p_openfile(int* rank,int* recproc,MPI_Comm * comm,char* filename, int *write_style,char *rw,
    int* ghostx, int* ghosty, int* ghostz,
    int* nk,int* kout,int* nj,int*jout,int* ni,int* iout,int* maxdim,
    int* kbnd,int* jbnd,int* ibnd,MPI_Datatype* datatype,int* write_step,
    int *nprec, int* prec,
    char* attr,char* attr_name,
    char* dimXname,char* dimYname,char* dimZname,char* dimTname,char* varname,
    MPI_Offset* start,MPI_Offset* count,MPI_Offset* start2,int* start2dim,int* ntime,int* xvarid,int* icount,int* ncid);
extern void c_p_writefile(int* rank,int* recproc,
        MPI_Datatype* datatype,int* write_step,int* nprec,int* prec,
        MPI_Offset* start,MPI_Offset* count,MPI_Offset* start2,int*start2dim,
        int *ntime,int* xvarid,int* icount,int* ncid,
        void*buffer,void* var,int *nx,int* ny,int* nz);
extern void c_p_readfile(int* rank,int *recproc,
        MPI_Datatype* datatype,int* write_step,int*nprec,int* prec,
        MPI_Offset* start,MPI_Offset* count,MPI_Offset* start2,int* start2dim,int*ntime,int* xvarid,int* icount,int*ncid,
        void*buffer,void*var,int*nx,int*ny,int*nz);
extern void c_p_closefile(int* ncid);
#endif /*PNETCDF*/

//ADIOS
#ifdef S_ADIOS
extern void c_a_initials(int* write_style,int* rank,int* nk,int* kout,int* nj,int* jout,int* ni,int*iout,int* kbnd,int* jbnd,int* ibnd,
        int* recproc,int* nprec,int* start2dim);
extern void c_a_openfile(int* rank,int* recproc,MPI_Comm * comm,char* filename, int *write_style,char *rw,
        int* ghostx, int* ghosty, int* ghostz,
        int* nk,int* kout,int* nj,int*jout,int* ni,int* iout,int* maxdim,
        int* kbnd,int* jbnd,int* ibnd,MPI_Datatype* datatype,int* write_step,
        int *nprec, int* prec,
        char* attr,char*attr_name, char*grpname,char* mdname,int *adios_buffer_size,
        char*dimXname,char*dimYname,char*dimZname,char*dimTname,
        MPI_Offset*start,MPI_Offset*count,MPI_Offset*start2,int* start2dim,int*ntime,char*varname,int*icount,
        int*iopen,long*m_adios_group,long*adios_handle);
extern void c_a_writefile(int* rank,int* recproc,
        MPI_Datatype* datatype,int* write_step,int* nprec,int* prec,
        MPI_Offset* start,MPI_Offset* count,MPI_Offset* start2,int* start2dim,int *ntime,char*varname,int*icount,
        int*iopen,long*adios_handle,
        void*buffer,void*var,int*nx,int*ny,int*nz);
extern void c_a_readfile(int* rank,int* recproc,
        MPI_Datatype* datatype,int* write_step,int* nprec,int* prec,
        MPI_Offset* start,MPI_Offset* count,MPI_Offset* start2,int* start2dim,int *ntime,char*varname,int*icount,
        int*iopen,long*adios_handle,
        void*buffer,void*var,int*nx,int*ny,int*nz);
extern void c_a_closefile(int64_t* adios_handle);
extern void c_a_read_closefile(int64_t* adios_handle);
#endif /*ADIOS*/





#endif	/* CIOLIBS_H */
