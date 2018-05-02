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
#include <stdint.h>
#include "mpi.h"

#include "seism_io_globals.h"
#include "seism_io_consts.h"
#include "seism_io.h"

#ifndef _CALCREC_ONLY
#ifndef _EXEC_SEISM_IO

#else   // _EXEC_SEISM_IO
#define MPI_Comm int
#define MPI_Barrier(comm) 0;
#endif  // end _EXEC_SEISM_IO

#ifdef S_ADIOS
#include "adios.h"
#include "seism_io_adios.h"
#endif
#ifdef S_PHDF5
#include "hdf5.h"
#include "seism_io_phdf5.h"
#endif
#ifdef S_PNETCDF
//    #include "netcdf.h"
#include "pnetcdf.h"
#include "seism_io_pnetcdf.h"
#endif

/* init function to initialize the SEISM lib. Reads in text input file
   and creates the initial group and variable objects.
   comm    : MPI_Comm pointer for initialization. ADIOS 1.3, 1.4 requires
        all processors to call init, not only the writers
   rank    : MPI_rank of the caller (int* because fortran interface passes
        the pointer)
   coords  : maxdim (3D) vector of 0-based processor indices in the virtual
        topology (0,0,0),(1,0,0)...(PX-1,PY-1,PZ-1)
   maxdim  : number of dimensions for the decomposition
   n[x,y,z]: simulation dimensions
   n[xyz]t : sub-domain dimensions
   ghost[] : number of ghost cells per side per dimension
      (if domain is handled by -1:nx+2, then ghostx=2)
   err     : error code
 */
void seism_init_(MPI_Comm *comm, int *rank, int *coords, int *maxdim,
        int *nx, int *ny, int *nz, int *nxt, int *nyt, int *nzt,
        int *ghostx, int *ghosty, int *ghostz,
        int *tpx, int *tpy, int *tpz,
        char *t_seism_method, int *err, unsigned int smlen){
    if(DEBUG>0 && *rank%400==0) printf("%d) Fortran interface is called\n",*rank);
    SEISM_APP_LANG = SEISM_APP_F;
    char *csm = (char*)malloc(sizeof(char)*(smlen+1));
    strncpy(csm,t_seism_method,smlen);
    strncpy(&csm[smlen],&EOS,1);
    if(DEBUG>0 && *rank%400==0) printf("%d) call seism_init with transport method=%s\n",*rank,csm);
    SEISM_APP_LANG = SEISM_APP_F;
    seism_init(comm, rank, coords, maxdim, nx, ny, nz, nxt, nyt, nzt,
           ghostx, ghosty, ghostz, tpx, tpy, tpz, csm, err);

}

/* creates a regular 3D grid to manipulate the variables for.
   nbg[x,y,z]    : int - beginning global index of the grid per dimension (1-based)
   ned[x,y,z]    : int - ending global index
   nskp[x,y,z]   : int - number of points to skip in each dimension (1=no skip)
   gridID      : the index of the PointSet (0-based, kept in the LinkedList)
 */
void seism_createregulargrid_(int *nbgx, int *nedx, int *nskpx,
        int *nbgy, int *nedy, int *nskpy, int *nbgz, int *nedz, int *nskpz,
        int *gridID, int *err)
{
    seism_createRegularGrid(nbgx,nedx,nskpx,nbgy,nedy,nskpy,nbgz,nedz,nskpz,gridID,err);
}

/* file open function fortran interface
   fname   : file name of the file
   fmode   : file mode = "w" "r" "a" as write-read-append
   write_step : write_step, the amount of aggregation in time
   psID    : PointSet ID (input)
   wpID    : WriteProfile ID (output)
   err     : error code pointer. 0 is no error.
   fnamelen: length of the file name
   fmodelen: length of the file mode
 */
void seism_file_open_(char *fname, char *fmode, int *write_step, char *fdata,
        int *psID, int *wpID, int *err,
        unsigned int fnamelen, unsigned int fmodelen, unsigned int fdatalen)
{

    // C version fname, fmode. They will be copied from
    //  Fortran names
    char *cfname, *cfmode, *cfdata;
    cfname = (char*)malloc(sizeof(char)*(fnamelen+1));
    cfmode = (char*)malloc(sizeof(char)*(fmodelen+1));
    cfdata = (char*)malloc(sizeof(char)*(fdatalen+1));
    strncpy(cfname,fname,fnamelen);
    strncpy(&cfname[fnamelen],&EOS,1);
    strncpy(cfmode,fmode,fmodelen);
    strncpy(&cfmode[fmodelen],&EOS,1);
    strncpy(cfdata,fdata,fdatalen);
    strncpy(&cfdata[fdatalen],&EOS,1);

    seism_file_open(cfname, cfmode, write_step, cfdata, psID, wpID, err);

    free(cfname);
    free(cfmode);
    free(cfdata);
}

/* readMesh() function reads in the binary media file
   fname:  media file name
   nvar:   number of variables to read
   mediaData:  multi dimentional matrix containing all the variables
   order:  if >0, first layer is top, if <=0, the first layer is bottom
   err:    returns error code
 */
void seism_readmesh_(char *fname, int *nvar, void *mediaData, int *order,
             int *err,
             unsigned int fnamelen)
{
    char *cfname;
    cfname = (char*)malloc(sizeof(char)*(fnamelen+1));
    strncpy(cfname,fname,fnamelen);
    strncpy(&cfname[fnamelen],&EOS,1);
    if(DEBUG>0 && mpi_rank%400==0) printf("fname:'%s'\ncfname:'%s'\nfnamelen=%d\n",fname,cfname,fnamelen);

    seism_readMesh(cfname, nvar, mediaData, order, err);

    free(cfname);
}

/* write function fortran interface
   seism_f : FileProfile ID integer. File has to be opened before.
   var   : pointer to multi-dimensional variable
   err     : error code, integer pointer. 0 is no error.
 */
void seism_write_(int *seism_f, void *var, int *err)
{

    seism_write(seism_f, var, err);
}

/* add by Hui 2016.08
   read function fortran interface
   seism_f : FileProfile ID integer. File has to be opened before.
   var   : pointer to multi-dimensional variable
   err     : error code, integer pointer. 0 is no error.
 */
void seism_read_(int *seism_f, void *var, int *err)
{

    seism_read(seism_f, var, err);
}

/* seism_file_close closes an already opened file
   seism_f : pointer to file
   err     : error code
 */
void seism_file_close_(int *seism_f, int *err)
{
    seism_file_close(seism_f, err);
}

/* seism_finalize finalizes the library
   err : error code
 */
void seism_finalize_(int *err)
{
    seism_finalize(err);
}
#endif
