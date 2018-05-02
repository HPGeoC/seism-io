/**
@section LICENSE
Copyright (c) 2013-2017, Regents of the University of California
All rights reserved.
Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/
#include "mpi.h"

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

#ifndef SEISM_IO_H
#define SEISM_IO_H
#define MAX_DIM 6
#define MAX_GROUP 20
#define MAX_CHAR_LEN 250
#define MAX_NUMS_DELIMS 75
#define DUMP(varname) fprintf(stderr,"%s = %x",# varname,varname);

extern int seism_method;

// MPI_rank of this processor
extern int mpi_rank;

extern int groupID;

extern int maxdim;
extern int nx, ny, nz, px, py, pz;
extern int *coords;
extern int nxt, nyt, nzt;
extern int ghostx, ghosty, ghostz;
extern int xbounds[2], ybounds[2], zbounds[2];

// MPI specific globals:
extern MPI_Comm *comm;

/* struct
   add by Hui Zhou, because cannot find any defination about struct Group.
   ID  : unique ID
   data: void * to the data which is carried by the linked list
   next/prev: next and previous pointers
 */
typedef struct Group {
    int id;
    char name[64];
    char method[64];
    void *val;
    struct Group *next;
    struct Group *prev;
    struct Group *firstVar;
} Group;

/* struct
   add by Hui Zhou, because cannot find any defination about struct Group.
   ID  : unique ID
   data: void * to the data which is carried by the linked list
   next/prev: next and previous pointers
 */
typedef struct Var {
    int id;
    char name[64];
    char type[64];
    struct Var *var;
    void *val;
    char *dims; //?
    char *globalDims;
    char *offsets;
    char *attrKey;
    char *attrVal;
    struct Group *next;
    struct Group *group;
    struct Group *prev;
    char *globalDimsNumeric; //?
} Var;

/* Basic Linked List
   ID  : unique ID
   data: void * to the data which is carried by the linked list
   next/prev: next and previous pointers
 */
typedef struct List {
    int ID;
    void *data;
    struct List *next;
    struct List *prev;
} List;

/* PointSet object
   isGrid    : true if this PointSet is a regular grid - boolean
   n[x,y,z]  : number of points in each direction if grid
        number of points if not grid - int
   loc[x,y,z]: location indices of the points (1-based)
   recproc   : =rank if this rank has a receiver
   nprec     : #receivers in this rank
   prec      : receiver locations (size = nprec*maxdim*sizeof(int))
   start2dim : int
 */
typedef struct PointSet {
    int isGrid;
    int nx;
    int ny;
    int nz;
    int *locx;
    int *locy;
    int *locz;
    int recproc;
    int nprec;
    int *prec;
    int start2dim;
    MPI_Offset *start2; //size:[4][start2dim];
} PointSet;

/* FileProfile object for (PointSet, variable, write_step)
   write_step  : write_step
   filename    : filename specific to this FileProfile
   filemode    : file mode 'r', 'w', 'a+', etc.
   pointSet    : pointer to the PointSet used in this profile
   datatype    : elementary (non-derived) MPI_Datatype of the variable
        to be written (real, float, int, char, byte, etc.)
   specificVars: pointer to the specificVars for the transport method used
   buffer      : the buffer to be used to aggregate data for writing
 */
typedef struct FileProfile {
    int write_step;
    char *filename;
    char *filemode;
    struct PointSet *pointSet;
    MPI_Datatype datatype;
    void *specificVars;
    void *buffer;
} FileProfile;

/* MPI_IO Specific variables: all are pointers.
   nrec  : int - number of records
   disp  : Offset - displacement
   filetype  : filetype
   icount    : counter in 0,..,write_step-1
   fh    : MPI file handler
 */
typedef struct SpecificVars_mpiio {
    //for Overflow issure modified by Hui Zhou 2016.11
    int64_t nrec;
    MPI_Offset disp;
    MPI_Datatype filetype;
    int icount;
    MPI_File fh;
} SpecificVars_mpiio;

#ifdef S_PHDF5
/* struct
   add by Hui Zhou 2016.08
   next/prev: next and previous pointers
 */
typedef struct SpecificVars_phdf5 {
    int nrec;
    MPI_Offset disp;
    int filetype;
    int icount;
    int fh;

    MPI_Offset start[4];
    MPI_Offset count[4];
    //int start2dim= 3;  //TODO ??

    int ntime;
    int iopen;
    hsize_t dimsf[4];

    hid_t file_id;
    hid_t dset_id;
    hid_t datatype;
    hid_t memspace;
    hid_t type_id, plist_id, dataspace;

    int id;
    char name[64];
    char type[64];
    struct Var *var;
    void *val;
    char *dims; //?
    char *globalDims;
    char *offsets;
    char *attrKey;
    char *attrVal;
    struct Group *next;
    struct Group *group;
    struct Group *prev;
    char *globalDimsNumeric; //?
} SpecificVars_phdf5;
#endif

#ifdef S_PNETCDF
/* struct
   add by Hui Zhou 2016.08
   next/prev: next and previous pointers
 */
typedef struct SpecificVars_pnetcdf {
    int nrec;
    MPI_Offset disp;
    int filetype;
    int icount;
    int fh;

    MPI_Offset start[4];
    MPI_Offset count[4];

    int ntime;
    int iopen;
    hsize_t dimsf[4];

    int ncid;
    int xvarid;

    int id;
    char name[64];
    char type[64];
    struct Var *var;
    void *val;
    char *dims; //?
    char *globalDims;
    char *offsets;
    char *attrKey;
    char *attrVal;
    struct Group *next;
    struct Group *group;
    struct Group *prev;
    char *globalDimsNumeric; //?
} SpecificVars_pnetcdf;
#endif

#ifdef S_ADIOS
/* struct
   add by Hui Zhou 2016.08
   next/prev: next and previous pointers
 */
typedef struct SpecificVars_adios {
    int nrec;
    MPI_Offset disp;
    int filetype;
    int icount;

    int64_t adios_handle;
    int64_t m_adios_group;
    int adios_buffer_size;

    MPI_Offset start[4];
    MPI_Offset count[4];

    int ntime;
    int iopen;
    hsize_t dimsf[4];

    int ncid;
    int xvarid;

    int id;
    char name[64];
    char type[64];
    struct Var *var;
    void *val;
    char *dims; //?
    char *globalDims;
    char *offsets;

} SpecificVars_adios;
#endif

#endif
