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
#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <unistd.h>

#include "seism_io_globals.h"
#include "seism_io_consts.h"
#include "seism_io_mpiio.h"
#include "seism_io_read_binary.h"
#include "seism_io_cuda.h"
#include "c_io_libs.h"

#ifndef _CALCREC_ONLY
#ifndef _EXEC_SEISM_IO
#include "mpi.h"
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

const int DEBUG = 2;
const char const *opt[] = {"group","var"};
const int n_opt = 2;

// constants
const char DEFAULT_TRANSPORT_METHOD[]="MPI";
const char EOS='\0';
const char intChar[]="integer";
int SEISM_APP_LANG = SEISM_APP_NOT_SET;

const char seismParamsFname[] = "./seism-params.txt";
int mpi_rank = -1;
int seism_method = SEISM_DEFAULT;

List *pointSets = NULL;
List *fileProfiles = NULL;

int nPointSets = 0;
int nFileProfiles = 0;
int maxdim;
int nx, ny, nz, nxt, nyt, nzt, ghostx, ghosty, ghostz;
int px, py, pz;
int *coords;
int xbounds[2], ybounds[2], zbounds[2];
int recproc, nprec, start2dim;

MPI_Comm *comm;




/* returns the PointSet with ID=ID, if does not exist, returns NULL
   ID  : int - PointSet ID
 */
PointSet* _getPointSet(int ID)
{
    List *tmp;
    tmp = pointSets;
    if(!tmp) return NULL;
    if(tmp->ID == ID) return (PointSet*)(tmp->data);
    while(tmp->next) {
        tmp = tmp->next;
        if(tmp->ID == ID) return (PointSet*)(tmp->data);
    }
    return NULL;
}

// The same as _getPointSet except that traverses FileProfiles
FileProfile* _getFileProfile(int ID)
{
    List *tmp;
    tmp = fileProfiles;
    if(!tmp) return NULL;
    if(tmp->ID == ID) return (FileProfile*)(tmp->data);
    while(tmp->next) {
        tmp = tmp->next;
        if(tmp->ID == ID) return (FileProfile*)(tmp->data);
    }
    return NULL;
}

// Free a PointSet
void _freePointSet(PointSet *tmp)
{
    if(tmp->locx) free(tmp->locx);
    if(tmp->locy) free(tmp->locy);
    if(tmp->locz) free(tmp->locz);
    if(tmp->prec) free(tmp->prec);
    //add Hui  2016.08
    if(tmp->start2) free(tmp->start2);
}

// Free a FileProfile (EXCLUDING PointSet it points to)
void _freeFileProfile(FileProfile *tmp)
{
    if(tmp->filename) free(tmp->filename);
    if(tmp->filemode) free(tmp->filemode);
    if(tmp->specificVars) free(tmp->specificVars);
    if(tmp->buffer) free(tmp->buffer);
}

// Free all the memory space (linked lists)
void _freeMemory()
{
    List *tmp, *tmpNext;
    tmp = pointSets;
    //printf("Free PointSets\n");
    while(tmp) {
        //printf("ID = %d\n",tmp->ID);
        if(tmp->data) {
            _freePointSet(tmp->data);
            free(tmp->data);
        }
        tmpNext = tmp->next;
        free(tmp);
        tmp = tmpNext;
    }
    pointSets = NULL;

    tmp = fileProfiles;
    //printf("Free FileProfiles\n");
    while(tmp) {
        //printf("ID = %d",tmp->ID);
        if(tmp->data) {
            _freeFileProfile(tmp->data);
            free(tmp->data);
        }
        tmpNext = tmp->next;
        free(tmp);
        tmp = tmpNext;
    }
    fileProfiles = NULL;
}

/* createRegularGrid sets a regular grid structure for use in read/write
 * This new PointSet is added to the beginning of the PointSet-LinkedList
   nbg[x,y,z]  : starting index of the grid in x,y,z directions (1-based)
   ned[x,y,z]  : end index of the grid in x,y,z directions
   nskp[x,y,z] : number of points to jump for consecutive points (1 means all points)
   gridID      : the index of the PointSet (0-based, kept in the LinkedList)
 */
void seism_createRegularGrid(int *nbgx, int *nedx, int *nskpx,
                 int *nbgy, int *nedy, int *nskpy, int *nbgz, int *nedz, int *nskpz,
                 int *gridID, int *err)
{
    int i,j,k;
    int cntx=0, cnty=0, cntz=0;
    //if(DEBUG>0 && mpi_rank%400==0) printf("%d) SEISM create regular grid\n",mpi_rank);
    if(!pointSets) nPointSets = 0;
    PointSet *newPS = (PointSet*)malloc(sizeof(PointSet));
    newPS->isGrid = 1;
    newPS->nx = (*nedx-*nbgx)/(*nskpx)+1;
    newPS->ny = (*nedy-*nbgy)/(*nskpy)+1;
    newPS->nz = (*nedz-*nbgz)/(*nskpz)+1;
    newPS->locx = (int*)malloc(sizeof(int)*(newPS->nx));
    newPS->locy = (int*)malloc(sizeof(int)*(newPS->ny));
    newPS->locz = (int*)malloc(sizeof(int)*(newPS->nz));
    // nbgx is 1-based, but locx is 0-based in C interface
    for(i=(*nbgx-1)+ghostx; i<*nedx+ghostx; i+=*nskpx)
        newPS->locx[cntx++] = i;
    for(j=(*nbgy-1)+ghosty; j<*nedy+ghosty; j+=*nskpy)
        newPS->locy[cnty++] = j;
    for(k=(*nbgz-1)+ghostz; k<*nedz+ghostz; k+=*nskpz)
        newPS->locz[cntz++] = k;
    newPS->recproc = -99;
    newPS->nprec = -1;
    newPS->prec = NULL;
    newPS->start2dim = -1;
    newPS->start2 = NULL;
    List *newItem = (List*)malloc(sizeof(List));
    newItem->ID = nPointSets;
    newItem->data = (void*)newPS;
    newItem->prev = NULL;
    if(!pointSets) newItem->next = NULL;
    else newItem->next = pointSets;
    pointSets = newItem;
    *gridID = nPointSets;
    nPointSets++;
    *err = 0;
    //if(DEBUG>0 && mpi_rank%400==0)
    //  printf("%d) SEISM regular grid is created with ID=%d\n",mpi_rank,*gridID);

    return;
}

/* creates a file profile and adds it to the beginning of the LinkedList
   fname   : file name of the file
   fmode   : file mode = "w" "r" "a" as write-read-append
   write_step : write_step, the amount of aggregation in time
   psID    : PointSet ID (input)
   fp      : FileProfile (output)
   err     : error code pointer. 0 is no error.
 */
void _createFileProfile(char *fname, char *fmode, int *write_step,
            char *fdata, int *psID, int *fpID, FileProfile **fp, int *err)
{
    *err = 0;

    //printf("creating new file profile\n");
    PointSet *ps = _getPointSet(*psID);
    if(!ps) {
        printf("ERROR! SEISM-IO: Invalid grid ID\n");
        *err = -1;
        return;
    }

    //check data type
    if(!strcmp(fdata,"real") && !strcmp(fdata,"float")
       && !strcmp(fdata,"int") && !strcmp(fdata,"double"))
    {
        printf("ERROR! SEISM_IO: Unrecognized data type\n");
        *err = -2;
        return;
    }

    //printf("point set is found\n");
    if(!fileProfiles) nFileProfiles = 0;
    FileProfile *newFP = (FileProfile*)malloc(sizeof(FileProfile));
    //printf("file profile is allocated\n");
    newFP->write_step = *write_step;
    newFP->filename = (char*)malloc(sizeof(char)*(strlen(fname)+1));
    strcpy(newFP->filename, fname);
    //printf("filename is copied\n");
    newFP->filemode = (char*)malloc(sizeof(char)*(strlen(fmode)+1));
    strcpy(newFP->filemode, fmode);
    newFP->pointSet = ps;

    if(!strcmp(fdata,"real") || !strcmp(fdata,"float"))
        newFP->datatype = MPI_FLOAT;
    else if(!strcmp(fdata,"int"))
        newFP->datatype = MPI_INT;
    else if(!strcmp(fdata,"double"))
        newFP->datatype = MPI_DOUBLE;
    //else if(!strcmp(fdata,"byte"))
    //  newFP->datatype = MPI_BYTE;

    newFP->specificVars = NULL;
    newFP->buffer = NULL;
    List *newItem = (List*)malloc(sizeof(List));
    //printf("List is allocated\n");
    newItem->ID = nFileProfiles;
    newItem->data = (void*)newFP;
    //printf("List->data=fp is set\n");
    newItem->prev = NULL;
    if(!fileProfiles) newItem->next = NULL;
    else newItem->next = fileProfiles;
    fileProfiles = newItem;
    *fp = newFP;
    *fpID = nFileProfiles;
    nFileProfiles++;

    //printf("err=%d\n",*err);
    return;
}

/* init function to initialize the SEISM lib. Reads in text input file
   and creates the initial group and variable objects.
   comm    : MPI_Comm pointer for initialization. Needs to be MPI_COMM_WORLD for ADIOS 1.3, 1.4 or readMesh()
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
void seism_init_new(MPI_Comm *tcomm, int *tmaxdim,
        int *tnx, int *tny, int *tnz,
        int *tghostx, int *tghosty, int *tghostz,
        int *tpx, int *tpy, int *tpz,
        char *t_seism_method, int *err)
{
   int trank;
   int tcoords[2];
   err       = MPI_Cart_coords(*comm, trank, 2, tcoords);
   int tnxt = (int)(*tnx) / (int)(*tpx);
   int tnyt = (int)(*tny) / (int)(*tpy);
   int tnzt = (int)(*tnz) / (int)(*tpz);

   seism_init(&tcomm, &trank, tcoords, (int*)&tmaxdim,
        &tnx, &tny, &tnz,
        &tnxt, &tnyt, &tnzt,
        &tghostx, &tghosty, &tghostz,
        &tpx, &tpy, &tpz,
        t_seism_method, &err);
}



void seism_init(MPI_Comm *tcomm, int *rank, int *tcoords, int *tmaxdim,
        int *tnx, int *tny, int *tnz,
        int *tnxt, int *tnyt, int *tnzt,
        int *tghostx, int *tghosty, int *tghostz,
        int *tpx, int *tpy, int *tpz,
        char *t_seism_method, int *err)
{

    int i;
    *err = 0;

    // if the fortran interface is not called, set app language to C
    if(SEISM_APP_LANG == SEISM_APP_NOT_SET)
        SEISM_APP_LANG = SEISM_APP_C;

    //printf("transport method=%s\tval=%d\n",t_seism_method,strncmp(t_seism_method,"mpi-io",6));
    if(!strncmp(t_seism_method,"mpiio",5)
       || !strncmp(t_seism_method,"mpi-io",6)
       || !strncmp(t_seism_method,"MPIIO",5)
       || !strncmp(t_seism_method,"MPI-IO",6))
        seism_method = SEISM_MPIIO;
    else if(!strncmp(t_seism_method,"adios",5)
        || !strncmp(t_seism_method,"ADIOS",5))
        seism_method = SEISM_ADIOS;
    else if(!strncmp(t_seism_method,"hdf5",4)
        || !strncmp(t_seism_method,"HDF5",4)
        || !strncmp(t_seism_method,"phdf5",5)
        || !strncmp(t_seism_method,"PHDF5",5))
        seism_method = SEISM_PHDF5;
    else if(!strncmp(t_seism_method,"pnetcdf",7)
        || !strncmp(t_seism_method,"PNETCDF",7)
        || !strncmp(t_seism_method,"PnetCDF",7))
        seism_method = SEISM_PNETC;
    else{
        printf("ERROR! SEISM_IO Init failed: Unrecognized transport method\n");
        *err = -1;
        return;
    }

    comm = tcomm;
    mpi_rank = *rank;
    nx = *tnx;
    ny = *tny;
    nz = *tnz;
    nxt = *tnxt;
    nyt = *tnyt;
    nzt = *tnzt;
    ghostx = *tghostx;
    ghosty = *tghosty;
    ghostz = *tghostz;
    px = *tpx;
    py = *tpy;
    pz = *tpz;
    maxdim = *tmaxdim;
    coords = (int*)malloc(sizeof(int)*maxdim);
    for(i=0; i<maxdim; i++) {
        coords[i] = tcoords[i];
    }
    xbounds[0] = nxt*coords[0]+ghostx;
    xbounds[1] = nxt*(coords[0]+1)-1+ghostx;
    ybounds[0] = nyt*coords[1]+ghosty;
    ybounds[1] = nyt*(coords[1]+1)-1+ghosty;
    int tmp = (maxdim>2 ? coords[2] : 0);
    zbounds[0] = nzt*tmp+ghostz;
    zbounds[1] = nzt*(tmp+1)-1+ghostz;
    int abort = 0;
    switch(seism_method) {
    case SEISM_ADIOS:
#ifdef S_ADIOS
        //seism_init_adios(tcomm, rank, tcoords, err);
#else
        printf("ERROR! ADIOS is not defined in compilation!\n");
        abort = 1;
#endif
        break;
    case SEISM_PHDF5:
#ifdef S_PHDF5
        //seism_init_phdf5(*comm);
#else
        printf("ERROR! PHDF5 is not defined in compilation!\n");
        abort = 2;
#endif
        break;
    case SEISM_MPIIO:
        // Do not do anything for now. We will initialize after the
        // write_style is defined
        break;
    case SEISM_PNETC:
#ifdef S_PNETCDF
        // seism_init_pnetc();
#else
        printf("ERROR! PNETCDF is not defined in compilation!\n");
        abort = 3;
#endif
        break;
    default:
        printf("ERROR! SEISM_IO method=%d is not known!\n",seism_method);
        abort = 4;
        break;
    }
    *err = abort;

    return;
}

/* file open function fortran interface
   fname   : file name of the file
   fmode   : file mode = "w" "r" "a" as write-read-append
   write_step : write_step, the amount of aggregation in time
   fdata   : file datatype ('real', 'int', 'double', 'byte', etc.)
   psID    : PointSet ID (input)
   fpID    : FileProfile ID (output)
   err     : error code pointer. 0 is no error.
 */
void seism_file_open(char *fname, char *fmode, int *write_step, char *fdata,
             int *psID, int *fpID, int *err)
{

    FileProfile *fp;
    PointSet *ps;
    *err = 0;
    _createFileProfile(fname, fmode, write_step, fdata, psID, fpID, &fp, err);
    if(*err !=0) {
        printf("ERROR! SEISM-IO: #=%d Cannot open file - error in creating File Profile\n",*err);
        return;
    }
    if(strcmp(fmode, "r") == 0)
    {
        if(access(fname, F_OK) != 0)
        {
            *err = -100;
            if(mpi_rank == 0) printf("ERROR! SEISM-IO: #=%d %s not exist.\n",*err, fname);
            return;
        }
    }
    if(strcmp(fmode, "w") == 0)
    {
        if(access(fname, F_OK) == 0)
        {
            if(mpi_rank ==0 && remove(fname) != 0)
            {
                *err = -101;
                printf("ERROR! SEISM-IO: #=%d %s exist.\n",*err, fname);
                return;
            }
        }
    }

    //printf("FileProfile is created\n");
#ifdef S_PHDF5
    hid_t dset_id;
    hid_t datatype;
#endif

    switch(seism_method) {
    case SEISM_ADIOS:
#ifdef S_ADIOS
        ps = _getPointSet(*psID);
        if(!ps) {
            printf("SEISM ERROR! Invalid grid ID\n");
            *err = -1;
            return;
        }
        seism_initials_adios(fp);
        seism_file_open_adios(fp, err);
#endif
        break;
    case SEISM_PHDF5:
#ifdef S_PHDF5
        ps = _getPointSet(*psID);
        if(!ps) {
            printf("SEISM ERROR! Invalid grid ID\n");
            *err = -1;
            return;
        }
        seism_initials_phdf5(ps);
        //seism_file_open_phdf5(*seism_f, cgname, cfname, comm, &dset_id, &datatype); //org
        seism_file_open_phdf5(fp, err);
#endif
        break;
    case SEISM_PNETC:
#ifdef S_PNETCDF
        ps = _getPointSet(*psID);
        if(!ps) {
            printf("SEISM ERROR! Invalid grid ID\n");
            *err = -1;
            return;
        }
        seism_initials_pnetcdf(ps);
        seism_file_open_pnetcdf(fp, err);
#endif
        break;
    case SEISM_MPIIO:
        ps = _getPointSet(*psID);
        if(!ps) {
            printf("SEISM ERROR! Invalid grid ID\n");
            *err = -1;
            return;
        }
        seism_initials_mpiio(ps);
        //printf("mpiio is initialized\n");
        seism_file_open_mpiio(fp, err);
        break;
    }
    //if(DEBUG>0 && mpi_rank%400==0)  printf("%d) File opened\n",mpi_rank);
    return;
}

/* reads in the binary media file
 * Since the file is binary, we only use MPI-IO.
   fname : file name
   nvar  : number of variables to read
   mediaData:  4D matrix of all the variables
   order:  if >0, first layer is top, if <=0, the first layer is bottom
   err   : error code, 0 means no error.
 */
void seism_readMesh(char *fname, int *nvar, void *mediaData, int *order, int *err)
{

    if(DEBUG>0 && mpi_rank%400==0) {
        printf("%d) Reading in mesh\n",mpi_rank);
    }
    *err = 0;
    seism_readMesh_binary(fname, nvar, mediaData, order, err);
    if(DEBUG>0 && mpi_rank%400==0) printf("%d) Mesh is read\n",mpi_rank);

    return;
}

/* write function fortran interface
   seism_f : FileProfile ID integer. File has to be opened before.
   var   : pointer to multi-dimensional variable
   err     : error code, integer pointer. 0 is no error.
 */
void seism_write(int *seism_f, void *var, int *err)
{

    if(DEBUG>0 && mpi_rank%400==0) {
        //printf("%d) Writing to file with ID=%d\n",mpi_rank,*seism_f);
    }
    *err = 0;
    FileProfile *fp = _getFileProfile(*seism_f);
    if(!fp) {
        printf("ERROR! SEISM-IO: File is not opened before write.\n");
        *err = -1;
        return;
    }
/*
    void *data = convertCudaPointer(fp, var, err);
    if(*err) {
        printf("ERROR! SEISM-IO: call convertCudaPointer error.\n");
        return;
    }
*/
    switch(seism_method) {
    case SEISM_ADIOS:
#ifdef S_ADIOS
        //seism_write_adios(seism_f, cvname, err);
        seism_write_adios(fp, var, err);
#endif
        break;
    case SEISM_PHDF5:
#ifdef S_PHDF5
        //seism_write_phdf5(cvname, err); //org
        seism_write_phdf5(fp, var, err);
#endif
        break;
    case SEISM_PNETC:
#ifdef S_PNETCDF
        seism_write_pnetcdf(fp, var, err);
#endif
        break;
    case SEISM_MPIIO:
        seism_write_mpiio(fp, var, err);
        break;
    }

    //freeCudaPointer(data);
    return;
}

/* Copy from above seism_write(). add by Hui 2016.08
   Read function C interface
   seism_f : FileProfile ID integer. File has to be opened before.
   var   : pointer to multi-dimensional variable
   err     : error code, integer pointer. 0 is no error.
 */
void seism_read(int *seism_f, void *var, int *err)
{

    *err = 0;
    FileProfile *fp = _getFileProfile(*seism_f);
    if(!fp) {
        printf("ERROR! SEISM-IO: File is not opened before read.\n");
        *err = -1;
        return;
    }
    switch(seism_method) {
    case SEISM_ADIOS:
#ifdef S_ADIOS
        //seism_read_adios(seism_f, cvname, err);
        seism_read_adios(fp, var, err);
#endif
        break;
    case SEISM_PNETC:
#ifdef S_PNETCDF
        seism_read_pnetcdf(fp, var, err);
#endif
        break;
    case SEISM_PHDF5:
#ifdef S_PHDF5
        //seism_read_phdf5(cvname, err); //org
        seism_read_phdf5(fp, var, err); //这里可能少一个cvname
#endif
        break;
    case SEISM_MPIIO:
        seism_read_mpiio(fp, var, err);
        break;
    }
    return;
}

/* seism_file_close closes an already opened file
   seism_f : pointer to file
   err     : error code
 */
void seism_file_close(int *seism_f, int *err)
{

    FileProfile *fp = _getFileProfile(*seism_f);

    switch(seism_method) {
    case SEISM_ADIOS:
#ifdef S_ADIOS
        //*err = adios_close(*seism_f);
        seism_file_close_adios(fp, err);
#endif
        break;
    case SEISM_PHDF5:
#ifdef S_PHDF5
        seism_file_close_phdf5(fp, err);
#endif
        break;
    case SEISM_PNETC:
#ifdef S_PNETCDF
        seism_file_close_pnetcdf(fp, err);
#endif
        break;
    case SEISM_MPIIO:
        seism_file_close_mpiio(fp, err);
        break;
    }
    return;
}

/* seism_finalize finalizes the library
   err : error code
 */
void seism_finalize(int *err)
{
    //printf("%d) Finalizing SEISM-IO...\n",mpi_rank);
    *err = 0;
    _freeMemory();
    switch(seism_method) {
    case SEISM_ADIOS:
#ifdef S_ADIOS
        //*err = adios_finalize(mpi_rank);
#endif
        break;
    case SEISM_MPIIO:
        // nothing here so far
        break;
    }
    //if(DEBUG>0 && mpi_rank%400==0) printf("%d) SEISM-IO is finalized\n",mpi_rank);
    return;
}
#endif

/* test function for receiver calculation only
 */
#ifdef _CALCREC_ONLY
int main(){
    int nxt=100;
    int nyt=120;
    int nzt=172;
    int nbgx=1;
    int nedx=522;
    int nskpx=7;
    int nbgy=133;
    int nedy=180;
    int nskpy=3;
    int nbgz=175;
    int nedz=175;
    int nskpz=2;
    int npbgx,npbgy,npbgz,npedx,npedy,npedz;
    int nprecx,nprecy,nprecz;
    int coord[]={2,1,1};

    seism_calc_rec_(&nxt,&nyt,&nzt,&nbgx,&nedx,&nskpx,
            &nbgy,&nedy,&nskpy,&nbgz,&nedz,&nskpz,
            &npbgx,&npedx,&npbgy,&npedy,&npbgz,&npedz,
            &nprecx,&nprecy,&nprecz,coord);

    printf("np: (%d,%d)\t(%d,%d)\t(%d,%d)\n",npbgx,npedx,npbgy,npedy,npbgz,npedz);
    printf("nprec: %d,%d,%d\n",nprecx,nprecy,nprecz);
    return 0;
}
#endif

/* test function for init, variable definitions, file open/close etc.
 */
#ifdef _EXEC_SEISM_IO
int main()
{

    size_t num;
    int err;
    int v1 = 100;
    int v2 = 1000;

    return 0;
}
#endif
