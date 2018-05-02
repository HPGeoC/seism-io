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

//define in f-mpiio-lib/readmesh.f90
extern void readmesh_(int *rank, MPI_Comm *comm,char* invel,int *coords,int *maxdim,
              int *nx,int *ny,int *nz,int *nxt,int *nyt,int *nzt,int *npx,int *npy,int *npz,
              int *nvar,int *artdeg,int *order,float *cube, int fnameLen);


void seism_readMesh_binary(char *fname, int *tnvar, void *mediaData, int *torder, int *err){

    int order = *torder; // if >0, first layer is top, if <=0, reverse order in Z
    int nvar = *tnvar;
    int *allCoords, *maxdims, *disps, size;
    int i, j, k, n;
    *err = 0;
    unsigned long long meshSize = nx*ny*nz*nvar;
    int submeshSize = nxt*nyt*nzt*nvar;
    float *buf;
    //if(meshSize < GB1){coords
    if(meshSize < 1) {
        // only master reads in serially and distributes according to coords
        if(!mpi_rank) {
            printf("Mesh size < 1GB, read in serially and distribute\n");
            size = px*py*pz;
            allCoords = (int*)malloc(sizeof(int)*maxdim*size);
            maxdims = (int*)malloc(sizeof(int)*size);
            disps = (int*)malloc(sizeof(int)*size);
            for(i=0; i<size; i++) {
                maxdims[i] = maxdim;
                disps[i] = i*maxdim;
            }
        }
        printf("\t%d) Before MPI_Gather for coords\n",mpi_rank);
        // everybody sends coords to the master
        MPI_Gatherv(coords, maxdim, MPI_INT, (void*)allCoords, maxdims, disps, MPI_INT, 0, *comm);
        printf("\t%d) coords are gathered\n",mpi_rank);
        free(maxdims);
        free(disps);
        // submeshSize is integer, because of the MPI_Isend count below!
        if(!mpi_rank) {
            FILE *f = fopen(fname, "rb");
            if(!f) {
                printf("SEISM ERROR! Media file cannot be opened: '%s'\n",fname);
                *err = -1;
                MPI_Abort(*comm, -1);
                MPI_Finalize();
                return;
            }
            buf = (float*)malloc(sizeof(float)*meshSize);
            size_t readBytes = fread(buf, sizeof(float), meshSize, f);
            if(readBytes != meshSize) {
                printf("SEISM ERROR! Error in reading media file serially: %s\n",fname);
                *err = -1;
                fclose(f);
                // WARNING! The application crashes here deliberately!
                MPI_Abort(*comm, -1);
                MPI_Finalize();
                return;
            }
            fclose(f);
            size_t blockSize = submeshSize*sizeof(float);
            size_t offset = blockSize;
            MPI_Request *req = (MPI_Request*)malloc(sizeof(MPI_Request)*(size-1));
            MPI_Status *stat = (MPI_Status*)malloc(sizeof(MPI_Status)*(size-1));
            int ipx, ipy, ipz;
            printf("\tSending data to every other rank\n");
            for(i=1; i<size; i++) {
                printf("\t\tIsend to %d\n",i);
                ipx = allCoords[i*maxdim];
                ipy = allCoords[i*maxdim+1];
                if(maxdim > 2) ipz = allCoords[i*maxdim+2];
                else ipz = 0;
                offset = nvar*nxt*ipx + nvar*nx*nyt*ipy;
                if(order>0)
                    offset += nvar*nx*ny*nzt*ipz;
                else
                    offset += nvar*nx*ny*nzt*(pz-ipz-1);
                MPI_Isend(&buf[offset], submeshSize, MPI_FLOAT, i, i, *comm, &req[i-1]);
            }
            MPI_Waitall(size-1, req, stat);
            free(req);
            free(stat);
        }
        else{
            MPI_Request req;
            MPI_Status stat;
            buf = (float*)malloc(sizeof(float)*submeshSize);
            printf("\t%d) Waiting for root to receive data\n",mpi_rank);
            MPI_Recv(buf, submeshSize, MPI_FLOAT, 0, mpi_rank, *comm, &stat);
            printf("\t%d) Mesh data is received\n",mpi_rank);
        }
    }
    else if(meshSize < GB10) {
        if(!mpi_rank) printf("Mesh size is medium/large, use 1-phase MPI-IO\n");
        buf = (float*)malloc(sizeof(float)*submeshSize);
        MPI_Datatype readtype;
        MPI_Status filestatus;
        MPI_File fh;
        int rmtype[3], rptype[3], roffset[3];
        rmtype[0]  = nz;
        rmtype[1]  = ny;
        rmtype[2]  = nx*nvar;
        rptype[0]  = nzt;
        rptype[1]  = nyt;
        rptype[2]  = nxt*nvar;
        if(maxdim>2) {
            if(order>0)
                roffset[0] = nzt*coords[2];
            else
                roffset[0] = nzt*(pz-coords[2]-1);
        }
        else
            roffset[0] = 0;
        roffset[1] = nyt*coords[1];
        roffset[2] = nxt*coords[0]*nvar;
        *err = MPI_Type_create_subarray(3, rmtype, rptype, roffset, MPI_ORDER_C, MPI_FLOAT, &readtype);
        *err = MPI_Type_commit(&readtype);
        *err = MPI_File_open(*comm,fname,MPI_MODE_RDONLY,MPI_INFO_NULL,&fh);
        *err = MPI_File_set_view(fh, 0, MPI_FLOAT, readtype, "native", MPI_INFO_NULL);
        *err = MPI_File_read_all(fh, buf, nvar*nxt*nyt*nzt, MPI_FLOAT, &filestatus);
        *err = MPI_File_close(&fh);
    }
    else{ // very large media file - read in using two-phase IO
        // partdeg must divide NY and PY and mpi_size
        // NX*NY*2*nvar/partdeg < MEMORY_SIZE
        // PY, NY > partdeg
        // mpi_size > NZ*partdeg
        buf = (float*)malloc(sizeof(float)*submeshSize);
        int partdeg;
        partdeg = (int)(meshSize/GB100)+1;
        if(!mpi_rank) printf("Mesh size is large, use MPI-IO with partdeg=%d\n",partdeg);

        // readmesh.f90 sets indices for FORTRAN + ORDER>0
        readmesh(&mpi_rank, comm, fname, coords, &maxdim,
                    &nx, &ny, &nz, &nxt, &nyt, &nzt, &px, &py, &pz, &nvar,
                    &partdeg, &order, buf, strlen(fname));
    }
    unsigned long long ind=0, _4dind = 0, tmp = 0;
    if(SEISM_APP_LANG == SEISM_APP_C) {
        for(i=0; i<nxt; i++) {
            tmp = nvar*i;
            for(j=0; j<nyt; j++) {
                tmp += nvar*nxt*j;
                _4dind = tmp;
                if(order>0) ind = _4dind;
                else ind = tmp + nvar*nxt*nyt*(nzt-1);
                for(k=0; k<nzt; k++) {
                    memcpy((void*)&((float*)mediaData)[_4dind], &buf[ind],
                           nvar*sizeof(float));
                    _4dind += nvar*nxt*nyt;
                    if(order>0) ind = _4dind;
                    else ind -= nvar*nxt*nyt;
                }
            }
        }
    }
    else{ // FORTRAN APP
        if(order>0)
            memcpy(mediaData, buf, submeshSize*sizeof(float));
        else{
            for(k=0; k<nzt; k++) {
                _4dind = k*nvar*nxt*nyt;
                ind = (nzt-k-1)*nvar*nxt*nyt;
                memcpy((void*)&((float*)mediaData)[_4dind], &buf[ind],
                       nvar*nxt*nyt*sizeof(float));
            }
        }
    }
    printf("%d) mediaData: %e, %e, %e\n",mpi_rank,
           ((float*)mediaData)[0],((float*)mediaData)[1],((float*)mediaData)[2]);

    free(buf);
    return;
}
