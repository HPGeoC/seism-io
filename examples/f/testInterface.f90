!**************************
!@section LICENSE
!Copyright (c) 2013-2017, Regents of the University of California
!All rights reserved.
!Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
!1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
!2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
!THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!**************************
program main
    use mpi
    implicit none

    integer :: Comm,MC1
    integer :: rank,rankmc1, size, err
    INTEGER :: n,nx,ny,nz,nxt,nyt,nzt
    INTEGER :: ghostx,ghosty,ghostz,npx,npy,npz
    INTEGER :: nbgx,nedx,nskpx,nbgy,nedy,nskpy,nbgz,nedz,nskpz
    INTEGER :: nbgx2,nedx2,nskpx2,nbgy2,nedy2,nskpy2,nbgz2,nedz2,nskpz2
    INTEGER :: seism_regGridID, seism_regGridID2, seism_regGridID3
    INTEGER :: seism_filex, seism_filex2, seism_filex3
    integer :: write_step

    integer,parameter ::maxdim1=3
    integer :: maxdim
    integer :: dims(maxdim1),coords(maxdim1)
    logical :: periodic(maxdim1),reorder

    CHARACTER (LEN=5) :: seism_method
    CHARACTER (LEN=1) :: fmode
    CHARACTER (LEN=5) :: fdata
    character (len=18) :: sxgro,sxgro2
    character (len=26) :: sxgro3

    REAL, dimension (:,:,:), allocatable :: u1
    integer ii,jj,kk,sum

    include "seism_io_consts.fh"

!   =========
!   START MPI
!   =========
    call MPI_INIT(err)
    call MPI_COMM_DUP(MPI_COMM_WORLD,Comm,err)
    call MPI_COMM_RANK(Comm,rank,err)
    call MPI_COMM_SIZE(Comm,size,err)
!!! 24000 Cores
!    nx=10800
!    ny=1600
!    nz=2048
!    nxt=180
!    nyt=64
!    nzt=128
!    npx=60
!    npy=25
!    npz=16

!!! 2048 Cores
!    nx=2880
!    ny=1024
!    nz=1024
!    nxt=180
!    nyt=64
!    nzt=128
!    npx=16
!    npy=16
!    npz=8

!!! 16 Cores
!    nx= 720
!    ny= 128
!    nz= 256
!    nxt=180
!    nyt=64
!    nzt=128
!    npx=4
!    npy=2
!    npz=2

!!! 8 Cores
    nx= 360
    ny= 128
    nz= 256
    nxt=180
    nyt=64
    nzt=128
    npx=2
    npy=2
    npz=2

!!! 4096 Cores
!    nx=2880
!    ny=1024
!    nz=2048
!    nxt=180
!    nyt=64
!    nzt=128
!    npx=16
!    npy=16
!    npz=16

!!! 1 Cores
!    nx= 180
!    ny= 64
!    nz= 128
!    nxt=180
!    nyt=64
!    nzt=128
!    npx=1
!    npy=1
!    npz=1

    IF (npx /= nx/nxt .OR. npy /= ny/nyt .OR.npz /= nz/nzt) THEN
        if(rank == 0) print *, "Size parameter nx,ny,nz,nxt,nyt,nzt,npx,npy,npz ERROR."
        CALL MPI_FINALIZE(err)
        STOP
    END IF

    npx = nx/nxt
    npy = ny/nyt
    npz = nz/nzt
    coords    = 0
    periodic  = .false.
    reorder   = .true.
    dims(0+1) = npx
    dims(1+1) = npy
    dims(2+1) = npz
    maxdim=3

    ghostx = 2
    ghosty = 2
    ghostz = 2

    seism_method = "mpiio"
    fmode = "w"
    fdata = "float"
    write_step = 1

    IF (size /= npx*npz*npy) THEN
        if(rank == 0) print *, "Size communicator (",size,") NOT equal Cartesian topology (",npx*npy*npz,")"
        CALL MPI_FINALIZE(err)
        STOP
    END IF

    !call MPI_CART_CREATE(MCW,maxdim,dims,periodic,reorder,MC1,err)
    !call MPI_CART_GET(MC1,maxdim,dims,periodic,coords,err)

    call MPI_CART_CREATE(Comm,maxdim,dims,periodic,reorder,MC1,err)
    if(mod(rank,200) == 0) write (6,*) 'err1 = ', err
    !call MPI_CART_COORDS(MC1,rank,maxdim,coords,err)
    call MPI_CART_GET(MC1, maxdim, dims, periodic, coords,err);
    if(mod(rank,200) == 0) write (6,*) "rank = ", rank, ", coords = ", coords
    call MPI_CART_RANK(MC1, coords, rankmc1,err)
    if(mod(rank,200) == 0) write (6,*) 'err3 = ', err

!   SEISM_IO Initialization
!   ==============================================
    if(mod(rank,200) == 0) print *, rank, ') initialize SEISM-IO'

    call seism_init(MC1,rankmc1,coords,maxdim,nx,ny,nz,nxt,nyt,nzt,2,2,2,npx,npy,npz, seism_method,err)

    if(err .ne. 0) then
        if(rank == 0) print *, 'SEISM ERROR! Init failed!'
        call MPI_ABORT(Comm,911,err)
        call MPI_FINALIZE(err)
    endif
    if(mod(rank,200) == 0) print *, 'SEISM: ', rank, ') after init'


    nbgx =1
    nedx =nx
    nskpx=1
    nbgy =1
    nedy =ny
    nskpy=1
    nbgz =1
    nedz =1
    nskpz=1
    !nbgz= nz-nbgz+1 ! nbgz and nedz are modified since NZ is the surface
    !nedz= nz-nedz+1
    call seism_createRegularGrid(nbgx, nedx, nskpx,nbgy, nedy, nskpy, nedz, nbgz, nskpz,seism_regGridID, err) !Oginal
    !call seism_createRegularGrid(nbgx, nedx, nskpx,nbgy, nedy, nskpy, nbgz, nedz, nskpz,seism_regGridID, err)


    nbgx2  =1
    nedx2  =nx
    nskpx2 =1
    nbgy2  =1
    nedy2  =ny
    nskpy2 =1
    nbgz2  =1
    nedz2  =nz
    nskpz2 =1
    !nedz= nz-nedz+1 ! nbgz and nedz are modified since NZ is the surface
    !nbgz= nz-nbgz+1
    call seism_createRegularGrid(nbgx2, nedx2, nskpx2, nbgy2, nedy2, nskpy2, nbgz2, nedz2, nskpz2,seism_regGridID2, err)
    call seism_createRegularGrid(nbgx2, nedx2, nskpx2, nbgy2, nedy2, nskpy2, nbgz2, nedz2, nskpz2,seism_regGridID3, err)

    sxgro = "output_sfc/CSX96PS"
    seism_filex = -1111
    call seism_file_open(sxgro, fmode, write_step, fdata, seism_regGridID, seism_filex, err)

    sxgro2 = "output_vlm/CSX96PV"
    seism_filex2 = -1111
    call seism_file_open(sxgro2, fmode, write_step, fdata, seism_regGridID2, seism_filex2, err)

    sxgro3 = "output_vlm/rewrite_CSX96PV"
    fmode = "w"
    call seism_file_open(sxgro3, fmode, write_step, fdata, seism_regGridID3, seism_filex3, err)

    allocate(u1(-1:nxt+2,-1:nyt+2,-1:nzt+2))

!   init the result matrix
    u1 = 0
    do kk =-1+ghostz, nzt+2-ghostz
        do jj =-1+ghosty, nyt+2-ghosty
            do ii=-1+ghostx, nxt+2-ghostx
                u1(ii,jj,kk) = -1*kk
            end do
        end do
    end do
    sum = 0
    do kk =-1+ghostz, nzt+2-ghostz
        do jj =-1+ghosty, nyt+2-ghosty
            do ii=-1+ghostx, nxt+2-ghostx
                u1(ii,jj,kk) = rank !*(20) !+sum
                !sum= sum + 1
            end do
        end do
    end do

!    IF(rankmc1==0)THEN
!    print *, 'Rank:',rankmc1,' ixjxk=',ii,'x',jj,'x',kk
!    DO kk=-1+ghostz, nzt+2-ghostz
!          print *, 'k:',kk,' level'
!          do jj=-1+ghosty, nyt+2-ghosty
!              !print *, 'j:',jj
!              write(*,'(I3)',advance='no'),jj
!              do ii=-1+ghostx, nxt+2-ghostx
!                  write(*,'(1x,I3)',advance='no') u1(ii,jj,kk)
!              enddo
!              print *,' '
!         enddo
!         print *,'================='
!    ENDDO
!    ENDIF

    call seism_write(seism_filex, u1, err)
    call MPI_BARRIER(Comm,err)
    if(mod(rank,400)== 0) write (6,*) rank,', Write 1 DONE!'
    call seism_write(seism_filex2, u1, err)
    call MPI_BARRIER(Comm,err)
    if(mod(rank,400)== 0) write (6,*) rank,', Write 2 DONE!'

!   init the result matrix
    u1 = -99999
!   read from vlm file
    call seism_read(seism_filex2, u1, err)
    if(mod(rank,400)== 0) write (6,*) rank,', Read 2 DONE!'

!   write to new vlm file
    call seism_write(seism_filex3, u1, err)
    call MPI_BARRIER(Comm,err)
    if(mod(rank,400)== 0) write (6,*) rank,', Write 3 DONE!'

    call seism_file_close(seism_filex, err)
    call seism_file_close(seism_filex2, err)
    call seism_file_close(seism_filex3, err)

    call MPI_BARRIER(Comm,err)
    if(rank == 0) write (6,*) '| DONE! rank = ', rank
    call seism_finalize(err)
    call MPI_BARRIER(Comm,err)
    call MPI_FINALIZE(err)
    stop
end
