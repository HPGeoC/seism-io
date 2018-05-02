!**************************
!@section LICENSE
!Copyright (c) 2013-2017, Regents of the University of California
!All rights reserved.
!Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
!1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
!2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
!THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!**************************
!!!!!!! Copied from AWP-ODC-v1.1.2 structure.f
!!!!!!! on Sep 17, 2014

!! Returns data specific to this rank in the order same as
!!  the media file. If order>0, then the distribution is
!!  done in the file order. Else, the distribution is in
!!  reversed order
SUBROUTINE readmesh(rank,comm,invelC,coords,maxdim, &
     nx,ny,nz,nxt,nyt,nzt,npx,npy,npz,nvar,partdeg,order,cube,invelLen) bind(C)
  use iso_c_binding, only: c_char, c_int
  IMPLICIT NONE
  INCLUDE 'mpif.h'

  INTEGER,PARAMETER :: i8 = selected_int_KIND(10)
  INTEGER nx, ny, nz, nxt, nyt, nzt, npx, npy, npz, nvar, comm, order
  INTEGER(i8) :: nxt_l,nyt_l,nzt_l,ny_l,nx_l,nvar_l
  INTEGER rank, zrank, est_partdeg, zpos, partdeg
  INTEGER(KIND=MPI_ADDRESS_KIND) :: partdeg_l
  REAL partdeg_f
  INTEGER, DIMENSION(:), ALLOCATABLE :: requests
  INTEGER, DIMENSION(:,:), ALLOCATABLE :: mpistatus
  INTEGER COMM_IO, yindex
  INTEGER(i8) i, j, k, l
  INTEGER x, y, z, mpirank
  INTEGER(kind = c_int) :: invelLen
  CHARACTER(kind = c_char), intent(in) :: invelC(*)
  CHARACTER(len=invelLen) :: invel
  INTEGER merr, mfh
  INTEGER :: PLANEUNIT, CUBEPLANEUNIT
  INTEGER(MPI_OFFSET_KIND) :: disp
  INTEGER :: maxdim, COMM_RECV
  INTEGER, DIMENSION(maxdim) :: coords
  REAL, DIMENSION (:,:,:,:,:), ALLOCATABLE :: xyplane
  REAL :: cube(nvar,nxt,nyt,nzt)
  INTEGER(KIND=MPI_ADDRESS_KIND),PARAMETER :: WSIZE = 4
  INTEGER,DIMENSION(MPI_STATUS_SIZE) :: mystatus

  ! copy c-filename to fortran
  FORALL(i = 1:invelLen) invel(i:i) = invelC(i)

  ! if (rank.eq.0) print *, "initmesh - CASE 2, MPI - for big files"

  nx_l=nx
  ny_l=ny
  nvar_l=nvar

  ALLOCATE(requests(1:npx*npy+nzt))
  ALLOCATE(mpistatus(MPI_STATUS_SIZE,1:npx*npy+nzt))

  zrank=rank   ! 0..(nz-1)

  ! maximum 1GB for temporal mesh buffers
  partdeg_f=REAL(nx_l*ny_l*(nvar_l+nvar)*4)/REAL(1024*1024*1024)
  IF (rank.EQ.0) PRINT *,partdeg_f,"GB"
  IF (partdeg_f.GT.1) THEN
     est_partdeg=FLOOR(LOG(partdeg_f)/LOG(2.))+1
     est_partdeg=2**est_partdeg
  ELSE
     est_partdeg=1
  END IF

  ! verify the ny and npy are divisible by partdeg and npx*npy*npz
  ! >= nz*partdeg
  partdeg_l=partdeg
  IF (rank.EQ.0) PRINT *,"Suggested amount of partitioning", &
       est_partdeg
  IF (rank.EQ.0) PRINT *,"Actual amount of partitioning", &
       partdeg_l

  IF (zrank<nz*partdeg) THEN
     CALL MPI_COMM_SPLIT(comm,1,0,COMM_IO,merr)
     CALL MPI_FILE_OPEN(COMM_IO,invel,MPI_MODE_RDONLY, &
          MPI_INFO_NULL,mfh,merr)

     ! new data type containing a contiguous data chunks
     CALL MPI_Type_contiguous(nx_l*(ny_l/partdeg_l)*nvar_l, &
          MPI_REAL,PLANEUNIT,merr)
     CALL MPI_Type_commit(PLANEUNIT,merr)

     ALLOCATE(xyplane(nvar,1:nxt,1:nyt,1:npx,1:npy/partdeg))

     ! t1=mpi_wtime()
     disp=INT(zrank,MPI_ADDRESS_KIND)*nx_l*(ny_l/partdeg_l)*nvar_l*WSIZE
     ! inverted z-axis, an absolute offset in bytes
     CALL MPI_File_set_view(mfh,disp,MPI_REAL,PLANEUNIT,"native", &
          MPI_INFO_NULL,merr)
     disp=0
     IF (rank.EQ.0) PRINT *, "read_at_all", rank
     CALL MPI_File_read_at_all(mfh,disp,xyplane,nx*(ny/partdeg)*nvar, &
          MPI_REAL,mystatus,merr)

     IF (rank.EQ.0) PRINT *, "read_at_all passed", rank
     ! t2=mpi_wtime()

     !call MPI_BARRIER(COMM_IO, merr)

     !if (rank.eq.0) then ! master
     !  !write (*,'("Reading time =",f10.3," sec")') t2-t1
     !  print *,"Reading time =",t2-t1," sec"
     !end if

     CALL MPI_FILE_CLOSE(mfh,merr)
     ! allocate(cube(nvar,1:nxt,1:nyt,1:nzt)) ! allocated in the app

     yindex=MOD(zrank,partdeg)

     ! t3=mpi_wtime()
     !if (rank.eq.0) then ! master
     !  !write (*,'("Restructuring time =",f10.3," sec")') t3-t2
     !  print *,"Reconstructing time =",t3-t2," sec"
     !end if

     !call MPI_BARRIER(COMM_IO, merr)

     ! DISTRIBUTING.....
     zpos=(zrank/partdeg)/nzt
     IF (order .LE. 0) THEN
        ! changing coordinate system: regular X-Y-Z -> MPI
        ! Z-Y-X
        z=npz-zpos-1   ! flip z-axis of the cubes
     ELSE
        z = zpos
     ENDIF
     DO i=yindex*(npy/partdeg)+1,(yindex+1)*npy/partdeg
        DO j=1,npx
           k=(i-yindex*(npy/partdeg)-1)*npx+j

           mpirank=(j-1)*npy*npz+(i-1)*npz + z

           PRINT *,rank,') j=',j,' i=',i,' zpos=',zpos,' z=',z,&
                ' dest=',mpirank,' tag=',MOD(zrank/partdeg,nzt)+1,' k=',k

           !if (mod(zrank,nzt).eq.0) print
           !*,"S",rank,mpirank,i,j,z

           CALL MPI_ISEND(xyplane(1,1,1,j,i-yindex*(npy/partdeg)), &
                nxt*nyt*nvar,MPI_REAL,mpirank, &
                MOD(zrank/partdeg,nzt)+1,comm,requests(k),merr)
        END DO
     END DO

     ! RECEIVING..... ! do we need flipping z-axis?
     x=rank/(npy*npz)
     y=(rank-x*npy*npz)/npz
     z=MOD(rank-x*npy*npz,npz)
     yindex=y/(npy/partdeg)

     IF (order .LE. 0) THEN
        z=npz-z-1   ! flip z-axis of the cubes
     ENDIF

     !print *,"SR",rank,sourcerank,x,y,z

     DO i=1,nzt
        PRINT *,'IRECV: ',rank,') x=',x,' y=',y,' z=',z, &
             ' yind=',yindex, &
             ' src=',(z*nzt+i-1)*partdeg+yindex, &
             ' tag=',i,' k=',i+npx*npy/partdeg
        CALL MPI_IRECV(cube(1,1,1,i),nxt*nyt*nvar,MPI_REAL, &
             (z*nzt+i-1)*partdeg+yindex,i,comm, &
             requests(i+npx*npy/partdeg),merr)
     END DO
     CALL MPI_WAITALL(npx*npy/partdeg+nzt,requests,mpistatus,merr)

  ELSE ! if(zrank<nz*partdeg)
     CALL MPI_COMM_SPLIT(comm,2,0,COMM_RECV,merr)
     ! allocate(cube(nvar,1:nxt,1:nyt,1:nzt))

     ! RECEIVING..... ! do wee need flipping z-axis?
     x=rank/(npy*npz)
     y=(rank-x*npy*npz)/npz
     z=MOD(rank-x*npy*npz,npz)
     yindex=y/(npy/partdeg)

     IF (order .LE. 0) THEN
        z=npz-z-1   ! flip z-axis of the cubes
     ENDIF

     !print *,"R", rank,sourcerank,x,y,z

     DO i=1,nzt
        CALL MPI_IRECV(cube(1,1,1,i),nxt*nyt*nvar,MPI_REAL, &
             (z*nzt+i-1)*partdeg+yindex,i,comm,requests(i),merr)
     END DO
     CALL MPI_WAITALL(nzt,requests,mpistatus,merr)
  END IF

  CALL MPI_BARRIER(comm, merr)
  !t4=mpi_wtime()
  !if (rank.eq.0) then ! master
  !  !write (*,'("Reception time =",f10.3," sec")') t4-t3
  !  print *,"Reception time =",t4-t3," sec"
  !end if

  IF (zrank<nz*partdeg) THEN
     DEALLOCATE(xyplane)
  END IF

  ! For test purpose only
  !            do iwrite=0,IOST
  !              if (mod(rank,IOST+1).eq.iwrite) then
  !                write(mediafile, '(a,i7.7,a)')
  !                'input_rst/meshpart/media', rank,'.bin'
  !                open(9, file=mediafile, form='unformatted',
  !                access='direct',
  !     $
  !     status='replace',recl=bsize*(nxt)*(nyt)*(nzt)*act_nvar)
  !                write(9, rec=1) cube
  !                close(9)
  !              end if
  !              call MPI_BARRIER(comm,merr)
  !            end do

  !t5=mpi_wtime()
  !if (rank.eq.0) then ! master
  !  !write (*,'("Writng file time =",f10.3," sec")') t5-t4
  !  print *,"Writing file time =",t5-t4," sec"
  !end if

  DEALLOCATE(requests)
  DEALLOCATE(mpistatus)

  RETURN
END SUBROUTINE readmesh
