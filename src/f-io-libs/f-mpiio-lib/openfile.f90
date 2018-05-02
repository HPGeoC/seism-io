!**************************
!@section LICENSE
!Copyright (c) 2013-2017, Regents of the University of California
!All rights reserved.
!Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
!1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
!2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
!THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!**************************
!!for each file these names must be unique:

!!FORTRAN
!!integer icount1,icount2,...
!!integer nprec1, nprec2...
!!integer, dimension(:,:),allocatable:: prec1,prec2,...
!!real, dimension (:),allocatable::buffer1,buffer2,...

!! need these in main code

!! 1) set nprec,recproc by calling initials
!! 2) allocate(prec(nprec,maxdim))
!! 3) allocate(buffer(nprec*write_step))
!! 4) declare integer(kind=MPI_OFFSET_KIND)::disp
!! 5) declare integer filetype, fh, icount etc

SUBROUTINE openfile(rank,recproc,comm,filenameC, write_style,rw, &
     ghostx,ghosty,ghostz, &
     nk,kout,nj,jout,ni,iout, maxdim, &
     kbnd,jbnd,ibnd,&
     datatype,write_step,nprec, prec, &
     nrec,disp,filetype,icount,fh,filenameLen) bind(C)
  use iso_c_binding, only: c_char, c_int
  IMPLICIT NONE
  INCLUDE 'mpif.h'

  INTEGER nk,nj,ni,comm,fh,recproc,nptotal,kbnd(2),jbnd(2),ibnd(2)
  INTEGER ghostx,ghosty,ghostz
  INTEGER kout(nk),jout(nj),iout(ni), prec(nprec,maxdim)
  INTEGER err,rank,maxdim,write_step,nprec,icount,filetype,read_write
  INTEGER(kind=MPI_OFFSET_KIND)::disp,nrec

  INTEGER datatype,comm2,onegroup,write_style,nbyte
  INTEGER(kind = c_int)                :: filenameLen
  CHARACTER(kind = c_char), intent(in) :: filenameC(*)
  CHARACTER(len = filenameLen)         :: filename
  INTEGER i
  CHARACTER*1 rw

  ! copy c-filename to fortran
  FORALL(i = 1:filenameLen) filename(i:i) = filenameC(i)

  disp=0

  onegroup=0
  IF (rank == recproc) onegroup = 1
  CALL MPI_COMM_SPLIT(comm,onegroup,0,comm2,err)

  IF(TRIM(rw)=='w')&
       CALL MPI_FILE_OPEN(comm2,TRIM(filename),MPI_MODE_WRONLY+MPI_MODE_CREATE,MPI_INFO_NULL,fh, err)

  IF(TRIM(rw)=='r') &
       CALL MPI_FILE_OPEN(comm2,TRIM(filename),MPI_MODE_RDONLY, MPI_INFO_NULL,fh,err)

  CALL MPI_BARRIER(comm,err)
  IF (rank .NE. recproc) RETURN

  !block region output
  read_write=write_step

  !need override
  IF(TRIM(rw)=='r') read_write=1

  CALL set_filetype_style(ghostx,ghosty,ghostz,write_style,nk,kout,nj,jout,ni,iout,rank, maxdim, &
       kbnd,jbnd,ibnd, datatype,read_write,&
       nprec, prec, nrec,filetype,recproc)

  if(rank==0) then
  !write(*,*)'After setup filetype:', rank,recproc,disp,datatype,nrec
  endif

  icount=0

  IF(TRIM(rw)=='r') THEN
     CALL MPI_TYPE_SIZE(datatype,nbyte,err)
     !nrec8=nrec
     disp=disp+(write_step-1)*nrec*nbyte
  ENDIF

  !write(*,*)' return from openfile'
  RETURN
END SUBROUTINE openfile
