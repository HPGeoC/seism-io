!**************************
!@section LICENSE
!Copyright (c) 2013-2017, Regents of the University of California
!All rights reserved.
!Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
!1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
!2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
!THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!**************************
SUBROUTINE h_readfile_real(rank,recproc,&
     datatype,write_step,nprec,prec, &
     start,count,start2,start2dim,ntime,icount,dimsf,iopen,type_id,&
     plist_id,dataspace,dset_id,memspace,file_id,&
     buffer,var,nx,ny,nz) bind(C)
  USE hdf5
  IMPLICIT NONE
  INCLUDE 'mpif.h'
  INTEGER nprec,write_step,nx,ny,nz

  REAL buffer(nprec*1),buffer_x(nprec*1),var(nx,ny,nz)

  INCLUDE "h_readfile_common.f90"

  RETURN
END SUBROUTINE h_readfile_real

SUBROUTINE h_readfile_integer(rank,recproc,&
     datatype,write_step,nprec,prec, &
     start,count,start2,start2dim,ntime,icount,dimsf,iopen,type_id,&
     plist_id,dataspace,dset_id,memspace,file_id,&
     buffer,var,nx,ny,nz) bind(C)
  USE hdf5
  IMPLICIT NONE
  INCLUDE 'mpif.h'
  INTEGER nprec,write_step,nx,ny,nz

  INTEGER buffer(nprec*1),buffer_x(nprec*1),var(nx,ny,nz)

  INCLUDE "h_readfile_common.f90"

  RETURN
END SUBROUTINE h_readfile_integer

SUBROUTINE h_readfile_double(rank,recproc,&
     datatype,write_step,nprec,prec, &
     start,count,start2,start2dim,ntime,icount,dimsf,iopen,type_id,&
     plist_id,dataspace,dset_id,memspace,file_id,&
     buffer,var,nx,ny,nz) bind(C)
  USE hdf5
  IMPLICIT NONE
  INCLUDE 'mpif.h'
  INTEGER nprec,write_step,nx,ny,nz

  DOUBLE PRECISION buffer(nprec*1),buffer_x(nprec*1),var(nx,ny,nz)

  INCLUDE "h_readfile_common.f90"

  RETURN
END SUBROUTINE h_readfile_double
