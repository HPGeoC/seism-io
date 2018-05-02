!**************************
!@section LICENSE
!Copyright (c) 2013-2017, Regents of the University of California
!All rights reserved.
!Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
!1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
!2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
!THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!**************************
SUBROUTINE a_readfile_real(rank,recproc,&
     datatype,write_step,nprec,prec, &
     start,count,start2,start2dim,ntime,varnameC,icount,iopen,adios_handle,&
     buffer,var,nx,ny,nz,varnameLen) bind(C)
  use iso_c_binding, only: c_char, c_int
  USE  adios_write_mod
  IMPLICIT NONE
  INCLUDE 'mpif.h'

  REAL var(nx,ny,nz),buffer(nprec)
  REAL,DIMENSION(:),ALLOCATABLE::buffer_block

  INCLUDE "a_readfile_common.f90"


  RETURN
END SUBROUTINE a_readfile_real


SUBROUTINE a_readfile_integer(rank,recproc,&
     datatype,write_step,nprec,prec, &
     start,count,start2,start2dim,ntime,varnameC,icount,iopen,adios_handle,&
     buffer,var,nx,ny,nz,varnameLen) bind(C)
  use iso_c_binding, only: c_char, c_int
  USE  adios_write_mod
  IMPLICIT NONE
  INCLUDE 'mpif.h'

  INTEGER var(nx,ny,nz),buffer(nprec)
  INTEGER,DIMENSION(:),ALLOCATABLE::buffer_block

  INCLUDE "a_readfile_common.f90"

  RETURN
END SUBROUTINE a_readfile_integer

SUBROUTINE a_readfile_double(rank,recproc,&
     datatype,write_step,nprec,prec, &
     start,count,start2,start2dim,ntime,varnameC,icount,iopen,adios_handle,&
     buffer,var,nx,ny,nz,varnameLen) bind(C)
     use iso_c_binding, only: c_char, c_int
  USE  adios_write_mod
  IMPLICIT NONE
  INCLUDE 'mpif.h'

  DOUBLE PRECISION var(nx,ny,nz),buffer(nprec)
  DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::buffer_block

  INCLUDE "a_readfile_common.f90"

  RETURN
END SUBROUTINE a_readfile_double


SUBROUTINE a_read_closefile(adios_handle) bind(C)
  IMPLICIT NONE
  INCLUDE 'mpif.h'

  INTEGER*8 adios_handle
  INTEGER adios_err

  CALL adios_read_close (adios_handle, adios_err)

  RETURN
END SUBROUTINE a_read_closefile
