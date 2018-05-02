!**************************
!@section LICENSE
!Copyright (c) 2013-2017, Regents of the University of California
!All rights reserved.
!Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
!1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
!2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
!THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!**************************
SUBROUTINE a_writefile_real(rank,recproc,&
     datatype,write_step,nprec,prec, &
     start,count,start2,start2dim,ntime,varnameC,icount,iopen,adios_handle,&
     buffer,var,nx,ny,nz,varnameLen) bind(C)
  use iso_c_binding, only: c_char, c_int
  USE  adios_write_mod
  IMPLICIT NONE
  INCLUDE 'mpif.h'

  INTEGER nprec,write_step,nx,ny,nz
  REAL buffer(nprec*write_step),var(nx,ny,nz)

  INCLUDE "a_writefile_common.f90"

  RETURN
END SUBROUTINE a_writefile_real


SUBROUTINE a_writefile_integer(rank,recproc,&
     datatype,write_step,nprec,prec, &
     start,count,start2,start2dim,ntime,varnameC,icount,iopen,adios_handle,&
     buffer,var,nx,ny,nz,varnameLen) bind(C)
  use iso_c_binding, only: c_char, c_int
  USE  adios_write_mod
  IMPLICIT NONE
  INCLUDE 'mpif.h'

  INTEGER nprec,write_step,nx,ny,nz
  INTEGER buffer(nprec*write_step),var(nx,ny,nz)

  INCLUDE "a_writefile_common.f90"

  RETURN
END SUBROUTINE a_writefile_integer

SUBROUTINE a_writefile_double(rank,recproc,&
     datatype,write_step,nprec,prec, &
     start,count,start2,start2dim,ntime,varnameC,icount,iopen,adios_handle,&
     buffer,var,nx,ny,nz,varnameLen) bind(C)
  use iso_c_binding, only: c_char, c_int
  USE  adios_write_mod
  IMPLICIT NONE
  INCLUDE 'mpif.h'

  INTEGER nprec,write_step,nx,ny,nz
  DOUBLE PRECISION buffer(nprec*write_step),var(nx,ny,nz)

  INCLUDE "a_writefile_common.f90"

  RETURN
END SUBROUTINE a_writefile_double



SUBROUTINE a_closefile(adios_handle) bind(C)
  !USE HDF5 ! This module contains all necessary modules

  IMPLICIT NONE
  INCLUDE 'mpif.h'

  INTEGER*8 adios_handle
  INTEGER adios_err

  CALL adios_close (adios_handle, adios_err)

  RETURN
END SUBROUTINE a_closefile
