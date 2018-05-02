!**************************
!@section LICENSE
!Copyright (c) 2013-2017, Regents of the University of California
!All rights reserved.
!Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
!1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
!2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
!THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!**************************
  INTEGER start2dim
  INTEGER(kind=MPI_OFFSET_KIND) start(4), COUNT(4),start2(4,start2dim)
  INTEGER ncid,xvarid,ntime,j,xn,yn,zn,icount,iopen,err
  INTEGER prec(nprec,3),rank,recproc,datatype
  INTEGER write_style
  INTEGER adios_err,istep,offx,offt
  INTEGER*8 adios_handle
  INTEGER                              :: i
  CHARACTER(kind = c_char), intent(in) :: varnameC(*)
  INTEGER(kind = c_int), intent(in) :: varnameLen
  CHARACTER(len = varnameLen) :: varname

  ! copy c-name to fortran
  FORALL(i = 1:varnameLen ) varname(i:i)  = varnameC(i)

  IF (rank .NE. recproc) RETURN

  !         write variables in parallel

  write_style=1
  IF ((start(1)+start(2)+start(3)) == 0 ) write_style = 2


  DO j=1,nprec
     xn=prec(j,1)
     yn=prec(j,2)
     zn=prec(j,3)
     IF (write_style==1 )THEN
        buffer(nprec*(icount-1)+j)=var(xn,yn,zn)
     ELSE
        buffer(icount+(j-1)*write_step)=var(xn,yn,zn)

     ENDIF
  ENDDO

  !        write(*,*)' IN write*:icount,write_step,ncid (3 reals)=', icount,write_step,ncid
  !        write(*,*)' rank buffer=',rank,buffer,write_step,icount
  !        write(*,*)' start : ',rank,start
  !        write(*,*)' count : ',rank,count
  !	write(*,*)' Using count(4)'

  IF (icount==write_step)THEN

     IF (write_style==1 )THEN


        IF (nprec>0)        CALL adios_write (adios_handle, TRIM(varname),buffer, adios_err)

     ELSE
	!! point-to-point write

        DO istep = 1,write_step
           offt = istep-1
           DO j = 1,nprec

              offx=start2(1,j)-1

              CALL adios_write (adios_handle, "offx",offx, adios_err)
              CALL adios_write (adios_handle, "offt",offt, adios_err)
              CALL adios_write (adios_handle, TRIM(varname),buffer(offt*nprec+j), adios_err)
           ENDDO
        ENDDO

     ENDIF

     !icount=0
     iopen=1  !mark the file is already opened

  ENDIF
