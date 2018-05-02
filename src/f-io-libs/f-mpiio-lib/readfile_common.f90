!**************************
!@section LICENSE
!Copyright (c) 2013-2017, Regents of the University of California
!All rights reserved.
!Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
!1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
!2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
!THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!**************************
  
  INTEGER icount,nptotal
  INTEGER rank,recproc,fh,datatype,filetype
  INTEGER j,xn,yn,zn,status(MPI_STATUS_SIZE),err,write_step
  INTEGER(kind=MPI_OFFSET_KIND)::disp,nrec
  INTEGER prec(nprec,3),nbyte

  IF (rank .NE. recproc) RETURN

  IF (write_step .NE. 1)THEN
     !	write(*,*)' Error in reading: write_step should equals 1. return'
     !	return
  ENDIF
  nptotal=nprec*1

  icount=0
  icount=icount+1


  !if(icount==write_step)then

  CALL MPI_FILE_SET_VIEW(fh,disp,datatype,filetype,'native',MPI_INFO_NULL,err)

  CALL MPI_FILE_READ_ALL(fh,buffer,nptotal,datatype,status,err)

  !endif

  DO j=1,nprec
     xn=prec(j,1)
     yn=prec(j,2)
     zn=prec(j,3)
     var(xn,yn,zn)= buffer(nprec*(icount-1)+j)
  ENDDO
