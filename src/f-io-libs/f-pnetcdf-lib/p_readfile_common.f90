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
  INTEGER(kind=MPI_OFFSET_KIND) start(4), COUNT(4),start2(4,start2dim),bufcount
  INTEGER ncid,xvarid,ntime,j,xn,yn,zn,icount,nprec,write_step,nx,ny,nz,err
  INTEGER prec(nprec,3),rank,recproc,datatype,write_style

  INTEGER, DIMENSION(:),ALLOCATABLE::request,st
  INTEGER status

  icount=1

  buffer=0

  !read one block: write_step is the cycle

  start(4)=write_step
  COUNT(4)=1
  IF (COUNT(1)==0) COUNT(4)=0

  !        write(*,*)' IN :', rank,ncid,xvarid,start,count, nprec*write_step,'ntime',ntime

  write_style=1
  IF ((start(1)+start(2)+start(3)) == 0 ) write_style = 2

  IF (write_style ==1)THEN
     bufcount=nprec

     err = nfmpi_get_vara_all(ncid, xvarid, start, count,buffer,bufcount,datatype)
  ELSE

     ALLOCATE(request(nprec))
     ALLOCATE(st(nprec))

     COUNT(1:3)=1
     DO j=1,nprec

        start2(4,j)=write_step
        bufcount=1
        err=nfmpi_iget_vara(ncid, xvarid, start2(1,j), &
             count,buffer((j-1)*1+1),bufcount,datatype,request(j))

        !       write(*,*)'rank,iput= ',rank,start2(1:3,j),count(1:3),buffer((j-1)*write_step+1)

     ENDDO

     err=nfmpi_wait_all(ncid,nprec,request,st)

     DEALLOCATE(request)
     DEALLOCATE(st)

  ENDIF
  IF (err .NE. NF_NOERR) THEN
     WRITE(6,*) "Error: get_vara  -", TRIM(nfmpi_strerror(err))
     WRITE(*,*)' rank=',rank,start,'count=',count
  ENDIF

  DO j=1,nprec
     xn=prec(j,1)
     yn=prec(j,2)
     zn=prec(j,3)
     var(xn,yn,zn)=buffer(nprec*(icount-1)+j)
  ENDDO
