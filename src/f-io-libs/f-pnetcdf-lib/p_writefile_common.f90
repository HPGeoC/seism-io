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
  INTEGER(kind=MPI_OFFSET_KIND) start(4), COUNT(4),start2(4,start2dim), bufcount
  INTEGER ncid,xvarid,ntime,j,xn,yn,zn,icount,err
  INTEGER prec(nprec,3),rank,recproc,datatype,datatype_tmp

  INTEGER write_style
  INTEGER status
  INTEGER, DIMENSION(:),ALLOCATABLE::request,st

  IF (rank .NE. recproc) RETURN

  !         write variables in parallel
  icount=icount+1

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

     start(4)=ntime
     COUNT(4)=write_step

     !        write(*,*)'W, IN :',rank, ncid,xvarid,start,count, nprec*write_step,'ntime',ntime

     IF (write_style==1 )THEN

        !	write(*,*)' check :',nprec,write_step,buffer,datatype==MPI_INTEGER,datatype==MPI_INT
        bufcount=nprec*write_step
        datatype_tmp=datatype

        err = nfmpi_put_vara_all(ncid, xvarid, start, count, buffer, bufcount, datatype_tmp)

        IF (err .NE. NF_NOERR .AND. rank == 0) THEN
           WRITE(6,*) "Error: put_vara  -", TRIM(nfmpi_strerror(err))
           !WRITE(*,*)' rank=',rank,start,'count=',count
           !WRITE(*,*)' ncid, xvarid =',ncid, xvarid,'buffer=',buffer
        ENDIF

     ELSE
        ALLOCATE(request(nprec))
        ALLOCATE(st(nprec))

        buffer_x(:)=buffer(:)

        COUNT(1:3)=1
        DO j=1,nprec

           start2(4,j)=ntime
           bufcount=write_step

           err=nfmpi_iput_vara(ncid, xvarid, start2(1,j), &
                count,buffer_x((j-1)*write_step+1),bufcount,datatype,request(j))

           !	write(*,*)'rank,iput= ',rank,start2(1:3,j),count(1:3),buffer((j-1)*write_step+1)
        ENDDO

        err=nfmpi_wait_all(ncid,nprec,request,st)

        DEALLOCATE(request)
        DEALLOCATE(st)

        IF (err .NE. NF_NOERR .AND. rank == 0) THEN
           WRITE(6,*) "Error: iput_vara_real  -", TRIM(nfmpi_strerror(err))
           WRITE(*,*)' rank=',rank,start,'count=',count
        ENDIF

     ENDIF

     icount=0
     ntime=ntime+write_step

  ENDIF
