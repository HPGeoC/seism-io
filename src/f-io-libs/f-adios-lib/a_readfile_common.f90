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
  INTEGER(kind=MPI_OFFSET_KIND) start(4),COUNT(4),start2(4,start2dim)
  INTEGER ncid,xvarid,ntime,j,xn,yn,zn,icount,iopen,nprec,write_step,nx,ny,nz,err
  INTEGER prec(nprec,3),rank,recproc,datatype
  INTEGER write_style
  INTEGER adios_err,istep,offx,offt
  INTEGER*8 adios_handle
  CHARACTER(kind = c_char), intent(in) :: varnameC(*)
  INTEGER(kind = c_int), intent(in) :: varnameLen
  CHARACTER(len = varnameLen) :: varname

  INTEGER*8               :: sel, npoints

  INTEGER write_step_read,iadvance_step,inside,i,lnx,lny,lnz,ndim,ntotal,index

  INTEGER*8,DIMENSION(:),ALLOCATABLE::points
  INTEGER from_step,nsteps

  ! copy c-name to fortran
  FORALL(i = 1:varnameLen ) varname(i:i)  = varnameC(i)

  IF (nprec <=0 ) RETURN

  write_style=1
  IF ((start(1)+start(2)+start(3)) == 0 ) write_style = 2

  ! select specific writer's data (using rank since each rank of
  ! the writer
  ! wrote one block of the data)

  CALL adios_selection_writeblock (sel, rank)
  from_step = 0
  nsteps = 1
  CALL adios_schedule_read (adios_handle,sel,'write_step',from_step,nsteps,write_step_read,adios_err)
  CALL adios_perform_reads (adios_handle, adios_err)
  ALLOCATE(buffer_block(nprec*write_step_read))

  IF (write_style==1 )THEN

     ! First get the scalars to calculate the size of the arrays.
     ! Note that we cannot use adios_get_scalar here because that
     !   retrieves the same NX for everyone (from writer rank 0).
     !       call adios_schedule_read (adios_handle, sel,
     !       trim(dimXname), 1, 1, nx, adios_err)
     !       call adios_schedule_read (adios_handle, sel,
     !       trim(dimYname), 1, 1, ny, adios_err)
     !       call adios_schedule_read (adios_handle, sel,
     !       trim(dimZname), 1, 1, nz, adios_err)

     CALL adios_schedule_read (adios_handle,sel,'lnx',from_step,nsteps,lnx,adios_err)
     CALL adios_schedule_read (adios_handle,sel,'lny',from_step,nsteps,lny,adios_err)
     CALL adios_schedule_read (adios_handle,sel,'lnz',from_step,nsteps,lnz,adios_err)
     CALL adios_perform_reads (adios_handle, adios_err)

     !                       write(*,*)' rank =',rank,lnx,lny,lnz,
     !                       write_step_read

     !read the arrays
     start(4)=1
     COUNT(4)=write_step_read

     start=start-1
     CALL adios_selection_boundingbox(sel,4,start,count)


  ELSE
     !!  by points --type 2

     npoints=nprec*write_step_read
     ndim =2
     ALLOCATE(points(nprec*write_step_read*ndim))

     CALL adios_schedule_read (adios_handle,sel,'Ntotal',from_step,nsteps,ntotal,adios_err)
     CALL adios_perform_reads(adios_handle,adios_err)

     !               write(*,*)' rank =',
     !               rank,ntotal,write_step_read,'nprec=',nprec

     index=0
     DO istep=1,write_step_read
        DO i=1,nprec
           index=index+1
           points(index)=start2(1,i)
           index=index+1
           points(index)=istep
        ENDDO
     ENDDO

     points=points-1

     CALL adios_selection_points(sel,ndim,npoints,points)

  ENDIF

  iadvance_step =(write_step-0.1)/write_step_read

  inside = write_step- write_step_read*iadvance_step

  DO i=1, iadvance_step
     CALL adios_advance_step(adios_handle,0,1.0,adios_err)

  ENDDO

  !               write(*,*) 'inside = ',inside,' write_step
  !               =',write_step,'write_step_read = ', write_step_read

  CALL adios_schedule_read(adios_handle,sel,TRIM(varname),from_step,nsteps,buffer_block,adios_err)
  CALL adios_perform_reads(adios_handle,adios_err)

  buffer(1:nprec)=buffer_block((inside-1)*nprec + 1:nprec)

  !               if (nprec >0) write(*,*)' cycle#: rank =',
  !               write_step,rank,nprec,buffer(1:nprec)

  IF (nprec >0)THEN
     DEALLOCATE(buffer_block)
     IF (write_style==2) DEALLOCATE(points)

  ENDIF

  icount=1
  DO j=1,nprec
     xn=prec(j,1)
     yn=prec(j,2)
     zn=prec(j,3)

     var(xn,yn,zn)= buffer(nprec*(icount-1)+j)

  ENDDO

  !CALL adios_read_close(adios_handle,adios_err)
