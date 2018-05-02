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
  INTEGER ncid,xvarid,ntime,j,xn,yn,zn,icount,err,i,k
  INTEGER prec(nprec,3),rank,recproc,datatype

  INTEGER write_style,iopen
  INTEGER status
  INTEGER, DIMENSION(:),ALLOCATABLE::request,st

  INTEGER(HID_T) :: file_id       ! File identifier
  INTEGER(HID_T) :: dset_id       ! Dataset identifier
  INTEGER(HID_T) :: dataspace     ! Dataspace identifier in file
  INTEGER(HID_T) :: memspace      ! Dataspace identifier in memory
  INTEGER(HID_T) :: plist_id      ! Property list identifier
  INTEGER(HSIZE_T), DIMENSION(4) :: dimsf  ! Dataset dimensions.
  INTEGER(HID_T), INTENT(IN) :: type_id  ! Datatype identifier

  INTEGER(HSSIZE_T),DIMENSION(4)::offset
  INTEGER(HSIZE_T), DIMENSION(4) ::size

  INTEGER(HSIZE_T), DIMENSION(1) :: dimsf_buffer  ! Dataset dimensions.

  INTEGER(SIZE_T) num_elements
  INTEGER :: data_rank=4  ! Dataset rank ----fix it
  INTEGER(HSIZE_T), DIMENSION(:,:),ALLOCATABLE ::coord

  dimsf_buffer(1)=nprec*write_step

  IF (rank .NE. recproc) RETURN

  ! write variables in parallel

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
        !buffer(nprec*(icount-1)+j)=var(xn,yn,zn)

     ENDIF
  ENDDO


  !write(*,*)rank,' IN write*:icount,write_step,ncid (3 reals)=', icount,write_step,ncid
  !write(*,*)rank,' rank buffer=',rank,buffer,write_step,icount
  !write(*,*)rank,' start : ',rank,start
  !write(*,*)rank,' count : ',rank,count
  !write(*,*)rank,' Using count(4)'

  IF (icount==write_step)THEN

     start(4)=ntime
     COUNT(4)=write_step

     !write(*,*)rank,' W, IN :',rank, ncid,xvarid,start,count, nprec*write_step,'ntime',ntime

     IF (write_style==1 )THEN

        ! Write the dataset collectively.
        !
        IF (iopen==0)THEN
           !! first time writing

           CALL h5dwrite_f(dset_id, type_id, buffer, dimsf_buffer, err, memspace,dataspace,plist_id)

           iopen=1
        ELSE
           ! extend the datasets

           SIZE(1:3)   = dimsf(1:3)
           SIZE(4) = start(4)-1+COUNT(4)

           CALL h5dset_extent_f(dset_id, size,err)
           IF(err <0) WRITE(*,*) rank,' Error in h5dset_extent_f error.'

           offset(1:4)=start(1:4)-1
           !write(*,*)rank,' SIZE:',size,offset
           CALL h5dget_space_f(dset_id, dataspace, err)
           IF(err <0) WRITE(*,*) rank,' Error in h5dget_space_f error.'
           !    h5sselect_hyperslab_f(space_id,  operator,         start,  count, hdferr, stride, block)
           CALL h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, offset, count, err)
           IF(err <0) WRITE(*,*) rank,' Error in h5sselect_hyperslab_f error.'

           !CALL H5dwrite_f(dset_id,type_id, buffer, dimsf_buffer, err, memspace, dataspace,plist_id)
           CALL H5dwrite_f(dset_id,type_id, buffer, dimsf_buffer, err, memspace, dataspace) ! Hui found there is no "plist_id" used
           IF(err <0) WRITE(*,*) rank,' Error in H5dwrite_f error.'

        ENDIF
	!
     ELSE
	! add scattered points
        IF (iopen==0)THEN
           !! first time writing

           !this call is bad somewhat

           !                CALL h5dget_space_f(dset_id, dataspace, err)

           IF (nprec >0) &
                CALL h5dwrite_f(dset_id, type_id, buffer,dimsf_buffer,err,memspace,dataspace,plist_id)
           iopen=1
        ELSE
           ! extend the datasets

           SIZE(1:3)   = dimsf(1:3)
           SIZE(4) = start(4)-1+COUNT(4)

           CALL h5dset_extent_f(dset_id, size,err)

           CALL h5dget_space_f(dset_id, dataspace, err)

           ALLOCATE(coord(4,nprec*write_step))
           DO j=1,write_step
              DO i =1,nprec
                 k=(j-1)*nprec + i
                 coord(1,k)=start2(1,i)
                 coord(2,k)=start2(2,i)
                 coord(3,k)=start2(3,i)
                 coord(4,k)=j+ntime-1
              ENDDO
           ENDDO

           num_elements=nprec*write_step

           IF (nprec >0) THEN
              CALL h5Sselect_elements_f(dataspace,H5S_SELECT_SET_F,data_rank,num_elements,coord,err)
              CALL h5dwrite_f(dset_id, type_id, buffer,dimsf_buffer,err,memspace,dataspace) !
           ENDIF

        ENDIF

     ENDIF
     ! done writing

     icount=0
     ntime=ntime+write_step

  ENDIF
