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

  !! one step is read
  icount=1

  dimsf_buffer(1)=nprec*1

  IF (rank .NE. recproc) RETURN

  !   write variables in parallel

  write_style=1
  IF ((start(1)+start(2)+start(3)) == 0 ) write_style = 2
  buffer = 0
  IF (write_style==1)THEN

     ! Select hyperslab in the file.
     offset(1:3)=start(1:3)-1

     !for only first time read
     offset(4)=write_step-1

     !for multiple read
     !start(4)=ntime
     !offset(4)=start(4)-1  !modified by Hui

     COUNT(4)=1

     CALL h5dget_space_f(dset_id, dataspace, err)
     CALL h5sselect_hyperslab_f (dataspace, H5S_SELECT_SET_F, offset, count, err)

     IF (nprec >0) &
          CALL h5dread_f(dset_id, type_id, buffer,dimsf_buffer,err,memspace,dataspace)

  ELSE
     !! write_stype==2 ---Green's points
     !! st step = write_step

     ALLOCATE(coord(4,nprec*1))

     DO i =1,nprec
        k= i
        coord(1,k)=start2(1,i)
        coord(2,k)=start2(2,i)
        coord(3,k)=start2(3,i)
        coord(4,k)=write_step
     ENDDO


     num_elements=nprec*1

     IF (nprec >0) THEN
        CALL h5dget_space_f(dset_id, dataspace, err)
        CALL h5Sselect_elements_f(dataspace,H5S_SELECT_SET_F,data_rank,num_elements,coord,err)
        CALL h5dread_f(dset_id, type_id, buffer,dimsf_buffer,err,memspace,dataspace)
     ENDIF

  ENDIF

  DO j=1,nprec
     xn=prec(j,1)
     yn=prec(j,2)
     zn=prec(j,3)
     var(xn,yn,zn)=  buffer(nprec*(icount-1)+j)
  ENDDO

  !for multiple read
  !ntime=ntime + 1
