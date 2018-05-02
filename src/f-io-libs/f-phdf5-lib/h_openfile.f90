!**************************
!@section LICENSE
!Copyright (c) 2013-2017, Regents of the University of California
!All rights reserved.
!Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
!1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
!2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
!THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!**************************
SUBROUTINE h_openfile(rank,recproc,comm,filenameC, write_style,rw, &
     ghostx,ghosty,ghostz, &
     nk,kout,nj,jout,ni,iout,maxdim,&
     kbnd,jbnd,ibnd,datatype,write_step, &
     nprec, prec,&
     dsetnameC,attrC,attr_nameC, &
     start,count,start2,start2dim,ntime,icount,dimsf,iopen,type_id,&
     plist_id,dataspace,dset_id,memspace,file_id,&
     filenameLen, dsetnameLen, attrLen, attr_nameLen) bind(C)
  use iso_c_binding, only: c_char, c_int
  USE hdf5

  IMPLICIT NONE
  INCLUDE 'mpif.h'

  INTEGER nk,nj,ni,comm,kbnd(2),jbnd(2),ibnd(2),start2dim
  INTEGER ghostx,ghosty,ghostz
  INTEGER kout(nk),jout(nj),iout(ni), nprec, prec(nprec,3)
  INTEGER err,write_step,rank,recproc,maxdim

  INTEGER datatype,write_style,icount,ntime,comm2,onegroup

  INTEGER ncid,dimID(4),xvarid
  INTEGER  info, mode, iopen,i,j,k
  CHARACTER(kind = c_char), intent(in) :: dsetnameC(*), attrC(*), attr_nameC(*), filenameC(*)
  INTEGER(kind = c_int), intent(in) :: filenameLen, dsetnameLen, attrLen, attr_nameLen
  CHARACTER(len = filenameLen ) :: filename
  CHARACTER(len = dsetnameLen ) :: dsetname
  CHARACTER(len = attrLen     ) :: attr
  CHARACTER(len = attr_nameLen) :: attr_name
  CHARACTER*1 rw

  INTEGER(kind=MPI_OFFSET_KIND) attlen
  INTEGER(kind=MPI_OFFSET_KIND) start(4), COUNT(4),nx,ny,nz,nt,start2(4,start2dim)

  INTEGER nx_id,ny_id,nz_id,nt_id

  INTEGER(HID_T) :: file_id       ! File identifier
  INTEGER(HID_T) :: dset_id       ! Dataset identifier
  INTEGER(HID_T) :: dataspace     ! Dataspace identifier in data
  INTEGER(HID_T) :: memspace      ! Dataspace identifier in memory
  INTEGER(HID_T) :: plist_id      ! Property list identifier
  INTEGER(HSIZE_T), DIMENSION(4) :: dimsf  ! Dataset dimensions.

  INTEGER :: data_rank=4  ! Dataset rank ----fix it

  INTEGER(HID_T) :: type_id  ! Datatype identifier
  INTEGER(HSIZE_T), DIMENSION(1:4) :: dimsc = (/5,5,5,5/)
 ! INTEGER(HSIZE_T), DIMENSION(1:4) :: dimsc = (/1024,4,1,1/)  !64
 ! INTEGER(HSIZE_T), DIMENSION(1:4) :: dimsc = (/16,16,16,1/)  !
  INTEGER(HSIZE_T), DIMENSION(1:4) :: maxdims,offset
  INTEGER(HSSIZE_T), DIMENSION(:),ALLOCATABLE ::pointDim
  INTEGER(HSIZE_T), DIMENSION(:,:),ALLOCATABLE ::coord
  INTEGER(SIZE_T) num_elements

  INTEGER(HSIZE_T), DIMENSION(1) :: dimsf_buffer  ! Dataset dimensions.

  INTEGER(HID_T) :: creation_id       ! creation identifier
  INTEGER layout,ndims_chunk
  INTEGER(HSIZE_T), DIMENSION(1:4)::dims_chunk
  
  ! copy c-strings to fortran
  FORALL(i = 1:filenameLen ) filename(i:i)  = filenameC(i)
  FORALL(i = 1:dsetnameLen ) dsetname(i:i)  = dsetnameC(i)
  FORALL(i = 1:attrLen     ) attr(i:i)      = attrC(i)
  FORALL(i = 1:attr_nameLen) attr_name(i:i) = attr_nameC(i)

  !Hui add, Need modify 
  !hbool_t actual_metadata_ops_collective


  dimsf_buffer(1)=nprec*write_step

  !!	real buffer(1)
  icount=0
  ntime=1
  iopen=0

  info = MPI_INFO_NULL

  !  calculate start and count, nprec,prec
  CALL setup(ghostx,ghosty,ghostz,nk,kout,nj,jout,ni,iout,kbnd,jbnd,ibnd,write_style,start,count,start2,start2dim,nprec,prec)
  nx=ni
  ny=nj
  nz=nk

  IF (write_style==2)THEN
     ny=1
     nz=1
  ENDIF

  ! Initialize FORTRAN interface.
  CALL h5open_f(err)

  IF (datatype==MPI_DOUBLE_PRECISION .OR. datatype==MPI_DOUBLE) type_id =H5T_NATIVE_DOUBLE
  IF (datatype==MPI_REAL .OR. datatype==MPI_FLOAT)  type_id = H5T_NATIVE_REAL
  IF (datatype==MPI_INTEGER .OR. datatype==MPI_INT)  type_id =H5T_NATIVE_INTEGER

  !! H5T_NATIVE_CHARACTER, H5T_NATIVE_HSIZE,H5T_NATIVE_HHSIZE

  onegroup=0
  IF (rank == recproc) onegroup = 1
  CALL MPI_COMM_SPLIT(comm,onegroup,0,comm2,err)

  !!part 1 ---write

  IF (TRIM(rw)=='w')THEN
     !add by Hui 2016.10
     IF (rank .NE. recproc) RETURN

     ! Initialize FORTRAN interface.
     !    	CALL h5open_f(err)

     ! Setup file access property list with parallel I/O access.
     !
     CALL h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, err)
     ! apply latest format according to HDF5 GROUP suggestion
     !H5Pset_libver_bounds(fapl, H5F_LIBVER_LATEST, H5F_LIBVER_LATEST); //c
     CALL H5Pset_libver_bounds_f(plist_id, H5F_LIBVER_LATEST_F, H5F_LIBVER_LATEST_F, err)
     CALL h5pset_fapl_mpio_f(plist_id, comm2, info, err)

     !Feature set collective metadata reads,this feature is only available in
     !1.10 and newer.
     !CALL h5pset_all_coll_metadata_ops_f(plist_id, .true. , err) 

     
     !
     ! Create the file collectively.
     !
     !CALL h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, err, H5P_DEFAULT_F, access_prp = plist_id)
     CALL h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, err, access_prp = plist_id)
     CALL h5pclose_f(plist_id, err)

     !
     ! Create the data space for the  dataset.
     !
     !----create data space

     !maxdims=(/H5S_UNLIMITED_F,H5S_UNLIMITED_F,H5S_UNLIMITED_F,H5S_UNLIMITED_F/)
     maxdims=(/nx,ny,nz,H5S_UNLIMITED_F/)
     dimsf(1)=nx
     dimsf(2)=ny
     dimsf(3)=nz
     dimsf(4)=write_step
     CALL h5screate_simple_f(data_rank, dimsf, dataspace, err, maxdims)
     !write(*,*)rank,') dimsf:', dimsf
     !write(*,*)rank,') count:', count
     !write(*,*)rank,') dimsc:', dimsc , data_rank

     !Modify dataset creation properties, i.e. enable chunking
     !
     CALL h5pcreate_f(H5P_DATASET_CREATE_F, plist_id, err)

     !Hui try to meet suggest from Gerd Heber HDF5 Group 2016.10
     IF(count(1)*count(2)*count(3)>0 .and. count(1)*count(2)*count(3)/(1024*1024)<=250) THEN
         dimsc = (/INT(count(1),HSIZE_T), INT(count(2),HSIZE_T), INT(count(3),HSIZE_T), INT(1,HSIZE_T)/)
         !dimsc = (/INT(1,HSIZE_T),INT(count(3),HSIZE_T), INT(count(2),HSIZE_T), INT(count(1),HSIZE_T)/)
     ELSE
         dimsc = (/nx, ny, INT(1,HSIZE_T), INT(1,HSIZE_T)/)
     ENDIF

     !Feature set collective metadata reads, this feature is only available in
     !1.10 and newer.
     !CALL h5pset_all_coll_metadata_ops_f(plist_id, .TRUE. , err) 

     !Hui add, Feature Note that this feature is only available in 1.10 and newer.
     ! verify that metadata ops actually performed collectively
     !h5pget_all_coll_metadata_ops_f(plist_id, actual_metadata_ops_collective, err) 

     CALL h5pset_chunk_f(plist_id, data_rank, dimsc, err) ! Hui found problem mightbe here the size "dimsc"
     
     !Feature Never fill
     CALL h5pset_fill_time_f(plist_id, H5D_FILL_TIME_NEVER_F, err) !Hui add, Never_flle feature 
     
     !Feature Early allocation
     CALL h5pset_alloc_time_f(plist_id, H5D_ALLOC_TIME_EARLY_F, err)!Hui add, Early allocation


     !
     ! Create the dataset with default properties.
     !
     CALL h5dcreate_f(file_id, dsetname, type_id, dataspace, dset_id, err ,plist_id)
     IF(err <0) WRITE(*,*) ' Error in h5dcreate_f error.'

     CALL h5sclose_f(dataspace, err)
     IF (write_style==1)THEN

        ! Each process defines dataset in memory and writes it to the hyperslab
        ! in the file.

        offset(1:3)=start(1:3)-1
        offset(4)=0

        COUNT(4)=write_step

        CALL h5screate_simple_f(data_rank, count, memspace, err)
        !
        ! Select hyperslab in the file.
        !
        CALL h5dget_space_f(dset_id, dataspace, err)
        CALL h5sselect_hyperslab_f (dataspace, H5S_SELECT_SET_F, offset, count, err)

     ELSE

        ALLOCATE(coord(4,nprec*write_step))
        ALLOCATE(pointDim(4))

        pointDim(1)=nprec*write_step
        pointDim(2)=1
        pointDim(3)=1
        pointDim(4)=1

        CALL h5screate_simple_f(4, pointDim, memspace, err)

        CALL h5dget_space_f(dset_id, dataspace, err)

        !The coord array is a two-dimensional array of size NUMP x RANK in C (RANK x NUMP in FORTRAN)

        DO j=1,write_step
           DO i =1,nprec
              k=(j-1)*nprec + i
              coord(1,k)=start2(1,i)
              coord(2,k)=start2(2,i)
              coord(3,k)=start2(3,i)
              coord(4,k)=j
           ENDDO
        ENDDO

        num_elements=nprec*write_step
        IF (nprec >0) &
             CALL h5Sselect_elements_f(dataspace,H5S_SELECT_SET_F,data_rank,num_elements,coord,err)

     ENDIF
     !
     ! Create property list for collective dataset write
     !
     CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, err)

     CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, err)

  ELSE
!!!!!!reader----------
     ! Setup file access property list with parallel I/O access.
     !
     CALL h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, err)
     CALL h5pset_fapl_mpio_f(plist_id, comm2, info, err)

     CALL H5Fopen_f(filename,H5F_ACC_RDONLY_F,file_id,err,plist_id)
     IF(err < 0) WRITE(*,*) ' Open file occur error.'

     CALL H5dopen_f(file_id,dsetname,dset_id,err);

     !  get dataset rank and dimension

     CALL H5Dget_space_f(dset_id,dataspace,err)
     CALL H5Sget_simple_extent_ndims_f(dataspace,data_rank,err)
     CALL H5Sget_simple_extent_dims_f(dataspace,dimsf,maxdims,err)

     !	write(*,*)' read: dimsf =',dimsf,'maxdims=',maxdims,'data_rank=',data_rank

     !  get creation properties list
     CALL H5Dget_create_plist_f(dset_id,creation_id,err)
     IF(err <0) WRITE(*,*) ' Error in getting creation property list id.'

     !  check if dataset is chunked

     IF (1>10)THEN

        ! this part crashed the c version for some reason

        CALL H5Pget_layout_f(creation_id,layout,err)
        IF(err <0) WRITE(*,*) ' Error in getting layout.'

        IF(H5D_CHUNKED_F==layout)THEN
           !get chunking information: rank and dimensions
           CALL H5Pget_chunk_f(creation_id,ndims_chunk,dims_chunk,err)
           IF(err <0) WRITE(*,*) ' Error in getting chunk.'
           !		write(*,*)' chunk dimension: ndims,dims= ',ndims_chunk,dims_chunk,'err=',err
        ENDIF
     ENDIF
     !  define memory space to read dataset ---use count instead of dimsf ---in parallel

     ! slice by slice in time

     !count(4)=1
     !CALL h5screate_simple_f(data_rank, count, memspace, err)

     ALLOCATE(pointDim(4))

     pointDim(1)=nprec*1
     pointDim(2)=1
     pointDim(3)=1
     pointDim(4)=1

     CALL h5screate_simple_f(4, pointDim, memspace, err)

     IF(err <0) WRITE(*,*) ' Error in create_simple.'

  ENDIF

  RETURN
END SUBROUTINE h_openfile
