!**************************
!@section LICENSE
!Copyright (c) 2013-2017, Regents of the University of California
!All rights reserved.
!Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
!1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
!2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
!THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!**************************
SUBROUTINE p_openfile(rank,recproc,comm,filenameC, write_style,rw, &
     ghostx,ghosty,ghostz, &
     nk,kout,nj,jout,ni,iout,maxdim,&
     kbnd,jbnd,ibnd,datatype,write_step, &
     nprec, prec,&
     attrC,attr_nameC, &
     dimXnameC,dimYnameC,dimZnameC,dimTnameC,varnameC,&
     start,count,start2,start2dim,ntime,xvarid,icount,ncid,&
     filenameLen, attrLen, attr_nameLen, dimXnameLen, dimYnameLen, dimZnameLen, dimTnameLen, varnameLen) bind(C)
  use iso_c_binding, only: c_char, c_int
  USE pnetcdf

  IMPLICIT NONE
  INCLUDE 'mpif.h'

  INTEGER nk,nj,ni,comm,kbnd(2),jbnd(2),ibnd(2),start2dim
  INTEGER ghostx,ghosty,ghostz
  INTEGER kout(nk),jout(nj),iout(ni), nprec, prec(nprec,3)
  INTEGER err,write_step,rank,recproc,maxdim

  INTEGER datatype,write_style,icount,ntime,comm2,onegroup

  INTEGER ncid,dimID(4),xvarid
  INTEGER  info, mode
  INTEGER i
  CHARACTER(kind = c_char), intent(in) :: filenameC(*), attrC(*), attr_nameC(*), dimXnameC(*), &
                                          dimYnameC(*), dimZnameC(*), dimTnameC(*), varnameC(*)
  INTEGER(kind = c_int), intent(in) :: filenameLen, attrLen, attr_nameLen, dimXnameLen, &
                                       dimYnameLen, dimZnameLen, dimTnameLen, varnameLen
  CHARACTER(len = filenameLen ) :: filename
  CHARACTER(len = attrLen     ) :: attr
  CHARACTER(len = attr_nameLen) :: attr_name
  CHARACTER(len = dimXnameLen ) :: dimXname
  CHARACTER(len = dimYnameLen ) :: dimYname
  CHARACTER(len = dimZnameLen ) :: dimZname
  CHARACTER(len = dimTnameLen ) :: dimTname
  CHARACTER(len = varnameLen )  :: varname
  CHARACTER*1 rw

  INTEGER(kind=MPI_OFFSET_KIND) attlen
  INTEGER(kind=MPI_OFFSET_KIND) start(4),COUNT(4),nx,ny,nz,nt,start2(4,start2dim)

  INTEGER nx_id,ny_id,nz_id,nt_id,nx_save,ny_save,nz_save

  ! copy c-strings to fortran
  FORALL(i = 1:filenameLen ) filename(i:i)  = filenameC(i)
  FORALL(i = 1:attrLen     ) attr(i:i)      = attrC(i)
  FORALL(i = 1:attr_nameLen) attr_name(i:i) = attr_nameC(i)
  FORALL(i = 1:dimXnameLen ) dimXname(i:i)  = dimXnameC(i)
  FORALL(i = 1:dimYnameLen ) dimYname(i:i)  = dimYnameC(i)
  FORALL(i = 1:dimZnameLen ) dimZname(i:i)  = dimZnameC(i)
  FORALL(i = 1:dimTnameLen ) dimTname(i:i)  = dimTnameC(i)
  FORALL(i = 1:varnameLen  ) varname(i:i)   = varnameC(i)

  !!	real buffer(1)
  icount=0
  ntime=1

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

  nx_save=nx
  ny_save=ny
  nz_save=nz

  onegroup=0
  IF (rank == recproc) onegroup = 1
  CALL MPI_COMM_SPLIT(comm,onegroup,0,comm2,err)

  IF (TRIM(rw)=='w')THEN
    !http://www.unidata.ucar.edu/software/netcdf/faq-lfs.html
     !mode = IOR(NF_CLOBBER, NF_64BIT_OFFSET)  !f77
     mode = or(nf90_clobber,nf90_64bit_offset) !f90

     err = nfmpi_create(comm2, TRIM(filename), mode, info, ncid)
     IF (err .NE. NF_NOERR .AND. rank == 0) THEN
        WRITE(6,*) "Error: create", TRIM(nfmpi_strerror(err))
     ENDIF

     !             define dimensions
     err = nfmpi_def_dim(ncid, TRIM(dimXname),   nx, dimID(1))
     IF (err .NE. NF_NOERR .AND. rank == 0) THEN
        WRITE(6,*) "Error: nx-", TRIM(nfmpi_strerror(err))
     ENDIF

     err = nfmpi_def_dim(ncid, TRIM(dimYname),   ny, dimID(2))
     IF (err .NE. NF_NOERR .AND. rank == 0) THEN
        WRITE(6,*) "Error: ny-", TRIM(nfmpi_strerror(err))
     ENDIF


     err = nfmpi_def_dim(ncid, TRIM(dimZname),   nz, dimID(3))
     IF (err .NE. NF_NOERR .AND. rank == 0) THEN
        WRITE(6,*) "Error: nz-", TRIM(nfmpi_strerror(err))
     ENDIF

     err = nfmpi_def_dim(ncid, TRIM(dimTname), NFMPI_UNLIMITED,dimID(4))
     IF (err .NE. NF_NOERR .AND. rank == 0) THEN
        WRITE(6,*) "Error: def dim-", TRIM(nfmpi_strerror(err))
     ENDIF

     !define variables
     !
     !NF_BYTE (NC_CHAR)
     IF (datatype==MPI_REAL .OR. datatype==MPI_FLOAT )THEN
        err = nfmpi_def_var(ncid, TRIM(varname), NF_REAL, 4, dimID, xvarid)
     ELSE
        IF (datatype==MPI_INTEGER .OR. datatype==MPI_INT)THEN
           err = nfmpi_def_var(ncid, TRIM(varname), NF_INT, 4, dimID, xvarid)
        ELSE
           IF (datatype==MPI_DOUBLE_PRECISION.OR.datatype==MPI_DOUBLE)THEN
              err = nfmpi_def_var(ncid, TRIM(varname), NF_DOUBLE, 4, dimID, xvarid)
           ENDIF
        ENDIF
     ENDIF
     !when 3K cores,occur error: One or more variable sizes violate format constraints
     !Acutually >512 cores, the .nc file became nknown file format, using ncdump -h read
     !Every var can not bigger than 4GiB, So every rank have 110M and X rank(400)> 4G
     !http://www.unidata.ucar.edu/software/netcdf/faq-lfs.html
     !http://www.unidata.ucar.edu/software/netcdf/netcdf-4/newdocs/netcdf/Limitations.html
     !http://www.unidata.ucar.edu/software/netcdf/netcdf-4/newdocs/netcdf/64-bit-Offset-Limitations.html#g_t64-bit-Offset-Limitations
     IF (err .NE. NF_NOERR .AND. rank == 0) THEN
        WRITE(6,*) "Error: def var ", TRIM(nfmpi_strerror(err))
     ENDIF

!!! add sample attribute for variable buffer
     ! attr = 'some text attribute annotation about var X'
     !attlen = len_TRIM(attr)

     !err = nfmpi_put_att_text(ncid, xvarid, attr_name, attlen, attr)
     !IF (err .NE. NF_NOERR .AND. rank == 0) THEN
     !   WRITE(6,*) "Error:att text", TRIM(nfmpi_strerror(err))
     !ENDIF

     !             end of define mode
     err = nfmpi_enddef(ncid)

  ELSE
!!!!!!reader----------

     mode= NF_NOWRITE
     err = nfmpi_open(comm2,TRIM(filename), mode, info, ncid)
     IF (err .NE. NF_NOERR) THEN
        WRITE(6,*) "Error: create", TRIM(nfmpi_strerror(err))
     ENDIF

     err= nfmpi_inq_dimid(ncid,TRIM(dimXname),nx_id)
     err= nfmpi_inq_dimid(ncid,TRIM(dimYname),ny_id)
     err= nfmpi_inq_dimid(ncid,TRIM(dimZname),nz_id)
     err= nfmpi_inq_dimid(ncid,TRIM(dimTname),nt_id)

     err=nfmpi_inq_dimlen(ncid,nx_id,nx)
     err=nfmpi_inq_dimlen(ncid,ny_id,ny)
     err=nfmpi_inq_dimlen(ncid,nz_id,nz)
     err=nfmpi_inq_dimlen(ncid,nt_id,nt)

     IF (nx.NE.nx_save .OR. ny.NE.ny_save .OR. nz.NE.nz_save)THEN
        WRITE(*,*)" Dimensions not match: read nx,ny,nz=", nx,ny,nz
        WRITE(*,*)" But nx_save,ny_save,nz_save =",nx_save,ny_save,nz_save
     ENDIF

     err=nfmpi_inq_varid(ncid,varname,xvarid)

  ENDIF


  RETURN
END SUBROUTINE p_openfile
