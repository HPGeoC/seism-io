!**************************
!@section LICENSE
!Copyright (c) 2013-2017, Regents of the University of California
!All rights reserved.
!Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
!1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
!2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
!THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!**************************
SUBROUTINE a_openfile(rank,recproc,comm,filenameC, write_style,rw, &
     ghostx,ghosty,ghostz, &
     nk,kout,nj,jout,ni,iout,maxdim,&
     kbnd,jbnd,ibnd,datatype,write_step, &
     nprec, prec,&
     attrC,attr_nameC, grpnameC,mdnameC,adios_buffer_size,&
     dimXnameC,dimYnameC,dimZnameC,dimTnameC,&
     start,count,start2,start2dim,ntime,varnameC,icount,iopen,m_adios_group,adios_handle,&
     filenameLen, attrLen, attr_nameLen, grpnameLen, mdnameLen,&
     dimXnameLen, dimYnameLen, dimZnameLen, dimTnameLen, varnameLen ) bind(C)
  use iso_c_binding, only: c_char, c_int
  USE adios_write_mod

  IMPLICIT NONE
  INCLUDE 'mpif.h'

  INTEGER nk,nj,ni,comm,kbnd(2),jbnd(2),ibnd(2),start2dim
  INTEGER ghostx,ghosty,ghostz
  INTEGER kout(nk),jout(nj),iout(ni), nprec, prec(nprec,3)
  INTEGER err,write_step,rank,recproc,maxdim

  INTEGER datatype,write_style,icount,iopen,ntime,comm2,onegroup

  !	integer ncid,dimID(4),xvarid
  INTEGER  info, mode,adios_buffer_size
  CHARACTER(kind = c_char), intent(in) :: filenameC(*), attrC(*), attr_nameC(*), grpnameC(*), mdnameC(*), &
                                          dimXnameC(*), dimYnameC(*), dimZnameC(*), dimTnameC(*), varnameC(*)
  INTEGER(kind = c_int), intent(in) :: filenameLen, attrLen, attr_nameLen, grpnameLen, mdnameLen, &
                                       dimXnameLen, dimYnameLen, dimZnameLen, dimTnameLen, varnameLen
  CHARACTER(len = filenameLen ) :: filename
  CHARACTER(len = attrLen     ) :: attr
  CHARACTER(len = attr_nameLen) :: attr_name
  CHARACTER(len = grpnameLen  ) :: grpname
  CHARACTER(len = mdnameLen   ) :: mdname
  CHARACTER(len = dimXnameLen ) :: dimXname
  CHARACTER(len = dimYnameLen ) :: dimYname
  CHARACTER(len = dimZnameLen ) :: dimZname
  CHARACTER(len = dimTnameLen ) :: dimTname
  CHARACTER(len = varnameLen  ) :: varname
  CHARACTER*1 rw,rwtemp

  INTEGER(kind=MPI_OFFSET_KIND) attlen
  INTEGER(kind=MPI_OFFSET_KIND) start(4), COUNT(4),nx,ny,nz,nt,start2(4,start2dim)

  INTEGER nx_id,ny_id,nz_id,nt_id

  INTEGER adios_err,i,istep,int_size
  INTEGER*8 varid,m_adios_group,adios_handle,adios_groupsize, adios_totalsize
  INTEGER is_integer,is_real,is_double,is_long,is_long_double,var_type,var_size
  INTEGER lnx,lny,lnz,offx,offy,offz,offt,nVar

  ! ADIOS variables declarations for reading

  INTEGER                 :: method = ADIOS_READ_METHOD_BP
  INTEGER*8               :: sel
  CHARACTER*250 parameters
  integer*8               :: f
  ! variables to read in
  integer                 :: Y, ierr

  ! copy c-strings to fortran
  FORALL(i = 1:filenameLen ) filename(i:i)  = filenameC(i)
  FORALL(i = 1:attrLen     ) attr(i:i)      = attrC(i)
  FORALL(i = 1:attr_nameLen) attr_name(i:i) = attr_nameC(i)
  FORALL(i = 1:grpnameLen  ) grpname(i:i)   = grpnameC(i)
  FORALL(i = 1:mdnameLen   ) mdname(i:i)    = mdnameC(i)
  FORALL(i = 1:dimXnameLen ) dimXname(i:i)  = dimXnameC(i)
  FORALL(i = 1:dimYnameLen ) dimYname(i:i)  = dimYnameC(i)
  FORALL(i = 1:dimZnameLen ) dimZname(i:i)  = dimZnameC(i)
  FORALL(i = 1:dimTnameLen ) dimTname(i:i)  = dimTnameC(i)
  FORALL(i = 1:varnameLen  ) varname(i:i)   = varnameC(i)

  !data types in adios--fortran
  is_integer=2
  is_long=4
  is_real=5
  is_double=6
  is_long_double=7

  CALL MPI_TYPE_SIZE(datatype,var_size,err)
  int_size = 4   ! fortran

  var_type= 2         !default integer
  IF (datatype==MPI_REAL .OR. datatype==MPI_FLOAT) THEN
     var_type= 5
  ELSE
     IF (datatype==MPI_DOUBLE_PRECISION .OR. datatype==MPI_DOUBLE) THEN
        var_type= 6
     ENDIF
  ENDIF

  !!	real buffer(1)
  ntime=1

  !  calculate start and count, nprec,prec
  CALL setup(ghostx,ghosty,ghostz,nk,kout,nj,jout,ni,iout,kbnd,jbnd,ibnd,write_style,start,count,start2,start2dim,nprec,prec)

  nx=ni
  ny=nj
  nz=nk

  IF (write_style==2)THEN
     IF ( (ni.NE.nj) .OR. (ni.NE.nk ))THEN
        WRITE(*,*)" Settings ni,nj,nk for write_stype=2 are incorrect !. Stop."
        STOP
     ENDIF

     ny=1
     nz=1
  ENDIF

  onegroup=0
  if (rank == recproc) onegroup = 1
  call MPI_COMM_SPLIT(comm,onegroup,0,comm2,err)

  IF (TRIM(rw)=='w')THEN

     icount=icount +1
     IF (icount .NE. write_step ) RETURN

     IF (iopen==0) THEN
        rwtemp='w'
     ELSE
        rwtemp='a'
     ENDIF

     !placed in main
     IF (TRIM(rw)=='x' )THEN
        CALL adios_init_noxml (comm2,adios_err);
        CALL adios_allocate_buffer (adios_buffer_size,adios_err);  !10MB
     ENDIF

     IF (iopen==0)THEN
        CALL adios_init_noxml (comm2, adios_err);
        CALL adios_allocate_buffer (adios_buffer_size, adios_err);

        CALL adios_declare_group (m_adios_group, TRIM(grpname), "", 1, adios_err);
        CALL adios_select_method (m_adios_group, TRIM(mdname), "", "", adios_err);
     ENDIF

     IF (rank .NE. recproc) RETURN

     IF (write_style==1)THEN
	!type 1--start-------------------------------------------------------------------------

        IF (iopen==0)THEN
           CALL adios_define_var (m_adios_group , "lnx"  ,"", is_integer,"","","",varid);
           CALL adios_define_var (m_adios_group , "lny"  ,"", is_integer,"","","",varid);
           CALL adios_define_var (m_adios_group , "lnz"  ,"", is_integer,"","","",varid);
           CALL adios_define_var (m_adios_group , "write_step"  ,"", is_integer,"","","",varid);

           CALL adios_define_var (m_adios_group , TRIM(dimXname),"", is_integer,"","","",varid);
           CALL adios_define_var (m_adios_group , TRIM(dimYname),"", is_integer,"","","",varid);
           CALL adios_define_var (m_adios_group , TRIM(dimZname),"", is_integer,"","","",varid);

           CALL adios_define_var (m_adios_group , "offx"  ,"", is_integer,"","","",varid);
           CALL adios_define_var (m_adios_group , "offy"  ,"", is_integer,"","","",varid);
           CALL adios_define_var (m_adios_group , "offz"  ,"", is_integer,"","","",varid);
           CALL adios_define_var (m_adios_group , "offt"  ,"", is_integer,"","","",varid);

           CALL adios_define_var (m_adios_group , TRIM(varname)  ,"", var_type, &
                "lnx,lny,lnz,write_step",&
                TRIM(dimXname)//','//TRIM(dimYname)//','//TRIM(dimZname)//','//"write_step", &
                "offx,offy,offz,offt",varid);
        ENDIF

        CALL  adios_open (adios_handle,TRIM(grpname),TRIM(filename),TRIM(rwtemp),comm2,adios_err)

        lnx=COUNT(1)
        lny=COUNT(2)
        lnz=COUNT(3)
        offx=start(1)-1
        offy=start(2)-1
        offz=start(3)-1
        offt=0

        adios_groupsize = 11 * int_size + var_size * (lnx*lny*lnz*write_step)

        CALL adios_group_size (adios_handle, adios_groupsize, adios_totalsize, adios_err)

        CALL adios_write (adios_handle, "lnx",lnx, adios_err)
        CALL adios_write (adios_handle, "lny",lny, adios_err)
        CALL adios_write (adios_handle, "lnz",lnz, adios_err)
        CALL adios_write (adios_handle, "write_step",write_step, adios_err)

        CALL adios_write (adios_handle,TRIM(dimXname),nx, adios_err)
        CALL adios_write (adios_handle,TRIM(dimYname),ny, adios_err)
        CALL adios_write (adios_handle,TRIM(dimZname),nz, adios_err)

        CALL adios_write (adios_handle, "offx",offx, adios_err)
        CALL adios_write (adios_handle, "offy",offy, adios_err)
        CALL adios_write (adios_handle, "offz",offz, adios_err)
        CALL adios_write (adios_handle, "offt",offt, adios_err)



	!type 1--end-------------------------------------------------------------------------

     ELSE

	!type 2--start ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        ! point-to-point write --nx--

        IF (iopen==0)THEN

           CALL adios_define_var (m_adios_group , "write_step"  ,"", is_integer,"","","",varid);
           CALL adios_define_var (m_adios_group , "Ntotal"  ,"", is_integer,"","","",varid);

           DO istep=1,write_step
              DO i=1,nprec
                 CALL adios_define_var (m_adios_group , "offx"  ,"", is_integer,"","","",varid);
                 CALL adios_define_var (m_adios_group , "offt"  ,"", is_integer,"","","",varid);

                 CALL adios_define_var (m_adios_group , TRIM(varname)  ,"", var_type, &
                      "1,1", "Ntotal,write_step",  "offx,offt",varid)
              ENDDO
           ENDDO

        ENDIF

        CALL  adios_open (adios_handle,TRIM(grpname),TRIM(filename),TRIM(rwtemp),comm2,adios_err)

        !! 2 integers: Ntotal + write_step; 2 offsets
        adios_groupsize = 2 * int_size &
             + (2 * int_size + var_size) * nprec * (write_step)

        CALL adios_group_size (adios_handle, adios_groupsize, adios_totalsize,adios_err)

        CALL adios_write (adios_handle, "Ntotal",nx, adios_err)
        CALL adios_write (adios_handle, "write_step",write_step, adios_err)


	!type 2--end ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

     ENDIF

     iopen=1

  ELSE
!!!!!!reader----------
     !parameters: a series of name=vales pairs separated by ";"

     method=ADIOS_READ_METHOD_BP
     parameters="verbose=3"
     CALL adios_read_init_method(method,comm2,parameters,adios_err)

     IF (adios_err .NE. 0) THEN
        WRITE(6,*) "Error: in a_opnefile() adios_read_init_method()"
     ENDIF

    ! CALL adios_read_open(f,TRIM(filename),method,comm,ADIOS_LOCKMODE_NONE, 1.0,adios_err)
     call adios_read_open_file (adios_handle, filename, method, comm, adios_err);
     IF (adios_err .NE. 0) THEN
        WRITE(6,*) "Error: in a_opnefile() adios_read()"
     ENDIF

     !CALL adios_selection_writeblock (sel, rank)

    ! CALL adios_schedule_read (adios_handle,sel,'write_step',0,1,Y,adios_err)
    ! CALL adios_perform_reads (adios_handle, adios_err)


  ENDIF


  RETURN
END SUBROUTINE a_openfile
