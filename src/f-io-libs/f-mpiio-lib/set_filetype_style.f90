!**************************
!@section LICENSE
!Copyright (c) 2013-2017, Regents of the University of California
!All rights reserved.
!Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
!1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
!2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
!THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!**************************
!set nprec,prec
SUBROUTINE set_filetype_style(ghostx,ghosty,ghostz,write_style,nk,kout,nj,jout,ni,iout,rank, maxdim, &
     kbnd,jbnd,ibnd, datatype,write_step,&
     nprec, prec, nrec,filetype,recproc)
  IMPLICIT NONE
  INCLUDE 'mpif.h'

  INTEGER,DIMENSION(:), ALLOCATABLE:: blocksize
  INTEGER(KIND=MPI_ADDRESS_KIND),DIMENSION(:), ALLOCATABLE::tpmap,pmap

  INTEGER err, filetype,write_step,nk,nj,ni,write_style
  INTEGER datatype,j1,j2,i1,i2
  INTEGER ghostx,ghosty,ghostz
  INTEGER k,j,i,kloc,jloc,iloc,nprec,maxdim,recproc,rank,nptotal
  INTEGER, DIMENSION(2):: kbnd,jbnd,ibnd
  INTEGER prec(nprec,3),kout(nk),jout(nj),iout(ni)
  INTEGER(KIND=MPI_OFFSET_KIND):: nrec
  INTEGER nbyte

  ! global values kout, kbnd(1) , etc
  IF (write_style ==1 )THEN
     nrec=nk*nj*INT(ni,MPI_OFFSET_KIND)
     !write(*,*)'In setup filetype 1:',MPI_OFFSET_KIND, nrec, nk,nj,ni
  ELSE
     IF (nk .NE.nj .OR. nk.NE.ni)THEN
        WRITE(*,*)' Wrong parameters. return'
        RETURN
     ENDIF
     nrec=nk
  ENDIF

  !
  !!just geto the number of receivers at each processor
  ALLOCATE(pmap(nprec))
  nprec=0
  DO k=1,nk
     kloc=kout(k)
     IF (kloc <kbnd(1) .OR. kloc >kbnd(2)) CYCLE
     IF (write_style==1)THEN
        j1=1
        j2=nj
     ELSE
        j1=k
        j2=k
     ENDIF
     DO j=j1,j2
        jloc=jout(j)
        IF (jloc <jbnd(1) .OR. jloc >jbnd(2)) CYCLE

        IF (write_style==1)THEN
           i1=1
           i2=ni
        ELSE
           i1=k
           i2=k
        ENDIF

        DO i=i1,i2
           iloc=iout(i)
           IF (iloc <ibnd(1) .OR. iloc >ibnd(2)) CYCLE

           nprec = nprec + 1

           !for overflow issure modified by Hui 2016.07
           IF (write_style==1)THEN
              pmap(nprec) = INT(nj,MPI_ADDRESS_KIND)*ni*(k-1) + INT(ni,MPI_ADDRESS_KIND)*(j-1) + INT((i-1),MPI_ADDRESS_KIND)
           ELSE
              pmap(nprec) =  i-1
           ENDIF

           ! MODIFIED TO SUPPORT GHOST CELLS
           prec(nprec,1) = iloc - ibnd(1) + ghostx +1
           prec(nprec,2) = jloc - jbnd(1) + ghosty +1
           prec(nprec,3) = kloc - kbnd(1) + ghostz +1
           !fixed 2 ghost cells
           !prec(nprec,1) = iloc - ibnd(1) +3
           !prec(nprec,2) = jloc - jbnd(1) +3
           !prec(nprec,3) = kloc - kbnd(1) +3
           !No ghost cells
           !prec(nprec,1) = iloc - ibnd(1) +1
           !prec(nprec,2) = jloc - jbnd(1) +1
           !prec(nprec,3) = kloc - kbnd(1) +1
        ENDDO !i
     ENDDO !j
  ENDDO !k

  !write(*,*)'nprec =', nprec, kbnd,jbnd,ibnd,rank,write_step

  !!
  !! set filetype
  nptotal=nprec*write_step
  ALLOCATE(tpmap(nptotal))
  ALLOCATE(blocksize(nptotal))

  DO i=1,write_step
     DO j=1,nprec
        tpmap((i-1)*nprec + j)= pmap(j) + INT((i-1)*nrec,MPI_ADDRESS_KIND)
     ENDDO
  ENDDO
  DEALLOCATE(pmap)

  blocksize=1

  CALL MPI_TYPE_SIZE(datatype,nbyte,err)
  tpmap=tpmap*nbyte         ! byte addressing instead of element addressing
  !write(*,*)' in type:', rank, tpmap,nbyte


  !write(*,*)' set type : rank ', rank, nprec, tpmap, blocksize,nptotal,datatype==MPI_REAL
  CALL MPI_TYPE_CREATE_HINDEXED(nptotal,blocksize,tpmap,datatype,filetype,err)
  CALL MPI_TYPE_COMMIT(filetype,err)

  DEALLOCATE(tpmap)
  DEALLOCATE(blocksize)

  RETURN
END SUBROUTINE set_filetype_style
