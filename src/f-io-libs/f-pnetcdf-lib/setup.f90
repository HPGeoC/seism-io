!**************************
!@section LICENSE
!Copyright (c) 2013-2017, Regents of the University of California
!All rights reserved.
!Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
!1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
!2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
!THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!**************************
SUBROUTINE setup(ghostx,ghosty,ghostz,nk,kout,nj,jout,ni,iout,kbnd,jbnd,ibnd,write_style,&
     start,count,start2,start2dim,nprec,prec) bind(C)

  IMPLICIT NONE
  INCLUDE 'mpif.h'

  INTEGER ghostx,ghosty,ghostz
  INTEGER nk,nj,ni,kbnd(2),jbnd(2),ibnd(2),start2dim
  INTEGER j1,j2,i1,i2,k,j,i,kloc,jloc,iloc

  INTEGER kout(nk),jout(nj),iout(ni)
  INTEGER(kind=MPI_OFFSET_KIND) start(4), COUNT(4),start2(4,start2dim)
  INTEGER nprec,prec(nprec,3),write_style,st0,ct0
  INTEGER,DIMENSION(:),ALLOCATABLE:: pmap

  ALLOCATE(pmap(nprec))

  nprec=0
  st0=-99
  ct0=0

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

           IF (write_style .NE.1)THEN
              pmap(nprec) =  i-1
           ENDIF

           ! MODIFIED TO SUPPORT GHOST CELLS
           prec(nprec,1) = iloc - ibnd(1) + ghostx +1
           prec(nprec,2) = jloc - jbnd(1) + ghosty +1
           prec(nprec,3) = kloc - kbnd(1) + ghostz +1

        ENDDO
     ENDDO
  ENDDO

  !	if (nprec ==0) return
  IF (write_style==1)THEN

     CALL find_start_count(kout,nk,kbnd,start(3),COUNT(3));
     CALL find_start_count(jout,nj,jbnd,start(2),COUNT(2));
     CALL find_start_count(iout,ni,ibnd,start(1),COUNT(1));

     IF (COUNT(1)*COUNT(2)*COUNT(3)==0) COUNT(1:3)=0

  ELSE
     !! flag to write as style #2
     start(1:3)=0
     IF (nprec >0) start2(:,:)=0
     DO j=1,nprec
        start2(1,j)=pmap(j)+1
        start2(2,j)=1
        start2(3,j)=1
        start2(4,j)=1

     ENDDO
     COUNT(1:3)=1



  ENDIF

  RETURN
END SUBROUTINE setup
