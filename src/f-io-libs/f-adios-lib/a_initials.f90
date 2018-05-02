!**************************
!@section LICENSE
!Copyright (c) 2013-2017, Regents of the University of California
!All rights reserved.
!Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
!1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
!2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
!THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!**************************
SUBROUTINE a_initials(write_style,rank,nk,kout,nj,jout,ni,iout,kbnd,jbnd,ibnd,recproc,nprec,start2dim) bind(C)
  IMPLICIT NONE

  INTEGER nk,nj,ni,recproc,nprec,kbnd(2),jbnd(2),ibnd(2),rank,start2dim
  INTEGER kout(nk),jout(nj),iout(ni)
  INTEGER k,j,i,kloc,jloc,iloc,write_style

  nprec=0
  recproc=-99

  IF (write_style==1)THEN

     DO k=1,nk
        kloc=kout(k)
        IF (kloc <kbnd(1) .OR. kloc >kbnd(2)) CYCLE

        DO j=1,nj
           jloc=jout(j)
           IF (jloc <jbnd(1) .OR. jloc >jbnd(2)) CYCLE

           DO i=1,ni
              iloc=iout(i)
              IF (iloc <ibnd(1) .OR. iloc >ibnd(2)) CYCLE
              nprec = nprec + 1
              recproc=rank

           ENDDO
        ENDDO
     ENDDO
  ELSE


     WRITE(*,*)' NOW implemented. '
     !	stop

     DO k=1,nk
        kloc=kout(k)
        IF (kloc <kbnd(1) .OR. kloc >kbnd(2)) CYCLE

        j=k
        jloc=jout(j)
        IF (jloc <jbnd(1) .OR. jloc >jbnd(2)) CYCLE
        i=k
        iloc=iout(i)
        IF (iloc <ibnd(1) .OR. iloc >ibnd(2)) CYCLE

        nprec = nprec + 1
        recproc=rank
     ENDDO

  ENDIF

  IF (write_style==1)start2dim=0
  IF (write_style==2)start2dim=nprec

  !! need this after the call for each file
  !! if (nprec.gt.0)then
  !!  allocate(prec(nprec,maxdim))
  !!      allocate(buffer(nprec*write_step))
  !!	allocate(start2(4,nprec))
  !!   declare start(4) etc.
  !!  endif

  RETURN
END SUBROUTINE a_initials
