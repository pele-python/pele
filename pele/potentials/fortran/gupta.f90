!   GMIN: A program for finding global minima
!   Copyright (C) 1999-2006 David J. Wales
!   This file is part of GMIN.
!   Loop structure recoded by J.A. Elliott 2009
!
!   GMIN is free software; you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation; either version 2 of the License, or
!   (at your option) any later version.
!
!   GMIN is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with this program; if not, write to the Free Software
!   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!
!*****************************************************************
!
!  Here we calculate the Gupta potential and gradient
!  Originally implemented by James Elliot in 2009 for GMIN.
!  Tidied up for PELE by Dmitri Schebarchov in 2014.
!  (Intended for clusters and not periodic systems...)
!  See Cleri and Rosato, PRB 48, 22 (1993), for parameters.
!
!*****************************************************************
!
SUBROUTINE GUPTA(NATOMS,X,V,p,q,GA,Gxi,PG)
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: NATOMS ! atom-count
  DOUBLE PRECISION, INTENT(IN) :: X(3*NATOMS) ! coordinates
  DOUBLE PRECISION, INTENT(IN) :: p, q, GA, Gxi ! model parameters
  DOUBLE PRECISION, INTENT(OUT) :: V(3*NATOMS), PG ! gradient and E
  !
  INTEGER :: J1, J2
  DOUBLE PRECISION :: DIST,PTEMP,RTEMP,twoq,RytoeV,dx,dy,dz, &
       GRHO(NATOMS),VTEMP,DUMMY,temp_dist(natoms,natoms), &
       VTEMP1(NATOMS),VTEMP2(NATOMS),VTEMP3(NATOMS),r0,GA2
  !  
  ! Put r0 in reduced units, and adjust some constant
  r0=1.0d0
  GA2=2.0d0*GA 
  ! Note that in orignal Cleri & Rosato paper, they define
  ! 2-body term as a sum over all atom pairs, but it is quicker to
  ! evaluate sum over just j2>j1 and double the constant A
  twoq=2.0d0*q
  !
  do j2=1,natoms ! initialise distance and density arrays
     do j1=1,natoms
        temp_dist(j1,j2)=0.0d0
     enddo
     grho(j2)=0.0d0
  enddo
  !
  rtemp = 0.0d0 ! initialise energy accumulators
  ptemp = 0.0d0 !
  !
  do j1=1,natoms-1 ! begin outer loop over all atoms except last
     do j2=j1+1,natoms ! begin first inner loop over all j2>j1
        !
        temp_dist(j2,j1)= & ! store distance between atoms j1,j2
             dsqrt(( X(3*(J2-1)+1)-X(3*(J1-1)+1) )**2 + &
             ( X(3*(J2-1)+2)-X(3*(J1-1)+2) )**2 + &
             ( X(3*(J2-1)+3)-X(3*(J1-1)+3) )**2)
        !
        dist=temp_dist(j2,j1) ! store the distance used in inner loop
        temp_dist(j1,j2)=dist ! impose symmetry on distance matrix
        !
        dist=dist/r0
        PTEMP=PTEMP+dexp(p*(1-dist)) ! accumulate two-body term
        RTEMP=dexp(twoq*(1-dist)) ! calculate many-body term
        !
        grho(j1)=grho(j1)+rtemp ! accumulate density
        grho(j2)=grho(j2)+rtemp ! accumulate density
        !
     enddo ! end inner loop over all j2>j1
     !
  enddo ! end outer loop over all atoms except last
  !
  ! Now, sum the potential energy over all atoms
  !
  pg=GA2*ptemp ! initialise potential energy
  !
  do j1=1,natoms
     grho(j1)=dsqrt(grho(j1)) ! square root density
     pg=pg-Gxi*grho(j1) ! accumulate potential energy
  enddo
  !
  ! Calculate gradient terms
  !
  do j1=1,natoms ! initialise total gradient terms
     V(3*(J1-1)+1)=0.d0
     V(3*(J1-1)+2)=0.d0
     V(3*(J1-1)+3)=0.d0
     vtemp1(J1)=0.0d0
     vtemp2(J1)=0.0d0
     vtemp3(J1)=0.0d0
  enddo
  !
  do j1=1,natoms-1 ! begin outer loop over all atoms except last
     !
     dummy=1.0d0/grho(j1) ! store reciprocal of density for atom j1
     !
     do j2=j1+1,natoms ! begin inner loop over all j2>j1
        !
        dist=temp_dist(j1,j2) ! recall distance from earlier loop
        dist=dist/r0
        !
        ! calculate gradient term
        VTEMP=(q*Gxi*(DUMMY+1.0D0/GRHO(J2))*dexp(twoq*(1-dist)) &
             -GA2*p*dexp(p*(1-dist)))/(r0*dist)
        !
        ! calculate Cartesian components of distance
        dx=(X(3*(J1-1)+1)-X(3*(J2-1)+1))         
        dy=(X(3*(J1-1)+2)-X(3*(J2-1)+2))
        dz=(X(3*(J1-1)+3)-X(3*(J2-1)+3))
        !
        ! accumulate primary gradient components
        vtemp1(j1)=vtemp1(j1)+vtemp*dx           
        vtemp2(j1)=vtemp2(j1)+vtemp*dy
        vtemp3(j1)=vtemp3(j1)+vtemp*dz
        !
        ! accumulate symmetric gradient components
        vtemp1(j2)=vtemp1(j2)-vtemp*dx
        vtemp2(j2)=vtemp2(j2)-vtemp*dy
        vtemp3(j2)=vtemp3(j2)-vtemp*dz
        !
     enddo
     !
  enddo
  !
  ! Finally, sum the gradient terms over all atoms
  !
  do j1=1,natoms
     V(3*(J1-1)+1)=V(3*(J1-1)+1)+vtemp1(j1)
     V(3*(J1-1)+2)=V(3*(J1-1)+2)+vtemp2(j1)
     V(3*(J1-1)+3)=V(3*(J1-1)+3)+vtemp3(j1)
  enddo
  !
  RETURN
  !
END SUBROUTINE GUPTA
