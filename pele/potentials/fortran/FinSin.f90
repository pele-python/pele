!   GMIN: A program for finding global minima
!   Copyright (C) 1999-2006 David J. Wales
!   This file is part of GMIN.
!   Finnis-Sinclair potential added by J.A. Elliott 2009
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
!**************************************************************
!
!  Here we calculate the Finnis-Sinclair potential and gradient
!  Originally implemented by James Elliot in 2009 for GMIN.
!  Tidied up for PELE by Dmitri Schebarchov in 2014.
!  (Intended for clusters and not periodic systems...)
!
!**************************************************************
!
SUBROUTINE FinSin(NATOMS,X,V,d,A,beta,c,c0,c1,c2,POT)
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: NATOMS
  DOUBLE PRECISION, INTENT(IN) :: X(3*NATOMS)
  DOUBLE PRECISION, INTENT(IN) :: d, A, beta, c, c0, c1, c2
  DOUBLE PRECISION, INTENT(OUT) :: V(3*NATOMS), POT
  !
  INTEGER :: J1, J2
  DOUBLE PRECISION :: DIST, PTEMP, fsrho(NATOMS), VTEMP, RTEMP, &
       DUMMY, VTEMP1(NATOMS), VTEMP2(NATOMS), VTEMP3(NATOMS), &
       dx, dy, dz, min_dist, cutp, cutr1, cutr2, temp_dist(natoms,natoms)
  !
  ! Finnis-Sinclar potential parameters
  ! Phil Mag A 50, 45 (1984)
  ! parameters for Fe subsequently modified in Erratum
  ! Phil Mag A 53, 161 (1986)
  ! (both sets for Fe are included)
  !
  ! bcc metals
  !
  ! d = density cut-off [Angstroms]
  ! A = binding energy [eV]
  ! beta = keeps phi within first-nearest-neighbor distance.
  ! c = two-body interaction cut-off [Angstroms]
  ! c0,c1,c2 = fitting parameters
  !
  ! If beta is non-zero, set a minimum distance cut-off to 
  ! avoid rtemp going negative
  if (beta.ne.0d0) then
     min_dist=d*(beta-1.0d0)/beta
  else
     min_dist=0.0d0
  endif
  !
  ! Calculate the potential and gradient terms
  !
  ! First, populate density and distance matrices
  !
  do j2=1,natoms ! initialise distance and density arrays
     do j1=1,natoms
        temp_dist(j1,j2)=0.0d0
     enddo
     fsrho(j2)=0.0d0
  enddo
  !
  rtemp = 0.0d0 ! initialise energy accumulators
  ptemp = 0.0d0
  !
  do j1=1,natoms-1 ! begin outer loop over all atoms except last
     do j2=j1+1,natoms ! begin first inner loop over all j2>j1
        !
        temp_dist(j2,j1)=  & ! distance between atoms j1,j2
             dsqrt( ( X(3*(J2-1)+1)-X(3*(J1-1)+1) )**2 + &
             ( X(3*(J2-1)+2)-X(3*(J1-1)+2) )**2 + &
             ( X(3*(J2-1)+3)-X(3*(J1-1)+3) )**2 )
        !
        dist=temp_dist(j2,j1) ! store for use in the inner loop
        temp_dist(j1,j2)=dist ! impose symmetry on distance matrix
        !
        cutp = dsign(0.5d0,c-dist)+0.5d0 ! cut-off for 2-body term
        ptemp = ptemp + & ! accumulate two-body potential term
             cutp*(((dist-c)**2)*(c0+dist*(c1+c2*dist)))           
        !
        cutr1=dsign(0.5d0,d-dist)+0.5d0 ! cut-off for many-body term
        cutr2=dsign(0.5d0,dist-min_dist)+0.5d0 ! min. distance cut-off
        rtemp=cutr1*cutr2*((dist-d)**2+(beta*(dist-d)**3/d))
        fsrho(j1)=fsrho(j1)+rtemp ! contribution to density of j1
        fsrho(j2)=fsrho(j2)+rtemp ! contribution to density of j2
        !
     enddo ! end inner loop over all j2>j
  enddo ! end outer loop over all atoms except last
  !
  ! Now, sum the potential energy over all atoms
  !
  POT=ptemp ! initialise potential energy
  do j1=1,natoms
     fsrho(j1)=dsqrt(fsrho(j1)) ! square root density
     POT=POT-a*fsrho(j1) ! accumulate potential energy
  enddo
  !
  ! Calculate gradient terms, if required
  !
  do j1=1,natoms ! initialise total gradient terms
     V(3*(J1-1)+1)=0.d0
     V(3*(J1-1)+2)=0.d0
     V(3*(J1-1)+3)=0.d0
     vtemp1(J1)=0.0d0
     vtemp2(J1)=0.0d0
     vtemp3(J1)=0.0d0
  enddo
  do j1=1,natoms-1 ! begin outer loop over all atoms except last
     !
     dummy=1.0d0/fsrho(j1) ! store reciprocal of density for j1
     !
     do j2=j1+1,natoms ! begin inner loop over all j2>j1
        !
        dist=temp_dist(j1,j2) ! recall distance from earlier loop
        cutp=dsign(0.5d0,c-dist)+0.5d0 ! 2-body cut-off
        vtemp= & ! accumulate 2-body gradient term
             cutp*(2.0d0*(dist-c)*(c0+dist*(c1+c2*dist)) &
             +((dist-c)**2)*(c1+2.0d0*c2*dist))/dist
        cutr1=dsign(0.5d0,d-dist)+0.5d0 ! many-body cut-off
        cutr2=dsign(0.5d0,dist-min_dist)+0.5d0 ! min. dist. cut-off
        vtemp=vtemp & ! accumulate many-body gradient term
             -cutr1*cutr2*((0.5d0*A)*(dummy+1.0d0/fsrho(j2))* &
             (2.0d0*(dist-d)+3.0d0*(beta/d)*(dist-d)**2))/dist
        !
        dx=(X(3*(J1-1)+1)-X(3*(J2-1)+1)) ! calculate Rij components
        dy=(X(3*(J1-1)+2)-X(3*(J2-1)+2))
        dz=(X(3*(J1-1)+3)-X(3*(J2-1)+3))
        !
        vtemp1(j1)=vtemp1(j1)+vtemp*dx ! primary gradient components
        vtemp2(j1)=vtemp2(j1)+vtemp*dy
        vtemp3(j1)=vtemp3(j1)+vtemp*dz
        !
        vtemp1(j2)=vtemp1(j2)-vtemp*dx ! symm. gradient components
        vtemp2(j2)=vtemp2(j2)-vtemp*dy
        vtemp3(j2)=vtemp3(j2)-vtemp*dz
        !
     enddo
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
END SUBROUTINE FinSin
