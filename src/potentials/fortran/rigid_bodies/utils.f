C   GMIN: A program for finding global minima
C   Copyright (C) 1999- David J. Wales
C   This file is part of GMIN.
C
C   GMIN is free software; you can redistribute it and/or modify
C   it under the terms of the GNU General Public License as published by
C   the Free Software Foundation; either version 2 of the License, or
C   (at your option) any later version.
C
C   GMIN is distributed in the hope that it will be useful,
C   but WITHOUT ANY WARRANTY; without even the implied warranty of
C   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C   GNU General Public License for more details.
C
C   You should have received a copy of the GNU General Public License
C   along with this program; if not, write to the Free Software
C   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
C

      INTEGER FUNCTION GETUNIT()
      IMPLICIT NONE
      LOGICAL :: INUSE
!
! start checking for available units > 10, to avoid system default units
! 
      INTEGER :: UNITNUM

      INUSE=.TRUE.
      UNITNUM=11

      DO WHILE (INUSE)
         INQUIRE(UNIT=UNITNUM,OPENED=INUSE)
         IF (.NOT.INUSE) THEN
            GETUNIT=UNITNUM 
         ELSE     
            UNITNUM=UNITNUM+1
         ENDIF
      ENDDO
      END FUNCTION GETUNIT


