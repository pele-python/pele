PROGRAM KEYWORD
IMPLICIT NONE
INTEGER NRESMIN , LUNIT,J1,NATOMS,J2
DOUBLE PRECISION, ALLOCATABLE :: EMIN(:), FVIBMIN(:), PFMIN(:)
DOUBLE PRECISION, ALLOCATABLE :: IXMIN(:), IYMIN(:), IZMIN(:), RESPOINTS(:,:)
DOUBLE PRECISION, ALLOCATABLE :: POINTSASCII(:,:)
INTEGER, ALLOCATABLE :: HORDERMIN(:)
LOGICAL YESNO
DOUBLE PRECISION DUMMY

CHARACTER(LEN=20) :: temp_str
INTEGER tempI , tempJ

print *, '===================================================================   '
print *, ' Program to convert '
print *, '     (a) binary points.min to ascii extractedmin, and '
print *, '     (b) binary points.ts to ascii extractedts '
print *, ' '
print *, '   Commented code can be used to read ascii files and write binary files, e.g. to export PyGmin database to '
print *, '    OPTIM/PATHSAMPLE format'
print *, ' '
print *, ' files: '
print *, '   min.data, points.min, ts.data, points.ts '
print *, '    (min.data and ts.data are just used to get the number of minima and ts) '
print *, ''
print *, ' Usage: '
print *, '  $ gfortram binpoints2ascii.f90 '
print *, '  $ ./a.out'
print *, ''
print *, ' Next, to create PyGMIN database: '
print *, '  $ python optim2sqlite.py  '
print *, '==================================================================='   
print *, ' '

! TOSET 
!NATOMS = 32 
PRINT *, 'Enter number of atoms ..'
READ  *, NATOMS 

LUNIT=23432 

! check if min.data and points.min exist 
         INQUIRE(FILE='min.data',EXIST=YESNO)
         IF (.NOT.YESNO) THEN
            PRINT '(A)','ERROR *** min.data not found'
            STOP
         ENDIF
         INQUIRE(FILE='points.min',EXIST=YESNO)
         IF (.NOT.YESNO) THEN
            PRINT '(A)','ERROR *** points.min not found'
            STOP
         ENDIF

! -- READING min.data  
         OPEN(UNIT=LUNIT,FILE='min.data',STATUS='OLD')
         ! find number of lines in min.data 
         NRESMIN=0
         DO
            READ(LUNIT,*,END=30) DUMMY
            NRESMIN=NRESMIN+1
         ENDDO
30       WRITE(*,'(A,I6)'), 'Number of lines in min.data = ',NRESMIN
         ! now read 
         ALLOCATE(EMIN(NRESMIN),FVIBMIN(NRESMIN),PFMIN(NRESMIN),IXMIN(NRESMIN),IYMIN(NRESMIN),IZMIN(NRESMIN),HORDERMIN(NRESMIN))
         ALLOCATE(RESPOINTS(3*NATOMS,NRESMIN))
         REWIND(LUNIT)
         DO J1=1,NRESMIN
            READ(LUNIT,*) EMIN(J1),FVIBMIN(J1),HORDERMIN(J1),IXMIN(J1),IYMIN(J1),IZMIN(J1)
         ENDDO
         CLOSE(LUNIT)

! -- READ binary points.min 
         OPEN(LUNIT,FILE='points.min',ACCESS='DIRECT',FORM='UNFORMATTED',STATUS='OLD',RECL=8*3*NATOMS)
         DO J1=1,NRESMIN
            READ(LUNIT,REC=J1) (RESPOINTS(J2,J1),J2=1,3*NATOMS)
         ENDDO
         CLOSE(LUNIT)
         WRITE(*,'(A,I6,A)') 'points for ',NRESMIN,' minima read from file points.min'

! -- WRITE points in ascii format with NRESMIN x NATOMS rows and three columns  
         OPEN(UNIT=23432,FILE='extractedmin',STATUS='UNKNOWN')
         WRITE(temp_str,'(A,I4,A)')'(',3,'G20.10)'  !  temp_str = (   3G20.10)        
         WRITE(23432 , temp_str ) RESPOINTS(:,:)
        ! WRITE(23432 , temp_str ) POINTSASCII(:,:)
         CLOSE(23432)

         PRINT *, 'extractedmin created.'
         PRINT *, ''
! ========== Transition States ==================
	DEALLOCATE(RESPOINTS) 

! check if ts.data and points.ts exist 
         INQUIRE(FILE='ts.data',EXIST=YESNO)
         IF (.NOT.YESNO) THEN
            PRINT '(A)','ERROR *** ts.data not found'
            STOP
         ENDIF
         INQUIRE(FILE='points.ts',EXIST=YESNO)
         IF (.NOT.YESNO) THEN
            PRINT '(A)','ERROR *** points.ts not found'
            STOP
         ENDIF

! -- READING min.data  
         OPEN(UNIT=LUNIT,FILE='ts.data',STATUS='OLD')
         ! find number of lines in min.data 
         NRESMIN=0
         DO
            READ(LUNIT,*,END=40) DUMMY
            NRESMIN=NRESMIN+1
         ENDDO
40       WRITE(*,'(A,I6)'), 'Number of lines in ts.data = ',NRESMIN
         ! now read 
!         ALLOCATE(EMIN(NRESMIN),FVIBMIN(NRESMIN),PFMIN(NRESMIN),IXMIN(NRESMIN),IYMIN(NRESMIN),IZMIN(NRESMIN),HORDERMIN(NRESMIN))
         ALLOCATE(RESPOINTS(3*NATOMS,NRESMIN))
!         REWIND(LUNIT)
!         DO J1=1,NRESMIN
!            READ(LUNIT,*) EMIN(J1),FVIBMIN(J1),HORDERMIN(J1),IXMIN(J1),IYMIN(J1),IZMIN(J1)
!         ENDDO
         CLOSE(LUNIT)

! -- READ binary points.min 
         OPEN(LUNIT,FILE='points.ts',ACCESS='DIRECT',FORM='UNFORMATTED',STATUS='OLD',RECL=8*3*NATOMS)
         DO J1=1,NRESMIN
            READ(LUNIT,REC=J1) (RESPOINTS(J2,J1),J2=1,3*NATOMS)
         ENDDO
         CLOSE(LUNIT)
         WRITE(*,'(A,I6,A)') 'points for ',NRESMIN,' ts read from file points.ts'


! -- WRITE points in ascii format with NRESMIN x NATOMS rows and three columns  
         OPEN(UNIT=23432,FILE='extractedts',STATUS='UNKNOWN')
         WRITE(temp_str,'(A,I4,A)')'(',3,'G20.10)'  !  temp_str = (   3G20.10)        
         WRITE(23432 , temp_str ) RESPOINTS(:,:)
        ! WRITE(23432 , temp_str ) POINTSASCII(:,:)
         CLOSE(23432)

         PRINT *, 'extractedts created.'
END PROGRAM 



! -- READ ascii points.min 
!         ALLOCATE(POINTSASCII(3*NATOMS,NRESMIN))
!         REWIND(LUNIT)
!         OPEN(UNIT=LUNIT,FILE='points.min.ascii',STATUS='OLD') 
!         DO J1=1,NRESMIN
!           READ(LUNIT,*) ( POINTSASCII(J2,J1),J2=1,3*NATOMS)  
!         ENDDO
!         CLOSE(LUNIT)
!         WRITE(*,'(A,I6,A)') 'points for ',NRESMIN,' minima read from file points.min.ascii'


! -- WRITE points in binary format  
!         PRINT '(A)', 'Writing points in binary'
!         LUNIT=23432
!         OPEN(UNIT=LUNIT,FILE='points.min.new.bin',ACCESS='DIRECT',FORM='UNFORMATTED',STATUS='UNKNOWN',RECL=8*3*NATOMS)
!         DO J1=1,NRESMIN
!            WRITE(LUNIT,REC=J1) (RESPOINTS(J2,J1),J2=1,3*NATOMS)
!         ENDDO 
!         CLOSE(23432)

! -- WRITE min.data.new
!         PRINT '(A)', 'Writing min.data.new'
!         LUNIT=23432
!         OPEN(UNIT=LUNIT,FILE='min.data.new',STATUS='UNKNOWN')
!         DO J1=1,NRESMIN
!            WRITE(LUNIT,'(2F20.10,I6,3F20.10)') EMIN(J1), FVIBMIN(J1), HORDERMIN(J1), IXMIN(J1), IYMIN(J1), IZMIN(J1)
!         ENDDO
!         CLOSE(23432)

