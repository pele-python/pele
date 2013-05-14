subroutine structure_read(file_name, k, n_atoms, read_in_coords)
  implicit none
  
  ! Returns the structure of minima k as a 3*n_atoms array
  ! Compile me with;
  ! $ f2py -c -m structure_read structure_read.f90
  
  integer :: j
  integer, intent(in) :: k, n_atoms
  character(len=120) :: file_name
  double precision, dimension(3*n_atoms), intent(out) :: read_in_coords
  
  
  ! =========== Obtain Structural Data =========== !
   
   open(unit=2,file=trim(file_name),access='direct',form='unformatted',status='old',recl=8*3*n_atoms)
   write(*,*) 'f', k, 8*3*n_atoms
   read(2,rec=k) (read_in_coords(j),j=1,3*n_atoms)
  
   
   close(2)
   

 end subroutine structure_read
