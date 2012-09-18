set(utils1 "_fortran_utils")
add_custom_command( 
  OUTPUT ${utils1}.so
  DEPENDS ${utils1}.f90
  COMMAND f2py -c -m ${utils1} --fcompiler=gfortran --link-lapack_opt ${utils1}.f90 > ${utils1}.setup
)
add_custom_target(${utils1} ALL DEPENDS ${utils1}.so)
