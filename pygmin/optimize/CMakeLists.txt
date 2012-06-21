set(utils2 "mylbfgs_updatestep")
add_custom_target(${utils2} ALL DEPENDS ${utils2}.so)
add_custom_command( 
  OUTPUT ${utils2}.so
  DEPENDS ${utils2}.f90
  COMMAND f2py -c -m ${utils2} --fcompiler=gfortran --link-lapack_opt ${utils2}.f90 > ${utils2}.setup
)

set(utils3 "mylbfgs_fort")
add_custom_target(${utils3} ALL DEPENDS ${utils3}.so)
add_custom_command( 
  OUTPUT ${utils3}.so
  DEPENDS ${utils3}.f90
  COMMAND f2py -c --fcompiler=gfortran --link-lapack_opt ${utils3}.pyf ${utils3}.f90 > ${utils3}.setup
)
