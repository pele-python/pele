set(utils2 "mylbfgs_updatestep")
add_custom_target(${utils2} ALL DEPENDS ${utils2}.so)
add_custom_command( 
  OUTPUT ${utils2}.so
  DEPENDS ${utils2}.f90
  COMMAND f2py -c -m ${utils2} --fcompiler=gfortran --link-lapack_opt ${utils2}.f90 > ${utils2}.setup
)
