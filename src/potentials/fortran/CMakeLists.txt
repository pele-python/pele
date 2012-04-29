set(utils "ljpshiftfort")

add_custom_command( 
  OUTPUT ${utils}.so
  DEPENDS ${utils}.f
  COMMAND f2py -c -m ${utils} --fcompiler=gfortran --link-lapack_opt ${utils}.f > ${utils}.setup
)

add_custom_target(${utils} ALL DEPENDS ${utils}.so)

set(utils2 "rmdrvt")

add_custom_command( 
  OUTPUT ${utils2}.so
  DEPENDS ${utils2}.f90
  COMMAND f2py -c -m ${utils2} --fcompiler=gfortran --link-lapack_opt ${utils2}.f90 > ${utils2}.setup
)

add_custom_target(${utils2} ALL DEPENDS ${utils2}.so)
