set(utils "overlap")

add_custom_command(
  OUTPUT ${utils}.so
  DEPENDS ${utils}.f90
  COMMAND f2py -c -m ${utils} --fcompiler=gfortran --link-lapack_opt ${utils}.f90 > ${utils}.setup
)

add_custom_target(${utils} ALL DEPENDS ${utils}.so)
