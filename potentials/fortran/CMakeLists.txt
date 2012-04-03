set(utils "ljpshiftfort")

add_custom_command(
  OUTPUT ${utils}.so
  DEPENDS ${utils}.f
  COMMAND f2py -c -m ${utils} --fcompiler=gfortran --link-lapack_opt ${utils}.f > ${utils}.setup
)

add_custom_target(${utils} ALL DEPENDS ${utils}.so)