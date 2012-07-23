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

set(utils3 "soft_sphere_pot")

add_custom_command( 
  OUTPUT ${utils3}.so
  DEPENDS ${utils3}.f90
  COMMAND f2py -c -m ${utils3} --fcompiler=gfortran --link-lapack_opt ${utils3}.f90 > ${utils3}.setup
)

add_custom_target(${utils3} ALL DEPENDS ${utils3}.so)

set(utils4 "lj")

add_custom_command( 
  OUTPUT ${utils4}.so
  DEPENDS ${utils4}.f90
  COMMAND f2py -c -m ${utils4} --fcompiler=gfortran --link-lapack_opt ${utils4}.f90 > ${utils4}.setup
)

add_custom_target(${utils4} ALL DEPENDS ${utils4}.so)

set(utils5 "AT")
add_custom_command( 
  OUTPUT ${utils5}.so
  DEPENDS ${utils5}.f
  COMMAND f2py -c -m ${utils5} --fcompiler=gfortran --link-lapack_opt ${utils5}.f > ${utils5}.setup
)
add_custom_target(${utils5} ALL DEPENDS ${utils5}.so)

set(utils6 "lj_hess")
add_custom_command( 
  OUTPUT ${utils6}.so
  DEPENDS ${utils6}.f
  COMMAND f2py -c -m ${utils6} --fcompiler=gfortran --f77flags="-ffixed-line-length-none" --link-lapack_opt ${utils6}.f > ${utils6}.setup
)
add_custom_target(${utils6} ALL DEPENDS ${utils6}.so)
