module sandbox_wrapper
    use sandbox_commons
    use sandbox_module

    implicit none

    contains



    subroutine input(fname)
    implicit none
    character(len=200) :: fname
    integer myun
    myun = 77
    OPEN(UNIT=myun,FILE=fname,STATUS='UNKNOWN')
    call sandbox_input(myun)
    close(myun)
    end subroutine input

end module sandbox_wrapper

