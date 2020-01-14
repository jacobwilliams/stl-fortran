!********************************************************************************
!>
!  Unit tests for STL module.

    program main

    use,intrinsic :: iso_fortran_env, only: wp => real64
    use stl_module

    implicit none

    type(stl_file) :: model
    integer :: istat

    ! add some spheres:
    call model%add_sphere([0.0_wp, 0.0_wp, 0.0_wp], 1.0_wp, 20, 40)
    call model%add_sphere([0.0_wp, 2.0_wp, 0.0_wp], 0.3_wp, 20, 40)
    call model%add_sphere([0.0_wp, 0.0_wp, 2.0_wp], 0.1_wp, 20, 40)
    call model%write_binary_stl_file('spheres.stl',istat)

    write(*,*) istat

    end program main
!********************************************************************************
