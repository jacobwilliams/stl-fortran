!********************************************************************************
!>
!  Unit tests for STL module.

    program main

    use,intrinsic :: iso_fortran_env, only: wp => real64
    use stl_module

    implicit none

    type(stl_file) :: model
    integer :: istat
    integer :: i
    real(wp),dimension(3) :: x,y,z

    ! add some spheres:
    call model%add_sphere([0.0_wp, 0.0_wp, 0.0_wp], 1.0_wp, 20, 40)
    call model%add_sphere([0.0_wp, 2.0_wp, 0.0_wp], 0.3_wp, 20, 40)
    call model%add_sphere([0.0_wp, 0.0_wp, 2.0_wp], 0.1_wp, 20, 40)
    call model%write_binary_stl_file('spheres.stl',istat)
    write(*,*) istat
    call model%destroy()

    ! spheres with different resolutions:
    do i = 1, 20
        call model%add_sphere([real(i,wp), 0.0_wp, 0.0_wp], 0.5_wp, i, i*2)
    end do
    call model%write_binary_stl_file('spheres2.stl',istat)
    write(*,*) istat
    call model%destroy()

    ! cylinders:
    call model%add_cylinder([0.0_wp,0.0_wp,0.0_wp],&
                            [0.0_wp,0.0_wp,1.0_wp],&
                             0.5_wp,10,initial_cap=.true.,final_cap=.true.)
    call model%add_cylinder([1.0_wp,0.0_wp,0.0_wp],&
                             [0.0_wp,0.0_wp,2.0_wp],&
                              0.4_wp,50,initial_cap=.true.,final_cap=.true.,&
                             initial_normal=[0.0_wp,0.0_wp,-1.0_wp])
     call model%add_cylinder([0.0_wp,3.0_wp,0.0_wp],&
                            [3.0_wp,3.0_wp,3.0_wp],&
                             0.2_wp,10,initial_cap=.true.,final_cap=.true.)
    call model%write_binary_stl_file('cylinders.stl',istat)
    write(*,*) istat
    call model%destroy()

    ! curves:
    x = [ 0.0_wp, 1.0_wp, 2.0_wp ]
    y = [ 0.0_wp, 1.0_wp, 3.0_wp ]
    z = [ 0.0_wp, 1.0_wp, 4.0_wp ]
    call model%add_curve(x,y,z,0.1_wp,10,&
                        initial_cap=.true.,final_cap=.true.)
    call model%write_binary_stl_file('curve.stl',istat)
    write(*,*) istat
    call model%destroy()


    end program main
!********************************************************************************
