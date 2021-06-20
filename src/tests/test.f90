!********************************************************************************
!>
!  Unit tests for STL module.

    program main

    use,intrinsic :: iso_fortran_env, only: wp => real64
    use stl_module

    implicit none

    type(stl_file) :: model
    integer :: istat !! status return code
    integer :: i !! counter
    integer :: j !! counter
    real(wp),dimension(:),allocatable :: x,y,z !! curve coordinates

    real(wp),parameter :: deg2rad = acos(-1.0_wp) / 180.0_wp !! degrees to radians

    write(*,*) 'Running tests...'

    ! add some spheres:
    call model%add_sphere([0.0_wp, 0.0_wp, 0.0_wp], 1.0_wp, 20, 40)
    call model%add_sphere([0.0_wp, 2.0_wp, 0.0_wp], 0.3_wp, 20, 40)
    call model%add_sphere([0.0_wp, 0.0_wp, 2.0_wp], 0.1_wp, 20, 40)
    call model%write_binary_stl_file('spheres.stl',istat)
    if (istat/=0) error stop 'error creating spheres'
    call model%destroy()

    ! spheres with different resolutions:
    do i = 1, 20
        call model%add_sphere([real(i,wp), 0.0_wp, 0.0_wp], 0.5_wp, i, i*2)
    end do
    call model%write_binary_stl_file('spheres2.stl',istat)
    if (istat/=0) error stop 'error creating spheres'
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
    if (istat/=0) error stop 'error creating cylinders'
    call model%destroy()

    ! curves:
    call model%add_curve(x = [ 0.0_wp, 1.0_wp, 2.0_wp, 5.0_wp ], &
                         y = [ 0.0_wp, 1.0_wp, 3.0_wp, 5.0_wp ], &
                         z = [ 0.0_wp, 1.0_wp, 4.0_wp, 5.0_wp ], &
                         radius = 0.1_wp,&
                         num_points = 10, &
                         initial_cap=.true.,final_cap=.true.)
    call model%write_binary_stl_file('curve.stl',istat)
    if (istat/=0) error stop 'error creating curve'
    call model%destroy()

    ! curves:
    allocate(x(360))
    allocate(y(360))
    allocate(z(360)); z = 0.0_wp
    do j = 1, 5
        do i = 1, 360
            x(i) = j*2.0_wp * cos(i*deg2rad) * sin(i/2.0_wp*deg2rad)
            y(i) = j*2.0_wp * sin(i*deg2rad) * sin(i/2.0_wp*deg2rad)
        end do
        call model%add_curve(x = x, &
                             y = y, &
                             z = z, &
                             radius = 0.1_wp,&
                             num_points = 10, &
                             initial_cap=.true.,final_cap=.true.)
    end do
    call model%add_axes(origin = [0.0_wp,0.0_wp,0.0_wp],&
                        vx = [4.0_wp,0.0_wp,0.0_wp],&
                        vy = [0.0_wp,1.0_wp,0.0_wp],&
                        vz = [0.0_wp,0.0_wp,1.0_wp],&
                        radius = 0.1_wp,&
                        num_points = 10,&
                        arrowhead_radius_factor = 2.0_wp,&
                        arrowhead_length_factor = 0.2_wp)
    call model%add_sphere([0.0_wp, 0.0_wp, 0.0_wp], 0.1_wp, 20, 40)
    call model%write_binary_stl_file('curve2.stl',istat)
    if (istat/=0) error stop 'error creating curves'
    call model%destroy()

    ! cone:
    call model%add_cone([0.0_wp,0.0_wp,0.0_wp],&
                        [0.0_wp,0.0_wp,1.0_wp],&
                        radius = 0.5_wp,&
                        num_points= 20, &
                        initial_cap = .true. )
    call model%add_cone([0.0_wp,0.0_wp,2.0_wp],&
                        [0.0_wp,0.0_wp,1.0_wp],&
                        radius = 0.5_wp,&
                        num_points= 20, &
                        initial_cap = .true. )
    call model%write_binary_stl_file('cones.stl',istat)
    if (istat/=0) error stop 'error creating cones'
    call model%destroy()

    ! read test:
    call model%read('cones.stl',istat)
    if (istat/=0) error stop 'error reading binary file'
    call model%write_binary_stl_file('cones_resaved.stl',istat)
    if (istat/=0) error stop 'error writing binary file'
    call model%destroy()

    write(*,*) 'Done.'

    end program main
!********************************************************************************
