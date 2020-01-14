!********************************************************************************
!>
!  STL (STereoLithography) file library.
!
!### Author
!  * Jacob Williams, Jan 12, 2020.
!
!### License
!  * BSD-3
!
!### Reference
!  * https://en.wikipedia.org/wiki/STL_(file_format)

    module stl_module

    use,intrinsic :: iso_c_binding
    use,intrinsic :: iso_fortran_env, only: wp => real64

    implicit none

    private

    ! constants:
    real(wp),parameter :: zero = 0.0_wp
    real(wp),parameter :: one  = 1.0_wp

    type :: plate
        !! a 3D triangular plate.
        !! [note that the order of the vertices defines the
        !! surface normal via the right-hand rule]
        real(wp),dimension(3) :: v1 = zero  !! first vertex
        real(wp),dimension(3) :: v2 = zero  !! second vertex
        real(wp),dimension(3) :: v3 = zero  !! third vertex
    end type plate

    type,public :: stl_file
        !! the main class for STL file I/O.
        private
        integer :: n_plates = 0         !! number of plates
        integer :: chunk_size = 1000    !! expand `plates` array in chunks of this size
        type(plate),dimension(:),allocatable :: plates !! the array of plates
        contains
        private
        procedure,public :: write_ascii_stl_file
        procedure,public :: write_binary_stl_file
        procedure,public :: add_plate
        procedure,public :: add_sphere
        procedure,public :: shift_mesh
        procedure,public :: set_chunk_size
        procedure,public :: destroy => destroy_stl_file
        procedure :: compute_vertex_scale
    end type stl_file

    contains
!********************************************************************************

!********************************************************************************
!>
!  Destroy an `stl_file`.

    subroutine destroy_stl_file(me)

    implicit none

    class(stl_file),intent(inout) :: me

    if (allocated(me%plates)) deallocate(me%plates)
    me%n_plates = 0

    end subroutine destroy_stl_file
!********************************************************************************

!********************************************************************************
!>
!  Set the chunk size in the class.

    subroutine set_chunk_size(me,chunk_size)

    class(stl_file),intent(inout) :: me
    integer,intent(in)            :: chunk_size !! must be >0

    me%chunk_size = max(1,chunk_size)

    end subroutine set_chunk_size
!********************************************************************************

!********************************************************************************
!>
!  Add a plate to the class.

    subroutine add_plate(me,v1,v2,v3)

    class(stl_file),intent(inout)    :: me
    real(wp),dimension(3),intent(in) :: v1 !! first vertex
    real(wp),dimension(3),intent(in) :: v2 !! second vertex
    real(wp),dimension(3),intent(in) :: v3 !! third vertex

    integer :: n !! actual size of `plates` array in the class
    type(plate),dimension(:),allocatable :: tmp !! for resizing the `plates` array

    if (allocated(me%plates)) then

        n = size(me%plates)

        if (me%n_plates == n) then
            ! have to add another chunk
            allocate(tmp(n+1))
            tmp(1:n) = me%plates
            tmp(n+1)%v1 = v1
            tmp(n+1)%v2 = v2
            tmp(n+1)%v3 = v3
            call move_alloc(tmp,me%plates)
            me%n_plates = me%n_plates + 1
            return
        else
            ! add to next element in the class
            me%n_plates = me%n_plates + 1
        end if

    else
        allocate(me%plates(me%chunk_size))
        me%n_plates = 1
    end if

    me%plates(me%n_plates)%v1 = v1
    me%plates(me%n_plates)%v2 = v2
    me%plates(me%n_plates)%v3 = v3

    end subroutine add_plate
!********************************************************************************

!********************************************************************************
!>
!  Generate a binary STL file.
!
!### Notes
!
!  The file format is:
!```
!  UINT8[80] – Header
!  UINT32 – Number of triangles
!  foreach triangle
!    REAL32[3] – Normal vector
!    REAL32[3] – Vertex 1
!    REAL32[3] – Vertex 2
!    REAL32[3] – Vertex 3
!    UINT16 – Attribute byte count
!  end
!```

    subroutine write_binary_stl_file(me,filename,istat,bounding_box)

    implicit none

    class(stl_file),intent(in)    :: me
    character(len=*),intent(in)   :: filename      !! STL file name
    integer,intent(out)           :: istat         !! `iostat` code
    real(wp),intent(in),optional  :: bounding_box  !! scale vertices so that model fits in a
                                                   !! box of this size (if <=0, no scaling is done)

    integer :: iunit                        !! file unit number
    integer :: i                            !! counter
    integer(c_int32_t) :: n_plates          !! number of plates [32 bits]
    real(c_float),dimension(3) :: n         !! normal vector [32 bits]
    real(c_float),dimension(3) :: v1,v2,v3  !! vertex vectors [32 bits]
    real(wp) :: scale                       !! scale factor

    integer(c_int16_t),parameter :: z = 0  !! Attribute byte count [16 bits]
    character(kind=c_char,len=80),parameter :: header = repeat(' ',80) !! [8 bits x 80]

    n_plates = int(me%n_plates,kind=c_int32_t)

    scale = me%compute_vertex_scale(bounding_box)

    ! open the binary file:
    open(newunit=iunit,&
         file    = filename,&
         action  = 'WRITE',&
         status  = 'REPLACE',&
         form    = 'UNFORMATTED', &
         access  = 'STREAM', &
         iostat  = istat)

    if (istat==0) then

        ! write the file:
        write(iunit) header,n_plates
        do i = 1, me%n_plates
            n = real(normal(me%plates(i)%v1,me%plates(i)%v2,me%plates(i)%v3), c_float)
            v1 = real(me%plates(i)%v1*scale, c_float)
            v2 = real(me%plates(i)%v2*scale, c_float)
            v3 = real(me%plates(i)%v3*scale, c_float)
            write(iunit) n,v1,v2,v3,z
        end do

        ! close the file:
        close(iunit)

    end if

    end subroutine write_binary_stl_file
!********************************************************************************

!********************************************************************************
!>
!  Generate an ascii STL file.

    subroutine write_ascii_stl_file(me,filename,modelname,istat,bounding_box)

    implicit none

    class(stl_file),intent(in)    :: me
    character(len=*),intent(in)   :: filename      !! STL file name
    character(len=*),intent(in)   :: modelname     !! the solid name (should not contain spaces)
    integer,intent(out)           :: istat         !! `iostat` code
    real(wp),intent(in),optional  :: bounding_box  !! scale vertices so that model fits in a
                                                   !! box of this size (if <=0, no scaling is done)

    integer  :: iunit       !! file unit number
    integer  :: i           !! counter
    real(wp) :: scale       !! scale factor

    character(len=*),parameter :: fmt = '(A,1X,E30.16,1X,E30.16,1X,E30.16)' !! format statement for vectors

    scale = me%compute_vertex_scale(bounding_box)

    ! open the text file:
    open(newunit=iunit, file=trim(filename), status='REPLACE', iostat=istat)

    if (istat==0) then

        ! write the file:
        write(iunit,'(A)') 'solid '//trim(modelname)
        do i = 1, me%n_plates
            write(iunit,fmt)    'facet normal', normal(me%plates(i)%v1,me%plates(i)%v2,me%plates(i)%v3)
            write(iunit,'(A)')  '    outer loop'
            write(iunit,fmt)    '        vertex', me%plates(i)%v1 * scale
            write(iunit,fmt)    '        vertex', me%plates(i)%v2 * scale
            write(iunit,fmt)    '        vertex', me%plates(i)%v3 * scale
            write(iunit,'(A)')  '    end loop'
            write(iunit,'(A)')  'end facet'
        end do
        write(iunit,'(A)') 'endsolid '//trim(modelname)

        ! close the file:
        close(iunit)

    end if

    end subroutine write_ascii_stl_file
!********************************************************************************

!********************************************************************************
!>
!  Compute the scale factor for the vertices (for writing to a file).

    pure function compute_vertex_scale(me,bounding_box) result(scale)

    implicit none

    class(stl_file),intent(in)    :: me
    real(wp),intent(in),optional  :: bounding_box  !! scale vertices so that model fits in a
                                                   !! box of this size (if <=0, no scaling is done)
    real(wp)                      :: scale         !! scale factor

    real(wp) :: max_value  !! largest absolute value of any vertex coordinate
    integer  :: i          !! counter

    scale = one
    if (present(bounding_box)) then
        if (bounding_box>zero) then
            max_value = -huge(one)
            do i = 1, size(me%plates)
                max_value = max(max_value, maxval(abs(me%plates(i)%v1)),&
                                           maxval(abs(me%plates(i)%v2)),&
                                           maxval(abs(me%plates(i)%v3)) )
            end do
            scale = bounding_box / max_value
        end if
    end if

    end function compute_vertex_scale
!********************************************************************************

!********************************************************************************
!>
!  Shift the vertex coordinates so that there are no non-positive components.

    subroutine shift_mesh(me)

    implicit none

    class(stl_file),intent(inout) :: me

    integer :: i !! counter
    integer :: j !! counter
    real(wp),dimension(3) :: offset !! offset vector for vertext coordinates [x,y,z]
    real(wp),dimension(3) :: mins   !! min values of vertex coordinates [x,y,z]

    real(wp),parameter :: tiny = 1.0e-4_wp !! small value to avoid zero

    ! first find the min value of each coordinate:
    mins = huge(one)
    do i = 1, size(me%plates)
        do concurrent (j = 1:3)
            mins(j) = min(mins(j), &
                            me%plates(i)%v1(j), &
                            me%plates(i)%v2(j), &
                            me%plates(i)%v3(j) )
        end do
    end do

    ! compute the offset vector:
    offset = zero
    do concurrent (j = 1:3)
        if (mins(j) <= zero) offset(j) = abs(mins(j)) + tiny
    end do

    if (any(offset/=zero)) then
        ! now add offset vector to each
        do i = 1, size(me%plates)
            me%plates(i)%v1 = me%plates(i)%v1 + offset
            me%plates(i)%v2 = me%plates(i)%v2 + offset
            me%plates(i)%v3 = me%plates(i)%v3 + offset
        end do
    end if

    end subroutine shift_mesh
!********************************************************************************

!********************************************************************************
!>
!  Add a sphere to an STL file.

    subroutine add_sphere(me,center,radius,num_lat_points,num_lon_points)

    implicit none

    class(stl_file),intent(inout)    :: me
    real(wp),dimension(3),intent(in) :: center         !! coordinates of sphere center [x,y,z]
    real(wp),intent(in)              :: radius         !! radius of the sphere
    integer,intent(in)               :: num_lat_points !! number of latitude points (not counting poles)
    integer,intent(in)               :: num_lon_points !! number of longitude points

    integer :: i  !! counter
    integer :: j  !! counter
    real(wp) :: delta_lat  !! step in latitude (deg)
    real(wp) :: delta_lon  !! step in longitude (deg)
    real(wp),dimension(:),allocatable :: lat !! array of latitude values (deg)
    real(wp),dimension(:),allocatable :: lon !! array of longitude value (deg)
    real(wp),dimension(3) :: v1,v2,v3,v4 !! vertices

    ! Example:
    !
    ! num_lat_points = 3
    ! num_lon_points = 5
    !
    !  90 -------------  North pole
    !     |  *  *  *  |
    !     |  *  *  *  |
    !     |  *  *  *  |
    ! -90 -------------  South pole
    !     0          360

    delta_lat = 180.0_wp / (1+num_lat_points)
    delta_lon = 360.0_wp / (1+num_lon_points)

    lat = -90.0_wp + [(delta_lat*(i-1), i = 1,num_lat_points+2)]
    lon =            [(delta_lon*(i-1), i = 1,num_lon_points+2)]

    ! generate all the plates on the sphere.
    ! start at bottom left and go right then up.
    ! each box is two triangular plates.
    do i = 1, num_lat_points+1
        do j = 1, num_lon_points+1

            !   3----2
            !   |  / |
            !   | /  |
            ! i 1----4
            !   j

            v1 = spherical_to_cartesian(radius,lon(j),  lat(i)  ) + center
            v2 = spherical_to_cartesian(radius,lon(j+1),lat(i+1)) + center
            v3 = spherical_to_cartesian(radius,lon(j),  lat(i+1)) + center
            v4 = spherical_to_cartesian(radius,lon(j+1),lat(i)  ) + center
            call me%add_plate(v1,v2,v3)
            call me%add_plate(v1,v4,v2)

        end do
    end do

    end subroutine add_sphere
!********************************************************************************

!********************************************************************************
!>
!  Normal vector for the plate (computed using right hand rule).

    pure function normal(v1,v2,v3) result(n)

    implicit none

    real(wp),dimension(3),intent(in) :: v1  !! first vertex of the triangle [x,y,z]
    real(wp),dimension(3),intent(in) :: v2  !! second vertex of the triangle [x,y,z]
    real(wp),dimension(3),intent(in) :: v3  !! third vertex of the triangle [x,y,z]
    real(wp),dimension(3) :: n  !! surface normal vector

    n = unit( cross( v2-v1, v3-v1 ) )

    end function normal
!********************************************************************************

!********************************************************************************
!>
!  3x1 Unit vector.

    pure function unit(r) result(rhat)

    implicit none

    real(wp),dimension(3)            :: rhat
    real(wp),dimension(3),intent(in) :: r

    real(wp) :: rmag

    rmag = norm2(r)

    if (rmag/=zero) then
        rhat = r/rmag
    else
        rhat = zero
    end if

    end function unit
!********************************************************************************

!********************************************************************************
!>
!  Vector cross product.

    pure function cross(a,b) result(axb)

    implicit none

    real(wp),dimension(3) :: axb
    real(wp),dimension(3),intent(in) :: a
    real(wp),dimension(3),intent(in) :: b

    axb(1) = a(2)*b(3) - a(3)*b(2)
    axb(2) = a(3)*b(1) - a(1)*b(3)
    axb(3) = a(1)*b(2) - a(2)*b(1)

    end function cross
!********************************************************************************

!********************************************************************************
!>
!  Convert spherical (r,alpha,beta) to Cartesian (x,y,z).

    pure function spherical_to_cartesian(r,alpha,beta) result(rvec)

    implicit none

    real(wp),intent(in)   :: r        !! magnitude
    real(wp),intent(in)   :: alpha    !! right ascension [deg]
    real(wp),intent(in)   :: beta     !! declination [deg]
    real(wp),dimension(3) :: rvec     !! [x,y,z] vector

    real(wp),parameter :: deg2rad = acos(-1.0_wp) / 180.0_wp !! degrees to radians

    rvec(1) = r * cos(alpha*deg2rad) * cos(beta*deg2rad)
    rvec(2) = r * sin(alpha*deg2rad) * cos(beta*deg2rad)
    rvec(3) = r * sin(beta*deg2rad)

    end function spherical_to_cartesian
!********************************************************************************

!********************************************************************************
    end module stl_module
!********************************************************************************