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
    use,intrinsic :: iso_fortran_env, only: wp => real64, error_unit

    implicit none

    private

    ! constants:
    real(wp),parameter :: zero = 0.0_wp
    real(wp),parameter :: one  = 1.0_wp
    real(wp),parameter :: deg2rad = acos(-1.0_wp) / 180.0_wp !! degrees to radians
    real(wp),dimension(3),parameter :: x_unit = [one,zero,zero] !! x-axis unit vector
    real(wp),dimension(3),parameter :: y_unit = [zero,one,zero] !! y-axis unit vector
    real(wp),dimension(3),parameter :: z_unit = [zero,zero,one] !! z-axis unit vector

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
        integer :: n_plates = 0       !! number of plates
        integer :: chunk_size = 1000  !! expand `plates` array in chunks of this size
        type(plate),dimension(:),allocatable :: plates !! the array of plates
        contains
        private
        procedure,public :: write_ascii_stl_file
        procedure,public :: write_binary_stl_file
        procedure,public :: read => read_binary_stl_file
        procedure,public :: read_tab_file
        procedure,public :: destroy => destroy_stl_file
        procedure,public :: add_plate
        procedure,public :: add_sphere
        procedure,public :: add_cylinder
        procedure,public :: add_curve
        procedure,public :: add_cone
        procedure,public :: add_arrow
        procedure,public :: add_axes
        procedure,public :: shift_mesh
        procedure,public :: set_chunk_size
        procedure :: generate_circle
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
            allocate(tmp(n+me%chunk_size))
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
    integer,intent(out)           :: istat         !! `iostat` code (=0 if no errors)
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
    open(newunit = iunit,&
         file    = filename,&
         action  = 'WRITE',&
         status  = 'REPLACE',&
         form    = 'UNFORMATTED', &
         access  = 'STREAM', &
         iostat  = istat)

    if (istat==0) then

        ! write the file:
        write(iunit,iostat=istat) header,n_plates
        if (istat==0) then
            do i = 1, me%n_plates
                n = real(normal(me%plates(i)%v1,me%plates(i)%v2,me%plates(i)%v3), c_float)
                v1 = real(me%plates(i)%v1*scale, c_float)
                v2 = real(me%plates(i)%v2*scale, c_float)
                v3 = real(me%plates(i)%v3*scale, c_float)
                write(iunit,iostat=istat) n,v1,v2,v3,z
                if (istat/=0) exit
            end do
        end if
        ! close the file:
        close(iunit)

    end if

    end subroutine write_binary_stl_file
!********************************************************************************

!********************************************************************************
!>
!  Read a binary STL file.

    subroutine read_binary_stl_file(me,filename,istat)

    implicit none

    class(stl_file),intent(out) :: me
    character(len=*),intent(in) :: filename !! STL file name
    integer,intent(out)         :: istat    !! `iostat` code (=0 if no errors)

    integer :: iunit                        !! file unit number
    integer :: i                            !! counter
    integer(c_int32_t) :: n_plates          !! number of plates [32 bits]
    real(c_float),dimension(3) :: n         !! normal vector [32 bits]
    real(c_float),dimension(3) :: v1,v2,v3  !! vertex vectors [32 bits]
    integer(c_int16_t) :: z                 !! Attribute byte count [16 bits]
    character(kind=c_char,len=80) :: header !! [8 bits x 80]

    call me%destroy()

    ! open the binary file:
    open(newunit = iunit,&
         file    = filename,&
         action  = 'READ',&
         status  = 'OLD',&
         form    = 'UNFORMATTED', &
         access  = 'STREAM', &
         iostat  = istat)

    if (istat==0) then

        ! read header:
        read(iunit,iostat=istat) header, n_plates

        if (istat==0) then

            ! size arrays:
            me%n_plates = int(n_plates)
            allocate(me%plates(me%n_plates))

            ! read the data from the file:
            do i = 1, me%n_plates
                read(iunit,iostat=istat) n,v1,v2,v3,z
                if (istat/=0) exit
                ! only need to save the plates:
                me%plates(i)%v1 = real(v1, wp)
                me%plates(i)%v2 = real(v2, wp)
                me%plates(i)%v3 = real(v3, wp)
            end do

        end if

        ! close the file:
        close(iunit, iostat=istat)

    end if

    end subroutine read_binary_stl_file
!********************************************************************************

!********************************************************************************
!>
!  Read a text vertex-facet file.
!
!### File Format
!
!  The file is a text file consisting of:
!
!  * Number of Vertices [int, %12d], Number of Triangular Plates [int, %12d]
!  * Vertex table:
!    * Vertex Number [int, %10d], x coordinate [real, %15.5f], y coordinate [real, %15.5f], z coordinate [real, %15.5f]
!  * Plate table:
!    * Plate Number [int, %10d], 1st Plate Vertex Number [int, %10d], 2nd Plate Vertex Number [int, %10d], 3rd Plate Vertex Number [int, %10d]
!
!### Example
!
!  * https://sbnarchive.psi.edu/pds4/non_mission/gaskell.phobos.shape-model/data/phobos_ver64q.tab
!
!```
! 25350        49152
! 1       -6.77444        6.26815        6.01149
! 2       -6.63342        6.34195        6.08444
! 3       -6.49302        6.41635        6.15759
! 4       -6.34883        6.48872        6.22619
!  ...
! 1         1        67         2
! 2         1        66        67
! 3        66       132        67
! 4        66       131       132
! 5       131       197       132
!  ...
!```

    subroutine read_tab_file(me,filename,istat)

    implicit none

    class(stl_file),intent(out) :: me
    character(len=*),intent(in) :: filename !! Vertex-facet file name
    integer,intent(out)         :: istat    !! `iostat` code (=0 if no errors)

    integer :: iunit  !! file unit
    integer :: number_of_vertices !! number of vertices in the file
    integer :: number_of_plates   !! number of plates defined in the file (three vertices)
    integer :: i !! counter
    integer :: ii !! vertex or plate index
    integer :: i1 !! 1st vertex index of plate
    integer :: i2 !! 2nd vertex index of plate
    integer :: i3 !! 3rd vertex index of plate
    real(wp),dimension(:,:),allocatable :: v !! vertices from the file
    real(wp) :: x !! x coordinate of vertex
    real(wp) :: y !! y coordinate of vertex
    real(wp) :: z !! z coordinate of vertex

    !initialize:
    call me%destroy()

    !open the file:
    open(newunit = iunit,&
         file    = filename,&
         action  = 'READ',&
         status  = 'OLD',&
         iostat  = istat)

    if (istat==0) then

        read(iunit,*,iostat=istat) number_of_vertices, number_of_plates

        if (istat==0) then

            me%n_plates = number_of_plates
            allocate(me%plates(me%n_plates))
            allocate(v(3,number_of_vertices))

            ! first accumulate the vertex coordinates:
            do i = 1, number_of_vertices

                read(iunit,*,iostat=istat) ii, x, y, z

                !get the data:
                if (istat==0) then
                    v(1,i) = x
                    v(2,i) = y
                    v(3,i) = z
                else
                    write(error_unit,'(A)') 'Error reading vertex from file: '//trim(filename)
                    call me%destroy()
                    close(iunit)
                    return
                end if

            end do

            ! now, read the plate vertex indices, and add them to the class
            do i = 1, number_of_plates

                read(iunit,*,iostat=istat) ii,i1,i2,i3

                if (istat==0) then
                    me%plates(i)%v1 = v(:,i1)
                    me%plates(i)%v2 = v(:,i2)
                    me%plates(i)%v3 = v(:,i3)
                else
                    write(error_unit,'(A)') 'Error reading plate from file: '//trim(filename)
                    call me%destroy()
                    close(iunit)
                    return
                end if

            end do

            ! close the file:
            close(iunit)
            deallocate(v)

        else
            write(error_unit,'(A)') 'Error reading first line from file: '//trim(filename)
            close(iunit)
        end if

    end if

    end subroutine read_tab_file
!********************************************************************************

!********************************************************************************
!>
!  Generate an ascii STL file.

    subroutine write_ascii_stl_file(me,filename,modelname,istat,bounding_box)

    implicit none

    class(stl_file),intent(in)    :: me
    character(len=*),intent(in)   :: filename      !! STL file name
    character(len=*),intent(in)   :: modelname     !! the solid name (should not contain spaces)
    integer,intent(out)           :: istat         !! `iostat` code (=0 if no errors)
    real(wp),intent(in),optional  :: bounding_box  !! scale vertices so that model fits in a
                                                   !! box of this size (if <=0, no scaling is done)

    integer  :: iunit  !! file unit number
    integer  :: i      !! counter
    real(wp) :: scale  !! scale factor

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
!  Add a cylinder to an STL file.
!
!  The cylinder is specified by the initial and final x,y,z coordinates. Optionally,
!  an initial and final normal vector can be specified (if not specified,
!  then a default one is constructed).

    subroutine add_cylinder(me,v1,v2,radius,num_points,initial_cap,final_cap,&
                            initial_normal,final_normal,final_normal_used,initial_vector,final_initial_vector_used)

    implicit none

    class(stl_file),intent(inout)              :: me
    real(wp),dimension(3),intent(in)           :: v1                        !! coordinates of initial point
    real(wp),dimension(3),intent(in)           :: v2                        !! coordinates of final point
    real(wp),intent(in)                        :: radius                    !! radius of the cylinder
    integer,intent(in)                         :: num_points                !! number of point on the circle (>=3)
    logical,intent(in)                         :: initial_cap               !! add a cap plate to the initial point
    logical,intent(in)                         :: final_cap                 !! add a cap plate to the final point
    real(wp),dimension(3),intent(in),optional  :: initial_normal            !! outward normal vector for initial circle
    real(wp),dimension(3),intent(in),optional  :: final_normal              !! outward normal vector for final circle
    real(wp),dimension(3),intent(out),optional :: final_normal_used         !! outward normal vector for final circle
                                                                            !! actually used
    real(wp),dimension(3),intent(in),optional  :: initial_vector            !! vector to use to generate the initial
                                                                            !! circle (x_unit by default)
    real(wp),dimension(3),intent(out),optional :: final_initial_vector_used !! the initial vector used for the final
                                                                            !! cap to generate the points

    integer :: i  !! counter
    integer :: nc !! number of points on the circle
    real(wp),dimension(3) :: n0 !! normal vector for initial circle
    real(wp),dimension(3) :: nf !! normal vector for final circle
    real(wp),dimension(:,:),allocatable :: n0_cap_points !! points for the initial cap
    real(wp),dimension(:,:),allocatable :: nf_cap_points !! points for the final cap

    nc = max(3, num_points)

    ! compute the end unit vectors
    !
    !        1 _________2
    !        |          |
    !  n0 <--*----------*--> nf
    !        |          |
    !         ----------

    if (present(initial_normal)) then
        n0 = unit(initial_normal)
    else
        n0 = unit(v1-v2)
    end if
    if (present(final_normal)) then
        nf = unit(final_normal)
    else
        nf = unit(v2-v1)
    end if
    if (present(final_normal_used)) final_normal_used = nf ! return if necessary

    ! create the points on the initial cap (optionally add the plate)
    call me%generate_circle(v1,radius,n0,nc,initial_cap,n0_cap_points,initial_vector=initial_vector)

    ! create the points on the final cap (optionally add the plate)
    ! [use the same initial vector to sure that the plate will form a good cylinder]
    call me%generate_circle(v2,radius,nf,nc,final_cap,nf_cap_points,&
                            initial_vector=unit(n0_cap_points(:,1)-v1),cw=.true.)
    if (present(final_initial_vector_used)) final_initial_vector_used = unit(nf_cap_points(:,1)-v2)

    ! now connect the points to form the cylinder:
    !   1----2  nf
    !   |  / |
    !   | /  |
    !   1----2  n0
    do i = 1, nc-1
        call me%add_plate(n0_cap_points(:,i),n0_cap_points(:,i+1),nf_cap_points(:,i+1))
        call me%add_plate(n0_cap_points(:,i),nf_cap_points(:,i+1),nf_cap_points(:,i))
    end do
    ! last one:
    !   n----1  nf
    !   |  / |
    !   | /  |
    !   n----1  n0
   call me%add_plate(n0_cap_points(:,nc),n0_cap_points(:,1),nf_cap_points(:,1))
   call me%add_plate(n0_cap_points(:,nc),nf_cap_points(:,1),nf_cap_points(:,nc))

    end subroutine add_cylinder
!********************************************************************************

!********************************************************************************
!>
!  Add a cone to an STL file.
!
!  The cylinder is specified by the initial and final x,y,z coordinates. Optionally,
!  an initial and final normal vector can be specified (if not specified,
!  then a default one is constructed).

    subroutine add_cone(me,v1,v2,radius,num_points,initial_cap,initial_normal)

    implicit none

    class(stl_file),intent(inout)             :: me
    real(wp),dimension(3),intent(in)          :: v1             !! coordinates of initial point (bottom of the cone)
    real(wp),dimension(3),intent(in)          :: v2             !! coordinates of final point (point of the cone)
    real(wp),intent(in)                       :: radius         !! radius of the cone (the bottom plate)
    integer,intent(in)                        :: num_points     !! number of point on the circle (>=3)
    logical,intent(in)                        :: initial_cap    !! add a cap plate to the initial point (bottom)
    real(wp),dimension(3),intent(in),optional :: initial_normal !! outward normal vector for initial plate (bottom)

    integer :: i  !! counter
    integer :: nc !! number of points on the circle
    real(wp),dimension(3) :: n0 !! normal vector for initial circle
    real(wp),dimension(:,:),allocatable :: n0_cap_points !! points for the initial cap

    nc = max(3, num_points)

    ! compute the end unit vector:
    if (present(initial_normal)) then
        n0 = unit(initial_normal)
    else
        n0 = unit(v1-v2)
    end if

    ! create the points on the initial cap (optionally add the plate)
    call me%generate_circle(v1,radius,n0,nc,initial_cap,n0_cap_points)

    ! draw the cone plates
    !      *     v2
    !     / \
    !    /   \
    !   1--*--2  v1
    do i = 1, nc-1
        call me%add_plate(n0_cap_points(:,i),n0_cap_points(:,i+1),v2)
    end do
    ! last one:
    call me%add_plate(n0_cap_points(:,nc),n0_cap_points(:,1),v2)

    end subroutine add_cone
!********************************************************************************

!********************************************************************************
!>
!  Add x,y,z axes to an STL file.

    subroutine add_axes(me,origin,vx,vy,vz,radius,num_points,&
                        arrowhead_radius_factor,arrowhead_length_factor)

    implicit none

    class(stl_file),intent(inout)    :: me
    real(wp),dimension(3),intent(in) :: origin                  !! coordinates of the origin of the axes
    real(wp),dimension(3),intent(in) :: vx                      !! x axis vector
    real(wp),dimension(3),intent(in) :: vy                      !! y axis vector
    real(wp),dimension(3),intent(in) :: vz                      !! z axis vector
    real(wp),intent(in)              :: radius                  !! radius of the cylinder
    integer,intent(in)               :: num_points              !! number of point on the circle (>=3)
    real(wp),intent(in)              :: arrowhead_radius_factor !! arrowhead cone radius factor
                                                                !! (multiple of cylinder radius)
    real(wp),intent(in)              :: arrowhead_length_factor !! arrowhead tip length factor
                                                                !! (multiple of vector length)

    call me%add_arrow(origin,vx,radius,num_points,arrowhead_radius_factor,arrowhead_length_factor)
    call me%add_arrow(origin,vy,radius,num_points,arrowhead_radius_factor,arrowhead_length_factor)
    call me%add_arrow(origin,vz,radius,num_points,arrowhead_radius_factor,arrowhead_length_factor)

    end subroutine add_axes
!********************************************************************************

!********************************************************************************
!>
!  Add an arrow to an STL file.

    subroutine add_arrow(me,origin,v,radius,num_points,&
                         arrowhead_radius_factor,arrowhead_length_factor)

    implicit none

    class(stl_file),intent(inout)    :: me
    real(wp),dimension(3),intent(in) :: origin                   !! coordinates of the origin of the axes
    real(wp),dimension(3),intent(in) :: v                        !! vector
    real(wp),intent(in)              :: radius                   !! radius of the cylinder
    integer,intent(in)               :: num_points               !! number of point on the circle (>=3)
    real(wp),intent(in)              :: arrowhead_radius_factor  !! arrowhead cone radius factor
                                                                 !! (multiple of cylinder radius)
    real(wp),intent(in)              :: arrowhead_length_factor  !! arrowhead tip length factor
                                                                 !! (multiple of vector length)

    call me%add_cylinder(origin,v,radius,num_points,&
                         initial_cap=.true.,final_cap=.true.)
    call me%add_cone(v,v+arrowhead_length_factor*v,&
                     arrowhead_radius_factor*radius,num_points,initial_cap=.true.)

    end subroutine add_arrow
!********************************************************************************

!********************************************************************************
!>
!  Generate the points in a circle, and optionally add it as a plate.

    subroutine generate_circle(me,c,radius,n,nc,add_circle,circle,initial_vector,cw)

    implicit none

    class(stl_file),intent(inout)                   :: me
    real(wp),dimension(3),intent(in)                :: c              !! center of the circle
    real(wp),intent(in)                             :: radius         !! radius of the cylinder
    real(wp),dimension(3),intent(in)                :: n              !! normal vector to the circle
    integer,intent(in)                              :: nc             !! number of points on the circle
                                                                      !! (must be at least 3)
    logical,intent(in)                              :: add_circle     !! to also add to the circle as a plate
    real(wp),dimension(:,:),allocatable,intent(out) :: circle         !! points on the circle
    real(wp),dimension(3),intent(in),optional       :: initial_vector !! vector to use to generate the initial
                                                                      !! circle (x_unit by default)
    logical,intent(in),optional                     :: cw             !! generate the points in the clockwise
                                                                      !! direction abound n (default is false)

    real(wp),dimension(3) :: v      !! initial vector for the circle
    integer               :: i      !! counter
    real(wp)              :: factor !! cw/ccw factor
    logical               :: compute_initial_vector !! if we need to compute an initial vector

    if (nc<3) error stop 'number of points on a circle must be at least 3'

    allocate(circle(3,nc))
    ! circle = -999

    factor = one
    if (present(cw)) then
        if (cw) factor = -one
    end if

    if (present(initial_vector)) then
        compute_initial_vector =  all(initial_vector==zero)
    else
        compute_initial_vector = .true.
    end if

    ! start with an initial vector on the circle (perpendicular to n0)
    ! [project x to circle (or y if x is parallel to n)]
    if (.not. compute_initial_vector) then
        v = unit(vector_projection_on_plane(initial_vector,n))
        if (.not. perpendicular(v, n) .or. all(v==zero)) then
            ! fall back to x or y axis
            v = unit(vector_projection_on_plane(x_unit,n))
            if (.not. perpendicular(v, n) .or. all(v==zero)) then
                v = unit(vector_projection_on_plane(y_unit,n))
            end if
        end if
    else
        v = unit(vector_projection_on_plane(x_unit,n))
        if (.not. perpendicular(v, n) .or. all(v==zero)) then
            v = unit(vector_projection_on_plane(y_unit,n))
        end if
    end if
    v = radius * unit(v)

    ! generate the points by rotating the initial vector around the circle:
    circle(:,1) = c + v
    do i = 2, nc
        circle(:,i) = c + axis_angle_rotation(v,n,(i-1)*factor*(360.0_wp/nc))
        if (add_circle) then
            ! draw the initial cap
            call me%add_plate(c,circle(:,i),circle(:,i-1))
        end if
    end do
    ! final plate that connects last to first
    if (add_circle) call me%add_plate(c,circle(:,1),circle(:,nc))

    end subroutine generate_circle
!********************************************************************************

!********************************************************************************
!>
!  Returns true if the two vectors are perpendicular.

    pure function perpendicular(v1, v2) result(is_parallel)

    implicit none

    real(wp),dimension(:),intent(in) :: v1
    real(wp),dimension(:),intent(in) :: v2
    logical :: is_parallel

    real(wp),parameter :: tol = 10.0_wp * epsilon(1.0_wp) !! tolerance

    is_parallel = abs(dot_product(unit(v1), unit(v2))) <= tol

    end function perpendicular
!********************************************************************************

!********************************************************************************
!>
!  Add a curve to an STL file.
!
!  A curve is a joined set of cylinders with no internal caps.

    subroutine add_curve(me,x,y,z,radius,num_points,&
                         initial_cap,initial_normal,final_cap,final_normal,initial_vector)

    implicit none

    class(stl_file),intent(inout)             :: me
    real(wp),dimension(:),intent(in)          :: x              !! x coordinate array
    real(wp),dimension(:),intent(in)          :: y              !! y coordinate array
    real(wp),dimension(:),intent(in)          :: z              !! z coordinate array
    real(wp),intent(in)                       :: radius         !! radius of the cylinder
    integer,intent(in)                        :: num_points     !! number of point on the cylinder perimeter
    logical,intent(in),optional               :: initial_cap    !! add a cap plate to the initial point
    real(wp),dimension(3),intent(in),optional :: initial_normal !! outward normal vector for initial circle
    logical,intent(in),optional               :: final_cap      !! add a cap plate to the final point
    real(wp),dimension(3),intent(in),optional :: final_normal   !! outward normal vector for final circle
    real(wp),dimension(3),intent(in),optional :: initial_vector !! vector to use to generate the first circle (x_unit by default)

    integer               :: i      !! counter
    integer               :: n      !! number of points
    real(wp),dimension(3) :: nv     !! for intermediate normal vectors
    real(wp),dimension(3) :: nv_tmp !! for intermediate normal vectors
    real(wp),dimension(3) :: v      !! for intermediate initial vectors

    n = min(size(x), size(y), size(z))
    if (n<2) error stop 'error: a curve must have more than one point'

    ! first cylinder [no final cap unless only two points]
    call me%add_cylinder([x(1),y(1),z(1)],&
                         [x(2),y(2),z(2)],&
                         radius,num_points,&
                         initial_cap=initial_cap,initial_normal=initial_normal,initial_vector=initial_vector,&
                         final_cap=n==2,final_normal_used=nv,final_initial_vector_used=v)

    if (n>3) then
        ! intermediate cylinders (the initial normal is the final normal from the previous cylinder)
        do i = 2, n-2
            call me%add_cylinder([x(i),y(i),z(i)],&
                                 [x(i+1),y(i+1),z(i+1)],&
                                 radius,num_points,&
                                 initial_cap=.false.,initial_normal=-nv,&
                                 final_cap=.false.,final_normal_used=nv_tmp,&
                                 initial_vector=v)
            nv = unit(nv_tmp)
        end do
    end if

    ! last cylinder [no initial cap]
    if (n>=3) then
        call me%add_cylinder([x(n-1),y(n-1),z(n-1)],&
                             [x(n),y(n),z(n)],&
                             radius,num_points,&
                             final_cap=final_cap,final_normal=final_normal,&
                             initial_normal=-nv,initial_cap=.false.,&
                             initial_vector=v)
    end if

    end subroutine add_curve
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

    rvec(1) = r * cos(alpha*deg2rad) * cos(beta*deg2rad)
    rvec(2) = r * sin(alpha*deg2rad) * cos(beta*deg2rad)
    rvec(3) = r * sin(beta*deg2rad)

    end function spherical_to_cartesian
!********************************************************************************

!********************************************************************************
!> author: Jacob Williams
!  date: 7/20/2014
!
!  Rotate a 3x1 vector in space, given an axis and angle of rotation.
!
!# Reference
!   * [Wikipedia](http://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula)

    pure function axis_angle_rotation(v,k,theta) result(vrot)

    implicit none

    real(wp),dimension(3),intent(in)  :: v      !! vector to rotate
    real(wp),dimension(3),intent(in)  :: k      !! rotation axis
    real(wp),intent(in)               :: theta  !! rotation angle [deg]
    real(wp),dimension(3)             :: vrot   !! result

    real(wp),dimension(3) :: khat
    real(wp) :: ct,st

    ct = cos(theta*deg2rad)
    st = sin(theta*deg2rad)
    khat = unit(k)   !rotation axis unit vector

    vrot = v*ct + cross(khat,v)*st + khat*dot_product(khat,v)*(one-ct)

    end function axis_angle_rotation
!********************************************************************************

!********************************************************************************
!>
!  The projection of one vector onto another vector.
!
!### Reference
!   * [Wikipedia](http://en.wikipedia.org/wiki/Gram-Schmidt_process)
!
!### History
!  * Jacob Williams : 7/21/2014
!  * JW : fixed a typo : 6/18/2021

    pure function vector_projection(a,b) result(c)

    implicit none

    real(wp),dimension(:),intent(in)       :: a  !! the original vector
    real(wp),dimension(size(a)),intent(in) :: b  !! the vector to project on to
    real(wp),dimension(size(a))            :: c  !! the projection of a onto b

    real(wp) :: bmag2

    bmag2 = dot_product(b,b)

    if (bmag2==zero) then
        c = zero
    else
        c = b * dot_product(a,b) / bmag2
    end if

    end function vector_projection
!********************************************************************************

!********************************************************************************
!>
!  Project a vector onto a plane.
!
!### Reference
!   * [Projection of a Vector onto a Plane](http://www.maplesoft.com/support/help/Maple/view.aspx?path=MathApps/ProjectionOfVectorOntoPlane)

    pure function vector_projection_on_plane(a,b) result(c)

    implicit none

    real(wp),dimension(3),intent(in)  :: a !! the original vector
    real(wp),dimension(3),intent(in)  :: b !! the plane to project on to (a normal vector)
    real(wp),dimension(3) :: c !! the projection of a onto the b plane

    c = a - vector_projection(a,b)

    end function vector_projection_on_plane
!********************************************************************************

!********************************************************************************
    end module stl_module
!********************************************************************************