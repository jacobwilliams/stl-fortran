!********************************************************************************
!>
!  STL (STereoLithography) file library.
!
!### Reference
!  * https://en.wikipedia.org/wiki/STL_(file_format)

    module stl_module

    use,intrinsic :: iso_c_binding
    use,intrinsic :: iso_fortran_env, only: wp => real64

    implicit none

    private

    type,public :: plate
        !! a 3D triangular plate.
        !! [note that the order of the vertices defines the
        !! surface normal via the right-hand rule]
        real(wp),dimension(3) :: v1 !! vertex
        real(wp),dimension(3) :: v2 !! vertex
        real(wp),dimension(3) :: v3 !! vertex
    end type plate

    public :: write_ascii_stl_file
    public :: write_binary_stl_file

    contains
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

    subroutine write_binary_stl_file(filename,plates,bounding_box,istat)

    implicit none

    character(len=*),intent(in)         :: filename      !! STL file name
    type(plate),dimension(:),intent(in) :: plates        !! array of triangular plates
    integer,intent(out)                 :: istat         !! `iostat` code
    real(wp),intent(in),optional        :: bounding_box  !! scale vertices so that model fits in a
                                                         !! box of this size (if <=0, no scaling is done)

    integer :: iunit                        !! file unit number
    integer :: i                            !! counter
    integer(c_int32_t) :: n_plates          !! number of plates [32 bits]
    real(c_float),dimension(3) :: n         !! normal vector [32 bits]
    real(c_float),dimension(3) :: v1,v2,v3  !! vertex vectors [32 bits]
    real(wp) :: scale                       !! scale factor

    integer(c_int16_t),parameter :: z = 0  !! Attribute byte count [16 bits]
    character(kind=c_char,len=80),parameter :: header = repeat(' ',80) !! [8 bits x 80]

    n_plates = size(plates,kind=c_int32_t)

    scale = compute_vertex_scale(plates,bounding_box)

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
        do i = 1, n_plates
            n = real(normal(plates(i)%v1,plates(i)%v2,plates(i)%v3), c_float)
            v1 = real(plates(i)%v1*scale, c_float)
            v2 = real(plates(i)%v2*scale, c_float)
            v3 = real(plates(i)%v3*scale, c_float)
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

    subroutine write_ascii_stl_file(filename,plates,modelname,istat,bounding_box)

    implicit none

    character(len=*),intent(in)          :: filename      !! STL file name
    type(plate),dimension(:),intent(in)  :: plates        !! array of triangular plates
    character(len=*),intent(in)          :: modelname     !! the solid name (should not contain spaces)
    integer,intent(out)                  :: istat         !! `iostat` code
    real(wp),intent(in),optional         :: bounding_box  !! scale vertices so that model fits in a
                                                          !! box of this size (if <=0, no scaling is done)
    integer :: iunit        !! file unit number
    integer  :: i           !! counter
    integer  :: n_plates    !! number of plates
    real(wp) :: scale       !! scale factor

    character(len=*),parameter :: fmt = '(A,1X,E30.16,1X,E30.16,1X,E30.16)' !! format statement for vectors

    n_plates = size(plates)
    scale = compute_vertex_scale(plates,bounding_box)

    ! open the text file:
    open(newunit=iunit, file=trim(filename), status='REPLACE', iostat=istat)

    if (istat==0) then

        ! write the file:
        write(iunit,'(A)') 'solid '//trim(modelname)
        do i = 1, n_plates
            write(iunit,fmt)    'facet normal', normal(plates(i)%v1,plates(i)%v2,plates(i)%v3)
            write(iunit,'(A)')  '    outer loop'
            write(iunit,fmt)    '        vertex', plates(i)%v1 * scale
            write(iunit,fmt)    '        vertex', plates(i)%v2 * scale
            write(iunit,fmt)    '        vertex', plates(i)%v3 * scale
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

    pure function compute_vertex_scale(plates,bounding_box) result(scale)

    implicit none

    type(plate),dimension(:),intent(in) :: plates        !! array of triangular plates
    real(wp),intent(in),optional        :: bounding_box  !! scale vertices so that model fits in a
                                                         !! box of this size (if <=0, no scaling is done)
    real(wp)                            :: scale         !! scale factor

    real(wp) :: max_value !! largest absolute value of any vertex coordinate
    integer :: i !! counter

    scale = 1.0_wp
    if (present(bounding_box)) then
        if (bounding_box>0.0_wp) then
            max_value = -huge(1.0_wp)
            do i = 1, size(plates)
                max_value = max(max_value, maxval(abs(plates(i)%v1)),&
                                           maxval(abs(plates(i)%v2)),&
                                           maxval(abs(plates(i)%v3)) )
            end do
            scale = bounding_box / max_value
        end if
    end if

    end function compute_vertex_scale
!********************************************************************************

!********************************************************************************
!>
!  Shift the vertex coordinates so that there are no non-positive components.

    subroutine shift_mesh(plates)

    implicit none

    type(plate),dimension(:),intent(inout) :: plates !! array of triangular plates

    integer :: i !! counter
    real(wp),dimension(3) :: offset !! offset vector for vertext coordinates [x,y,z]
    real(wp),dimension(3) :: mins   !! min values of vertex coordinates [x,y,z]

    real(wp),parameter :: tiny = 1.0e-4_wp !! small value to avoid zero

    ! first find the min value of each coordinate:
    mins = huge(1.0_wp)
    do i = 1, size(plates)
        mins(1) = min(mins(1), plates(i)%v1(1), plates(i)%v2(1), plates(i)%v3(1) )
        mins(2) = min(mins(2), plates(i)%v1(2), plates(i)%v2(2), plates(i)%v3(2) )
        mins(3) = min(mins(3), plates(i)%v1(3), plates(i)%v2(3), plates(i)%v3(3) )
    end do

    ! compute the offset vector:
    offset = 0.0_wp
    if (mins(1) <= 0.0_wp) offset(1) = abs(mins(1)) + tiny
    if (mins(2) <= 0.0_wp) offset(2) = abs(mins(2)) + tiny
    if (mins(3) <= 0.0_wp) offset(3) = abs(mins(3)) + tiny

    if (any(offset/=0.0_wp)) then
        ! now add offset vector to each
        do i = 1, size(plates)
            plates(i)%v1 = plates(i)%v1 + offset(1)
            plates(i)%v2 = plates(i)%v2 + offset(2)
            plates(i)%v3 = plates(i)%v3 + offset(3)
        end do
    end if

    end subroutine shift_mesh
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

    if (rmag/=0.0_wp) then
        rhat = r/rmag
    else
        rhat = 0.0_wp
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
    end module stl_module
!********************************************************************************