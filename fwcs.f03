!   This file is part of gaialaxy
!
!   Copyright (C) 2022 C. Ringeval
!   
!   gaialaxy is free software: you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation, either version 3 of the License, or
!   (at your option) any later version.
!
!   gaialaxy is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with gaialaxy.  If not, see <https://www.gnu.org/licenses/>.


module fwcs
  use, intrinsic :: iso_c_binding
  implicit none

  private

  integer, parameter :: lenheader = 80
  integer, parameter :: lenkeyval = 80 - 8
  

  integer, parameter :: cd = C_DOUBLE
  integer, parameter :: ci = C_INT

  logical, parameter :: display = .false.
  
  interface
     
     subroutine ini_wcsprm(naxis, crpix, cdelt, crval, cptr_pcij, cptr_cunit, cptr_ctype, &
          dateref ) bind(C)
       import C_INT, C_DOUBLE, C_PTR, C_CHAR, lenkeyval
       integer(C_INT), value :: naxis
       real(C_DOUBLE), dimension(*) :: crpix, cdelt, crval
       type(C_PTR), value :: cptr_pcij, cptr_cunit, cptr_ctype
       character(C_CHAR), dimension(lenkeyval) :: dateref
     end subroutine ini_wcsprm

     subroutine destroy_wcsprm() bind(C)

     end subroutine destroy_wcsprm

     subroutine size_wcsprm(sizes) bind(C)
       import C_INT
       integer(C_INT), dimension(2) :: sizes
     end subroutine size_wcsprm
     
     subroutine fix_wcsprm() bind(C)

     end subroutine fix_wcsprm

     subroutine hdo_alloc_wcsprm(cptr_header,n) bind(C)
       import C_CHAR, C_PTR, C_INT
       type(C_PTR), intent(out) :: cptr_header
       integer(C_INT) :: n
     end subroutine hdo_alloc_wcsprm

     subroutine hdo_free_wcsprm(hptr) bind(C)
       import C_PTR
       type(C_PTR), value :: hptr
     end subroutine hdo_free_wcsprm

     subroutine p2s_wcsprm(ncoord, nelem, pixcrd, imgcrd, phi, theta, world) bind(C)
       import C_INT, C_DOUBLE
       integer(C_INT), value :: ncoord, nelem
       real(C_DOUBLE), dimension(nelem,ncoord), intent(in) :: pixcrd
       real(C_DOUBLE), dimension(nelem,ncoord) :: imgcrd
       real(C_DOUBLE), dimension(ncoord) :: phi, theta
       real(C_DOUBLE), dimension(nelem,ncoord) :: world
     end subroutine p2s_wcsprm

     subroutine s2p_wcsprm(ncoord, nelem, world, phi, theta, imgcrd, pixcrd) bind(C)
       import C_INT, C_DOUBLE
       integer(C_INT), value :: ncoord, nelem
       real(C_DOUBLE), dimension(nelem,ncoord), intent(in) :: world
       real(C_DOUBLE), dimension(ncoord) :: phi, theta
       real(C_DOUBLE), dimension(nelem,ncoord) :: imgcrd, pixcrd
     end subroutine s2p_wcsprm

     subroutine ccs_wcsprm(lng2p1, lat2p1, lng1p2, clng, clat, radesys, &
          equinox, alt) bind(C)
       import C_DOUBLE, C_PTR, C_CHAR
       real(C_DOUBLE), value :: lng2p1, lat2p1, lng1p2
       character(C_CHAR), dimension(*), intent(in) :: clng, clat, radesys
       real(C_DOUBLE), value :: equinox
       character(C_CHAR), dimension(*), intent(in) :: alt
     end subroutine ccs_wcsprm
       
  end interface


  public cd, lenheader
  public create_platecarree_wcsheader, get_wcsheader, print_wcsheader
  public world_to_pixel_wcsheader, pixel_to_world_wcsheader, galactic_to_equatorial_wcsheader
  public free_wcsheader, check_wcsheader
  
contains


!fortran to c string converter
  function f_c_string(fname)
    implicit none
    character(len=*), intent(in) :: fname
    character(kind=C_CHAR, len=len(fname)+1) :: f_c_string
    integer :: i,n

    n = len(fname)
    do i=1,n
       f_c_string(i:i) = fname(i:i)
    enddo
    f_c_string(n+1:n+1)=C_NULL_CHAR

  end function f_c_string


  
  subroutine create_platecarree_wcsheader(nx,ny,pixscales)
    implicit none
    integer, intent(in) :: nx,ny
    integer(ci), parameter :: naxis = 2
    real(cd), dimension(naxis), intent(out), optional :: pixscales
        
    
    real(cd), dimension(naxis) :: crpix, cdelt, crval

    character(kind=C_CHAR, len=lenkeyval), dimension(:), pointer :: ctype, cunit
    real(C_DOUBLE), dimension(:,:), pointer :: pcij

    character(kind=C_CHAR, len=lenkeyval) :: dateref = '2022-01-01'//C_NULL_CHAR
    
    crpix = (/ 0.5_cd*nx + 0.5_cd, 0.5_cd*ny + 0.5_cd /)
    cdelt = (/ 360._cd/nx, 180._cd/ny /)
    crval = (/ 0._cd, 0._cd /)

    allocate(pcij(naxis,naxis))
    pcij = reshape( (/ 1.0_cd, 0.0_cd, 0.0_cd, 1.0_cd /), shape=(/ naxis, naxis /), order = (/ 2, 1 /) )
    
    allocate(ctype(naxis), cunit(naxis))
    
    cunit(1) = f_c_string('deg')
    cunit(2) = f_c_string('deg')
    
    ctype(1) = f_c_string('GLON-CAR')
    ctype(2) = f_c_string('GLAT-CAR')
    
    call ini_wcsprm(naxis, crpix, cdelt, crval, C_LOC(pcij), C_LOC(cunit), &
         C_LOC(ctype), dateref)
    
    call fix_wcsprm()

    if (present(pixscales)) then
       pixscales = cdelt
    end if
    
  end subroutine create_platecarree_wcsheader


  function check_wcsheader()
    implicit none
    logical :: check_wcsheader

    integer(C_INT), dimension(2) :: sizes

    check_wcsheader = .false.
    
    call size_wcsprm(sizes)

    if (sizes(2).ne.0) then
       check_wcsheader = .true.
       write(*,*)'check_wcsheader: ',sizes
    endif
    
  end function check_wcsheader
  

  subroutine free_wcsheader()
    implicit none
    
    call destroy_wcsprm()

  end subroutine free_wcsheader

  
  subroutine print_wcsheader()
    type(C_PTR) :: cptr_header
    character(kind=C_CHAR, len=lenheader), dimension(:), pointer :: fptr_header

    integer(ci) :: nkeys

    integer :: i
    
    call hdo_alloc_wcsprm(cptr_header,nkeys)
    call C_F_POINTER(CPTR=cptr_header , FPTR=fptr_header , shape=[nkeys])

    do i=1,nkeys
       write(*,*)fptr_header(i)
    enddo
       
    call hdo_free_wcsprm(C_LOC(fptr_header))
    
  end subroutine print_wcsheader


  subroutine get_wcsheader(header)
    character(len=lenheader), dimension(:), allocatable :: header

    type(C_PTR) :: cptr_header
    character(kind=C_CHAR, len=lenheader), dimension(:), pointer :: fptr_header
    integer(ci) :: nkeys

    integer :: i

    if (allocated(header)) stop 'get_wcsheader: header already set!'
    
    call hdo_alloc_wcsprm(cptr_header,nkeys)
    call C_F_POINTER(CPTR=cptr_header , FPTR=fptr_header , shape=[nkeys])

    allocate(header(nkeys))
    
    do i=1,nkeys
!       write(*,*)fptr_header(i)
       header(i) = fptr_header(i)
    enddo

    call hdo_free_wcsprm(C_LOC(fptr_header))
    
  end subroutine get_wcsheader
  
  
 


  subroutine world_to_pixel_wcsheader(world, pixel)
    implicit none
    real(C_DOUBLE), dimension(:,:), intent(in) :: world
    real(C_DOUBLE), dimension(size(world,1),size(world,1)), intent(out) :: pixel

    real(C_DOUBLE), dimension(size(world,2)) :: phi, theta
    real(C_DOUBLE), dimension(size(world,1),size(world,2)) :: imgcrd

    integer(C_INT) :: nelem, ncoord
    
    nelem = int(size(world,1),C_INT)
    ncoord = int(size(world,2),C_INT)
    
    call s2p_wcsprm(ncoord, nelem, world, phi, theta, imgcrd, pixel)

    if (display) then
       write(*,*)'phi= theta= ',phi,theta
       write(*,*)'imgcrd=     ',imgcrd
    end if
    
  end subroutine world_to_pixel_wcsheader


  subroutine pixel_to_world_wcsheader(pixel, world)
    implicit none
    real(C_DOUBLE), dimension(:,:), intent(in) :: pixel
    real(C_DOUBLE), dimension(size(pixel,1),size(pixel,2)), intent(out) :: world

    real(C_DOUBLE), dimension(size(pixel,1)) :: phi, theta
    real(C_DOUBLE), dimension(size(pixel,1),size(pixel,2)) :: imgcrd

    integer(C_INT) :: nelem, ncoord
    
    nelem = int(size(pixel,1),C_INT)
    ncoord = int(size(pixel,2),C_INT)
    
    call p2s_wcsprm(ncoord, nelem, pixel, imgcrd, phi, theta, world)

    if (display) then
       write(*,*)'imgcrd=     ',imgcrd
       write(*,*)'phi= theta= ',phi,theta
    end if
    
  end subroutine pixel_to_world_wcsheader


  subroutine galactic_to_equatorial_wcsheader()
    implicit none

!J2000 equatorial coordinates of the galactic pole + longitude (PA)    
    real(C_DOUBLE), parameter :: lng2p1 = 122.9319
    real(C_DOUBLE), parameter :: lat2p1 = 27.1283
    real(C_DOUBLE), parameter :: lng1p2 = 12.860114

!to equatorial J2000    
    character(kind=C_CHAR, len=*), parameter :: clng = 'RA--'//C_NULL_CHAR
    character(kind=C_CHAR, len=*), parameter :: clat = 'DEC-'//C_NULL_CHAR
    character(kind=C_CHAR, len=*), parameter :: radesys = 'ICRS'
    real(C_DOUBLE), parameter :: equinox = 2000.0
    character(C_CHAR), parameter :: alt = C_NULL_CHAR
    
    call ccs_wcsprm(lng2p1, lat2p1, lng1p2, clng, clat, radesys, equinox, alt)
    

  end subroutine galactic_to_equatorial_wcsheader
  
  
  
end module fwcs
