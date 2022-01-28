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


program gaialaxy

  use precision, only : fdp, pidp
  use fwcs, only : cd, lenheader
  use fwcs, only : create_platecarree_wcsheader
  use fwcs, only : print_wcsheader, get_wcsheader
  use fwcs, only : world_to_pixel_wcsheader

  use iogaia, only : s2p
  use iogaia, only : lentablename
  use iogaia, only : iglon, iglat, ngaiacolors
  use iogaia, only : set_gridding_method, set_flux_units
  use iogaia, only : accumulate_gaia_votable

  use gaiaconst, only : GAIA2sRGB
  
  use iofits, only : write_wcsimage_fits
  
  implicit none

  
  integer :: nx,ny
  real(fdp), dimension(:,:), allocatable :: image

  real(cd), dimension(2) :: pixangles
  real(fdp) :: pixsterad
  
  integer :: icolor
  integer, parameter :: ncolors = 3
  real(fdp), dimension(ncolors,ngaiacolors) :: colormatrix

  character(len=40) :: outfitsname
  character(len=2), dimension(ncolors), parameter :: suffix = (/'_R','_G','_B'/)
    
  character(len=lentablename), dimension(:), allocatable :: tablenames
  procedure(s2p), pointer :: ptr_s2p => null()

  character(len=lenheader), dimension(:), allocatable :: header


  logical, parameter :: ask4input = .true.
  
    
!image size  
  nx = 32768
  ny = 16384

! color channel R(1), G(2), B(3)  
  icolor = 3

  if (ask4input) call get_args()

!search data files in that directory (as "gaiadata_01.fits,...")
  call get_votable_files('data/')

  
!image output in linear sRGB color space  
  colormatrix = GAIA2sRGB
!native GAIA filters
  !colormatrix = I3
  
  
  write(*,*)
  write(*,*)'creating WCS galactic header...'
  call create_platecarree_wcsheader(nx,ny,pixangles)
  call print_wcsheader()
  call get_wcsheader(header)

  write(*,*)
  write(*,*)'setting stacking method and units...'
  allocate(image(nx,ny))
  image = 0._fdp
  
  ptr_s2p => s2p_converter

  call set_gridding_method('cic')

  call set_flux_units('wavelength')

  call accumulate_gaia_votable(tablenames,colormatrix,iglon,iglat,ptr_s2p,icolor,image)

  write(*,*)
  write(*,*)'normalizing flux in sr^-1 and dumping image...'

  outfitsname = 'gaialaxy'//suffix(icolor)//'.fits'
  pixsterad = product(pixangles)*(pidp/180.0_fdp)**2
  write(*,*)'pixel solid angle (sr)= ',pixsterad

  call write_wcsimage_fits(trim(outfitsname),image/pixsterad,header)

  deallocate(tablenames)
  deallocate(header)
  deallocate(image)

contains

  
    
  subroutine get_args()
    implicit none

    write(*,*)'image dimensions: (nx x ny)'
    write(*,*)'nx=? ny=? '
    read(*,*)nx,ny
    write(*,*)'color channel: (1,2,3)'
    write(*,*)'icolor=? '
    read(*,*)icolor

  end subroutine get_args
  
  

  subroutine get_votable_files(cdir)
    implicit none
    character(len=*), intent(in) :: cdir

    integer, parameter :: nmaxfiles = 99

    integer, parameter :: lennum = 2
    character(len=lennum) :: cnum

    character(len=*), parameter :: prefix = 'gaiadata_'

    character(len=lentablename) :: filename
    character(len=lentablename), dimension(nmaxfiles) :: foundnames

    logical :: here
    integer :: i,nfound

    logical, parameter :: display = .false.

    nfound = 0
    do i=1,nmaxfiles

       call int2char_zero(i,lennum,cnum)

       filename = cdir//prefix//cnum//'.fits'

       inquire(file=filename,exist=here)

       if (.not.here) cycle

       nfound = nfound + 1
       foundnames(nfound) = filename

    enddo

    allocate(tablenames(nfound))
    tablenames = foundnames(1:nfound)

    if (display) then
       write(*,*)'get_votable_files:'
       write(*,*)'data files found are:'
       do i=1,nfound
          write(*,*)trim(tablenames(i))
       end do
    end if
       
  end subroutine get_votable_files



  subroutine int2char_zero(icount,clen,ccount)
    implicit none
    integer, intent(in) :: icount
    integer, intent(in) :: clen
    character(len=clen), intent(out) :: ccount

    integer, parameter :: leniomax = 64
    character(len=leniomax) :: numToStrg
    character(len=clen) :: strg

    if (clen.gt.leniomax) then
       stop 'int2char_zero: clen > lenIoMax!'
    endif
    
    write(numToStrg,*) icount
    strg = trim(adjustl(numToStrg)) 
    strg = adjustr(strg)
    call replace_char(strg,' ','0')

    ccount(1:clen) = strg(1:clen)

  end subroutine int2char_zero



  subroutine replace_char(strg,charold,charnew)
    implicit none
    character(len=*), intent(inout) :: strg
    character, intent(in) :: charold, charnew

    integer :: position
    position = 1

    do 
       position = scan(strg,charold)
       if (position.ne.0) then
          strg(position:position) = charnew
       else
          exit
       endif
    enddo

  end subroutine replace_char

      

  
  subroutine s2p_converter(world,pixel)
    implicit none
    real(fdp), dimension(:,:), intent(in) :: world
    real(fdp), dimension(size(world,1),size(world,2)), intent(out) :: pixel

    real(cd), dimension(:,:), allocatable :: cworld
    real(cd), dimension(:,:), allocatable :: cpixel


    if (cd == fdp) then
       call world_to_pixel_wcsheader(world, pixel)
       return
    endif

    allocate(cworld(size(world,1),size(world,2)))
    allocate(cpixel(size(world,1),size(world,2)))

    cworld = world
    call world_to_pixel_wcsheader(world, cpixel)
    pixel = cpixel

    deallocate(cworld,cpixel)        
    
  end subroutine s2p_converter


  subroutine s2p_test(world,pixel)
    implicit none
    real(fdp), dimension(:,:), intent(in) :: world
    real(fdp), dimension(size(world,1),size(world,2)), intent(out) :: pixel

    
    pixel(1,:) = world(1,:)*nx/360._fdp
    pixel(2,:) = (world(2,:)+90.0)*ny/180._fdp
    
  end subroutine s2p_test
  
  
end program gaialaxy
