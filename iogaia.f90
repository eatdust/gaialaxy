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


module iogaia
  use precision, only : fdp
  use gaiaconst, only : clambdaR, clambdaG, clambdaB
  use gaiaconst, only : cnuR, cnuG, cnuB
  use gaiaconst, only : cphotR, cphotG, cphotB
  use gaiaconst, only : ngaiacolors
  
  
  implicit none

  private

  
  abstract interface
     subroutine s2p(world,pixel)
       import fdp
       implicit none
       real(fdp), dimension(:,:), intent(in) :: world
       real(fdp), dimension(size(world,1),size(world,2)), intent(out) :: pixel
     end subroutine s2p

     subroutine gridding(nx,ny,pixel,flux,image)
       import fdp
       implicit none
       integer, intent(in) :: nx, ny
       real(fdp), dimension(2), intent(in) :: pixel
       real(fdp), intent(in) :: flux
       real(fdp), dimension(nx,ny), intent(inout) :: image
     end subroutine gridding

  end interface

  
  integer, parameter :: lentablename = 72
  integer, parameter :: lenttype = 72
  integer, parameter :: lentunit = 72

!max number of columns in the votables
  integer, parameter :: nfieldmax = 6

!tables indices for galactic coordinates and fluxes (in photo-electron/s)
  integer, parameter :: iglon = 2
  integer, parameter :: iglat = 3
  integer, parameter :: igaiaR = 5
  integer, parameter :: igaiaG = 4
  integer, parameter :: igaiaB = 6



!default flux scaling coefficients
  real(fdp), save :: cR = 1._fdp
  real(fdp), save :: cG = 1._fdp
  real(fdp), save :: cB = 1._fdp
  
  
!which default gridding method to use
  procedure(gridding), pointer :: ptr_gridding_func => nearest_neighbour


  integer, parameter :: debuglevel = 0
  logical, parameter :: display = .true.
  
  
 
  public s2p
  public lentablename
  public iglon, iglat, igaiaG, igaiaR, igaiaB, ngaiacolors

  public set_gridding_method, set_flux_units
  public accumulate_gaia_votable
  
  
contains

  subroutine set_gridding_method(cmethod)
    implicit none
    character(len=*), intent(in) :: cmethod

    select case (cmethod)

    case ('cic')
       ptr_gridding_func => cloud_in_cell

    case('nearest')
       ptr_gridding_func => nearest_neighbour

    case default
       stop 'set_gridding_method: not implemented!'
    end select
       
  end subroutine set_gridding_method


  
  subroutine set_flux_units(cunits)
    implicit none
    character(len=*), intent(in) :: cunits

    select case (cunits)

    case('raw')
       cR = 1._fdp
       cG = 1._fdp
       cB = 1._fdp
       if(display) write(*,*)'flux units: raw (e-/s)'
       
    case('wavelength')
       cR = clambdaR
       cG = clambdaG
       cB = clambdaB
       if(display) write(*,*)'flux units: wavelength (W/m^2/nm)'
       
    case('frequency')
       cR = cnuR
       cG = cnuG
       cB = cnuB
       if(display) write(*,*)'flux units: frequency (W/m^2/Hz)'
       
    case('photon')
       cR = cphotR
       cG = cphotG
       cB = cphotB
       if(display) write(*,*)'flux units: photon (photo/s/m^2)'
       
    case default
       stop 'set_flux_units: not found!'
    end select

  end subroutine set_flux_units
  
  
!given a list of fits table files, a color transformation matrix
!(applied to gaiaRGB after flux scaling), the table indices to use as
!world coordinates, a world to pixel converter, a transformed color
!index to use, this routine stacks the resulting flux into an already
!allocated image array (this array *is not* reset on entry and should
!be filled with 0 at first call)
  subroutine accumulate_gaia_votable(tablenames, colormatrix, ix, iy, ptr_s2p_func, icolor, image)
    use iofits, only : print_table_hdr
    use iofits, only : open_table_fits, stream_table_fits, close_table_fits
    implicit none
    character(len=lentablename), dimension(:), intent(in) :: tablenames
    real(fdp), dimension(:,:), intent(in) :: colormatrix

    integer, intent(in) :: ix, iy
    procedure(s2p), pointer :: ptr_s2p_func

    integer, intent(in) :: icolor
    real(fdp), dimension(:,:), intent(inout) :: image

    
    integer :: nx, ny, nfiles
    character(len=lentablename) :: tablename        
    character(len=lenttype), dimension(nfieldmax) :: ttype
    character(len=lentunit), dimension(nfieldmax) :: tunit

    real(fdp), dimension(ngaiacolors) :: gaiaRGB
    real(fdp), dimension(:,:), allocatable :: buffer

    real(fdp), dimension(size(colormatrix,1)) :: colorfluxes
    real(fdp), dimension(:,:), allocatable :: world, pixel

    integer :: i,j,k
    integer :: unit, ncut, nmod, nread
    integer :: nrows,frow, nopt, tfields

    integer :: nstacked
    
    nfiles = size(tablenames,1)
    nx = size(image,1)
    ny = size(image,2)

    if ((size(colormatrix,2)).ne.ngaiacolors) then
       write(*,*)'accumulate_gaia_votable: colormatrix rows= ',size(colormatrix,2)
       stop
    endif    

    nstacked = 0
    
    do k=1,nfiles
       tablename = trim(adjustl(tablenames(k)))

       if (display) then
          write(*,*)
          write(*,*)'accumulate_gaia_votable:'
          write(*,*)'reading table:   ',tablename
          if (debuglevel==3) call print_table_hdr(tablename)
       end if
          
       call open_table_fits(tablename,unit,nrows,nopt,tfields,ttype,tunit)

       ncut =  nrows / nopt
!       ncut = nrows

       nmod = modulo(nrows,nopt)

       nstacked = nstacked + nrows
       
       if (display) then
          write(*,*)'nrows=        ',nrows
          write(*,*)'n_*  =        ',nstacked          
       end if
       
       if (debuglevel == 1) then
          write(*,*)'unit=            ',unit
          write(*,*)'tfields=         ',tfields      
          write(*,*)'nopt=            ',nopt
          write(*,*)'ncut=            ',ncut
          write(*,*)'nmod=            ',nmod
          write(*,*)
       end if

       
!$omp parallel &       
!$omp default(shared) &
!$omp private(buffer,world,pixel)

       allocate(buffer(nopt,tfields))
       allocate(world(2,nopt),pixel(2,nopt))

!$omp do &
!$omp private(i,j,frow,nread) &
!$omp private(gaiaRGB, colorfluxes)
       do j=1,ncut+1

          frow = 1 + nopt*(j-1)

          nread = nopt

          if (j.eq.ncut+1) then
             if (nmod.eq.0) cycle
             nread = nmod
          endif
          
          if (debuglevel==2) then
             write(*,*)'icut=             ',j
             write(*,*)'first row=        ',frow
             write(*,*)'nopt= nmod= nread= ',nopt,nmod,nread
             write(*,*)
          endif

!$omp critical(iolockfits)          
          call stream_table_fits(unit,frow,buffer(1:nread,:))
!$omp end critical(iolockfits)
          
          do i=1,nread
             world(1,i) = buffer(i,ix)
             world(2,i) = buffer(i,iy)
          enddo
          
          call ptr_s2p_func(world,pixel)

          do i=1,nread
             
             gaiaRGB = (/cR*buffer(i,igaiaR),cG*buffer(i,igaiaG),cB*buffer(i,igaiaB)/)
             colorfluxes = matmul(colormatrix,gaiaRGB)

             !          call nearest_neighbour(nx,ny,pixel,colorfluxes(icolor),image)
             call ptr_gridding_func(nx,ny,pixel(:,i),colorfluxes(icolor),image)

          enddo

       enddo
!$omp end do
       
       deallocate(buffer, world, pixel)
!$omp end parallel       
      
       
       call close_table_fits(unit)

    enddo

  end subroutine accumulate_gaia_votable


  subroutine nearest_neighbour(nx,ny,pixel,flux,image)
    implicit none
    integer, intent(in) :: nx, ny
    real(fdp), dimension(2), intent(in) :: pixel
    real(fdp), intent(in) :: flux
    real(fdp), dimension(nx,ny), intent(inout) :: image
    
    integer :: i,j    

!assumes that first pixel *center* has value 1. For pixels < 1.5, we
!get 1, for pixel >=1.5 we get 2
    i = int(pixel(1)+0.5_fdp)
    j = int(pixel(2)+0.5_fdp)
       
!periodize the celestial sphere (should not be needed)
    if (i.eq.0) then
       if (debuglevel >= 1) write(*,*)'nearest_neighbour: i= j= ',i,j,pixel(:)
       i = nx
    end if
    if (j.eq.0) then
       if (debuglevel >= 1) write(*,*)'nearest_neighbour: i= j= ',i,j,pixel(:)
       j = ny
    end if

    if ((i.gt.nx).or.(j.gt.ny)) then
       write(*,*)'i= j= pixel= ',i,j,pixel(:)
       stop
    end if
       
    image(i,j) = image(i,j) + flux
    
  end subroutine nearest_neighbour

  
  subroutine cloud_in_cell(nx,ny,pixel,flux,image)
    implicit none
    integer, intent(in) :: nx, ny
    real(fdp), dimension(2), intent(in) :: pixel
    real(fdp), intent(in) :: flux
    real(fdp), dimension(nx,ny), intent(inout) :: image
    
    integer :: ic,jc
    integer :: i,ip1,im1,j,jp1,jm1
    integer :: in,jn

    real(fdp) :: tx, ty, dx, dy
    
   
    ic = int(pixel(1)+0.5_fdp)
    jc = int(pixel(2)+0.5_fdp)

    dx = pixel(1) - real(ic,fdp)
    dy = pixel(2) - real(jc,fdp)

    tx = 1._fdp - abs(dx)
    ty = 1._fdp - abs(dy)
    
    if ((tx.lt.0._fdp).or.(ty.lt.0._fdp)) then
       stop 'cloud_in_cell: improper weights!'
    endif

    if ((tx.gt.1._fdp).or.(ty.gt.1._fdp)) then
       stop 'cloud_in_cell: improper weights!'
    endif
    
!find the closest neighbours on a periodic grid
    i=ic       
    ip1 = modulo(ic + 1,nx)
    im1 = modulo(ic - 1,nx)
    j=jc
    jp1 = modulo(jc + 1,ny)
    jm1 = modulo(jc - 1,ny)

!periodize
    if (i.eq.0) i = nx
    if (ip1.eq.0) ip1 = nx
    if (im1.eq.0) im1 = nx
    if (j.eq.0) j = ny
    if (jp1.eq.0) jp1 = ny
    if (jm1.eq.0) jm1 = ny

!determine the 4 closest neighbour indices
    if (dx.gt.0._fdp) then
       in = ip1
    else
       in = im1
    endif

    if (dy.gt.0._fdp) then
       jn = jp1
    else
       jn = jm1
    endif

    image(i,j) = image(i,j) + tx*ty*flux

    image(i,jn) = image(i,jn) + tx*abs(dy)*flux

    image(in,j) = image(in,j) + abs(dx)*ty*flux

    image(in,jn) = image(in,jn) + abs(dx)*abs(dy)*flux

  end subroutine cloud_in_cell
    
end module iogaia
