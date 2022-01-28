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


module iofits
  use precision


  private

  logical :: verbose = .true.
  integer, parameter :: lenerr = 30

  integer, parameter :: lenrec = 80

  interface read_bintable_key
     module procedure read_floattable_key
     module procedure read_doubletable_key
  end interface read_bintable_key
  
  interface append_table_fits
     module procedure append_doubletable_fits
     module procedure append_floattable_fits
  end interface append_table_fits

  interface read_table_fits
     module procedure read_doubletable_fits
     module procedure read_floattable_fits
  end interface read_table_fits

  interface close_table_fits
     module procedure close_doubletable_fits
     module procedure close_floattable_fits
     module procedure close_anytable_fits
  end interface close_table_fits

  interface open_table_fits
     module procedure open_bintable_fits
  end interface open_table_fits

  interface stream_table_fits
     module procedure stream_doubletable_fits
  end interface stream_table_fits

  
  public write_twod_fits, check_fits_header, read_twod_fits

  public write_wcsimage_fits

  public read_xtension_fits, read_bintable_key, read_asctable_fits 

  public print_table_hdr
  public create_table_fits, close_table_fits
  public append_table_fits, read_table_fits

  public open_table_fits, stream_table_fits

    
  
contains  

  subroutine delete_fitsfile(filename,status)

    implicit none
    integer :: status,unit,blocksize
    character(len=*) :: filename

    if (status .gt. 0)return

    call ftgiou(unit,status)
    call ftopen(unit,filename,1,blocksize,status)
    
    if (status .eq. 0) then
       call ftdelt(unit,status)
    else if (status .eq. 103) then
       status=0
       call ftcmsg
    else
       status=0
       call ftcmsg
       call ftdelt(unit,status)
    end if
    
    call ftfiou(unit, status)
  end subroutine delete_fitsfile


  
  subroutine check_fits_header(filename,npixtot,nx,ny)
    implicit none
    character(len=*) :: filename    
    integer :: npixtot,nx,ny

    integer, dimension(2) :: naxes
    integer :: status,unit,readwrite,blocksize,nfound
 
    status=0
    call ftgiou(unit,status)
    readwrite=0
    call ftopen(unit,filename,readwrite,blocksize,status)
    call ftgknj(unit,'NAXIS',1,2,naxes,nfound,status)
    
    if (nfound .ne. 2) then
       write(*,*) 'ERROR in check_fits_header :'
       write(*,*) 'Failure in reading the NAXIS keywords'
       STOP
    endif
    
    nx=naxes(1)
    ny=naxes(2)
    npixtot=nx*ny
    if (verbose) write(*,*) 'Dimensions of the image detected ',nx,ny    
    call ftclos(unit, status)
    call ftfiou(unit, status)
    
  end subroutine check_fits_header


  subroutine read_twod_fits(filename,npixtot,nx,ny,map)
    implicit none
    character(len=*) :: filename
    integer :: npixtot,nx,ny
    real(fsp), dimension(0:npixtot-1) :: map

    integer, dimension(2) :: naxes
    integer :: status,unit,readwrite,blocksize,nfound,nbuffer
    integer :: group,firstpix,naxis1,naxis2
    logical :: anynull
    real(fsp) :: rnullval
    integer :: i,j,inc
    character(len=lenerr) :: errtext
    
    status=0
    call ftgiou(unit,status)
    readwrite=0
    
    call ftopen(unit,filename,readwrite,blocksize,status)
    call ftgknj(unit,'NAXIS',1,2,naxes,nfound,status)
    
    naxis1 = naxes(1)
    naxis2 = naxes(2)
    
    nbuffer=naxis1*naxis2
    group=1
    firstpix=1
    rnullval=1e30
    
    call ftgpve(unit,group,firstpix,nbuffer,rnullval,map,anynull,status)
    
    call ftgerr(status,errtext)
    call ftclos(unit, status)
    call ftfiou(unit, status)
    if (status > 0) then
       write(*,*) 'ERROR in read_twod_fits :',errtext
       STOP
    endif
    
  end subroutine read_twod_fits



  subroutine write_twod_fits(filename,map,npixtot,nx,ny)
    implicit none
    character(len=*) :: filename    
    integer, intent(in) :: npixtot,nx,ny
    real(fsp),  dimension(0:npixtot-1), intent(in) :: map
    
    integer, dimension(2) :: naxes
    integer :: status,unit,blocksize,bitpix,naxis
    integer :: i,j,group,fpixel,nelements
    logical :: simple,extend
    character(len=lenerr) :: errtext
    
    status=0
    call ftgiou(unit,status)
    blocksize=1
    call ftinit(unit,filename,blocksize,status)
    simple=.true.
    bitpix=-32
    naxis=2
    naxes(1)=nx
    naxes(2)=ny
    extend=.true.
    call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
    group=1
    fpixel=1
    nelements=npixtot
    call ftppre(unit,group,fpixel,nelements,map,status)
    call ftclos(unit, status)
    call ftfiou(unit, status)
    if (status > 0) then
       write(*,*) 'ERROR in write_twod_fits :',errtext
       STOP
    endif
    
  end subroutine write_twod_fits



  subroutine write_wcsimage_fits(filename,image,wcsheader)
    implicit none
    character(len=*), intent(in) :: filename
    real(fdp), dimension(:,:), intent(in) :: image
    character(len=*), dimension(:), intent(in), optional :: wcsheader
    
    integer :: unit, status
    logical :: simple, extend
    
    integer, parameter :: naxis = 2
    integer, dimension(naxis) :: naxes
    integer ::  blocksize, bitpix
    
    integer(idp) :: nelements, group, fpixel
    

    integer :: i, nkeys
       
    status=0
    call ftgiou(unit,status)
    blocksize=1
    call ftinit(unit,filename,blocksize,status)

    if (status.ne.0) then
       write(*,*) 'write_wcs_image: cannot open file: ',filename
       stop
    end if
    
    simple=.true.
    bitpix=-32
    naxes(1)=size(image,1)
    naxes(2)=size(image,2)
    extend=.true.

    call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)

    if (present(wcsheader)) then
       nkeys = size(wcsheader,1)
       do i=1,nkeys
          call ftprec(unit,wcsheader(i),status)
          if (status.ne.0) then
             write(*,*)'issue writing card= ',wcsheader(i)
          end if
       end do
    end if

    group=1
    fpixel=1

    nelements = naxes(1)*naxes(2)

    !this punk overlflows for 64k images due to nelements not being
    !    long long...
    !call ftpprd(unit,group,fpixel,nelements,image,status)
    
    call ftp2dd(unit,group,naxes(1),naxes(1),naxes(2),image,status)
    
    call ftclos(unit,status)
    call ftfiou(unit,status)
    
  end subroutine write_wcsimage_fits
  
  
!open the ieme extension image of filename
  subroutine read_xtension_fits(filename,naxes,image)
    implicit none

    character(len=*), intent(in) :: filename
    integer, dimension(:), allocatable :: naxes
    real(fsp), dimension(:), allocatable :: image

    integer :: nhdu, status, readwrite
    character(len=lenerr) :: errtext

    integer :: naxis, unit
    integer :: group, firstpix, nelements
    logical :: anynull

    real(fsp) :: rnullval
    
    if (allocated(image)) then
       stop 'read_xtension_fits: image already associated!'
    endif

    if (allocated(naxes)) then
       stop 'read_xtension_fits: naxes already associated!'
    endif
    
    status=0
    call ftgiou(unit,status)
    readwrite=0
    
    call ftnopn(unit,filename,readwrite,status)
    call ftghdn(unit,nhdu)

    write(*,*)'read_xtension_image:'
    write(*,*)'opened hdu number: ',nhdu

    call ftgidm(unit,naxis,status)

    write(*,*)'naxis= ',naxis
    if (naxis.eq.0) stop 'naxis null!'
        
    allocate(naxes(naxis))

    call ftgisz(unit,naxis,naxes,status)


    nelements = product(naxes)
    write(*,*)'naxes= nelements= ',naxes,nelements

    allocate(image(nelements))
    
    group=1
    firstpix=1
    rnullval=1e30
    
    call ftgpve(unit,group,firstpix,nelements,rnullval,image,anynull,status)
    call ftgerr(status,errtext)
    call ftclos(unit, status)
    call ftfiou(unit, status)
    if (status > 0) then
       write(*,*) 'ERROR in read_xtension_fits :',errtext
       STOP
    endif
    
  end subroutine read_xtension_fits



  subroutine read_floattable_key(filename,keyname,keyval,nhdu)
    implicit none
    character(len=*), intent(in) :: filename
    character(len=*), intent(in) :: keyname
    real(fsp), intent(out) :: keyval
    integer, intent(in), optional :: nhdu
    
    character(len=lenrec) :: comment

    integer :: status, unit
    integer :: readwrite, hdutype


    status=0
    call ftgiou(unit,status)

    readwrite=0
    call fttopn(unit,filename,readwrite,status)

    if (present(nhdu)) then
       call ftmahd(unit,nhdu, hdutype,status)
       write(*,*)'opening nhdu= of type= ',nhdu,hdutype
    endif

    call ftgkye(unit,keyname,keyval,comment,status)

    call ftclos(unit, status)
    call ftfiou(unit, status)

    write(*,*)'read_floattable_key: ',comment

    if (status > 0) then
       write(*,*) 'ERROR in read_floattable_fits :',status
       STOP
    endif


  end subroutine read_floattable_key


  subroutine read_doubletable_key(filename,keyname,keyval,nhdu)
    implicit none
    character(len=*), intent(in) :: filename
    character(len=*), intent(in) :: keyname
    real(fdp), intent(out) :: keyval
    integer, intent(in), optional :: nhdu
    
    character(len=lenrec) :: comment

    integer :: status, unit
    integer :: readwrite, hdutype


    status=0
    call ftgiou(unit,status)

    readwrite=0
    call fttopn(unit,filename,readwrite,status)

    if (present(nhdu)) then
       call ftmahd(unit,nhdu, hdutype,status)
       write(*,*)'opening nhdu= of type= ',nhdu,hdutype
    endif

    call ftgkyd(unit,keyname,keyval,comment,status)

    call ftclos(unit, status)
    call ftfiou(unit, status)

    write(*,*)'read_doubletable_key: ',comment

    if (status > 0) then
       write(*,*) 'ERROR in read_doubletable_key :',status
       STOP
    endif


  end subroutine read_doubletable_key


  


  subroutine read_asctable_fits(filename,ttype,tunit,sdata)
    implicit none
    character(len=*), intent(in) :: filename

    character(len=*), dimension(:), intent(out) :: ttype
    character(len=*), dimension(:), intent(out) :: tunit

    character(len=3), dimension(size(ttype,1)) :: tform
    character(len=*), dimension(:,:), intent(inout) :: sdata
    
    real(fsp) :: buffer
    integer :: nullval
    logical :: anynull

    integer :: icol
    character(len=18) :: comment, extname

    character :: nullstrg
    
    integer :: maxdim, nrow, tfields
    integer :: status, unit
    integer :: readwrite, varidat
    integer :: dispw

    integer :: frow, felem


    maxdim = size(ttype,1)
    if (size(tunit,1).ne.maxdim) stop 'read_asctable_fits: tunit/ttype size missmatch!'
    
    status=0
    call ftgiou(unit,status)

    readwrite=0
    call ftnopn(unit,filename,readwrite,status)
      
    call ftghbn(unit,maxdim,nrow,tfields,ttype,tform,tunit,extname,varidat,status)

    if ((maxdim.ne.tfields)) then
       write(*,*)'tfields= ',tfields,maxdim
       stop 'read_asctable_fits: TFIELDS missmatch!'
    endif    

    do icol=1,tfields
       write(*,*)'icol= tform= tunit= ',icol,tform(icol)
    enddo

    nullstrg = ' '
    frow = 1
    felem = 1

    anynull = .false.

    do icol=1,tfields
       call ftgcdw(unit,icol,dispw,status)

       if (dispw.gt.len(sdata(1,1))) then
          write(*,*)'display w= ',icol, dispw
          stop
       endif

       call ftgcvs(unit,icol,frow,felem,nrow,nullstrg,sdata(:,icol),anynull,status)

       if (anynull) then
          write(*,*)'icol= ',icol
          stop 'read_ascrable_fits: undefined value read!'
       endif
    enddo
   
    call ftclos(unit, status)
    call ftfiou(unit, status)
    
    if (status > 0) then
       write(*,*) 'ERROR in read_asctable_fits :',status
       STOP
    endif


  end subroutine read_asctable_fits



  function create_doubletable_fits(filename,ttype,tunit,unit)
    implicit none
    integer :: create_doubletable_fits
    character(len=*), intent(in) :: filename
    character(len=*), dimension(:), intent(in) :: ttype
    character(len=*), dimension(size(ttype,1)), intent(in) :: tunit
    integer, intent(in), optional :: unit
    
    integer :: unitf
    integer :: blocksize ,status
    
    if (.not.present(unit)) then

       call ftgiou(unitf,status)
       blocksize=1
       call ftinit(unitf,filename,blocksize,status)

    else

       unitf=unit

    endif

    
    create_doubletable_fits = create_kindtable_fits(filename,ttype,tunit,unitf,'1D')

  end function create_doubletable_fits



  function create_floattable_fits(filename,ttype,tunit,unit)
    implicit none
    integer :: create_floattable_fits
    character(len=*), intent(in) :: filename
    character(len=*), dimension(:), intent(in) :: ttype
    character(len=*), dimension(size(ttype,1)), intent(in) :: tunit
    integer, intent(in), optional :: unit

    integer :: unitf
    integer :: blocksize ,status
    
    
    if (.not.present(unit)) then

       call ftgiou(unitf,status)
       blocksize=1
       call ftinit(unitf,filename,blocksize,status)

    else

       unitf=unit

    endif

       
    create_floattable_fits = create_kindtable_fits(filename,ttype,tunit,unitf,'1M')

  end function create_floattable_fits
  


  function create_table_fits(filename,ttype,tunit,tkind,unit)
    implicit none
    integer :: create_table_fits
    character(len=*), intent(in) :: filename
    character(len=*), dimension(:), intent(in) :: ttype
    character(len=*), dimension(size(ttype,1)), intent(in) :: tunit
    integer, intent(in) :: tkind
    integer, intent(in), optional :: unit

    integer :: unitf
    integer :: blocksize ,status
    character(len=2) :: charkind
    

    
    if (tkind == kind(1.0_4)) then
       charkind = '1M'
    elseif (tkind == kind(1.0_8)) then
       charkind = '1D'
    else
       write(*,*)'kind= ',tkind,kind(1.0_4),kind(1.0_8)
       stop 'create_table_fits: kind not implemented'
    endif

    if (.not.present(unit)) then

       call ftgiou(unitf,status)
       blocksize=1
       call ftinit(unitf,filename,blocksize,status)

    else

       unitf=unit

    endif

       
    create_table_fits = create_kindtable_fits(filename,ttype,tunit,unitf,charkind)

  end function create_table_fits
  
  
  
  
  function create_kindtable_fits(filename,ttype,tunit,unitf,charform)
    implicit none
    integer :: create_kindtable_fits
    character(len=*), intent(in) :: filename
    character(len=*), dimension(:), intent(in) :: ttype
    character(len=*), dimension(size(ttype,1)), intent(in) :: tunit
    integer, intent(in) :: unitf
    character(len=2), intent(in) :: charform
    
    character(len=2), dimension(size(ttype,1)) :: tform

    integer :: nrow, tfields, nelements
    integer :: status, blocksize

    integer, parameter :: varidat= 0
    character, parameter :: extname= ''

    if ((charform.ne.'1D').and.(charform.ne.'1M')) then
       write(*,*)'charform= ',charform
       stop 'create_kindtable_fits: tform unknown!'
    endif

    nrow=0
    tfields=size(ttype,1)
    
    tform(:)=charform

       
    status=0
        
    !put binary into HDU=1
    call ftibin(unitf,nrow,tfields,ttype,tform,tunit,extname,varidat,status)

    create_kindtable_fits = int(unitf)

    if (status > 0) then
       write(*,*) 'ERROR in create_table_fits :',status
       stop
    endif

  end function create_kindtable_fits



  subroutine append_floattable_fits(unit,data)
    implicit none
    integer, intent(in) :: unit
    real(fsp), dimension(:,:), intent(in) :: data

    integer :: tfields, status
    integer :: frow, nrows, nelems
    integer, parameter :: felem = 1
    integer :: icol, unitf

    tfields=size(data,2)
    unitf = unit
    nelems = size(data,1)

    status = 0

    call ftgnrw(unitf,nrows,status)

    frow = nrows + 1

    do icol=1,tfields
       call ftpcle(unitf,icol,frow,felem,nelems,data(:,icol),status)
    enddo


    if (status > 0) then
       write(*,*) 'ERROR in append_doubletable_fits :',status
       stop
    endif

  end subroutine append_floattable_fits


  
  subroutine append_doubletable_fits(unit,data)
    implicit none
    integer, intent(in) :: unit
    real(fdp), dimension(:,:), intent(in) :: data

    integer :: tfields, status
    integer :: frow, nrows, nelems
    integer, parameter :: felem = 1
    integer :: icol, unitf

    tfields=size(data,2)
    unitf = unit
    nelems = size(data,1)

    status = 0

    call ftgnrw(unitf,nrows,status)

    frow = nrows + 1

    do icol=1,tfields
       call ftpcld(unitf,icol,frow,felem,nelems,data(:,icol),status)
    enddo

    if (status > 0) then
       write(*,*) 'ERROR in append_doubletable_fits :',status
       stop
    endif

  end subroutine append_doubletable_fits

  
  

  subroutine close_floattable_fits(unit,keyname,keyval,comment)
    implicit none
    integer, intent(in) :: unit
    character(len=8), dimension(:), intent(in) :: keyname
    real(fsp), dimension(:), intent(in) :: keyval
    character(len=*), dimension(:), intent(in) :: comment

    integer :: status, unitf
    integer :: nkeys, ikey

    integer, parameter :: decimals=-16
    logical, parameter :: addchksum = .true.
    
    status = 0
    unitf = unit

    
    nkeys = size(keyname,1)
    if (size(keyval,1).ne.nkeys) then
       stop 'close_floattable_fits: keysize missmatch!'
    endif

    do ikey=1,nkeys
       call ftpkye(unitf,keyname(ikey),keyval(ikey),decimals,comment(ikey),status)
    enddo

    if (addchksum) call ftpcks(unitf, status)
    call close_anytable_fits(unitf)
    
  end subroutine close_floattable_fits


  


  subroutine close_doubletable_fits(unit,keyname,keyval,comment)
    implicit none
    integer, intent(in) :: unit
    character(len=8), dimension(:), intent(in) :: keyname
    real(fdp), dimension(:), intent(in) :: keyval
    character(len=*), dimension(:), intent(in) :: comment

    integer :: status, unitf
    integer :: nkeys, ikey

    integer, parameter :: decimals=-16
    logical, parameter :: addchksum = .true.
    
    
    status = 0
    unitf = unit

    
    nkeys = size(keyname,1)
    if (size(keyval,1).ne.nkeys) then
       stop 'close_doubletable_fits: keysize missmatch!'
    endif

    do ikey=1,nkeys
       call ftpkyd(unitf,keyname(ikey),keyval(ikey),decimals,comment(ikey),status)
    enddo

    if (addchksum) call ftpcks(unitf, status)

    call close_anytable_fits(unitf)
    
  end subroutine close_doubletable_fits


  

  subroutine close_anytable_fits(unitf)
    implicit none
    integer, intent(in) :: unitf

    integer :: status

    status = 0    

    call ftclos(unitf, status)
    call ftfiou(unitf, status)

    if (status > 0) then
       write(*,*) 'ERROR in close_table_fits :',status
       STOP
    endif

  end subroutine close_anytable_fits


  

  
  
  subroutine print_table_hdr(filename,nhdu)
    implicit none
    character(len=*), intent(in) :: filename
    integer, intent(in), optional :: nhdu
    
    character(len=lenrec) :: record
    

    integer :: status, unit
    integer :: readwrite
    integer :: hdutype
    
    integer :: ikey
    integer :: nkeys

   
    status=0
    call ftgiou(unit,status)

    readwrite=0
    call fttopn(unit,filename,readwrite,status)

    if (present(nhdu)) then
       call ftmahd(unit,nhdu, hdutype,status)
       write(*,*)'opening nhdu= of type= ',nhdu,hdutype
    endif
       
    call ftghps(unit,nkeys,ikey,status)

    do ikey=1,nkeys
       call ftgrec(unit,ikey,record,status)
       write(*,*) record
    end do

    call ftclos(unit, status)
    call ftfiou(unit, status)

    if (status > 0) then
       write(*,*) 'ERROR in print_table_fits :',status
       STOP
    endif
    
  end subroutine print_table_hdr


  
  subroutine read_floattable_fits(filename,ttype,tunit,data,keyname,keyval,nhdu)
    implicit none
    character(len=*), intent(in) :: filename

    character(len=*), dimension(:), intent(out) :: ttype
    character(len=*), dimension(:), intent(out) :: tunit
    real(fsp), dimension(:,:), intent(out) :: data

    character(len=2), dimension(size(ttype,1)) :: tform

    character(len=8), dimension(:), intent(in), optional :: keyname
    real(fsp), dimension(:), intent(out), optional :: keyval

    integer, intent(in), optional :: nhdu
    
    real :: buffer
    real, parameter :: nullval = -1
    logical :: anynull

    integer :: nkeys, ikey, icol
    character(len=lenrec) :: comment, extname

    integer :: maxdim, nrow, tfields
    integer :: status, unit
    integer :: readwrite, varidat, hdutype

    integer, parameter :: frow=1, felem=1    


    logical, parameter :: display = .true.

    maxdim = size(ttype,1)
    if (size(tunit,1).ne.maxdim) stop 'read_floattable_fits: tunit/ttype size missmatch!'
    
    status=0
    call ftgiou(unit,status)

    readwrite=0
    call fttopn(unit,filename,readwrite,status)

    if (present(nhdu)) then
       call ftmahd(unit,nhdu, hdutype,status)
       write(*,*)'opening nhdu= of type= ',nhdu,hdutype
    endif
    
    call ftghbn(unit,maxdim,nrow,tfields,ttype,tform,tunit,extname,varidat,status)

    if ((maxdim.ne.tfields).or.(size(data,2).ne.tfields)) then
       stop 'read_floattable_fits: TFIELDS missmatch!'
    endif

    if (size(data,1).ne.nrow) stop 'read_doubletable_fits: NROW missmatch!'

    do icol=1,tfields
       if (tform(icol).ne.'1M') then
          write(*,*)'icol= tform= ',icol,tform(icol)
          stop 'read_floattable_fits: data type missmatch!'
       endif
    enddo
    
    anynull = .false.

   
    do icol=1,tfields
       call ftgcve(unit,icol,frow,felem,nrow,nullval,data(:,icol),anynull,status)
       if (anynull) then
          write(*,*)'icol= ',icol
          stop 'read_floattable_fits: undefined value read!'
       endif
    enddo
    
    if (present(keyname)) then
       if (.not.present(keyval)) then
          stop 'read_floattable_fits: missing value!'
       endif
       nkeys = size(keyname,1)
       if (size(keyval,1).ne.nkeys) then
          stop 'read_floattable_fits: keysize missmatch!'
       endif

       do ikey=1,nkeys
          call ftgkye(unit,keyname(ikey),keyval(ikey),comment,status)
          if (display) write(*,*)'keyname= keyval= ', keyname,keyval
       enddo
       
    end if

    call ftclos(unit, status)
    call ftfiou(unit, status)
    
    if (status > 0) then
       write(*,*) 'ERROR in read_floattable_fits :',status
       STOP
    endif


  end subroutine read_floattable_fits

  

  subroutine read_doubletable_fits(filename,ttype,tunit,data,keyname,keyval,nhdu)
    implicit none
    character(len=*), intent(in) :: filename

    character(len=*), dimension(:), intent(out) :: ttype
    character(len=*), dimension(:), intent(out) :: tunit
    real(fdp), dimension(:,:), intent(out) :: data

    character(len=2), dimension(size(ttype,1)) :: tform

    character(len=8), dimension(:), intent(in), optional :: keyname
    real(fdp), dimension(:), intent(out), optional :: keyval

    integer, intent(in), optional :: nhdu
    
    real :: buffer
    real, parameter :: nullval = -1
    logical :: anynull

    integer :: nkeys, ikey, icol
    character(len=lenrec) :: comment, extname

    integer :: maxdim, nrow, tfields
    integer :: status, unit
    integer :: readwrite, varidat, hdutype

    integer, parameter :: frow=1, felem=1    


    logical, parameter :: display = .true.

    maxdim = size(ttype,1)
    if (size(tunit,1).ne.maxdim) stop 'read_doubletable_fits: tunit/ttype size missmatch!'
    
    status=0
    call ftgiou(unit,status)

    readwrite=0
    call fttopn(unit,filename,readwrite,status)

    if (present(nhdu)) then
       call ftmahd(unit,nhdu, hdutype,status)
       write(*,*)'opening nhdu= of type= ',nhdu,hdutype
    endif
    
    call ftghbn(unit,maxdim,nrow,tfields,ttype,tform,tunit,extname,varidat,status)

    if ((maxdim.ne.tfields).or.(size(data,2).ne.tfields)) then
       stop 'read_doubletable_fits: TFIELDS missmatch!'
    endif

    if (size(data,1).ne.nrow) stop 'read_doubletable_fits: NROW missmatch!'

    do icol=1,tfields
       if (tform(icol).ne.'1D') then
          write(*,*)'icol= tform= ',icol,tform(icol)
          stop 'read_doubletable_fits: data type missmatch!'
       endif
    enddo
    
    anynull = .false.

    
    do icol=1,tfields
       call ftgcvd(unit,icol,frow,felem,nrow,nullval,data(:,icol),anynull,status)
       if (anynull) then
          write(*,*)'icol= ',icol
          stop 'read_doubletable_fits: undefined value read!'
       endif
    enddo
    
    
    if (present(keyname)) then
       if (.not.present(keyval)) then
          stop 'read_doubletable_fits: missing value!'
       endif
       nkeys = size(keyname,1)
       if (size(keyval,1).ne.nkeys) then
          stop 'read_doubletable_fits: keysize missmatch!'
       endif

       do ikey=1,nkeys
          call ftgkyd(unit,keyname(ikey),keyval(ikey),comment,status)
          if (display) write(*,*)'keyname= keyval= ', keyname,keyval
       enddo
       
    end if

    call ftclos(unit, status)
    call ftfiou(unit, status)
    
    if (status > 0) then
       write(*,*) 'ERROR in read_doubletable_fits :',status
       STOP
    endif


  end subroutine read_doubletable_fits

  

  
  subroutine open_bintable_fits(filename,unit,nrows,noptrows,tfields,ttype,tunit,nhdu)
    implicit none
    character(len=*), intent(in) :: filename

    integer, intent(out) :: unit,nrows,noptrows,tfields
    
    character(len=*), dimension(:), intent(out) :: ttype
    character(len=*), dimension(:), intent(out) :: tunit

    character(len=2), dimension(size(ttype,1)) :: tform

    integer, intent(in), optional :: nhdu
    
    character(len=lenrec) :: extname

    integer :: maxdim
    integer :: status
    integer :: readwrite, varidat, hdutype

    logical, parameter :: display = .true.

    maxdim = size(ttype,1)
    if (size(tunit,1).ne.maxdim) stop 'open_bintable_fits: tunit/ttype size missmatch!'
    
    status=0
    call ftgiou(unit,status)

    readwrite=0
    call fttopn(unit,filename,readwrite,status)

    if (present(nhdu)) then
       call ftmahd(unit,nhdu, hdutype,status)
       write(*,*)'opening nhdu= of type= ',nhdu,hdutype
    endif

    call ftgrsz(unit,noptrows,status)
    
    call ftghbn(unit,maxdim,nrows,tfields,ttype,tform,tunit,extname,varidat,status)

    if (maxdim.ne.tfields) then
       stop 'read_doubletable_fits: TFIELDS missmatch!'
    endif

    if (status > 0) then
       write(*,*) 'ERROR in open_bintable_fits :',status
       STOP
    endif


  end subroutine open_bintable_fits


  subroutine stream_doubletable_fits(unit,frow,data)
    implicit none
    integer, intent(in) :: unit,frow
    real(fdp), dimension(:,:), intent(out) :: data

    logical :: anynull
    integer :: icol, tfields, status, nelements, felem
    real, parameter :: nullval = -1.0
    
    nelements = size(data,1)
    tfields = size(data,2)
    
    felem = 1
    status = 0
    anynull = .false.
    
    do icol=1,tfields
       call ftgcvd(unit,icol,frow,felem,nelements,nullval,data(:,icol),anynull,status)
       if (anynull) then
          write(*,*)'icol= ',icol
          stop 'stream_doubletable_fits: undefined value read!'
       endif
       if (status.ne.0) then
          write(*,*)'stream_doubletable_fits: error= ',status
          stop
       end if
    enddo


  end subroutine stream_doubletable_fits

  
end module iofits
