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


module precision
  implicit none
  
  public
  integer, parameter :: fsp=kind(1.0_4)
  integer, parameter :: fdp=kind(1.0_8)
  integer, parameter :: fqp=kind(1.0_16)
  
  integer, parameter :: isp=4
  integer, parameter :: idp=8

  real(fsp), parameter :: pisp = 3.141592653589793238462643383279502884197169399375105820974944592_fsp
  real(fdp), parameter :: pidp = 3.141592653589793238462643383279502884197169399375105820974944592_fdp

end module precision
