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


module gaiaconst
  use precision, only : fdp

  public

!number of colors   
  integer, parameter :: ngaiacolors = 3
  
!conversion factors from GAIA fluxes to SI units (from
!https://www.cosmos.esa.int/web/gaia-users/archive/gedr3-documentation-pdf)
!to W/m^2/nm
  real(fdp), parameter :: clambdaR = 1.638483d-21
  real(fdp), parameter :: clambdaG = 1.346109d-21
  real(fdp), parameter :: clambdaB = 3.009167d-21

!to W/m^2/Hz  
  real(fdp), parameter :: cnuR = 3.298815d-33
  real(fdp), parameter :: cnuG = 1.736011d-33
  real(fdp), parameter :: cnuB = 2.620707d-33

!to photons/s/m^2
  real(fdp), parameter :: cphotR = 0.0064544055_fdp
  real(fdp), parameter :: cphotG = 0.0043306082_fdp
  real(fdp), parameter :: cphotB = 0.0078508254_fdp

  integer, parameter :: nsRGB = 3
!color transformation matrix from GAIA (RP G BP) flux (in physical
!units) to sRGB (linear)
  real(fdp), dimension(nsRGB,ngaiacolors) :: GAIA2sRGB = reshape( (/ &
       0.20739144_fdp,-0.01605617_fdp,-0.00222121_fdp, &
       0.72273442_fdp,0.67452044_fdp,0.55640397_fdp, &
       0.59382999_fdp,0.62803003_fdp,0.61312493_fdp &
       /), shape = (/nsRGB,ngaiacolors/) )

!identity matrix, map GAIA RP,G,BP into R,G,B
  real(fdp), dimension(ngaiacolors,ngaiacolors) :: I3 = reshape( (/ &
       1._fdp,0._fdp,0._fdp, &
       0._fdp,1._fdp,0._fdp, &
       0._fdp,0._fdp,1._fdp &
       /), shape = (/ngaiacolors,ngaiacolors/) )
  
end module gaiaconst
