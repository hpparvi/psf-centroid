!!=== Lorentzian profile module ===
!!
!! Routines to calculate an undersampled Lorentzian profile (where FWHM ~ 1 pixel)
!! accurately by integrating the profile analytically over each pixel.
!!
!! -GPL-
!!
!! Copyright (C) 2013--2016  Hannu Parviainen
!!
!! This program is free software: you can redistribute it and/or modify
!! it under the terms of the GNU General Public License as published by
!! the Free Software Foundation, either version 3 of the License, or
!! (at your option) any later version.
!!
!! This program is distributed in the hope that it will be useful,
!! but WITHOUT ANY WARRANTY; without even the implied warranty of
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!! GNU General Public License for more details.
!!
!! You should have received a copy of the GNU General Public License
!! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!! -GPL-
!!
!! Author
!!  Hannu Parviainen <hannu.parviainen@physics.ox.ac.uk>
!!
!! Date 
!!  4.07.2016
!!
module lorentzian
  use omp_lib
  implicit none
  
  real(8), parameter :: PI = 3.14159265359d0
  real(8), parameter :: H_SQRT_PI = 0.5d0*sqrt(PI)
  real(8), parameter :: LLC = 0.5d0*log(2.0d0*PI)
  real(8), parameter :: FWHM_TO_SIGMA = 1.0d0/(2.0d0*sqrt(2.0d0*log(2.0d0)))

contains
  !! One-dimensional Lorentzian
  !! --------------------------
  subroutine lorentzian1d(center, amplitude, fwhm, npx, flux)
    implicit none
    real(8), intent(in)  :: center, amplitude, fwhm
    integer, intent(in)  :: npx
    real(8), intent(out) :: flux(npx)

    real(8) :: i1, i2, hfwhm
    integer :: i, wstart, wwidth

    flux  = 0.0d0
    hfwhm = 0.5d0*fwhm
    wwidth = min(npx, 10*ceiling(fwhm))
    wstart = max(0, floor(center) - wwidth/2)
    wwidth = min(npx-wstart, wwidth)

    i1 = hfwhm * atan((wstart-center-0.5d0)/hfwhm)
    do i = 1, wwidth
       i2 = hfwhm * atan((wstart+i-center-0.5d0)/hfwhm)
       flux(wstart+i) = amplitude*(i2-i1)
       i1 = i2
    end do
  end subroutine lorentzian1d

  !! Multithreaded one-dimensional Gaussian 
  !! --------------------------------------
  subroutine lorentzian1dmt(center, amplitude, fwhm, npx, nthr, flux)
    implicit none
    real(8), intent(in)  :: center, amplitude, fwhm
    integer, intent(in)  :: npx, nthr
    real(8), intent(out) :: flux(npx)

    real(8), allocatable :: ai(:)
    real(8) :: hfwhm
    integer :: i, wstart, wwidth

    !$ call omp_set_num_threads(nthr)

    flux  = 0.0d0
    hfwhm = 0.5d0*fwhm
    wwidth = min(npx, 10*ceiling(fwhm))
    wstart = max(0, floor(center) - wwidth/2)
    wwidth = min(npx-wstart, wwidth)

    allocate(ai(wwidth+1))
    !$omp parallel default(none) shared(flux,ai,hfwhm,center,amplitude,wstart,wwidth) private(i)
    !$omp do
    do i = 0, wwidth
       ai(i+1) = hfwhm * atan((wstart+i-center-0.5d0)/hfwhm) 
    end do
    !$omp end do
    !$omp do
    do i = 1, wwidth
       flux(wstart+i) = amplitude*(ai(i+1)-ai(i))
    end do
    !$omp end do
    !$omp end parallel 
    deallocate(ai)
  end subroutine lorentzian1dmt

  !! Multiple one-dimensional Gaussian profiles
  !! ------------------------------------------
  subroutine lorentzians1d(centers, amplitudes, fwhm, npx, nlines, flux)
    implicit none
    integer, intent(in)  :: npx, nlines
    real(8), intent(in)  :: centers(nlines), amplitudes(nlines), fwhm
    real(8), intent(out) :: flux(npx)

    real(8) :: e1, e2, hfwhm
    integer :: i, iline, wstart, wwidth

    flux  = 0.0d0
    hfwhm = 0.5d0 * fwhm
    wwidth = min(npx, 10*ceiling(fwhm))

    do iline = 1, nlines
       wstart = max(0, floor(centers(iline)) - wwidth/2)
    
       e1 = hfwhm * atan((wstart-centers(iline)-0.5d0)/hfwhm)
       do i = 1, min(npx-wstart, wwidth)
          e2 = hfwhm * atan((wstart+i-centers(iline)-0.5d0)/hfwhm)
          flux(wstart+i) = flux(wstart+i) + amplitudes(iline)*(e2-e1)
          e1 = e2
       end do
    end do
  end subroutine lorentzians1d


  real(8) function lnlike_lorentzian1d(center, amplitude, fwhm, error, sky, npx, fobs)
    implicit none
    integer, intent(in)  :: npx
    real(8), intent(in)  :: center,amplitude,fwhm,error,sky,fobs(npx)
    real(8) :: fmod(npx)

    call lorentzian1d(center, amplitude, fwhm, npx, fmod)
    lnlike_lorentzian1d = -LLC*npx - 0.5d0*npx*log(error**2) -0.5d0*sum((fobs - (sky+fmod))**2/error**2)
  end function lnlike_lorentzian1d

!!$  real(8) function lnlike_gaussians1d(centers, amplitudes, fwhm, error, sky, npx, nlines, fobs)
!!$    implicit none
!!$    integer, intent(in)  :: npx, nlines
!!$    real(8), intent(in)  :: centers(nlines), amplitudes(nlines), fwhm, error, sky, fobs(npx)
!!$    real(8) :: fmod(npx)
!!$
!!$    call gaussians1d(centers, amplitudes, fwhm, npx, nlines, fmod)
!!$    lnlike_gaussians1d = -LLC*npx - 0.5d0*npx*log(error**2) -0.5d0*sum((fobs - (sky+fmod))**2/error**2)
!!$  end function lnlike_gaussians1d

end module lorentzian
