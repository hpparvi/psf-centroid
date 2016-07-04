!!=== Gaussian profile module ===
!!
!! Routines to calculate an undersampled Gaussian profile (where sigma ~ 1 pixel)
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
module gaussian
  use omp_lib
  implicit none
  
  real(8), parameter :: PI = 3.14159265359d0
  real(8), parameter :: H_SQRT_PI = 0.5d0*sqrt(PI)
  real(8), parameter :: LLC = 0.5d0*log(2.0d0*PI)
  real(8), parameter :: FWHM_TO_SIGMA = 1.0d0/(2.0d0*sqrt(2.0d0*log(2.0d0)))

contains
  subroutine psf_g1d(center, amplitude, fwhm, npx, flux)
    implicit none
    real(8), intent(in)  :: center, amplitude, fwhm
    integer, intent(in)  :: npx
    real(8), intent(out) :: flux(npx)

    real(8) :: aerf(npx+1)
    integer :: i

    aerf = erf(([(i, i=0,npx)]-center-0.5d0)/(fwhm*FWHM_TO_SIGMA))
    flux = amplitude*fwhm*FWHM_TO_SIGMA*H_SQRT_PI*(aerf(2:)-aerf(1:npx))
  end subroutine psf_g1d

  !! One-dimensional Gaussian 
  !! ------------------------
  subroutine gaussian1d(center, amplitude, fwhm, npx, flux)
    implicit none
    real(8), intent(in)  :: center, amplitude, fwhm
    integer, intent(in)  :: npx
    real(8), intent(out) :: flux(npx)

    real(8) :: e1, e2, sigma
    integer :: i, wstart, wwidth

    flux  = 0.0d0
    sigma = fwhm * FWHM_TO_SIGMA
    wwidth = min(npx, 4*ceiling(fwhm))
    wstart = max(0, floor(center) - wwidth/2)
    wwidth = min(npx-wstart, wwidth)

    e1 = erf((wstart-center-0.5d0) / sigma)
    do i = 1, wwidth
       e2 = erf((wstart+i-center-0.5d0) / sigma)
       flux(wstart+i) = amplitude*sigma*H_SQRT_PI*(e2-e1)
       e1 = e2
    end do
  end subroutine gaussian1d

  !! Multithreaded one-dimensional Gaussian 
  !! --------------------------------------
  subroutine gaussian1dmt(center, amplitude, fwhm, npx, nthr, flux)
    implicit none
    real(8), intent(in)  :: center, amplitude, fwhm
    integer, intent(in)  :: npx, nthr
    real(8), intent(out) :: flux(npx)

    real(8), allocatable :: aerf(:)
    real(8) :: sigma
    integer :: i, wstart, wwidth

    !$ call omp_set_num_threads(nthr)

    flux  = 0.0d0
    sigma = fwhm * FWHM_TO_SIGMA
    wwidth = min(npx, 4*ceiling(fwhm))
    wstart = max(0, floor(center) - wwidth/2)
    wwidth = min(npx-wstart, wwidth)

    allocate(aerf(wwidth+1))
    !$omp parallel default(none) shared(flux,aerf,sigma,center,amplitude,wstart,wwidth) private(i)
    !$omp do
    do i = 0, wwidth
       aerf(i+1) = erf((wstart+i-center-0.5d0) / sigma)
    end do
    !$omp end do
    !$omp do
    do i = 1, wwidth
       flux(wstart+i) = amplitude*sigma*H_SQRT_PI*(aerf(i+1)-aerf(i))
    end do
    !$omp end do
    !$omp end parallel 
    deallocate(aerf)
  end subroutine gaussian1dmt

  !! Multiple one-dimensional Gaussian profiles
  !! ------------------------------------------
  subroutine gaussians1d(centers, amplitudes, fwhm, npx, nlines, flux)
    implicit none
    integer, intent(in)  :: npx, nlines
    real(8), intent(in)  :: centers(nlines), amplitudes(nlines), fwhm
    real(8), intent(out) :: flux(npx)

    real(8) :: e1, e2, sigma
    integer :: i, iline, wstart, wwidth

    flux  = 0.0d0
    sigma = fwhm * FWHM_TO_SIGMA
    wwidth = min(npx, 4*ceiling(fwhm))

    do iline = 1, nlines
       wstart = max(0, floor(centers(iline)) - wwidth/2)
    
       e1 = erf((wstart-centers(iline)-0.5d0) / sigma)
       do i = 1, min(npx-wstart, wwidth)
          e2 = erf((wstart+i-centers(iline)-0.5d0) / sigma)
          flux(wstart+i) = flux(wstart+i) + amplitudes(iline)*sigma*H_SQRT_PI*(e2-e1)
          e1 = e2
       end do
    end do
  end subroutine gaussians1d


  real(8) function logl_g1d(center, amplitude, fwhm, error, sky, npx, fobs)
    implicit none
    integer, intent(in)  :: npx
    real(8), intent(in)  :: center,amplitude,fwhm,error,sky,fobs(npx)
    real(8) :: fmod(npx)

    call gaussian1d(center, amplitude, fwhm, npx, fmod)
    logl_g1d = -LLC*npx - 0.5d0*npx*log(error**2) -0.5d0*sum((fobs - (sky+fmod))**2/error**2)
  end function logl_g1d

  real(8) function lnlike_gaussian1d(center, amplitude, fwhm, error, sky, npx, fobs)
    implicit none
    integer, intent(in)  :: npx
    real(8), intent(in)  :: center,amplitude,fwhm,error,sky,fobs(npx)
    real(8) :: fmod(npx)

    call gaussian1d(center, amplitude, fwhm, npx, fmod)
    lnlike_gaussian1d = -LLC*npx - 0.5d0*npx*log(error**2) -0.5d0*sum((fobs - (sky+fmod))**2/error**2)
  end function lnlike_gaussian1d

  real(8) function lnlike_gaussians1d(centers, amplitudes, fwhm, error, sky, npx, nlines, fobs)
    implicit none
    integer, intent(in)  :: npx, nlines
    real(8), intent(in)  :: centers(nlines), amplitudes(nlines), fwhm, error, sky, fobs(npx)
    real(8) :: fmod(npx)

    call gaussians1d(centers, amplitudes, fwhm, npx, nlines, fmod)
    lnlike_gaussians1d = -LLC*npx - 0.5d0*npx*log(error**2) -0.5d0*sum((fobs - (sky+fmod))**2/error**2)
  end function lnlike_gaussians1d

end module gaussian
