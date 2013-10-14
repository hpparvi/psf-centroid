module cntr
  implicit none
  
  real(8), parameter :: PI = 3.14159265359d0
  real(8), parameter :: H_SQRT_PI = 0.5d0*sqrt(PI)
  real(8), parameter :: LLC = 0.5d0*log(2.0d0*PI)
  real(8), parameter :: FWHM_TO_SIGMA = 1.0d0/(2.0d0*sqrt(2.0d0*log(2.0d0)))

contains
  subroutine psf_g1d(center,amplitude,fwhm,npx,flux)
    implicit none
    real(8), intent(in)  :: center, amplitude, fwhm
    integer, intent(in)  :: npx
    real(8), intent(out) :: flux(npx)

    real(8) :: aerf(npx+1)
    integer :: i

    aerf = erf(([(i, i=0,npx)]-center)/(fwhm*FWHM_TO_SIGMA))
    flux = amplitude*fwhm*FWHM_TO_SIGMA*H_SQRT_PI*(aerf(2:)-aerf(1:npx))
  end subroutine psf_g1d


  real(8) function logl_g1d(center,amplitude,fwhm,error,sky,npx,fobs)
    implicit none
    integer, intent(in)  :: npx
    real(8), intent(in)  :: center,amplitude,fwhm,error,sky,fobs(npx)

    real(8) :: fmod(npx)

    call psf_g1d(center,amplitude,fwhm,npx,fmod)
    logl_g1d = -LLC*npx - 0.5d0*npx*log(error**2) -0.5d0*sum((fobs - (sky+fmod))**2/error**2)
  end function logl_g1d
end module cntr
