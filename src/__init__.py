"""Fast PSF modelling and centroiding in python

   Calculates an accurate one-dimensional undersampled Gaussian profile by
   integrating the profile over each pixel analytically. 

   Contains
   --------
     gaussian1d   - 1D Gaussian profile
     gaussian1dmt - 1D Gaussian profile (multithreaded)
     gaussians1d  - Multiple 1D Gaussian profiles

   Notes
   -----
     - The routines use the center of the first pixel as the origo. Use 
       an appropriate offset when calculating over a subarray. 

"""
import numpy as np
from scipy.optimize import fmin
from numpy.random import multivariate_normal
from cntrf import cntr

try:
    from emcee import EnsembleSampler
    with_emcee = True
except ImportError:
    with_emcee = False

__all__ = ['psf_g1d', 'logl_g1d', 'centroid', 'gaussian1d', 'gaussian1dmt', 'gaussians1d']

gaussian1d   = cntr.gaussian1d
gaussian1dmt = cntr.gaussian1dmt
gaussians1d  = cntr.gaussians1d

logl_gaussian1d = cntr.logl_g1d

psf_g1d  = cntr.psf_g1d
logl_g1d = cntr.logl_g1d

def centroid(img, x0, y0, fwhm_x=8., fwhm_y=8., verbose=False, **kwargs):
    def prior_bounds(pv):
        return -1e18 if not ((0 < pv[0] < img.shape[1]) | (0 < pv[1] < img.shape[0])) else 0

    estimate_errors = kwargs.get('estimate_errors', True)
    return_chains = kwargs.get('return_chains', False)
    operation   = kwargs.get('operation', 'mean')
    maxiter     = kwargs.get('maxiter',     5000)
    maxfun      = kwargs.get('maxfun',      5000)
    mc_threads  = kwargs.get('mc_threads',     1)
    mc_nwalkers = kwargs.get('mc_nwalkers',   50)
    mc_niter    = kwargs.get('mc_niter',     300)
    mc_thinning = kwargs.get('mc_thinning',  300)
    mc_burn     = kwargs.get('mc_burn',      100)

    if operation == 'mean':
        x, y = img.mean(axis=0), img.mean(axis=1)
    elif operation == 'max':
        x, y = img.max(axis=0), img.max(axis=1)
    else:
        raise TypeError

    vmin, vmax = 0.5*(x.min()+y.min()), 0.5*(x.max()+y.max())
    
    pv0 = np.array([x0, y0, fwhm_x, fwhm_y, vmax-vmin, 1e-2*(vmax-vmin), vmin])

    lpfun = lambda pv: (  logl_g1d(pv[0], pv[4], pv[2], pv[5], pv[6], x)
                        + logl_g1d(pv[1], pv[4], pv[3], pv[5], pv[6], y)
                        + prior_bounds(pv))

    pv = fmin(lambda pv:-lpfun(pv), pv0, disp=verbose, maxfun=maxfun, maxiter=maxiter)

    if not (with_emcee and estimate_errors):
        return pv, -np.ones(pv.size)
    else:
        sampler = EnsembleSampler(mc_nwalkers, pv.size, lpfun, threads=1)
        sampler.run_mcmc(multivariate_normal(pv, 5e-3*np.eye(pv.size), size=mc_nwalkers), mc_niter);

        fc = sampler.chain[:,mc_burn::mc_thinning,:].reshape([-1,pv.size])
        pc = np.array(np.percentile(fc, [50,16,84], axis=0))
        
        if return_chains:
            pc[0,:], np.mean(np.abs(pc[1:,:]-pc[0,:]), axis=0), fc
        else:
            return pc[0,:], np.mean(np.abs(pc[1:,:]-pc[0,:]), axis=0)
