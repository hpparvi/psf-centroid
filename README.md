psf-centroid
============

Calculates an accurate one-dimensional undersampled Gaussian profile by
integrating the profile over each pixel analytically.

That is, instead of evaluating

\[
f(a,c,fwhm) = \exp \left ( - \left (\frac{x-c}{s}\right )^2 \right )
\]
