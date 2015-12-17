# zaspe
A Code to Measure Stellar Atmospheric Parameters and their Covariance from Spectra

Author: Rafael Brahm (rbrahm@astro.puc.cl)

# About the code
ZASPE computes the atmospheric stellar parameters (Teff, log(g), [Fe/H] and vsin(i)) from echelle spectra via least squares minimization with a pre-computed library of synthetic spectra. The minimization is performed only in the most sensitive spectral zones to changes in the atmospheric parameters. The uncertainities and covariances computed by ZASPE assume that the principal source of error is the systematic missmatch between the observed spectrum and the sythetic one that produces the best fit. A detailed description of this code is available in Brahm et al. (2015). If you use this code for your research, please cite that publication.

In order to run ZASPE you will require a grid of synthetic spectra. ZASPE supports 3 public available libraries, namely Coelho et al. (2005), Husser et al. (2013) and Brahm et al. (2015). For a better performance we recommend to use the latter one, which is available in the following link. After minor modifications to the code, ZASPE is able to use any pre-computed library.

ZASPE accepts two formats of input spectra. The first one corresponds to the standard output of the CERES pipelines (Brahm et al. 2016 in prep) and the second one to a text file with 3 columns, where the first column refers to the echelle order, the second column to the wavelength and the third one to the flux. Each echelle order must have the same number of data points.




