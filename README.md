# EFTCAMB-mods

This repo contains the modified files for [EFTCAMB](https://github.com/EFTCAMB) needed for my projects.

## Dark Energy Denisty Reconstruction project
We needed to study the implications of the non-trivial observed Dark Energy density (soon a paper coming on arXiv) for the generalized Brans-Dicke (GBD) models. 


## Dark Energy Fractional Density modules
The reconstruction of a GBD theory is implemented inside EFTCAMB in the files
- `09p3_Designer_GBD.f90` using the DE eos wDE(a). 
- `09p4_Designer_GBD_mod.f90` using the DE fractional density X(a).
 The first reconstructs the GBD theory provided the equation of state of the effective dark energy (wDE). For theoretical reasons, wDE might not be enough general so I created `09p4_Designer_GBD_mod.f90` to use the fractional density of Dark Energy.

There is also a module that implements the fractional Dark Energy density instead of the equation of state for the pure EFT models,
 - `07p2_Pure_EFT_mode.f90`
 
 The DE density can assume several forms, we are mostly focused on a non-Bayesian reconstruction of the DE density. Hence I introduced an interpolation function from a table a, X(a) in
 - `04p15_interpolated_function_1D.f90`.
 
 ## Sampling the X(a)
 To study the implications of the reconstructed Dark Energy (fractional) density X(a) it is useful to sample through the allowed X(a) distribution. 
 The reconstruction of the dark energy density is provided in two files (see `params_EFT.ini`)
 - `xDE_means_filename` which contains the  means of the reconstructed Dark Energy density
 - `xDE_covmat_filename` which contains the  covariance matrix of the reconstructed Dark Energy density
 
 We assume that the DE reconstruction is performed on a number of bins specified by `xDE_n_bins` to be specified in `params_EFT.ini`.
 
 Drawing the sample of the DE density his is handled by the module contained in
 - `12_EFT_sampler.f90` 
 
 and is driven by the sampling drivers
 - `xDE_sampler.f90` - the standard CAMB version
 - `xDE_sampler_sources.f90` - for CAMB sources
 
 
### Installing (need to be updated)

Download [EFTCAMB sources](https://github.com/EFTCAMB/EFTCAMB_sources)

I assume that your EFTCAMB sources directory is `EFTCAMB/`

Put the following files in the folder `EFTCAMB/` (replace if asked to):
- `params_EFT.ini` 
- `equations_EFT_new.f90`   (in the `Makefile_main` file set EQUATIONS ?= equations_EFT_new )


Put the folloing files in the folder `EFTCAMB/eftcamb`
- `02_GBD_utils.f90`                                            -> in `EFTCAMB/eftcamb/`
- `02p10_gl10.f90`                                              -> in `EFTCAMB/eftcamb/02_bundled_algorithms`
- `07p1_Pure_EFT_std.f90`                                       -> in `EFTCAMB/eftcamb/`
- `09p1_Designer_fR.f90`                                        -> in `EFTCAMB/eftcamb/`
- `09p3_Designer_GBD.f90`                                       -> in `EFTCAMB/eftcamb/`
- `09p4_Designer_GBD_mod.f90`                                   -> in `EFTCAMB/eftcamb/`
- `11_EFTCAMB_main.f90`                                         -> in `EFTCAMB/eftcamb/`
                                             


Also remember, when usign the PureEFT approach, to provide two files for Omega(a)  and wDE(a) as said above. Please, make sure that the number of points in the files matches the numbet_of_points variable declared at the beginning of the file `04p9_interpolated_parametrizations_1D.f90`.

Then run 
  - `make clean`
  - `make eftcamb_dep`
  - `make eftcamb`
  
 ### Run the program
 Explore the updated `params_EFT.ini` and choose the model you want to run.
 Then run EFTCAMB with `./camb params.ini` as always.
  

