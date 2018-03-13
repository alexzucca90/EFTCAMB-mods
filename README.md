# EFTCAMB-mods
A few mods for EFTCAMB (https://github.com/EFTCAMB) that I need for my projects.

## Generalized-Brans-Dicke project
This mod allows to implement in the Pure_EFT case $\Omega(a)$ and $w_{\rm DE}(a)$ from a table of the form

a f(a) f'(a), \int_a^1 dx (1+f(x))/x

The integral is needed for wDE(a).

NEW: The reconstruction of a GBD theory is implemented inside EFTCAMB in the files
- 09p3_Designer_GBD.f90
- 09p4_Designer_GBD_mod.f90
 The first reconstructs the GBD theory provided the equation of state of the effective dark energy (wDE). For theoretical reasons, wDE might not be enough general so I created 09p4_Designer_GBD_mod.f90 to use the fractional density of Dark Energy.


### Installing

I assume that your EFTCAMB directory is EFTCAMB/ 

Put the files:
- 02_GBD_utils.f90                                            -> in EFTCAMB/eftcamb/
- 02p10_gl10.f90                                              -> in EFTCAMB/eftcamb/02_bundled_algorithms
- 04p9_interpolated_parametrizations_1D.f90                   -> in EFTCAMB/eftcamb/04f_parametrizations_1D/
- 04p10_reconstructed_fit_parametrizations_1D.f90             -> in EFTCAMB/eftcamb/04f_parametrizations_1D/
- 04p11_reconstructed_fit_DE_parametrizations_1D.f90          -> in EFTCAMB/eftcamb/04f_parametrizations_1D/
- 04p12_reconstructed_fit_DE_tracking_parametrizations_1D.f90 -> in EFTCAMB/eftcamb/04f_parametrizations_1D/
- 04p13_power_law_DE_parametrizations.f90                     -> in EFTCAMB/eftcamb/04f_parametrizations_1D/
- 07p1_Pure_EFT_std.f90                                       -> in EFTCAMB/eftcamb/
- 09p1_Designer_fR.f90                                        -> in EFTCAMB/eftcamb/
- 09p3_Designer_GBD.f90                                       -> in EFTCAMB/eftcamb/
- 09p4_Designer_GBD_mod.f90                                   -> in EFTCAMB/eftcamb/
- 11_EFTCAMB_main.f90                                         -> in EFTCAMB/eftcamb/
- params_EFT.ini                                              -> in EFTCAMB/


Also remember, when usign the PureEFT approach, to provide two files for Omega(a)  and wDE(a) as said above. Please, make sure that the number of points in the files matches the numbet_of_points variable declared at the beginning of the file 04p9_interpolated_parametrizations_1D.f90.

Then run 
  - make clean
  - make eftcamb_dep
  - make eftcamb
  
 ### Run the program
 Explore the updated params_EFT.ini and choose the model you want to run.
 Then run EFTCAMB with ./camb params.ini as always.
  

