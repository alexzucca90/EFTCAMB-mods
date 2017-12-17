# EFTCAMB-mods
A few mods for EFTCAMB (https://github.com/EFTCAMB) that I need for my projects.

## Generalized-Brans-Dicke project
This mod allows to implement in the Pure_EFT case $\Omega(a)$ and $w_{\rm DE}(a)$ from a tableof the form

a f(a) f'(a), \int_a^1 dx (1+f(x))/x

The integral is needed for wDE(a).

### Installation

I assume that your EFTCAMB directory is EFTCAMB/ 

Put the files:
- GBD_utils.f90                              -> in EFTCAMB/eftcamb/
- 04p9_interpolated_parametrizations_1D.f90  -> in EFTCAMB/eftcamb/04f_parametrizations_1D/
- 07p1_Pure_EFT_std.f90                      -> in EFTCAMB/eftcamb/
- params_EFT.ini                             -> in EFTCAMB/

Also remember to provide two files for Omega(a)  and wDE(a) as said above.

Then run 
  - make clean
  - make eftcamb_dep
  - make eftcamb
  
 ### Run the program
 
 Simply with ./camb params.ini as always.
  

