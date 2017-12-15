# EFTCAMB-mods
A few mods for EFTCAMB that I need for my projects.

## Generalized-Brans-Dicke project
This mod allows to implement in the Pure_EFT case $\Omega(a)$ and $w_{\rm DE}(a)$ from a tableof the form

a f(a) f'(a)

### Installation

Put the files:
- GBD_utils.f90                              -> in /eftcamb/
- 04p9_interpolated_parametrizations_1D.f90  -> in /eftcamb/04f_parametrizations_1D/
- 07p1_Pure_EFT_std.f90                      -> in /eftcamb/

Run 
  - make clean
  - make eftcamb_dep
  - make eftcamb
  
  

