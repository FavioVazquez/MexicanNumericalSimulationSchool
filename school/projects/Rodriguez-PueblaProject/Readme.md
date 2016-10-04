# Aldo Rodriguez-Puebla Project

The following link is for a galaxy mock catalogue based on the technique described today during my lecture. It contains 575,884 galaxies with stellar mass above log (M_* / Msol) > 9 when using Mvir as the main halo property.  The simulation box has 250 Mpc/h of length. Note that every halo has 8 different stellar masses calculated using different halo properties. The file is organized as follows:

https://www.dropbox.com/s/eyhx3opluipvh9h/mock_%40_0.00000.fits?dl=0


  - col 1 (cs?): -1 if this is a distinct halo otherwise a subhalo
  - col 2-9 (logMs_P): Stellar mass calculated using the property P of the halo, described below.
  - col 10 (X): Position X in Mpc
  - col 11 (Y): Position Y in Mpc
  - col12 (Z): Position Z in Mpc
  - col13 (logMvir): Present virial mass for halo, virial mass at infall for subhalos
  - col14 (logMpeak): Maximum mass of the (sub)halo has ever had in its history
  - col15 (logM200b):  Present mass of the halo enclosing 200 times the mean matter density. Same for subhalos but at the mass at infall. 
  - col16 (logM200c):  Present mass of the halo enclosing 200 times the critical density. Same for subhalos but at the mass at infall. 
  - col17 (logM500c):  Present mass of the halo enclosing 500 times the critical density. Same for subhalos but at the mass at infall. 
  - col18 (logM2500c):  Present mass of the halo enclosing 2500 times the critical density. Same for subhalos but at the mass at infall. 
  - col19 (vmax):  Present maximum circular velocity of the halo. Same for subhalos but at  infall. 
  - col20 (vpeak):  The highest maximum circular velocity the (sub)halo has ever had in its history
  - col21(Rvir): Virial radius of the halo in kpc. 

Below are some suggested exercises:

  1. You can start by calculating the galaxy stellar mass function using any stellar mass definition that you want. 

  2. Calculate the fraction of satellites as a function of stellar mass for all the stellar mass definitions. Do they produce similar results? Which gives stellar mass definition gives the highest and the lowest fraction?

  3. Plot the mean logMs_vmax as a function of logMvir for all, centrals and satellites separately. What did you notice? Do it on the other way around, i.e., logMs_vir as a function vmax  for all, centrals and satellites separately. Based on the observed fractions on ii) what would you conclude?
