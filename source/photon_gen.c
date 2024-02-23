
/***********************************************************/
/** @file  photon_gen.c
 * @author ksl
 * @date   May, 2018
 *
 * @brief  Primary routines for creating photons for use in
 * the radiative transfer calculation.
 *
 ***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "atomic.h"
#include "python.h"

// These are external variables that are used to determine whether one needs to reinitialize
// by running xdefine_phot
double f1_old = 0;
double f2_old = 0;
int iwind_old = 0;

#define PRINT_OFF 0
#define PRINT_ON  1


/**********************************************************/
/**
 * @brief
 * the controlling routine for creating the underlying photon distribution.
 *
 * @param [out] PhotPtr  p   The structure where all photons are stored
 * @param [in] double  f1   The mininum frequency
 * @param [in] double  f2   The maximum frequency if a uniform distribution
 * @param [in] long  nphot_tot   The total number of photons that need to be generated to reach the total
 * luminosity, not necessarilly the number of photons which will be generated by the call to define_phot,
 * which instead is defined by NPHOT
 * @param [in] int  ioniz_or_extract   CYCLE_IONIZ -> this is for the wind ionization calculation,
 * CYCLE_EXTRACT-> it is for the final spectrum calculation
 * @param [in] int  iwind   A variable (see below) that controls the generation of phtons for the wind
 * @param [in] int  freq_sampling   If 0, generate photons over the full frequency range without banding,
 * else use banding
 * @return     Always returns 0
 *
 * @details
 * This routine oversees the generation of photons for the radiative transfer calculation.  Most of the actual
 * work is done in other routines, particularly xdefine_phot).  There can be multiple
 * radiation sources (e.g the star and a disk)
 *
 * There are two basic approaches to generating photons.  One can sample the entire wavelength interval in which
 * case all photons have the same weight, but photon frequencies are selected based on the the specific luminiosty.
 * In the banded approach, one defines bands in which photons are generated.  Specific numbers of photons are
 * are generated in each band, and within a band the photons are selected according the the specific luminosity,
 * but the weights of photons in different bands differ.
 * still used for detailed spectrum calculation. Which of this choices to use is controlled by freq_sampling
 * (The weights are established here)
 *
 * iwind is a variable that determines how or whether to create photons from the wind:
 * * -1-> Do not consider wind photons under any circumstances
 * * 0  ->Consider wind photons.  There is no need to recalculate the
 *   relative contribution to number of photons from the wind
 *   unless the frequencies have changed
 * * 1 -> Force a recalculation of the fraction of photons from the wind.
 *   This is necessary because after each ionization recalculation
 *   one needs to recalculate the flux fraction from the wind.
 *
 * ### Notes ###
 * @bug Is this correct, have subcycles been removed.
 *
 **********************************************************/

int
define_phot (p, f1, f2, nphot_tot, ioniz_or_extract, iwind, freq_sampling)
     PhotPtr p;
     double f1, f2;
     long nphot_tot;
     int ioniz_or_extract;
     int iwind;
     int freq_sampling;         // 0 --> old uniform approach, 1 --> minimum fractions ins various bins

{
  double natural_weight, weight;
  double ftot;
  int n;
  int iphot_start, nphot_rad, nphot_k;
  long nphot_tot_rad, nphot_tot_k;
  nphot_k = nphot_tot_k = natural_weight = iphot_start = 0;     // Initialize to avoid compiler warnings

  /* if we are generating nonradiative kpackets, then we need to subtract 
     off the fraction reserved for k-packets */
  if (geo.nonthermal && (geo.rt_mode == RT_MODE_MACRO) && (ioniz_or_extract == CYCLE_IONIZ))
  {
    nphot_k = (geo.frac_extra_kpkts * NPHOT);
    nphot_rad = NPHOT - nphot_k;
    nphot_tot_k = (geo.frac_extra_kpkts * nphot_tot);
    nphot_tot_rad = nphot_tot - nphot_tot_k;
  }
  else
  {
    nphot_rad = NPHOT;
    nphot_tot_rad = nphot_tot;
  }

  if (freq_sampling == 0)
  {
    /* Original approach, uniform sampling of entire wavelength interval,
       still used for detailed spectrum calculation */

    if (f1 != f1_old || f2 != f2_old || iwind != iwind_old)
    {                           // The reinitialization is required
      xdefine_phot (f1, f2, ioniz_or_extract, iwind, PRINT_ON, 1);
    }
    /* The weight of each photon is designed so that all of the photons add up to the
       luminosity of the photosphere.  This implies that photons must be generated in such
       a way that it mimics the energy distribution of the star. */

    geo.weight = (weight) = (geo.f_tot) / (nphot_tot_rad);

    for (n = 0; n < NPHOT; n++)
      p[n].path = -1.0;         /* SWM - Zero photon paths */

    xmake_phot (p, f1, f2, ioniz_or_extract, iwind, weight, 0, nphot_rad);
  }
  else
  {
    /* Use banding, create photons with different weights in different wavelength
       bands.  This is used for the for ionization calculation where one wants to assure
       that you have "enough" photons at high energy */

    ftot = populate_bands (ioniz_or_extract, iwind, &xband);

    for (n = 0; n < NPHOT; n++)
      p[n].path = -1.0;         /* SWM - Zero photon paths */

    /* Now generate the photons */

    iphot_start = 0;

    for (n = 0; n < xband.nbands; n++)
    {

      if (xband.nphot[n] > 0)
      {
        /* Reinitialization is required here always because we are changing
         * the frequencies around all the time */
        Log ("Defining photons for band %d...\n", n);
        if (n == 0)
          xdefine_phot (xband.f1[n], xband.f2[n], ioniz_or_extract, iwind, PRINT_ON, 1);
        else
          xdefine_phot (xband.f1[n], xband.f2[n], ioniz_or_extract, iwind, PRINT_ON, 0);
        /* The weight of each photon is designed so that all of the photons add up to the
           luminosity of the photosphere.  This implies that photons must be generated in such
           a way that it mimics the energy distribution of the star. */
        geo.weight = (natural_weight) = (ftot) / (nphot_tot_rad);
        xband.weight[n] = weight = natural_weight * xband.nat_fraction[n] / xband.used_fraction[n];
        xmake_phot (p, xband.f1[n], xband.f2[n], ioniz_or_extract, iwind, weight, iphot_start, xband.nphot[n]);

        iphot_start += xband.nphot[n];
      }

    }

  }

  /* deal with k-packets generated from nonradiative heating */
  if (geo.nonthermal && (geo.rt_mode == RT_MODE_MACRO) && (ioniz_or_extract == 0))
  {
    /* calculate the non-radiative kpkt luminosity throughout the wind */
    geo.f_kpkt = get_kpkt_heating_f ();

    /* get the number of photons we have reserved in the photon structure */
    //nphot_k = geo.frac_extra_kpkts * NPHOT; 
    weight = (geo.f_kpkt) / (nphot_tot_k);

    /* throw an error if the k-packet weight is too high or low */
    if (weight > (100.0 * natural_weight) || weight < (0.01 * natural_weight))
    {
      Error ("define_phot: kpkt weight is %8.4e compared to characteristic photon weight %8.4e\n", weight, natural_weight);
    }
    if (sane_check (weight))
      Error ("define_phot: kpkt weight is %8.4e!\n", weight);

    Log ("!! xdefine_phot: total & banded kpkt luminosity due to non-radiative heating: %8.2e %8.2e \n", geo.heat_shock, geo.f_kpkt);


    /* generate the actual photons produced by the k-packets */
    photo_gen_kpkt (p, weight, iphot_start, nphot_k);
  }


  for (n = 0; n < NPHOT; n++)
  {
    p[n].w_orig = p[n].w;
    p[n].freq_orig = p[n].freq;
    p[n].origin_orig = p[n].origin;
    p[n].np = n;
    p[n].ds = 0;
    p[n].line_res = NRES_NOT_SET;
    p[n].frame = F_OBSERVER;
    if (geo.reverb != REV_NONE && p[n].path < 0.0)      // SWM - Set path lengths for disk, star etc.
      simple_paths_gen_phot (&p[n]);
  }

  return (0);

}





/**********************************************************/
/**
 * @brief      Determine how many photons to allocate to each band
 *
 * @param [in] int  ioniz_or_extract   indicate whether this is for an ioniz_or_extract
 * @param [in] int  iwind   indicates where wind photons are created
 * @param [in,out] struct xbands *  band   is a pointer to the band sturctue where everythin is
 * @return     ftot which is the sum of the band limited luminosities
 *
 * @details
 * This routine determines how many photons to allocate to each
 * band.  Inputs determine a minimum fraction to allocate to each
 * band, but these need not sum to 1, and so the rest of the photons
 * are allocated naturally (that is to say) by the luminosities within
 * the bands.
 *
 * Much of the actual work is carried out in xdefine_phot.
 *
 * ### Notes ###
 *
 **********************************************************/

double
populate_bands (ioniz_or_extract, iwind, band)
     int ioniz_or_extract;
     int iwind;
     struct xbands *band;

{
  double ftot, frac_used, z;
  int n, nphot, most, nphot_rad;

  /* Get all of the band limited luminosities */
  ftot = 0.0;

  /* this is the number of photons minus the number reserved for k-packets */
  if (geo.nonthermal && (geo.rt_mode == RT_MODE_MACRO))
    nphot_rad = NPHOT - (geo.frac_extra_kpkts * NPHOT);
  else
    nphot_rad = NPHOT;

  for (n = 0; n < band->nbands; n++)    // Now get the band limited luminosities
  {
    if (band->f1[n] < band->f2[n])
    {
      if (n == 0)
        xdefine_phot (band->f1[n], band->f2[n], ioniz_or_extract, iwind, PRINT_OFF, 1); //We only need to compute the total wind lum for the first time thruogh
      else
        xdefine_phot (band->f1[n], band->f2[n], ioniz_or_extract, iwind, PRINT_OFF, 0);

      ftot += band->flux[n] = geo.f_tot;
    }
    else
      band->flux[n] = 0.0;

    if (band->flux[n] == 0.0)
      band->min_fraction[n] = 0;
  }


/* So now we can distribute the photons among the different bands */

  frac_used = 0;

  for (n = 0; n < band->nbands; n++)
  {
    band->nat_fraction[n] = band->flux[n] / ftot;
    frac_used += band->min_fraction[n];
  }

  nphot = 0;
  z = 0;
  most = 0;


  for (n = 0; n < band->nbands; n++)
  {
    band->used_fraction[n] = band->min_fraction[n] + (1 - frac_used) * band->nat_fraction[n];
    nphot += band->nphot[n] = nphot_rad * band->used_fraction[n];
    if (band->used_fraction[n] > z)
    {
      z = band->used_fraction[n];
      most = n;
    }
  }

  /* Because of roundoff errors nphot may not sum to the desired value, namely NPHOT less kpackets.  
     So add a few more photons to the band with most photons already. It should only be a few, at most
     one photon for each band. */

  if (nphot < nphot_rad)
  {
    band->nphot[most] += (nphot_rad - nphot);
  }

  return (ftot);

}



/**********************************************************/
/**
 * @brief      calculates the band limited fluxes & luminosities, and prepares
 * for the producting of photons (by reinitializing the disk, among other things).
 *
 * @param [in] double  f1   The minimum freqeuncy
 * @param [in] double  f2   The maximum frequency
 * @param [in] int  ioniz_or_extract   0 -> this is for the wind ionization calculation;
 * 1-> it is for the final spectrum calculation
 * @param [in] int  iwind   if 0, include wind photons; if 1 include wind photons and force a recalcuation of
 * ion denisities, if -1, ignore the possibility of wind photons
 * @param [in] int print_mode Detemines whether certain summary messages are printed out or not
 * @param [in] int tot_flag, we only need to compute the total wind luminosity for the first band
 * @return     Always returns 0
 *
 * @details
 * This is a routine that initializes variables that are used to determine how many
 * photons from a particular source in a particular wavelength range are to be generated
 * It does not generate photons itself.
 *
 * The routine calls various other routines to calcuate the band limited luminosities
 * of various radiation sources (which is used in allocating how many photons to create
 * from each).  
 * 
 * The results are stored in elements of the geo structure. (By convention
 * elements of the geo structure that begin with f, such a geo.f_star refer
 * to the band limited luminosity of a particular source, in theis case the
 * star (or central object).
 *
 * ### Notes ###
 * 
 *
 * @bug This routine is something of a kluge.  Information is passed back to the calling
 * routine through the geo structure, rather than a more obvious method.  The routine was
 * created when a banded approach was defined, but the whole section might be more obvious.
 *
 **********************************************************/

int
xdefine_phot (f1, f2, ioniz_or_extract, iwind, print_mode, tot_flag)
     double f1, f2;
     int ioniz_or_extract;
     int iwind;
     int print_mode;
     int tot_flag;
{
  /* First determine if you need to reinitialize because the frequency boundaries are
     different than previously */

  geo.lum_star = geo.lum_disk = 0;      // bl luminosity is an input so it is not zeroed
  geo.f_star = geo.f_disk = geo.f_bl = 0;
  geo.f_kpkt = geo.f_matom = 0; //SS - kpkt and macro atom luminosities set to zero

  if (geo.star_radiation)
  {
    star_init (f1, f2, ioniz_or_extract, &geo.f_star);
  }
  if (geo.disk_radiation)
  {

/* Note  -- disk_init not only calculates fluxes and luminosity for the disk.  It
calculates the boundaries of the various disk annulae depending on f1 and f2 */

    disk_init (geo.disk_rad_min, geo.disk_rad_max, geo.mstar, geo.disk_mdot, f1, f2, ioniz_or_extract, &geo.f_disk);
  }
  if (geo.bl_radiation)
  {
    bl_init (geo.lum_bl, geo.t_bl, f1, f2, ioniz_or_extract, &geo.f_bl);
  }
  if (geo.agn_radiation)
  {
    agn_init (geo.rstar, geo.lum_agn, geo.alpha_agn, f1, f2, ioniz_or_extract, &geo.f_agn);
  }

/* The choices associated with iwind are
iwind = -1 	Don't generate any wind photons at all
         1  Create wind photons and force a reinitialization of the wind
         0  Create wind photons but remain open to the question of whether
		the wind needs to be reinitialized.  Initialization is forced
		in that case by init
*/

  if (iwind == -1)
    geo.f_wind = geo.lum_wind = 0.0;

  if (iwind == 1 || (iwind == 0))
  {
    /* Find the luminosity of the wind in the CMF, and the energy emerging in the simulation time step */
    if (tot_flag == 1)
    {
      geo.lum_wind = wind_luminosity (0.0, VERY_BIG, MODE_CMF_TIME);
    }

//    xxxpdfwind = 1;             // Turn on the portion of the line luminosity routine which creates pdfs
    geo.f_wind = wind_luminosity (f1, f2, MODE_OBSERVER_FRAME_TIME);
//    xxxpdfwind = 0;             // Turn off the portion of the line luminosity routine which creates pdfs

  }

  /* Handle the initialization of emission via k-packets and macro atoms. SS */

  if (geo.matom_radiation)
  {
    /* Only calculate macro atom emissivity during ionization cycles ant
       at the beginning of the spectral cycles.  Otherwise we can 
       can use the saved emissivities.  The routine  returns the specific luminosity
       in the spectral band of interest */

    if (!modes.jumps_for_detailed_spectra)
    {
      Log ("Using accelerated calculation of emissivities\n");
      if (geo.pcycle == 0)
      {
        geo.f_matom = get_matom_f_accelerate (CALCULATE_MATOM_EMISSIVITIES);
      }
      else
        geo.f_matom = get_matom_f_accelerate (USE_STORED_MATOM_EMISSIVITIES);
    }
    else
    {
      Log ("Using old slow calculation of emissivities\n");
      if (geo.pcycle == 0)
      {
        geo.f_matom = get_matom_f (CALCULATE_MATOM_EMISSIVITIES);
      }
      else
        geo.f_matom = get_matom_f (USE_STORED_MATOM_EMISSIVITIES);
    }



    geo.f_kpkt = get_kpkt_f (); /* This returns the specific luminosity
                                   in the spectral band of interest */

    matom_emiss_report ();
  }

  geo.f_tot = geo.f_star + geo.f_disk + geo.f_bl + geo.f_wind + geo.f_kpkt + geo.f_matom + geo.f_agn;
  geo.lum_tot = geo.lum_star + geo.lum_disk + geo.lum_bl + geo.lum_agn + geo.lum_wind;
  /* Store the 3 variables that have to remain the same to avoid reinitialization */

  geo.f1 = f1_old = f1;
  geo.f2 = f2_old = f2;
  iwind_old = iwind;

  if (print_mode == PRINT_ON)
  {
    phot_status ();
  }


  return (0);

}



/**********************************************************/
/**
 * @brief    Log information about total and band limited
 * luminosities
 *
 * @return     Always returns 0
 *
 * @details
 * This routines logs information from that exists in the
 * geo structure at the time it is called.
 *
 * ### Notes ###
 *
 **********************************************************/

int
phot_status ()
{

  Log
    ("!! xdefine_phot: lum_tot %8.2e lum_star %8.2e lum_bl %8.2e lum_bh %8.2e lum_disk %8.2e lum_wind %8.2e\n",
     geo.lum_tot, geo.lum_star, geo.lum_bl, geo.lum_agn, geo.lum_disk, geo.lum_wind);

  Log
    ("!! xdefine_phot:   f_tot %8.2e   f_star %8.2e   f_bl %8.2e   f_bh %8.2e   f_disk %8.2e   f_wind %8.2e   f_matom %8.2e   f_kpkt %8.2e \n",
     geo.f_tot, geo.f_star, geo.f_bl, geo.f_agn, geo.f_disk, geo.f_wind, geo.f_matom, geo.f_kpkt);

  Log
    ("!! xdefine_phot: wind ff %8.2e       fb %8.2e   lines  %8.2e  for freq %8.2e %8.2e\n",
     geo.lum_ff, geo.lum_rr, geo.lum_lines, geo.f1, geo.f2);
  if (geo.lum_star > 0)
  {
    Log
      ("!! xdefine_phot: star  tstar  %8.2e   %8.2e   lum_star %8.2e %8.2e  %8.2e \n",
       geo.tstar, geo.tstar_init, geo.lum_star, geo.lum_star_init, geo.lum_star_back);
  }
  if (geo.lum_disk > 0)
  {
    Log
      ("!! xdefine_phot: disk                               lum_disk %8.2e %8.2e  %8.2e \n",
       geo.lum_disk, geo.lum_disk_init, geo.lum_disk_back);
  }
  if (geo.adiabatic)
    Log ("!! xdefine_phot: heating & cooling  due to adiabatic processes:         %8.2e %8.2e \n", geo.heat_adiabatic, geo.cool_adiabatic);

  return (0);
}


/**********************************************************/
/**
 * @brief      just makes photons (in a particular wavelength range for
 * all radiation sources)
 *
 * @param [out] PhotPtr  p   The entire photon structure
 * @param [in] double  f1   The minimum frequency for generating photons
 * @param [in] double  f2   The maximum frequency for generating photons
 * @param [in] int  ioniz_or_extract   A flag indicating whether this if for an ionization or detailed
 * spectral cycle (Used to determine what underlying spectrum, e.g bb or detailed models) to sample
 * @param [in] int  iwind   A flag indicating whether or not to generate any wind photons.
 * @param [in] double  weight   The weight of photons to generate
 * @param [in] int  iphot_start   The position in the photon structure to start storing photons
 * @param [in] int  nphotons   The number of photons to generate
 * @return     Always returns 0
 *
 * @details
 * make_phot controls the actual generation of photons.  All of the initializations should
 * have been done previously (xdefine_phot).  xmake_phot cycles through the various possible
 * sources of the wind, including for example, the disk, the central object, and the
 * wind, and creates photons for each, using the ratio of the band limited luminosites to the
 * total band limited luminosity to determine how many photons to select from each source.
 *
 * ### Notes ###
 *
 **********************************************************/

int
xmake_phot (p, f1, f2, ioniz_or_extract, iwind, weight, iphot_start, nphotons)
     PhotPtr p;
     double f1, f2;
     int ioniz_or_extract;
     int iwind;
     double weight;
     int iphot_start;           //The place to begin putting photons in the photon structure in this call
     int nphotons;              //The total number of photons to generate in this call
{

  int nphot, nn;
  int nstar, nbl, nwind, ndisk, nmatom, nagn, nkpkt;
  double agn_f1;

  /* For the diagnostic searchlight mode
     we intercept the normal procedure for generating
     photons and substitute searchlight mode.  
   */

  if (ioniz_or_extract == CYCLE_EXTRACT && modes.searchlight)
  {

    weight = 1;

    if (modes.searchlight == 1)
    {
      if (geo.system_type == SYSTEM_TYPE_STAR || geo.system_type == SYSTEM_TYPE_CV)
      {
        photo_gen_star (p, geo.rstar, geo.tstar, weight, f1, f2, geo.star_spectype, iphot_start, nphotons);
      }
      else if (geo.system_type == SYSTEM_TYPE_BH || geo.system_type == SYSTEM_TYPE_AGN)
      {
        photo_gen_agn (p, geo.rstar, geo.alpha_agn, weight, f1, f2, geo.agn_spectype, iphot_start, nphotons);
      }
      else
      {
        Error ("xmake_phot: Uknown system type (%d) for central object searchlight option\n", geo.system_type);
        Exit (1);
      }
    }
    else if (modes.searchlight == 2)
    {
      photo_gen_disk (p, weight, f1, f2, geo.disk_spectype, iphot_start, nphotons);
    }
    else
    {
      Error ("xmake_phot: Unknown searchlight mode %d\n", modes.searchlight);
      Exit (1);
    }
  }


/* End of generation of photons via the seaach light mode */

  nstar = nbl = nwind = ndisk = 0;
  nagn = nkpkt = nmatom = 0;

  if (geo.star_radiation)
  {
    nstar = geo.f_star / geo.f_tot * nphotons;
  }
  if (geo.bl_radiation)
  {
    nbl = geo.f_bl / geo.f_tot * nphotons;
  }
  if (iwind >= 0)
  {
    nwind = geo.f_wind / geo.f_tot * nphotons;
  }
  if (geo.disk_radiation)
  {
    ndisk = geo.f_disk / geo.f_tot * nphotons;  /* Ensure that nphot photons are created */
  }
  if (geo.agn_radiation)
  {
    nagn = geo.f_agn / geo.f_tot * nphotons;    /* Ensure that nphot photons are created */
  }
  if (geo.matom_radiation || geo.nonthermal)
  {
    nkpkt = geo.f_kpkt / geo.f_tot * nphotons;

    if (geo.matom_radiation)
      nmatom = geo.f_matom / geo.f_tot * nphotons;
  }


  nphot = ndisk + nwind + nbl + nstar + nagn + nkpkt + nmatom;

  if (nphot < nphotons)
  {
    if (ndisk > 0)
      ndisk += (nphotons - nphot);
    else if (nwind > 0)
      nwind += (nphotons - nphot);
    else if (nbl > 0)
      nbl += (nphotons - nphot);
    else if (nagn > 0)
      nagn += (nphotons - nphot);
    else
      nstar += (nphotons - nphot);
  }


  Log
    ("photon_gen: band %6.2e to %6.2e weight %6.2e nphotons %8d ndisk %7d nwind %7d nstar %7d npow %d \n",
     f1, f2, weight, nphotons, ndisk, nwind, nstar, nagn);

  /* Generate photons from the star, the bl, the wind and then from the disk */
  /* Now adding generation from kpkts and macro atoms too (SS June 04) */


  if (geo.star_radiation)
  {
    nphot = nstar;
    if (nphot > 0)
    {
      if (ioniz_or_extract == CYCLE_EXTRACT)
        photo_gen_star (p, geo.rstar, geo.tstar, weight, f1, f2, geo.star_spectype, iphot_start, nphot);
      else
        photo_gen_star (p, geo.rstar, geo.tstar, weight, f1, f2, geo.star_ion_spectype, iphot_start, nphot);
    }
    iphot_start += nphot;
  }
  if (geo.bl_radiation)
  {
    nphot = nbl;

    if (nphot > 0)
    {
      if (ioniz_or_extract == CYCLE_EXTRACT)
        photo_gen_star (p, geo.rstar, geo.t_bl, weight, f1, f2, geo.bl_spectype, iphot_start, nphot);
      else
        photo_gen_star (p, geo.rstar, geo.t_bl, weight, f1, f2, geo.bl_ion_spectype, iphot_start, nphot);
/* Reassign the photon type since we are actually using the same routine as for generating
stellar photons */
      nn = 0;
      while (nn < nphot)
      {
        p[iphot_start + nn].origin = PTYPE_BL;
        nn++;
      }
    }
    iphot_start += nphot;
  }

/* Generate the wind photons */

  if (iwind >= 0)
  {
    nphot = nwind;
    if (nphot > 0)
      photo_gen_wind (p, weight, f1, f2, iphot_start, nphot);
    iphot_start += nphot;
  }


/* Generate the disk photons */

  if (geo.disk_radiation)
  {
    nphot = ndisk;
    if (nphot > 0)
    {
      if (ioniz_or_extract == CYCLE_EXTRACT)
        photo_gen_disk (p, weight, f1, f2, geo.disk_spectype, iphot_start, nphot);
      else
        photo_gen_disk (p, weight, f1, f2, geo.disk_ion_spectype, iphot_start, nphot);
    }
    iphot_start += nphot;
  }

  /* Generate the agn photons */

  if (geo.agn_radiation)
  {
    nphot = nagn;
    if (nphot > 0)
    {
      /* JM 1502 -- lines to add a low frequency power law cutoff. accessible
         only in advanced mode */
      if (geo.pl_low_cutoff != 0.0 && geo.pl_low_cutoff > f1)
        agn_f1 = geo.pl_low_cutoff;

      /* error condition if user specifies power law cutoff below that hardwired in
         ionization cycles */
      else if (geo.pl_low_cutoff > f1 && ioniz_or_extract == CYCLE_IONIZ)
      {
        Error ("photo_gen_agn: power_law low f cutoff (%8.4e) is lower than hardwired minimum frequency (%8.4e)\n", geo.pl_low_cutoff, f1);
        agn_f1 = f1;
      }
      else
        agn_f1 = f1;


      if (ioniz_or_extract == CYCLE_EXTRACT)
        photo_gen_agn (p, geo.rstar, geo.alpha_agn, weight, agn_f1, f2, geo.agn_spectype, iphot_start, nphot);
      else
        photo_gen_agn (p, geo.rstar, geo.alpha_agn, weight, agn_f1, f2, geo.agn_ion_spectype, iphot_start, nphot);
    }
    iphot_start += nphot;
  }

  /* Now do macro atoms and k-packets. SS June 04 */

  if (geo.matom_radiation)
  {
    nphot = nkpkt;
    if (nphot > 0)
    {
      if (ioniz_or_extract == CYCLE_IONIZ)
      {
        Error ("xmake_phot: generating photons by k-packets when performing ionization cycle without shock heating. Abort.\n");
        Exit (1);               //The code shouldn't be doing this - something has gone wrong somewhere. (SS June 04)
      }
      else
      {
        photo_gen_kpkt (p, weight, iphot_start, nphot);
      }
    }
    iphot_start += nphot;

    nphot = nmatom;
    if (nphot > 0)
    {
      if (ioniz_or_extract == CYCLE_IONIZ)
      {
        Error ("xmake_phot: generating photons by macro atoms when performing ionization cycle. Abort.\n");
        Exit (0);               //The code shouldn't be doing this - something has gone wrong somewhere. (SS June 04)
      }
      else
      {
        photo_gen_matom (p, weight, iphot_start, nphot);
      }
    }
    iphot_start += nphot;
  }

  return (0);
}






/**********************************************************/
/**
 * @brief      (r, tstar, freqmin, freqmax, ioniz_or_extract, f)
 *
 * @param [in] double  freqmin   The minimum freqency for the band
 * @param [in] double  freqmax   The maximum freqency for the band
 * @param [in out] int  ioniz_or_extract  A flag indicating whether this is for
 * an ionization cycle or a detailed spectrum cycle
 * @param [out] double *  f   The band limited luminosity of the star
 * @return     Always returns 0
 *
 * @details
 * This routine calculates the luminosity of the star and the luminosity within
 * the frequency boundaries.
 *
 * ### Notes ###
 * The routine allows for backscattering of light onto the star
 * assuming the user has opted to include this in the determination
 * of the effective temperature of the star
 *
 **********************************************************/

int
star_init (freqmin, freqmax, ioniz_or_extract, f)
     double freqmin, freqmax, *f;
     int ioniz_or_extract;
{
  double r, tstar, log_g;
  double emit, emittance_bb (), emittance_continuum ();
  int spectype;

  log_g = geo.gstar = log10 (GRAV * geo.mstar / (geo.rstar * geo.rstar));
  r = geo.rstar;

  tstar = geo.tstar = geo.tstar_init;
  geo.lum_star = geo.lum_star_init;

  if (geo.absorb_reflect == BACK_RAD_ABSORB_AND_HEAT && geo.lum_star_back > 0)
  {
    geo.lum_star = geo.lum_star + geo.lum_star_back;
    tstar = geo.tstar = pow (geo.lum_star / (4 * PI * STEFAN_BOLTZMANN * r * r), 0.25);
  }


  if (ioniz_or_extract == CYCLE_EXTRACT)
    spectype = geo.star_spectype;       /* type for final spectrum */
  else
    spectype = geo.star_ion_spectype;   /*type for ionization calculation */

  if (spectype >= 0)
  {
    emit = emittance_continuum (spectype, freqmin, freqmax, tstar, log_g);
  }
  else
  {
    emit = emittance_bb (freqmin, freqmax, tstar);
  }

  *f = emit;                    // Calculate the surface flux between freqmin and freqmax
  *f *= (4. * PI * r * r);


  return (0);

}



/**********************************************************/
/**
 * @brief
 * Generate nphot photons from the star in the frequency interval f1 to f2
 *
 * @param [out] PhotPtr  p   The entire photon structure
 * @param [in] double  r   The radiius at which to generate photons
 * @param [in] double  t   The effective temperature of the stat
 * @param [in] double  weight   The weight of photons to generate
 * @param [in] double  f1   The minimum frequencey
 * @param [in] double  f2   The maximum frequencey
 * @param [in] int  spectype   The type of spectrum (bb or spectral model) to use in generating photons
 * @param [in] int  istart   The positions in the photon structure where the first photons will be stored
 * @param [in] int  nphot   The nubmer of photons to generate
 * @return     Always returns 0
 *
 * @details
 * The routine generates photons emergying from a sphere.  It uses
 * randvcos to determine the direction of the photon.
 *
 * ### Notes ###
 * This routine is also used in generating photons from a boundary layer
 *
 * The routine allows for the possibility that part of the star will be
 * covered by a vertically extended disk and avoids generating photons there.
 *
 **********************************************************/

int
photo_gen_star (p, r, t, weight, f1, f2, spectype, istart, nphot)
     PhotPtr p;
     double r, t, weight;
     double f1, f2;             /* The freqency mininimum and maximum if a uniform distribution is selected */
     int spectype;              /*The spectrum type to generate: 0 is bb, 1 (or in fact anything but 0)
                                   is uniform in frequency space */
     int istart, nphot;         /* Respecitively the starting point in p and the number of photons to generate */
{
  double freqmin, freqmax;
  int i, iend;
  if ((iend = istart + nphot) > NPHOT)
  {
    Error ("photo_gen_star: iend %d > NPHOT %d\n", iend, NPHOT);
    Exit (0);
  }
  if (f2 < f1)
  {
    Error ("photo_gen_star: Cannot generate photons if freqmax %g < freqmin %g\n", f2, f1);
  }
  Log_silent ("photo_gen_star creates nphot %5d photons from %5d to %5d \n", nphot, istart, iend);
  freqmin = f1;
  freqmax = f2;
  r = (1. + EPSILON) * r;       /* Generate photons just outside the photosphere */
  for (i = istart; i < iend; i++)
  {
    p[i].origin = PTYPE_STAR;   // For BL photons this is corrected in photon_gen
    p[i].frame = F_OBSERVER;    // Stellar photons are not redshifted
    p[i].w = weight;
    p[i].istat = p[i].nscat = p[i].nrscat = p[i].nmacro = 0;
    p[i].grid = 0;
    p[i].tau = 0.0;
    p[i].nres = p[i].line_res = -1;     // It's a continuum photon
    p[i].nnscat = 1;

    if (spectype == SPECTYPE_BB)
    {
      p[i].freq = planck (t, freqmin, freqmax);
    }
    else if (spectype == SPECTYPE_UNIFORM)
    {                           /* Kurucz spectrum */
      /*Produce a uniform distribution of frequencies */
      p[i].freq = random_number (freqmin, freqmax);     //Generate a random frequency - this will exclude freqmin,freqmax.

    }
    else if (spectype == SPECTYPE_MONO)
    {
      p[i].w = 1. / geo.pcycles;
      p[i].freq = geo.mono_freq;
    }

    else
    {
      p[i].freq = one_continuum (spectype, t, geo.gstar, freqmin, freqmax);
    }

    if (p[i].freq < freqmin || freqmax < p[i].freq)
    {
      Error_silent ("photo_gen_star: phot no. %d freq %g out of range %g %g\n", i, p[i].freq, freqmin, freqmax);
    }

    if (modes.searchlight && geo.ioniz_or_extract == CYCLE_EXTRACT)
    {
      stuff_v (geo.searchlight_x, p[i].x);
      stuff_v (geo.searchlight_lmn, p[i].lmn);
    }
    else
    {
      randvec (p[i].x, r);

      if (geo.disk_type == DISK_VERTICALLY_EXTENDED)
      {
        while (fabs (p[i].x[2]) < zdisk (r))
        {
          randvec (p[i].x, r);
        }

      }

      randvcos (p[i].lmn, p[i].x);

    }

    /* This is set up for looking at photons in spectral cycles at present */
    //if (modes.save_photons && geo.ioniz_or_extract == CYCLE_EXTRACT)
    //  save_photons (&p[i], "STAR");

  }
  return (0);
}

/* THESE ROUTINES ARE FOR THE BOUNDARY LAYER */




/**********************************************************/
/**
 * @brief      calculate parameters need to intialize the boundary layer (namely
 * the band-liminted luminoisity)
 *
 * @param [in out] double  lum_bl   The desired luminosity for the boundary layaer
 * @param [in out] double  t_bl   The temperature
 * @param [in out] double  freqmin   The minimum freqency in the band
 * @param [in out] double  freqmax   The maximum freqency in the band
 * @param [in out] int  ioniz_or_extract   The spectral type to use in the calucation
 * (NOT USED) but see below
 * @param [in out] double *  f   The band limited luminosity
 * @return   The total luminosity
 *
 * The only thing that is actually calculated here is f, the luminosity
 * 	within the frequency range that is specified.
 *
 * @details
 * This routine calculates the  luminosity of the bl within the frequency boundaries.
 * BB functions are assumed.  It was derived from the same routine for a star,
 * but here we have assumed that the temperature and the luminosity are known
 *
 * ### Notes ###
 * @bug  At present bl_init assumes a BB regardless of the spectrum.
 * This is not really correct, and is different for what is done in initializing
 * the star
 *
 * 0703 - ksl - This is rather an odd little routine.  As noted all that is
 * calculated is f.  ioniz_or_extract is not used, and lum_bl which
 * is returned is only the luminosity that was passed.
 *
 **********************************************************/

double
bl_init (lum_bl, t_bl, freqmin, freqmax, ioniz_or_extract, f)
     double lum_bl, t_bl, freqmin, freqmax, *f;
     int ioniz_or_extract;
{
//OLD  double q1;
  double integ_planck_d ();
//OLD  double alphamin, alphamax;

//OLD  q1 = 2. * PI * (BOLTZMANN * BOLTZMANN * BOLTZMANN * BOLTZMANN) / (PLANCK * PLANCK * PLANCK * VLIGHT * VLIGHT);
//OLD  alphamin = PLANCK * freqmin / (BOLTZMANN * t_bl);
//OLD  alphamax = PLANCK * freqmax / (BOLTZMANN * t_bl);
//  *f = q1 * integ_planck_d (alphamin, alphamax) * lum_bl / STEFAN_BOLTZMANN;

  *f = emittance_bb (freqmin, freqmax, t_bl) * lum_bl / (t_bl * t_bl * t_bl * t_bl * STEFAN_BOLTZMANN);



  return (lum_bl);
}




/**********************************************************/
/**
 * @brief
 * Perform some simple checks on the photon distribution just produced.
 *
 * @param [in] PhotPtr  p  The photon structure
 * @param [in] double  freqmin   The minimum fequency that was used to generate the photons
 * @param [in] double  freqmax   The maximum fequency that was used to generate the photons
 * @param [in] char *  comment   A comment that accompanies this particular call
 * @return     Always returns 0
 *
 * @details
 * The routine checks photons to see if they are "reasonable".  The checks mostly have to
 * do with the frequency of the photons.   
 *
 * The program will exit is too many photons fail the
 * checks. The number of photons that are permitted to "fail" depends on the
 * total number of photons in a flight of photons
 *
 * ### Notes ###
 *
 * In checking for the reasonableness of the frequencies, their is a
 * fairly generous allowance for doppler shifts in the wind.
 *
 * The frequency limits are not enforced on photons that have excited
 * macro-atoms.
 *
 * The routine also determines the calculates some features of the
 * photon distribution, specifically having to do with the ionizing
 * photons.  It is not entirely clear why this is where this is done
 *
 * 181009 - ksl - Previously, this routine caused Python to exit 
 * if photon_checks produced more than a small number of errors. I
 * have removed this extreme measure but that does not mean that
 * photon checks should be ignored.
 *
 **********************************************************/

int
photon_checks (p, freqmin, freqmax, comment)
     char *comment;
     PhotPtr p;
     double freqmin, freqmax;
{
  int nnn, nn;
  int nlabel;
  geo.n_ioniz = 0;
  geo.cool_tot_ioniz = 0.0;
  nnn = 0;
  nlabel = 0;


  /* Next two lines are to allow for fact that photons generated in
   * a frequency range may be Doppler shifted out of that range, especially
   * if they are disk photons generated right up against one of the frequency
   * limits
   * 04aug--ksl-increased limit from 0.02 to 0.03, e.g from 6000 km/s to 9000 km/s
   * 11apr--NSH-decreased freqmin to 0.4, to take account of double redshifted photons.
   * shift.
   */

  Debug ("photon_checks: %s\n", comment);

  freqmax *= (1.8);
  freqmin *= (0.6);

  for (nn = 0; nn < NPHOT; nn++)
  {
    if (PLANCK * p[nn].freq > ion[0].ip)
    {
      geo.cool_tot_ioniz += p[nn].w;
      geo.n_ioniz += p[nn].w / (PLANCK * p[nn].freq);
    }

    if (sane_check (p[nn].freq) != 0 || sane_check (p[nn].w))
    {
      if (nlabel == 0)
      {
        Error ("photon_checks:   nphot  origin  freq     freqmin    freqmax\n");
        nlabel++;
      }
      Error ("photon_checks: %5d %5d %5d %10.4e %10.4e %10.4e freq or weight are not sane\n", nn, p[nn].origin, p[nn].nres, p[nn].freq,
             freqmin, freqmax);
      p[nn].freq = freqmax;
      nnn++;
    }
    if (p[nn].origin < 10 && (p[nn].freq < freqmin || freqmax < p[nn].freq))
    {
      if (nlabel == 0)
      {
        Error ("photon_checks:   nphot  origin  freq     freqmin    freqmax\n");
        nlabel++;
      }
      Error ("photon_checks: %5d %5d %5d %10.4e %10.4e %10.4e freq out of range\n", nn, p[nn].origin, p[nn].nres, p[nn].freq, freqmin,
             freqmax);
      p[nn].freq = freqmax;
      nnn++;
    }

    if (length (p[nn].lmn) < 0.9999999 || length (p[nn].lmn) > 1.0000001)
    {
      Error ("photon_checks: %5d %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e  origin %d length lmn out of range %10.6f\n",
             nn, p[nn].x[0], p[nn].x[1], p[nn].x[2], p[nn].lmn[0], p[nn].lmn[1], p[nn].lmn[2], p[nn].origin, length (p[nn].lmn));
    }

    if (check_frame (&p[nn], F_OBSERVER, "photon_checks: all photons shouuld be in OBSERVER frame\n"))
    {
      nnn++;
    }
  }

  if (nnn == 0)
    Debug ("photon_checks: All photons passed checks successfully\n");
  else
  {
    Log ("photon_checks: %d of %d or %e per cent of photons failed checks\n", nnn, NPHOT, nnn * 100. / NPHOT);
  }

  return (0);
}
