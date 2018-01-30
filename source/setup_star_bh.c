#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "atomic.h"
#include "python.h"


/***********************************************************
             University of Southampton

Synopsis:
  get_stellar_params sets rstar, mstar, tstar as well
  as secondary parameters based on user inputs

Arguments:

Returns:


Description:

Notes:

History:
	1502  JM 	Moved here from main()

**************************************************************/

double
get_stellar_params ()
{

    rdpar_comment ("Parameters for the Central Object");

  /* Describe the basic binary star system */

  geo.mstar /= MSOL;		// Convert to MSOL for ease of data entry
  rddoub ("Central_object.mass(msol)", &geo.mstar);
  geo.mstar *= MSOL;

  /* If a BH we want geo.rstar to be at least as large as the last stable orbit for
   * a non-rotating BH
   */

  if (geo.system_type == SYSTEM_TYPE_AGN)
    {
      geo.rstar = 6. * G * geo.mstar / (C * C);	//correction - ISCO is 6x Rg NSH 121025
    }

  rddoub ("Central_object.radius(cm)", &geo.rstar);


  geo.r_agn = geo.rstar;	/* At present just set geo.r_agn to geo.rstar */
  geo.rstar_sq = geo.rstar * geo.rstar;
  if (geo.system_type != SYSTEM_TYPE_AGN)
    {
      rdint ("Central_object.radiation(y=1)", &geo.star_radiation);
      get_spectype (geo.star_radiation,
		    //"Rad_type_for_star(0=bb,1=models)_to_make_wind",
		    "Central_object.rad_type_to_make_wind(0=bb,1=models)",
		    &geo.star_ion_spectype);

      if (geo.star_radiation)
	rddoub ("Central_object.temp", &geo.tstar_init);
    }
  else
    {
      geo.star_radiation = 0;
      geo.tstar_init = 0;
    }

  /* tstar_init and lum_star_init refer to values without the effects of backscattering */

  geo.tstar = geo.tstar_init;

  geo.lum_star = geo.lum_star_init =
    4 * PI * geo.rstar * geo.rstar * STEFAN_BOLTZMANN * pow (geo.tstar, 4.);


  /* Describe the secondary if that is required */

  if (geo.system_type == SYSTEM_TYPE_BINARY)	/* It's a binary system */
    {

      geo.m_sec /= MSOL;	// Convert units for ease of data entry
      rddoub ("msec(msol)", &geo.m_sec);
      geo.m_sec *= MSOL;

      geo.period /= 3600.;	// Convert units to hours for easy of data entry
      rddoub ("period(hr)", &geo.period);
      geo.period *= 3600.;	// Put back to cgs immediately
    }

  return (geo.lum_star_init);
}




/***********************************************************
             University of Southampton

Synopsis:
  get_bl_and_agn_params sets up the boundary layer and agn power law parameters
  based on user input and system type

Arguments:
  lstar     double
            star luminosity as calculated by get_stellar_params
Returns:

Description:

Notes:

History:
	1502  JM 	Moved here from main()

**************************************************************/

int
get_bl_and_agn_params (lstar)
     double lstar;
{
  double xbl;
  double temp_const_agn;

  rdpar_comment ("Parameters for BL or AGN");

  if (geo.system_type == SYSTEM_TYPE_AGN)	/* If it is an AGN */
    {
      geo.star_radiation = 0;	// 70b - AGN do not have a star at the center */
      geo.bl_radiation = 0;
      geo.agn_radiation = 1;
      rdint ("QSO_BH_radiation(y=1)", &geo.agn_radiation);
    }
  else
    {
      rdint ("Boundary_layer.radiation(y=1)", &geo.bl_radiation);
      geo.agn_radiation = 0;	// So far at least, our star systems don't have a BH
    }


  get_spectype (geo.bl_radiation,
		//"Rad_type_for_bl(0=bb,1=models,3=pow)_to_make_wind",
		"Boundary_layer.rad_type_to_make_wind(0=bb,1=models,3=pow)",
		&geo.bl_ion_spectype);
  get_spectype (geo.agn_radiation,
		"Rad_type_for_agn(0=bb,1=models,3=power_law,4=cloudy_table,5=bremsstrahlung)_to_make_wind",
		&geo.agn_ion_spectype);

  /* 130621 - ksl - This is a kluge to add a power law to stellar systems.  What id done
     is to remove the bl emission, which we always assume to some kind of temperature
     driven source, and replace it with a power law source

     Note that the next 3 or 4 lines just tell you that there is supposed to be a power
     law source.  They don't teel you what the parameters are.
   */

  if (geo.bl_ion_spectype == SPECTYPE_POW)
    {
      geo.agn_radiation = 1;
      geo.agn_ion_spectype = SPECTYPE_POW;
      geo.bl_radiation = 0;
      Log ("Trying to make a start with a power law boundary layer\n");
    }
  else
    {
      Log ("Not Trying to make a start with a power law boundary layer %d\n",
	   geo.bl_ion_spectype);
    }


  /* Describe the boundary layer */

  if (geo.bl_radiation && geo.bl_ion_spectype != SPECTYPE_POW)
    {
      xbl = geo.lum_bl = 0.5 * G * geo.mstar * geo.disk_mdot / geo.rstar;

      rddoub ("Boundary_layer.luminosity(ergs/s)", &geo.lum_bl);
      Log ("OK, the bl lum will be about %.2e the disk lum\n",
	   geo.lum_bl / xbl);
      rddoub ("Boundary_layer.temp(K)", &geo.t_bl);
    }
  else
    {
      geo.lum_bl = 0;
      geo.t_bl = 0;
    }

  /* Describe the agn */

  if (geo.agn_radiation && geo.system_type == SYSTEM_TYPE_AGN)	/* This peculiar line is to enamble us to add a star with a power law component */
    {
      xbl = geo.lum_agn = 0.5 * G * geo.mstar * geo.disk_mdot / geo.r_agn;

      /* If there is no disk, initilize geo.lum to the luminosity of a star */
      if (geo.disk_type == DISK_NONE)
	{
	  geo.lum_agn = lstar;
	}

      // At present we have set geo.r_agn = geo.rstar, and encouraged the user
      // set the default for the radius of the BH to be 6 R_Schwartschild.
      // rddoub("R_agn(cm)",&geo.r_agn);

      /* if we have a "blackbody agn" the luminosity is set by Stefan Boltzmann law
         once the AGN blackbody temp is read in, otherwise set by user */
      if (geo.agn_ion_spectype != SPECTYPE_BB)
	rddoub ("lum_agn(ergs/s)", &geo.lum_agn);


      Log ("OK, the agn lum will be about %.2e the disk lum\n",
	   geo.lum_agn / xbl);
      if (geo.agn_ion_spectype == SPECTYPE_POW
	  || geo.agn_ion_spectype == SPECTYPE_CL_TAB)
	{
	  geo.alpha_agn = (-1.5);
	  rddoub ("agn_power_law_index", &geo.alpha_agn);

	  if (geo.alpha_agn == -1.0)	//deal with the pathological case
	    {
	      geo.const_agn = geo.lum_agn / (log (2.42e18) - log (4.84e17));
	    }
	  else
	    {
	      geo.const_agn =
		geo.lum_agn /
		(((pow (2.42e18, geo.alpha_agn + 1.)) -
		  pow (4.84e17,
		       geo.alpha_agn + 1.0)) / (geo.alpha_agn + 1.0));
	    }
	  Log ("AGN Input parameters give a power law constant of %e\n",
	       geo.const_agn);
	}
      else if (geo.agn_ion_spectype == SPECTYPE_BREM)
	{

	  geo.brem_temp = 1.16e8;	//10kev
	  geo.brem_alpha = -0.2;	//This is the cloudy form of bremstrahlung
	  geo.const_agn = 1.0;
	  rddoub ("agn_bremsstrahlung_temp(K)", &geo.brem_temp);
	  rddoub ("agn_bremsstrahlung_alpha", &geo.brem_alpha);
	  temp_const_agn =
	    geo.lum_agn / qromb (integ_brem, 4.84e17, 2.42e18, 1e-4);
	  geo.const_agn = temp_const_agn;
	  Log ("AGN Input parameters give a Bremsstrahlung constant of %e\n",
	       temp_const_agn);

	}
      else if (geo.agn_ion_spectype == SPECTYPE_BB)
	{
	  /* note that alpha_agn holds the temperature in the case of "blackbody agn" */
	  rddoub ("agn_blackbody_temp(K)", &geo.alpha_agn);
	  geo.lum_agn =
	    4 * PI * geo.r_agn * geo.r_agn * STEFAN_BOLTZMANN *
	    pow (geo.alpha_agn, 4.);
	}

      /* JM 1502 -- lines to add a low frequency power law cutoff. accessible
         only in advanced mode and for non broken power law.
         default is zero which is checked before we call photo_gen_agn */
      geo.pl_low_cutoff = 0.0;
      if (modes.iadvanced && (geo.agn_ion_spectype == SPECTYPE_POW))
	rddoub ("@agn_power_law_cutoff", &geo.pl_low_cutoff);

      rdint ("geometry_for_pl_source(0=sphere,1=lamp_post)",
	     &geo.pl_geometry);

      if (geo.pl_geometry == PL_GEOMETRY_LAMP_POST)
	{
	  rddoub ("lamp_post.height(r_g)", &geo.lamp_post_height);
	  geo.lamp_post_height *= G * geo.mstar / C / C;	//get it in CGS units
	  Log ("lamp_post_height is cm is %g\n", geo.lamp_post_height);
	}
      else if (geo.pl_geometry != PL_GEOMETRY_SPHERE)	// only two options at the moment
	{
	  Error ("Did not understand power law geometry %i. Fatal.\n",
		 geo.pl_geometry);
	  exit (0);
	}



      /* Computes the constant for the power law spectrum from the input alpha and 2-10 luminosity.
         This is only used in the sim correction factor for the first time through.
         Afterwards, the photons are used to compute the sim parameters. */



      if (geo.agn_ion_spectype == SPECTYPE_CL_TAB)	/*NSH 0412 - option added to allow direct comparison with cloudy power law table option */
	{
	  geo.agn_cltab_low = 1.0;
	  geo.agn_cltab_hi = 10000;
	  rddoub ("low_energy_break(ev)", &geo.agn_cltab_low);	/*lo frequency break - in ev */
	  rddoub ("high_energy_break(ev)", &geo.agn_cltab_hi);
	  geo.agn_cltab_low_alpha = 2.5;	//this is the default value in cloudy
	  geo.agn_cltab_hi_alpha = -2.0;	//this is the default value in cloudy
	}
    }
  else if (geo.agn_radiation)	/* We want to add a power law to something other than an AGN */
    {
      xbl = geo.lum_agn = 0.5 * G * geo.mstar * geo.disk_mdot / geo.r_agn;

      // At present we have set geo.r_agn = geo.rstar, and encouraged the user
      // set the default for the radius of the BH to be 6 R_Schwartschild.
      // rddoub("R_agn(cm)",&geo.r_agn);

      rddoub ("lum_agn(ergs/s)", &geo.lum_agn);
      Log ("OK, the agn lum will be about %.2e the disk lum\n",
	   geo.lum_agn / xbl);
      geo.alpha_agn = (-1.5);
      rddoub ("agn_power_law_index", &geo.alpha_agn);

      /* JM 1502 -- lines to add a low frequency power law cutoff. accessible
         only in advanced mode. default is zero which is checked before we call photo_gen_agn */
      geo.pl_low_cutoff = 0.0;
      if (modes.iadvanced)
	rddoub ("@agn_power_law_cutoff", &geo.pl_low_cutoff);


      /* Computes the constant for the power law spectrum from the input alpha and 2-10 luminosity.
         This is only used in the sim correction factor for the first time through.
         Afterwards, the photons are used to compute the sim parameters. */


      if (geo.alpha_agn == -1.0)	//deal with the pathological case
	{
	  geo.const_agn = geo.lum_agn / (log (2.42e18) - log (4.84e17));
	}
      else
	{
	  geo.const_agn =
	    geo.lum_agn /
	    (((pow (2.42e18, geo.alpha_agn + 1.)) -
	      pow (4.84e17, geo.alpha_agn + 1.0)) / (geo.alpha_agn + 1.0));
	}


      Log ("AGN Input parameters give a power law constant of %e\n",
	   geo.const_agn);

      if (geo.agn_ion_spectype == SPECTYPE_CL_TAB)	/*NSH 0412 - option added to allow direct comparison with cloudy power law table option */
	{
	  geo.agn_cltab_low = 1.0;
	  geo.agn_cltab_hi = 10000;
	  rddoub ("low_energy_break(ev)", &geo.agn_cltab_low);	/*lo frequency break - in ev */
	  rddoub ("high_energy_break(ev)", &geo.agn_cltab_hi);
	  geo.agn_cltab_low_alpha = 2.5;	//this is the default value in cloudy
	  geo.agn_cltab_hi_alpha = -2.0;	//this is the default value in cloudy
	}
    }

  else
    {
      geo.r_agn = 0.0;
      geo.lum_agn = 0.0;
      geo.alpha_agn = 0.0;
      geo.const_agn = 0.0;
    }
  return (0);
}




/***********************************************************
             University of Southampton

Synopsis:
  get_standard_care_factors provides more control over how the program is
  run

Arguments:

Returns:

Description:

Notes:
    ksl - It is not obvious that much recent thought has been
    given to the choices that are here.  The fractional distance
    that a photon travel is intended to make sure the velocity
    along the line of sight can be approximated linearly.  If 
    a photon travels too far in an azimuthal direction the sense
    of the velocity can change and this prevensts this

    The lowest ion density contributing to photoionization is used
    to determine what ions one has to calculate the photoionzation
    xsection for.  The lower this density; the more that have to be
    calculated, and as a result the slower the program.

    Keeping photoionizaion during final spectrum allows one to
    check the contribution of photoabsorption.

History:
  1502  JM  Moved here from main()

**************************************************************/
int
get_standard_care_factors ()
{
  int istandard;
  istandard = 1;
  SMAX_FRAC = 0.5;
  DENSITY_PHOT_MIN = 1.e-10;

  /* 141116 - ksl - Made care factors and advanced command as this is clearly somethng that is diagnostic */

  if (modes.iadvanced)
    {
      rdint ("@Use.standard.care.factors(1=yes)", &istandard);

      if (!istandard)
	{
	  rddoub ("@Fractional.distance.photon.may.travel", &SMAX_FRAC);
	  rddoub ("@Lowest.ion.density.contributing.to.photoabsorption",
		  &DENSITY_PHOT_MIN);
	  rdint ("@Keep.photoabs.during.final.spectrum(1=yes)",
		 &modes.keep_photoabs);
	}
    }
  return (0);
}