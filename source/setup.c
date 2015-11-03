#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "atomic.h"
#include "python.h"

/***********************************************************
             University of Southampton

Synopsis: 
  parse_command_line parses the command line and communicates stuff
  to the loggin routines (e.g. verbosity).
   
Arguments:   
  argc            command line arguments 

Returns:
  restart_start   1 if restarting
 
 
Description:  

Notes:

	The switches should be described in the introduction to
	python.c, instead of here

History:
  1502  JM  Moved here from main()

**************************************************************/


int
parse_command_line (argc, argv)
     int argc;
     char *argv[];
{
  int restart_stat, verbosity, time_to_quit, i;
  char dummy[LINELENGTH];
  int mkdir ();
  double time_max;

  restart_stat = 0;

  if (argc == 1)
    {
      printf ("Input file (interactive=stdin):");
      fgets (dummy, LINELENGTH, stdin);
      get_root (files.root, dummy);
      strcpy (files.diag, files.root);
      strcat (files.diag, ".diag");
    }
  else
    {

      for (i = 1; i < argc; i++)
	{
	  if (strcmp (argv[i], "-h") == 0)
	    {
	      help ();
	    }
	  else if (strcmp (argv[i], "-r") == 0)
	    {
	      Log ("Restarting %s\n", files.root);
	      restart_stat = 1;
	    }
	  else if (strcmp (argv[i], "-t") == 0)
	    {
	      if (sscanf (argv[i + 1], "%lf", &time_max) != 1)
		{
		  Error ("python: Expected time after -t switch\n");
		  exit (0);
		}
	      i++;

	    }
	  else if (strcmp (argv[i], "-v") == 0)
	    {
	      if (sscanf (argv[i + 1], "%d", &verbosity) != 1)
		{
		  Error ("python: Expected verbosity after -v switch\n");
		  exit (0);
		}
	      Log_set_verbosity (verbosity);
	      i++;

	    }
	  else if (strcmp (argv[i], "-e") == 0)
	    {
	      if (sscanf (argv[i + 1], "%d", &time_to_quit) != 1)
		{
		  Error ("python: Expected max errors after -e switch\n");
		  exit (0);
		}
	      Log_quit_after_n_errors (time_to_quit);
	      i++;

	    }
	  else if (strcmp (argv[i], "-d") == 0)
	    {
	      modes.iadvanced = 1;
	      i++;
	    }
	  else if (strcmp (argv[i], "-f") == 0)
	    {
	      modes.fixed_temp = 1;
	      i++;
	    }
    else if (strcmp (argv[i], "-z") == 0)
      {
        modes.zeus_connect = 1;
        Log ("setting zeus_connect to %i\n",modes.zeus_connect);
        i++;
      }
	  else if (strcmp (argv[i], "-i") == 0)
	    {
	      modes.quit_after_inputs = 1;
	    }

    else if (strcmp (argv[i], "--version") == 0)
      {
        /* give information about the pyhon version, such as commit hash */
        Log ("Python Version %s \n", VERSION);  //54f -- ksl -- Now read from version.h
        Log ("Built from git commit hash %s\n", GIT_COMMIT_HASH);
        /* warn the user if there are uncommited changes */
        int git_diff_status = GIT_DIFF_STATUS;
        if (git_diff_status > 0)
          Log("This version was compiled with %i files with uncommitted changes.\n",
              git_diff_status);
        exit(0);
      }

	  else if (strncmp (argv[i], "-", 1) == 0)
	    {
	      Error ("python: Unknown switch %s\n", argv[i]);
	      help ();
	    }
	}

      /* The last command line variable is always the .pf file */

      strcpy (dummy, argv[argc - 1]);
      get_root (files.root, dummy);

      /* This completes the parsing of the command line */

      /* JM130722 we now store diag files in a subdirectory if in parallel */
      /* ksl - I believe this is created in all cases, and that is what we want */

      sprintf (files.diagfolder, "diag_%s/", files.root);
      mkdir (files.diagfolder, 0777);
      strcpy (files.diag, files.diagfolder);
      sprintf (dummy, "_%d.diag", rank_global);
      strcat (files.diag, files.root);
      strcat (files.diag, dummy);


    }

  return (restart_stat);
}



/***********************************************************
             University of Southampton

Synopsis: 
  init_files 
   
Arguments:   
  argc            command line arguments 

Returns:
  restart_start   1 if restarting
 
 
Description:  

Notes:

History:
  1502  JM  Moved here from main()

**************************************************************/

int
init_log_and_windsave (restart_stat)
     int restart_stat;
{
  FILE *fopen (), *qptr;

  if (restart_stat == 0)
    {				// Then we are simply running from a new model
      xsignal_rm (files.root);	// Any old signal file
      xsignal (files.root, "%-20s %s \n", "START", files.root);
      Log_init (files.diag);
    }
  else
    {
      /* Note that alghough we chekc that we dan open the windsave file, it is not read here.   */

      strcpy (files.windsave, files.root);
      strcat (files.windsave, ".wind_save");
      qptr = fopen (files.windsave, "r");

      if (qptr != NULL)
	{
	  /* Then the file does exist and we can restart */
	  fclose (qptr);
	  xsignal (files.root, "%-20s %s\n", "RESTART", files.root);
	  Log_append (files.diag);
	}
      else
	{
	  /* It does not exist and so we start from scratch */
	  restart_stat = 0;
	  xsignal_rm (files.root);	// Any old signal file
	  xsignal (files.root, "%-20s %s \n", "START", files.root);
	  Log_init (files.diag);
	}
    }

  return (0);
}

/***********************************************************
             University of Southampton

Synopsis: 
  get_grid_params reads information on the coordinate system
  and grid dimensions and sets the corresponding variables
  in the geo structure
   
Arguments:    

Returns:
 
 
Description:  

Notes:

History:
  1502  JM  Moved here from main()
  1508	ksl	Updated for domains

**************************************************************/

int
get_grid_params (ndom)
     int ndom;
{
  int input_int;

  // XXX - If we keep this error, then we need to assure that geo.ndomain is incremented
  // before this statement

  if (ndom >= geo.ndomain)
    Error ("Trying to get grid params for a non-existent domain!\n");

  input_int=1;

  /* ksl - The if statement seems superflous.  Why are we entering this routine if we are continuing and earlier calculation? */
  if (geo.run_type != SYSTEM_TYPE_PREVIOUS)
    {
      /* Define the coordinate system for the grid and allocate memory for the wind structure */
      rdint
	("Coord.system(0=spherical,1=cylindrical,2=spherical_polar,3=cyl_var)",
	 &input_int);
      switch (input_int)
	{
	case 0:
	  zdom[ndom].coord_type = SPHERICAL;
	  break;
	case 1:
	  zdom[ndom].coord_type = CYLIND;
	  break;
	case 2:
	  zdom[ndom].coord_type = RTHETA;
	  break;
	case 3:
	  zdom[ndom].coord_type = CYLVAR;
	  break;
	default:
	  Error
	    ("Invalid parameter supplied for 'Coord_system'. Valid coordinate types are: \n\
          0 = Spherical, 1 = Cylindrical, 2 = Spherical polar, 3 = Cylindrical (varying Z)");
	}

      rdint ("Wind.dim.in.x_or_r.direction", &zdom[ndom].ndim);
      if (zdom[ndom].coord_type)
	{
	  rdint ("Wind.dim.in.z_or_theta.direction", &zdom[ndom].mdim);
	  if (zdom[ndom].mdim < 4)
	    {
	      Error
		("python: domain mdim must be at least 4 to allow for boundaries\n");
	      exit (0);
	    }
	}
      else
	zdom[ndom].mdim = 1;

    }

/* 130405 ksl - Check that NDIM_MAX is greater than NDIM and MDIM.  */

  if ((zdom[ndom].ndim > NDIM_MAX) || (zdom[ndom].mdim > NDIM_MAX))
    {
      Error
	("NDIM_MAX %d is less than NDIM %d or MDIM %d. Fix in python.h and recompile\n",
	 NDIM_MAX, zdom[ndom].ndim, zdom[ndom].mdim);
      exit (0);
    }


  /* If we are in advanced then allow the user to modify scale lengths */
  if (modes.iadvanced)
    {
      rdint ("adjust_grid(0=no,1=yes)", &modes.adjust_grid);

      if (modes.adjust_grid)
	{
	  Log ("You have opted to adjust the grid scale lengths\n");
	  rddoub ("geo.xlog_scale", &zdom[ndom].xlog_scale);
	  if (&zdom[ndom].coord_type)
	    rddoub ("geo.zlog_scale", &zdom[ndom].zlog_scale);
	}
    }

  zdom[ndom].ndim2 = zdom[ndom].ndim * zdom[ndom].mdim;


  return (0);
}



/***********************************************************
             University of Southampton

Synopsis: 
  get_line_transfer_mode reads in the variable geo.line_mode
  and sets the variables geo.line_mode, geo.scatter_mode,
  geo.rt_mode and geo.macro_simple accordingly
   
Arguments:		

Returns:
 
 
Description:	

Notes:

History:
	1502  JM 	Moved here from main()

**************************************************************/


int
get_line_transfer_mode ()
{
  rdint
    ("Line_transfer(0=pure.abs,1=pure.scat,2=sing.scat,3=escape.prob,6=macro_atoms,7=macro_atoms+aniso.scattering)",
     &geo.line_mode);


  /* ?? ksl Next section seems rather a kluge.  Why don't we specifty the underlying variables explicitly 
     It also seems likely that we have mixed usage of some things, e.g geo.rt_mode and geo.macro_simple */

  /* JM 1406 -- geo.rt_mode and geo.macro_simple control different things. geo.rt_mode controls the radiative
     transfer and whether or not you are going to use the indivisible packet constraint, so you can have all simple 
     ions, all macro-atoms or a mix of the two. geo.macro_simple just means one can turn off the full macro atom 
     treatment and treat everything as 2-level simple ions inside the macro atom formalism */

  /* For now handle scattering as part of a hidden line transfermode ?? */
  if (geo.line_mode == 4)
    {
      geo.scatter_mode = 1;	// Turn on anisotropic scattering
      geo.line_mode = 3;	// Drop back to escape probabilities
      geo.rt_mode = 1;		// Not macro atom (SS)
    }
  else if (geo.line_mode == 5)
    {
      geo.scatter_mode = 2;	// Thermal trapping model
      geo.line_mode = 3;	// Single scattering model is best for this mode
      geo.rt_mode = 1;		// Not macro atom (SS) 
    }
  else if (geo.line_mode == 6)
    {
      geo.scatter_mode = 0;	// isotropic
      geo.line_mode = 3;	// Single scattering
      geo.rt_mode = 2;		// Identify macro atom treatment (SS)
      geo.macro_simple = 0;	// We don't want the all simple case (SS)
    }
  else if (geo.line_mode == 7)
    {
      geo.scatter_mode = 2;	// thermal trapping
      geo.line_mode = 3;	// Single scattering
      geo.rt_mode = 2;		// Identify macro atom treatment (SS)
      geo.macro_simple = 0;	// We don't want the all simple case (SS)
    }
  else if (geo.line_mode == 8)
    {
      geo.scatter_mode = 0;	// isotropic
      geo.line_mode = 3;	// Single scattering
      geo.rt_mode = 2;		// Identify macro atom treatment i.e. indivisible packets
      geo.macro_simple = 1;	// This is for test runs with all simple ions (SS)
    }
  else if (geo.line_mode == 9)	// JM 1406 -- new mode, as mode 7, but scatter mode is 1
    {
      geo.scatter_mode = 1;	// anisotropic scatter mode 1
      geo.line_mode = 3;	// Single scattering
      geo.rt_mode = 2;		// Identify macro atom treatment 
      geo.macro_simple = 0;	// We don't want the all simple case 
    }
  else
    {
      geo.scatter_mode = 0;	// isotropic
      geo.rt_mode = 1;		// Not macro atom (SS)
    }

  return (0);
}



/***********************************************************
             University of Southampton

Synopsis: 
  get_radiation_sources
   
Arguments:    

Returns:
 
 
Description:  

Notes:

History:
  1502  JM  Moved here from main()

**************************************************************/

int
get_radiation_sources ()
{
  if (geo.system_type != SYSTEM_TYPE_AGN)
    {				/* If is a stellar system */
      rdint ("Star_radiation(y=1)", &geo.star_radiation);
      if (geo.disk_type!=DISK_NONE){
      rdint ("Disk_radiation(y=1)", &geo.disk_radiation);
      }
      else {
	      geo.disk_radiation=0;
      }
      rdint ("Boundary_layer_radiation(y=1)", &geo.bl_radiation);
      rdint ("Wind_radiation(y=1)", &geo.wind_radiation);
      geo.agn_radiation = 0;	// So far at least, our star systems don't have a BH
    }
  else				/* If it is an AGN */
    {
      geo.star_radiation = 0;	// 70b - AGN do not have a star at the center */
      rdint ("Disk_radiation(y=1)", &geo.disk_radiation);
      geo.bl_radiation = 0;
      rdint ("Wind_radiation(y=1)", &geo.wind_radiation);
      geo.agn_radiation = 1;
      rdint ("QSO_BH_radiation(y=1)", &geo.agn_radiation);
    }

  if (!geo.star_radiation && !geo.disk_radiation && !geo.bl_radiation
      && !geo.bl_radiation && !geo.agn_radiation)
    {
      Error ("python: No radiation sources so nothing to do but quit!\n");
      exit (0);
    }

  /* 
     With the macro atom approach we won't want to generate photon 
     bundles in the wind so switch it off here. (SS)
   */

  if (geo.rt_mode == 2)
    {
      Log
	("python: Using Macro Atom method so switching off wind radiation.\n");
      geo.wind_radiation = 0;
    }


  /* 080517 - ksl - Reassigning bb to -1, etc is to make room for reading in model
     grids, but complicates what happens if one tries to restart a model.  This needs
     to be updated so one can re-read the geo file, proabbly by defining variaables 
     BB etc, and then by checking whether or not the type is assigned to BB or read
     in as 0.  Also need to store each of these model list names in geo structure.
   */

  get_spectype (geo.star_radiation,
		"Rad_type_for_star(0=bb,1=models)_to_make_wind",
		&geo.star_ion_spectype);

  get_spectype (geo.disk_radiation,
		"Rad_type_for_disk(0=bb,1=models)_to_make_wind",
		&geo.disk_ion_spectype);

  get_spectype (geo.bl_radiation,
		"Rad_type_for_bl(0=bb,1=models,3=pow)_to_make_wind",
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

  return (0);
}



/***********************************************************
             University of Southampton

Synopsis: 
  get_wind_params calls the relevant subroutine to get wind parameters
  according to the wind type specified 
   
Arguments:		

Returns:
 
 
Description:	

Notes:

History:
	1502  JM 	Moved here from main()

**************************************************************/

int
get_wind_params (ndom)
     int ndom;
{

  if (geo.system_type == SYSTEM_TYPE_AGN)
    {
      geo.rmax = 50. * geo.r_agn;
    }

  // XXX These need to be initalized sensibly and 
  // it is not obvious that is happenning

  zdom[ndom].rmax = 1e12;
  zdom[ndom].twind = 1e5;

  rddoub ("wind.radmax(cm)", &zdom[ndom].rmax);
  
  geo.rmax_sq=geo.rmax*geo.rmax;    
  
  
  /* The next lines assure that geo.rmax_sq really does define the edge of the grid */
  if (zdom[ndom].rmax>geo.rmax){
	  geo.rmax=zdom[ndom].rmax;
	  geo.rmax_sq=geo.rmax*geo.rmax;
	  Log("XXXX  size  %e  %e\n",geo.rmax,geo.rmax_sq);
  }

  /* ksl XXX - There is something of a philosophicla problem that needs to be worked
   * out with geo.rmax and zdom[ndom].rmax for the general case of winds.  Suppose
   * we wish to create, say a spherical outflow with two domains one going from 
   * r1 to r2 and the other going from r2 to r3.  Then we want to keep geo.rmax which is 
   * intended to be the distance beyond which photons are moving through free space separate
   * from the values in the wind zones.  Right now we are setting the outer limit of each
   * wind to be geo.rmax regardless, in routines like get_stellar_wind_params and get_sv_wind
   * This is not what we want.  What should happen is that for each componetn where it is
   * relevant we should ask for the outer edge of the domain and then at the end we should determine
   * what geo.rmax should be set to.  There are some cases, e.g. get_hydor_wind where one should not
   * need to ask the question about rmax, but others where it is necessary
   */

  rddoub ("wind.t.init", &geo.twind);



  /* Now get parameters that are specific to a given wind model

     Note: When one adds a new model, the only things that should be read in and modified
     are parameters in geo.  This is in order to preserve the ability to continue a calculation
     with the same basic wind geometry, without reading in all of the input parameters.  
   */

  if (zdom[ndom].wind_type == 1)
    {
      get_stellar_wind_params (ndom);
    }
  else if (zdom[ndom].wind_type == 0)
    {
      get_sv_wind_params (ndom);
    }
  else if (zdom[ndom].wind_type == 3)
    {
      get_hydro_wind_params (ndom);
    }
  else if (zdom[ndom].wind_type == 4)
    {
      get_corona_params (ndom);
    }
  else if (zdom[ndom].wind_type == 5)
    {
      get_knigge_wind_params (ndom);
    }
  else if (zdom[ndom].wind_type == 6)
    {
      get_homologous_params (ndom);
    }
  else if (zdom[ndom].wind_type == 7)
    {
      get_yso_wind_params (ndom);
    }
  else if (zdom[ndom].wind_type == 8)
    {
      get_elvis_wind_params (ndom);
    }
  else if (zdom[ndom].wind_type == 9)	//NSH 18/2/11 This is a new wind type to produce a thin shell.
    {
      get_shell_wind_params (ndom);
/*NSH 121219 moved	  dfudge = (geo.wind_rmax - geo.wind_rmin) / 1000.0;	Stop photons getting pushed out of the cell 
Modified again in python 71b to take account of change in parametrisation of shell wind 
	  DFUDGE = dfudge; */
    }
  else if (zdom[ndom].wind_type != 2)
    {
      /* XXX this is part of the new problem with a previous wind model */
      Error ("python: Unknown wind type %d\n", zdom[ndom].wind_type);
      exit (0);
    }

  /* Get the filling factor of the wind */
  // XXX  This is not in the right place.  It provides a golal filling factor to our
  // models

	geo.fill = 1.;
  if (geo.wind_type != 3) //At present, we wont ask this question if we have a read in hydro model.
  	rddoub ("wind.filling_factor(1=smooth,<1=clumped)", &geo.fill);


  return (0);
}


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
  double lstar;

  /* Describe the basic binary star system */

  geo.mstar /= MSOL;		// Convert to MSOL for ease of data entry
  rddoub ("mstar(msol)", &geo.mstar);
  geo.mstar *= MSOL;

  /* If a BH we want geo.rstar to be at least as large as the last stable orbit for
   * a non-rotating BH
   */

  if (geo.system_type == SYSTEM_TYPE_AGN)
    {
      geo.rstar = 6. * G * geo.mstar / (C * C);	//correction - ISCO is 6x Rg NSH 121025
    }

  rddoub ("rstar(cm)", &geo.rstar);


  geo.r_agn = geo.rstar;	/* At present just set geo.r_agn to geo.rstar */
  geo.rstar_sq = geo.rstar * geo.rstar;
  if (geo.star_radiation)
    rddoub ("tstar", &geo.tstar);

  lstar =
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

  return (lstar);
}

/***********************************************************
             University of Southampton

Synopsis: 
  get_disk_params sets up the disk parameters according to user inputs, 
  e.g. the temperature profile, accretion rate etc.
   
Arguments:		

Returns:
  disk_illum - this is used by python.c and so needs to be returned
 
Description:	

Notes:

History:
	1502  JM 	Moved here from main()
	1510	ksl	Modified to restore illumination
			options, which were brokedn

**************************************************************/


double
get_disk_params ()
{
	// XXX Commenting lines out instead of fixing a problem is not good practice
//        if (geo.disk_radiation) /*NSH 130906 - Commented out this if loop. It was causing problems with restart - bug #44
//          {
  geo.disk_mdot /= (MSOL / YR);	// Convert to msol/yr to simplify input
  rddoub ("disk.mdot(msol/yr)", &geo.disk_mdot);
  geo.disk_mdot *= (MSOL / YR);
  rdint
    ("Disk.illumination.treatment(0=no.rerad,1=high.albedo,2=thermalized.rerad,3=analytic)",
     &geo.disk_illum);
  rdint ("Disk.temperature.profile(0=standard;1=readin)", &geo.disk_tprofile);
  if (geo.disk_tprofile == 1)
    {
      rdstr ("T_profile_file", files.tprofile);
    }
//          }
//        else
//          {
//            geo.disk_mdot = 0;
//            disk_illum = 0;
//          }



  /* Set a default for diskrad for an AGN */
  if (geo.system_type == SYSTEM_TYPE_AGN)
    {
      geo.diskrad = 100. * geo.r_agn;
    }

  rddoub ("disk.radmax(cm)", &geo.diskrad);
  Log ("geo.diskrad  %e\n", geo.diskrad);

  geo.diskrad_sq = geo.diskrad * geo.diskrad;

/* If diskrad <= geo.rstar set geo.disk_type = DISK_NONE to make any disk transparent anyway. */

  if (geo.diskrad < geo.rstar)
    {
      Log ("Disk radius is less than star radius, so assuming no disk)\n");
      geo.disk_type = DISK_NONE;
    }

  if (geo.disk_type == DISK_VERTICALLY_EXTENDED)
    {				/* Get the additional variables need to describe a vertically extended disk */
      rddoub ("disk.z0(fractional.height.at.diskrad)", &geo.disk_z0);
      rddoub ("disk.z1(powerlaw.index)", &geo.disk_z1);
    }
  return (0);
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

  /* Describe the boundary layer */

  //OLD 130622      if (geo.bl_radiation )    Change made to allow a power law boundary layer
  if (geo.bl_radiation && geo.bl_ion_spectype != SPECTYPE_POW)
    {
      xbl = geo.lum_bl = 0.5 * G * geo.mstar * geo.disk_mdot / geo.rstar;

      rddoub ("lum_bl(ergs/s)", &geo.lum_bl);
      Log ("OK, the bl lum will be about %.2e the disk lum\n",
	   geo.lum_bl / xbl);
      rddoub ("t_bl", &geo.t_bl);
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

      rddoub ("lum_agn(ergs/s)", &geo.lum_agn);
      Log ("OK, the agn lum will be about %.2e the disk lum\n",
	   geo.lum_agn / xbl);
		if (geo.agn_ion_spectype == SPECTYPE_POW || geo.agn_ion_spectype == SPECTYPE_CL_TAB)
		{
      geo.alpha_agn = (-1.5);
      rddoub ("agn_power_law_index", &geo.alpha_agn);
      geo.const_agn =
	geo.lum_agn /
	(((pow (2.42e18, geo.alpha_agn + 1.)) -
	  pow (4.84e17, geo.alpha_agn + 1.0)) / (geo.alpha_agn + 1.0));
      Log ("AGN Input parameters give a power law constant of %e\n",
	   geo.const_agn);
	}
	else if (geo.agn_ion_spectype == SPECTYPE_BREM)
	{
		geo.brem_temp=1.16e8; //10kev
		geo.const_agn=1.0;
		rddoub ("agn_bremsstrahung_temp(K)",&geo.brem_temp);
		temp_const_agn = geo.lum_agn / qromb(integ_brem,4.84e17,2.42e18,1e-4);
		geo.const_agn=temp_const_agn;
      Log ("AGN Input parameters give a Bremsstrahlung constant of %e\n", temp_const_agn);
		
	}

      /* JM 1502 -- lines to add a low frequency power law cutoff. accessible
         only in advanced mode. default is zero which is checked before we call photo_gen_agn */
      geo.pl_low_cutoff = 0.0;
      if (modes.iadvanced)
	rddoub ("agn_power_law_cutoff", &geo.pl_low_cutoff);


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
	rddoub ("agn_power_law_cutoff", &geo.pl_low_cutoff);


      /* Computes the constant for the power law spectrum from the input alpha and 2-10 luminosity. 
         This is only used in the sim correction factor for the first time through. 
         Afterwards, the photons are used to compute the sim parameters. */

      geo.const_agn =
	geo.lum_agn /
	(((pow (2.42e18, geo.alpha_agn + 1.)) -
	  pow (4.84e17, geo.alpha_agn + 1.0)) / (geo.alpha_agn + 1.0));
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
  get_meta_params reads in data pertaining to simulation meta-
  properties like reverberation mapping settings and variance
  reduction techniques.
   
Arguments:    
Returns:
 
Description:  
Notes:
History:
  1504  SWM   Added
**************************************************************/

int
get_meta_params (void)
{
  int meta_param;

  rdint ("reverb.type", &meta_param);
  switch (meta_param)
    {
    case 0:
      geo.reverb = REV_NONE;
      break;
    case 1:
      geo.reverb = REV_PHOTON;
      break;
    case 2:
      geo.reverb = REV_WIND;
      break;
    default:
      geo.reverb = REV_NONE;
      break;
    }

  if (geo.reverb == REV_WIND)
    {
      geo.reverb_path_bins = 30;
      geo.reverb_theta_bins = 30;
      rdint ("reverb.path_bins", &geo.reverb_path_bins);
      rdint ("reverb.theta_bins", &geo.reverb_theta_bins);
    }

  return (0);
}


/***********************************************************
             University of Southampton

Synopsis: 
  setup_dfudge works out dfudge and returns it to the user.
  the global variable DFUDGE is not altered here.
   
Arguments:		

Returns:
  dfudge 	double	
  			the push through distance
Description:	

Notes:

History:
	1502  JM 	Moved here from main()

**************************************************************/

double
setup_dfudge ()
{
  double dfudge;

  /* 121219 NSH Set up DFUDGE to be a value that makes some kind of sense
     given the scale of the wind. Up till py74b2 it was set to be fixed at
     1e5, so we ensure that this is a minimum, so any winds of CV type scale
     will keep the old dfudge, and hopefully look the same. We also need to
     set dfudge slightly differently for the shell wind. */

  if (geo.wind_type == 9)
    {
      dfudge = (geo.rmax - geo.rmin) / 1000.0;
    }
  else
    {
      if (geo.rmax / 1.e10 < 1e5)
	{
	  dfudge = 1e5;
	  Log ("DFUDGE set to minimum value of %e\n", dfudge);
	}
      else
	{
	  dfudge = geo.rmax / 1.e10;
	  Log ("DFUDGE set to %e based on geo.rmax\n", dfudge);
	}
    }

  return (dfudge);		
}



/***********************************************************
             University of Southampton

Synopsis: 
  setup_windcone sets up the windcone 
   
Arguments:		

Returns:
  dfudge 	double	
  			the push through distance
Description:	

Notes:
  The angles thetamin and
  thetamax are all defined from the z axis, so that an angle of 0
  is a flow that is perpeindicular to to the disk and one that is
  close to 90 degrees will be parallel to the plane of the disk
  geo.wind_thetamin and max are defined in the routines that initialize
  the various wind models, e. g. get_sv_wind_parameters. These
  have been called at this point.  

  z is the place where the windcone intercepts the z axis
  dzdr is the slope 

  111124 fixed notes on this - ksl
History:
	1502  JM 	Moved here from main()
	1508	ksl	Modified to construct wind cones  fo
			all domains

**************************************************************/

int
setup_windcone ()
{
  int ndom;

  for (ndom = 0; ndom < geo.ndomain; ndom++)
    {

      if (zdom[ndom].wind_thetamin > 0.0)
	{
	  zdom[ndom].windcone[0].dzdr = 1. / tan (zdom[ndom].wind_thetamin);
	  zdom[ndom].windcone[0].z =
	    (-zdom[ndom].wind_rho_min / tan (zdom[ndom].wind_thetamin));
	}
      else
	{
	  zdom[ndom].windcone[0].dzdr = VERY_BIG;
	  zdom[ndom].windcone[0].z = -VERY_BIG;;
	}


      if (zdom[ndom].wind_thetamax > 0.0)
	{
	  zdom[ndom].windcone[1].dzdr = 1. / tan (zdom[ndom].wind_thetamax);
	  zdom[ndom].windcone[1].z =
	    (-zdom[ndom].wind_rho_max / tan (zdom[ndom].wind_thetamax));
	}
      else
	{
	  zdom[ndom].windcone[1].dzdr = VERY_BIG;
	  zdom[ndom].windcone[1].z = -VERY_BIG;;
	}
    }
  return (0);
}



/***********************************************************
             University of Southampton

Synopsis: 
  get_spectrum_params 
   
Arguments:    

Returns:
  dfudge  double  
        the push through distance
Description:  

Notes:
  The angles thetamin and
  thetamax are all defined from the z axis, so that an angle of 0
  is a flow that is perpeindicular to to the disk and one that is
  close to 90 degrees will be parallel to the plane of the disk
  geo.wind_thetamin and max are defined in the routines that initialize
  the various wind models, e. g. get_sv_wind_parameters. These
  have been called at this point.  

  z is the place where the windcone intercepts the z axis
  dzdr is the slope 

  111124 fixed notes on this - ksl
History:
  1502  JM  Moved here from main()

**************************************************************/









/***********************************************************
             University of Southampton

Synopsis: 
  setup_created_files 
   
Arguments:    

Returns:

Description:  

Notes:

History:
  1502  JM  Moved here from main()

**************************************************************/

int
setup_created_files ()
{
  int opar_stat;

  opar_stat = 0;		/* 59a - ksl - 08aug - Initialize opar_stat to indicate that if we do not open a rdpar file, 
				   the assumption is that we are reading from the command line */

  if (strncmp (files.root, "dummy", 5) == 0)
    {
      Log
	("Proceeding to create rdpar file in dummy.pf, but will not run prog\n");
    }

  else if (strncmp (files.root, "stdin", 5) == 0
	   || strncmp (files.root, "rdpar", 5) == 0 || files.root[0] == ' '
	   || strlen (files.root) == 0)
    {
      strcpy (files.root, "mod");
      Log
	("Proceeding in interactive mode\n Output files will have rootname mod\n");
    }

  else
    {
      strcpy (files.input, files.root);
      strcat (files.input, ".pf");

      if ((opar_stat = opar (files.input)) == 2)
	{
	  Log ("Reading data from file %s\n", files.input);
	}
      else
	{
	  Log ("Creating a new parameter file %s\n", files.input);
	}

    }


  /* Now create the names of all the files which will be written.  Note that some files
     have the same root as the input file, while others have a generic name of python.
     This is intended so that files which you really want to keep have unique names, while
     those which are for short-term diagnostics are overwritten.  ksl 97aug. */

  strcpy (basename, files.root);	//56d -- ksl --Added so filenames could be created by routines as necessary

  strcpy (files.wspec, files.root);
  strcpy (files.lspec, files.root);

  strcpy (files.spec, files.root);
  strcpy (files.windrad, "python");
  strcpy (files.windsave, files.root);
  strcpy (files.specsave, files.root);

  /* 130722 JM we now save python.phot and disk.diag files under diag_root folder */
  strcpy (files.phot, files.diagfolder);
  strcpy (files.disk, files.diagfolder);
  strcat (files.phot, "python");
  strcat (files.disk, files.root);

  strcat (files.wspec, ".spec_tot");
  strcat (files.lspec, ".log_spec_tot");
  strcat (files.spec, ".spec");
  strcat (files.windrad, ".wind_rad");
  strcat (files.windsave, ".wind_save");
  strcat (files.specsave, ".spec_save");
  strcat (files.phot, ".phot");
  strcat (files.disk, ".disk.diag");


  return (opar_stat);
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
      rdint ("Use.standard.care.factors(1=yes)", &istandard);

      if (!istandard)
	{
	  rddoub ("Fractional.distance.photon.may.travel", &SMAX_FRAC);
	  rddoub ("Lowest.ion.density.contributing.to.photoabsorption",
		  &DENSITY_PHOT_MIN);
	  rdint ("Keep.photoabs.during.final.spectrum(1=yes)",
		 &modes.keep_photoabs);
	}
    }
  return (0);
}
