
/***********************************************************/
/** @file  inspect_wind.c
 * @author ksl
 * @date   October, 2021
 *
 * @brief  Routines to inspect variables in a wind structure 
 *
 *###Notes###
 * This is intended just as a diagnostic routine 
 * so that one can print out whatever variables in
 * a windstrucutre one wants in order to diagnose
 * a problem.  It was written so that we could inspect
 * some of the macro atom variables in paralell mode
 * in diagnosing issue #898 and #910, but anyon 
 * should change it so that other problems might 
 * be addressed.
 *
 *
 ***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "atomic.h"
#include "sirocco.h"


char inroot[LINELENGTH], outroot[LINELENGTH], model_file[LINELENGTH], folder[LINELENGTH];
int model_flag, ksl_flag, cmf2obs_flag, obs2cmf_flag;

#define MAX_HR1985_LINES 20
int hr1985_flag, n_hr1985;
double hr1985_lambda[MAX_HR1985_LINES];

double line_matom_lum_single (double lum[], PlasmaPtr xplasma, int uplvl);
int line_matom_lum (int uplvl);
int create_matom_level_map ();
int hr1985_dump ();

/**********************************************************/
/**
 * @brief      parses the command line options
 *
 * @param [in]  int  argc   the number of command line arguments
 * @param [in]  char *  argv[]   The command line arguments
 *
 *
 * ###Notes###
 *
 * The general purpose of each of the command line options
 * should be fairly obvious from reading the code.
 *
 *
 * Although this routine uses the standard Log and Error commands
 * the diag files have not been created yet and so this information
 * is really simply written to the terminal.
 *
 **********************************************************/

int
xparse_command_line (argc, argv)
     int argc;
     char *argv[];
{
  int j = 0;
  int i;
  char dummy[LINELENGTH];
  int mkdir ();
  char *fgets_rc;


  sprintf (outroot, "%s", "new");

  model_flag = ksl_flag = obs2cmf_flag = cmf2obs_flag = 0;
  hr1985_flag = n_hr1985 = 0;

  if (argc == 1)
  {
    printf ("Parameter file name (e.g. my_model.pf, or just my_model):");
    fgets_rc = fgets (dummy, LINELENGTH, stdin);
    if (!fgets_rc)
    {
      printf ("Input rootname is NULL or invalid\n");
      exit (1);
    }
    get_root (inroot, dummy);
  }
  else
  {

    for (i = 1; i < argc; i++)
    {
      if (strcmp (argv[i], "-out_root") == 0)
      {
        if (sscanf (argv[i + 1], "%s", dummy) != 1)
        {
          printf ("sirocco: Expected out_root after -out_root switch\n");
          exit (0);
        }

        get_root (outroot, dummy);
        i++;
        j = i;

      }
      if (strcmp (argv[i], "-model_file") == 0)
      {
        if (sscanf (argv[i + 1], "%s", dummy) != 1)
        {
          printf ("sirocco: Expected a model file containing density, velocity and temperature after -model_file switch\n");
          exit (0);
        }
        get_root (model_file, dummy);
        i++;
        j = i;
        printf ("got a model file %s\n", model_file);
        model_flag = 1;
      }
      else if (strcmp (argv[i], "-ksl") == 0)
      {
        printf ("Carrying out a simple hard wired ion modification\n");
        ksl_flag = 1;
      }
      else if (strcmp (argv[i], "--dry-run") == 0)
      {
        modes.quit_after_inputs = 1;
        j = i;
      }
      else if (strcmp (argv[i], "-cmf") == 0)
      {
        obs2cmf_flag = 1;
      }
      else if (strcmp (argv[i], "-obs") == 0)
      {
        cmf2obs_flag = 1;
      }
      else if (strcmp (argv[i], "-hr1985") == 0)
      {
        /* comma separated list of wavelengths in Angstroms, e.g. -hr1985 1215.67,1548.19,6562.8 */
        char *tok;
        if (i + 1 >= argc || argv[i + 1][0] == '-')
        {
          printf ("sirocco: Expected a comma separated list of wavelengths (A) after -hr1985 switch\n");
          exit (0);
        }
        strcpy (dummy, argv[i + 1]);
        n_hr1985 = 0;
        for (tok = strtok (dummy, ","); tok != NULL && n_hr1985 < MAX_HR1985_LINES; tok = strtok (NULL, ","))
        {
          hr1985_lambda[n_hr1985++] = atof (tok);
        }
        hr1985_flag = 1;
        i++;
        j = i;
      }
      else if (strncmp (argv[i], "-", 1) == 0)
      {
        printf ("sirocco: Unknown switch %s\n", argv[i]);
        exit (0);
      }
    }

    /* The last command line variable is always the windsave file */

    if (j + 1 == argc)
    {
      printf ("All of the command line has been consumed without specifying a parameter file name, so exiting\n");
      exit (1);
    }
    strcpy (dummy, argv[argc - 1]);
    get_root (inroot, dummy);

  }

  return (0);
}



/**********************************************************/
/**
 * @brief      the main routine which carries out the effort
 *
 * @param [in]  int  argc   the number of command line arguments
 * @param [in]  char *  argv[]   The command line arguments
 *
 *
 * ###Notes###
 *
 * This routine oversees the effort.  The basic steps are
 *
 * - parse the command line to get the names of files
 * - read the old windsave file
 * - read the densities from in this case H
 * - modify the densities
 * - write out the new windsave file
 *
 *
 **********************************************************/


int
main (argc, argv)
     int argc;
     char *argv[];
{

  char infile[LINELENGTH], outfile[LINELENGTH];
  int n, i;
  FILE *fptr, *fopen ();
  int ii, jj, ndom, nnwind;
  int mkdir ();


  xparse_command_line (argc, argv);

  sprintf (infile, "%.150s.wind_save", inroot);
  sprintf (outfile, "%.150s.txt", inroot);

  printf ("Reading %s and writing to %s\n", infile, outfile);

  zdom = calloc (MAX_DOM, sizeof (domain_dummy));
  if (zdom == NULL)
  {
    printf ("Unable to allocate memory for domain\n");
    return EXIT_FAILURE;
  }

  wind_read (infile);

  if (hr1985_flag)
  {
    hr1985_dump ();
    exit (0);
  }

  if (nlevels_macro == 0)
  {
    printf ("Currently this routine only looks at macro atom values, and this is a simple atom file\n");
    exit (0);
  }



  fptr = fopen (outfile, "w");

  fprintf (fptr, "# Results for %s\n", infile);

  fprintf (fptr, "# Size  %d  %d %d \n", size_Jbar_est, size_gamma_est, size_alpha_est);

  fprintf (fptr, "%15s %4s %4s %4s ", "Variable", "np", "i", "j");
  for (i = 0; i < nlevels_macro; i++)
  {
    fprintf (fptr, "MacLev%02d ", i);
  }
  fprintf (fptr, "\n");

  fprintf (fptr, "%15s %4s %4s %4s ", "---------------", "----", "----", "----");
  for (i = 0; i < nlevels_macro; i++)
  {
    fprintf (fptr, "-------- ");
  }
  fprintf (fptr, "\n");

  for (n = 0; n < NPLASMA; n++)
  {

    fprintf (fptr, "%15s", "jbar");

    nnwind = plasmamain[n].nwind;
    ndom = wmain[nnwind].ndom;
    wind_n_to_ij (ndom, nnwind, &ii, &jj);
    fprintf (fptr, " %4d %4d %4d ", n, ii, jj);

    for (i = 0; i < nlevels_macro; i++)
    {
      fprintf (fptr, "%8.2e ", macromain[n].jbar[i]);
    }
    fprintf (fptr, "\n");
  }



  for (n = 0; n < NPLASMA; n++)
  {
    fprintf (fptr, "%15s", "jbar_old");

    nnwind = plasmamain[n].nwind;
    ndom = wmain[nnwind].ndom;
    wind_n_to_ij (ndom, nnwind, &ii, &jj);
    fprintf (fptr, " %4d %4d %4d ", n, ii, jj);

    for (i = 0; i < nlevels_macro; i++)
    {
      fprintf (fptr, "%8.2e ", macromain[n].jbar_old[i]);
    }
    fprintf (fptr, "\n");
  }



  for (n = 0; n < NPLASMA; n++)
  {
    fprintf (fptr, "%15s", "gamma");

    nnwind = plasmamain[n].nwind;
    ndom = wmain[nnwind].ndom;
    wind_n_to_ij (ndom, nnwind, &ii, &jj);
    fprintf (fptr, " %4d %4d %4d ", n, ii, jj);

    for (i = 0; i < nlevels_macro; i++)
    {
      fprintf (fptr, "%8.2e ", macromain[n].gamma[i]);
    }
    fprintf (fptr, "\n");
  }




  for (n = 0; n < NPLASMA; n++)
  {
    fprintf (fptr, "%15s", "alpha_st");

    nnwind = plasmamain[n].nwind;
    ndom = wmain[nnwind].ndom;
    wind_n_to_ij (ndom, nnwind, &ii, &jj);
    fprintf (fptr, " %4d %4d %4d ", n, ii, jj);

    for (i = 0; i < nlevels_macro; i++)
    {
      fprintf (fptr, "%8.2e ", macromain[n].alpha_st[i]);
    }
    fprintf (fptr, "\n");
  }




  for (n = 0; n < NPLASMA; n++)
  {
    fprintf (fptr, "%15s", "recomb_sp");

    nnwind = plasmamain[n].nwind;
    ndom = wmain[nnwind].ndom;
    wind_n_to_ij (ndom, nnwind, &ii, &jj);
    fprintf (fptr, " %4d %4d %4d ", n, ii, jj);

    for (i = 0; i < nlevels_macro; i++)
    {
      fprintf (fptr, "%8.2e ", macromain[n].recomb_sp[i]);
    }
    fprintf (fptr, "\n");
  }




  for (n = 0; n < NPLASMA; n++)
  {
    fprintf (fptr, "%15s", "matom_abs");

    nnwind = plasmamain[n].nwind;
    ndom = wmain[nnwind].ndom;
    wind_n_to_ij (ndom, nnwind, &ii, &jj);
    fprintf (fptr, " %4d %4d %4d ", n, ii, jj);

    for (i = 0; i < nlevels_macro; i++)
    {
      fprintf (fptr, "%8.2e ", macromain[n].matom_abs[i]);
    }
    fprintf (fptr, "\n");
  }




  for (n = 0; n < NPLASMA; n++)
  {
    fprintf (fptr, "%15s", "matom_emiss");

    nnwind = plasmamain[n].nwind;
    ndom = wmain[nnwind].ndom;
    wind_n_to_ij (ndom, nnwind, &ii, &jj);
    fprintf (fptr, " %4d %4d %4d ", n, ii, jj);

    for (i = 0; i < nlevels_macro; i++)
    {
      fprintf (fptr, "%8.2e ", macromain[n].matom_emiss[i]);
    }
    fprintf (fptr, "\n");
  }
  fclose (fptr);


  /* calculate all the line luminosities in macro-atom mode */
  /* these are stored in a folder with a separate file for each upper level */
  sprintf (folder, "matom_linelum_%.200s", inroot);
  mkdir (folder, 0777);
  printf ("Calculating macro-atom line luminosities for all lines in wavelength range %.1f to %.1f Angstroms...\n",
          VLIGHT / geo.sfmax / ANGSTROM, VLIGHT / geo.sfmin / ANGSTROM);

  /* create a file that tells the user the wavelengths of every macro-atom line and the 
     corresponding upper and lower levels */
  create_matom_level_map ();

  /* loop over all macro atom upper levels */
  for (i = 0; i < nlevels_macro; i++)
  {
    line_matom_lum (i);
  }
  printf ("Done for %d levels.\n", nlevels_macro);

  exit (0);

}

/**********************************************************/
/**
 * @brief print out which line wavelengths correspond to which upper and lower levels 
 *
 **********************************************************/

int
create_matom_level_map ()
{
  int uplvl, nbbd, n;
  char outfile[LINELENGTH];
  FILE *fptr, *fopen ();

  /* open a file in the folder where we store the matom line luminosities */
  sprintf (outfile, "%.200s/line_map.txt", folder);

  /* print some header information to the file */
  fptr = fopen (outfile, "w");
  fprintf (fptr, "# model name %s\n", inroot);
  fprintf (fptr, "# which line wavelengths (Angstroms) correspond to which upper and lower levels\n");
  fprintf (fptr, "upper lower wavelength\n");

  for (uplvl = 0; uplvl < nlevels_macro; uplvl++)
  {
    nbbd = xconfig[uplvl].n_bbd_jump;
    for (n = 0; n < nbbd; n++)
    {
      fprintf (fptr, "%d %d %12.2f\n", uplvl, line[xconfig[uplvl].bbd_jump[n]].nconfigl,
               VLIGHT / line[xconfig[uplvl].bbd_jump[n]].freq / ANGSTROM);
    }
  }
  fclose (fptr);
  return (0);
}


/**********************************************************/
/**
 * @brief calculate line luminosities for all cells for a given upper level and save to file
 *
 **********************************************************/

int
line_matom_lum (uplvl)
     int uplvl;
{
  int n, nbbd, i, ii, jj, nnwind, ndom, inwind;
  double lum[NBBJUMPS];
  char outfile[LINELENGTH];
  FILE *fptr, *fopen ();

  nbbd = xconfig[uplvl].n_bbd_jump;

  /* open a file in the folder where we store the matom line luminosities */
  sprintf (outfile, "%.100s/linelums_%.100s_lev%d.txt", folder, inroot, uplvl);

  /* print some header information to the file */
  fptr = fopen (outfile, "w");
  fprintf (fptr, "# model name %s\n", inroot);
  fprintf (fptr, "# Line luminosities from upper level %d from %.2f to %.2f Angstroms\n", uplvl, VLIGHT / geo.sfmax / ANGSTROM,
           VLIGHT / geo.sfmin / ANGSTROM);
  fprintf (fptr, "# Lower Levels:     ");
  for (n = 0; n < nbbd; n++)
  {
    fprintf (fptr, "%12d", line[xconfig[uplvl].bbd_jump[n]].nconfigl);
  }
  fprintf (fptr, "\n# Line Wavelengths:");
  for (n = 0; n < nbbd; n++)
  {
    fprintf (fptr, " %12.2f", VLIGHT / line[xconfig[uplvl].bbd_jump[n]].freq / ANGSTROM);
  }
  fprintf (fptr, "\n");

  /* columns of wind cells */
  fprintf (fptr, "%12s %12s %12s %4s %4s %12s %12s %12s", "nwind", "x", "z", "i", "j", "inwind", "wind_vol", "plasma_vol");

  /* header info for each lower level */
  for (n = 0; n < nbbd; n++)
  {
    fprintf (fptr, "  LowerLev%03d  ", line[xconfig[uplvl].bbd_jump[n]].nconfigl);
  }
  fprintf (fptr, "\n");
  // fprintf (fptr, "%12s %12s %12s %4s %4s %12s %12s %12s", "------------", "------------", "------------", "----", "----", "------------",
  //          "------------", "------------");
  // for (n = 0; n < nbbd; n++)
  // {
  //   fprintf (fptr, "  -----------  ");
  // }
  // fprintf (fptr, "\n");


  /* loop over each plasma cell and do the calculation */
  for (nnwind = 0; nnwind < NDIM2; nnwind++)
  {
    ndom = wmain[nnwind].ndom;
    inwind = wmain[nnwind].inwind;
    wind_n_to_ij (ndom, nnwind, &ii, &jj);
    fprintf (fptr, "%12d %12.4e %12.4e %4d %4d %12d %12.4e", nnwind, wmain[nnwind].xcen[0], wmain[nnwind].xcen[2], ii, jj, inwind,
             wmain[nnwind].vol);
    if (wmain[nnwind].inwind >= 0)
    {
      n = wmain[nnwind].nplasma;
      line_matom_lum_single (lum, &plasmamain[n], uplvl);
      /* print the filled volume */
      fprintf (fptr, " %13.4e", plasmamain[n].vol);
      for (i = 0; i < nbbd; i++)
        fprintf (fptr, " %13.4e", lum[i]);
    }
    else
    {
      /* print a 0.0 instead of the filled volume and line lums outside the wind */
      fprintf (fptr, " %13.4e", 0.0);
      for (i = 0; i < nbbd; i++)
        fprintf (fptr, " %13.4e", 0.0);
    }

    fprintf (fptr, "\n");
  }
  fclose (fptr);
  return (0);
}

/**********************************************************/
/**
 * @brief calculate line luminosities for a single cell for a given upper level 
 *
 **********************************************************/

double
line_matom_lum_single (lum, xplasma, uplvl)
     double lum[];
     PlasmaPtr xplasma;
     int uplvl;
{
  int n, nbbd, m;
  double penorm, bb_cont;
  double freq_min, freq_max, lum_tot;
  struct lines *line_ptr;
  double eprbs[NBBJUMPS];
  freq_min = geo.sfmin;
  freq_max = geo.sfmax;
  /* identify number of bb downward jumps */
  nbbd = xconfig[uplvl].n_bbd_jump;
  penorm = 0.0;
  m = 0;
  lum_tot = 0.0;
  /* work out how often we come out in each of these locations */
  for (n = 0; n < nbbd; n++)
  {
    line_ptr = &line[xconfig[uplvl].bbd_jump[n]];
    if ((line_ptr->freq > freq_min) && (line_ptr->freq < freq_max))     // correct range
    {
      // bb_cont = (a21 (line_ptr) * p_escape (line_ptr, xplasma));
      bb_cont = (a21 (line_ptr) * 1.0);
      eprbs[m] = bb_cont * (xconfig[uplvl].ex - xconfig[line[xconfig[uplvl].bbd_jump[n]].nconfigl].ex); //energy difference
      penorm += eprbs[m];
    }
    else
    {
      eprbs[m] = 0.0;
    }
    m++;
  }

  /* correctly normalise the probabilities */
  for (n = 0; n < nbbd; n++)
  {
    if (penorm == 0)
    {
      lum[n] = 0.0;
    }
    else
    {
      eprbs[n] = eprbs[n] / penorm;
      lum[n] = eprbs[n] * macromain[xplasma->nplasma].matom_emiss[uplvl];
    }
    lum_tot += lum[n];
  }

  return (lum_tot);
}


/**********************************************************/
/**
 * @brief The HR1985 two-sided continuum loss probability, for diagnostics
 *
 * @param [in] double beta   continuum-to-line opacity ratio k_C/k_L
 * @param [in] double gamma  1/tau_Sobolev
 * @param [in] double P      the Sobolev escape probability
 * @param [in] int use_table TRUE to use the table read with the atomic data,
 * FALSE for the analytic fit
 * @return  A = beta F(beta) f(beta,gamma), limited so that P + A <= 1, or -1
 * if the table is requested but none was read
 *
 * ### Notes ###
 * Hummer & Rybicki (1985) eqs. 2.30 and 2.38, using hr1985_g_table or
 * hr1985_g_fit from lines.c. A is zero for gamma >= 1 (tau <= 1), as in
 * hr1985_es_loss.
 **********************************************************/

double
hr1985_A (beta, gamma, P, use_table)
     double beta, gamma, P;
     int use_table;
{
  double A;

  if (use_table && hr1985_npts == 0)
    return (-1.0);
  if (beta <= 0.0 || gamma >= 1.0)
    return (0.0);

  A = 2.0 * (use_table ? hr1985_g_table (beta, gamma) : hr1985_g_fit (beta, gamma));
  if (P + A > 1.0)
    A = 1.0 - P;
  return (A);
}


/**********************************************************/
/**
 * @brief Write out the quantities needed for the Hummer & Rybicki (1985)
 * continuum-absorption correction to the escape probability
 *
 * @return Always returns 0
 *
 * @details
 * For each wavelength given with -hr1985, the nearest line in the atomic
 * data is used and, for every plasma cell, we write tau_Sobolev, the
 * thermal speed of the ion, the continuum opacities at the line frequency,
 * beta = kappa_c v_th / (dvds tau), gamma = 1/tau, the Sobolev escape
 * probability P and the continuum absorption probability A from the analytic
 * fit and (if one was read with the atomic data) the HR1985 table. kappa_c is
 * true absorption only (bf + ff). beta_es and A_es repeat the calculation with
 * electron scattering alone as the continuum, in which case A_es is the
 * probability that a trapped line photon is electron scattered out of the
 * resonance zone rather than destroyed. P_sirocco is what p_escape actually
 * returns, which for Lyman alpha includes A_es when a table was read.
 *
 * HR1985 say the continuum term dominates where gamma < 15 beta. The counts
 * printed to the screen are for the electron scattering case, using the table
 * if there is one and the fit otherwise.
 **********************************************************/

int
hr1985_dump ()
{
  int k, n, nline, nplasma, nwind, ndom, ii, jj;
  int ncells, nthick, ndominant, nten;
  double lambda, dlam, dlam_best;
  double tau, dvds, vth, kap_bf_c, kap_ff_c, kap_es_c, beta, beta_es, gamma, P, P_sirocco;
  double A_fit, A_es_fit, A_tab, A_es_tab, A_es;
  double max_ratio;
  char outfile[LINELENGTH];
  struct lines *line_ptr;
  PlasmaPtr xplasma;
  WindPtr one;
  FILE *fptr, *fopen ();

  sprintf (outfile, "%.150s.hr1985.txt", inroot);
  fptr = fopen (outfile, "w");
  printf ("Writing HR1985 diagnostics to %s\n", outfile);

  fprintf (fptr, "# HR1985 continuum absorption diagnostics for %s\n", inroot);
  fprintf (fptr, "# beta = (kappa_bf+kappa_ff) v_th / (dvds tau); beta_es = kappa_es v_th / (dvds tau); gamma = 1/tau\n");
  fprintf (fptr, "# P = Sobolev escape probability, A = beta F(beta) f(beta,gamma) from the analytic fit (_fit) or the table (_tab)\n");
  fprintf (fptr, "# A_es uses beta_es; _tab columns are -1 if no HR1985 table was read; P_sirocco is what p_escape returns\n");
  fprintf (fptr, "# HR1985 table read: %s\n", hr1985_npts > 0 ? "yes" : "no");
  fprintf (fptr,
           "%9s %3s %3s %6s %6s %4s %4s %11s %11s %9s %9s %11s %11s %11s %11s %11s %11s %11s %11s %11s %11s %11s %11s %11s %11s %11s\n",
           "lambda", "z", "ion", "nplasm", "nwind", "i", "j", "x", "zz", "ne", "t_e", "dvds", "tau", "v_th", "kappa_bf", "kappa_ff",
           "kappa_es", "beta", "beta_es", "gamma", "P", "A_fit", "A_es_fit", "A_tab", "A_es_tab", "P_sirocco");

  for (k = 0; k < n_hr1985; k++)
  {
    /* find the line nearest to the requested wavelength */
    nline = -1;
    dlam_best = VERY_BIG;
    for (n = 0; n < nlines; n++)
    {
      dlam = fabs (VLIGHT / line[n].freq / ANGSTROM - hr1985_lambda[k]);
      if (dlam < dlam_best)
      {
        dlam_best = dlam;
        nline = n;
      }
    }
    if (nline < 0)
    {
      printf ("hr1985_dump: no line found near %.2f A\n", hr1985_lambda[k]);
      continue;
    }

    line_ptr = &line[nline];
    lambda = VLIGHT / line_ptr->freq / ANGSTROM;
    ncells = nthick = ndominant = nten = 0;
    max_ratio = 0.0;

    for (nplasma = 0; nplasma < NPLASMA; nplasma++)
    {
      xplasma = &plasmamain[nplasma];
      nwind = xplasma->nwind;
      one = &wmain[nwind];
      ndom = one->ndom;
      wind_n_to_ij (ndom, nwind, &ii, &jj);
      dvds = one->dvds_ave;

      if (dvds <= 0.0)
        continue;

      /* same tau as p_escape uses */
      tau = sobolev (one, one->x, xplasma->density[line_ptr->nion], line_ptr, dvds);
      P = p_escape_from_tau (tau);      /* the Sobolev part of p_escape, without any HR1985 correction */
      P_sirocco = p_escape (line_ptr, xplasma);

      vth = sqrt (2. * BOLTZMANN * xplasma->t_e / (MPROT * ele[ion[line_ptr->nion].nelem].atomic_weight));
      kap_bf_c = kappa_bf (xplasma, line_ptr->freq, 0);
      kap_ff_c = kappa_ff (xplasma, line_ptr->freq);
      kap_es_c = xplasma->ne * THOMPSON * zdom[ndom].fill;

      if (tau > 0.0)
      {
        beta = (kap_bf_c + kap_ff_c) * vth / (dvds * tau);
        beta_es = kap_es_c * vth / (dvds * tau);
        gamma = 1.0 / tau;
        A_fit = hr1985_A (beta, gamma, P, FALSE);
        A_es_fit = hr1985_A (beta_es, gamma, P, FALSE);
        A_tab = hr1985_A (beta, gamma, P, TRUE);
        A_es_tab = hr1985_A (beta_es, gamma, P, TRUE);
      }
      else
      {
        beta = beta_es = gamma = -1.0;
        A_fit = A_es_fit = 0.0;
        A_tab = A_es_tab = hr1985_npts > 0 ? 0.0 : -1.0;
      }
      A_es = hr1985_npts > 0 ? A_es_tab : A_es_fit;

      ncells++;
      if (tau > 1.0)
        nthick++;
      if (tau > 0.0 && gamma < 15. * beta_es)
        ndominant++;
      if (P > 0.0 && A_es / P > 0.1)
        nten++;
      if (P > 0.0 && A_es / P > max_ratio)
        max_ratio = A_es / P;

      fprintf (fptr,
               "%9.2f %3d %3d %6d %6d %4d %4d %11.4e %11.4e %9.3e %9.3e %11.4e %11.4e %11.4e %11.4e %11.4e %11.4e %11.4e %11.4e %11.4e %11.4e %11.4e %11.4e %11.4e %11.4e %11.4e\n",
               lambda, line_ptr->z, line_ptr->istate, nplasma, nwind, ii, jj, one->xcen[0], one->xcen[2], xplasma->ne, xplasma->t_e,
               dvds, tau, vth, kap_bf_c, kap_ff_c, kap_es_c, beta, beta_es, gamma, P, A_fit, A_es_fit, A_tab, A_es_tab, P_sirocco);
    }

    printf
      ("%9.2f A (z=%d ion=%d, requested %.2f): %d cells, %d with tau>1, %d with gamma<15beta_es, %d with A_es/P>0.1, max A_es/P %.3e (%s)\n",
       lambda, line_ptr->z, line_ptr->istate, hr1985_lambda[k], ncells, nthick, ndominant, nten, max_ratio,
       hr1985_npts > 0 ? "table" : "fit");
  }

  fclose (fptr);
  return (0);
}
