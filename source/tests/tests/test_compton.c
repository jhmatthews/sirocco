/** ********************************************************************************************************************
 *
 *  @file test_compton.c
 *  @author Edward J. Parkinson (e.parkinson@soton.ac.uk)
 *  @date August 2023
 *
 *  @brief
 *
 * ****************************************************************************************************************** */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <CUnit/CUnit.h>

#include "../../atomic.h"
#include "../../sirocco.h"
#include "../assert.h"

/** *******************************************************************************************************************
 *
 * @brief Test the Compton formula, which calculates the fractional energy change
 *
 * @details
 *
 * In this test, a (photon) frequency of 5e16 Hz is used. Testing this function is a bit trickier, as we need access
 * to some (static) global variables.
 *
 * ****************************************************************************************************************** */

void
test_compton_func (void)
{
  const double test_frequency = 5e18;
  const double cross_section_test = 0.5;
  const double energy_ratio = (PLANCK * test_frequency) / (MELEC * VLIGHT * VLIGHT);
  const double f_min = 1.0;
  const double f_max = 1 + (2 * energy_ratio);
  const double mid_f_max = 0.5 * (f_min + f_max);
  const double three_quarters_f_max = 0.75 * (f_min + f_max);
  const double cross_section_max = sigma_compton_partial (f_max, energy_ratio);
  set_comp_func_values (cross_section_test, cross_section_max, energy_ratio);

  CU_ASSERT_DOUBLE_EQUAL_FATAL (compton_func (f_min, NULL), -cross_section_test, EPSILON);
  CU_ASSERT_DOUBLE_EQUAL_FATAL (compton_func (mid_f_max, NULL), 2.130244e-02, EPSILON);
  CU_ASSERT_DOUBLE_EQUAL_FATAL (compton_func (three_quarters_f_max, NULL), 1.457703e+02, EPSILON);
}

/** *******************************************************************************************************************
 *
 * @brief Test the Compton cooling cross-section formula
 *
 * @details
 *
 * This contains two tests. One test where the frequency is below a threshold, where the cross section is equal to 1.0.
 * The other test has a cross-section < 1.
 *
 * ****************************************************************************************************************** */

void
test_compton_beta (void)
{
  CU_ASSERT_DOUBLE_EQUAL_FATAL (compton_beta (1e16), 1.0, EPSILON);
  CU_ASSERT_DOUBLE_EQUAL_FATAL (compton_beta (1e17), 0.9991066689339109, EPSILON);
}

/** *******************************************************************************************************************
 *
 * @brief Test the Compton heating cross-section formula
 *
 * @details
 *
 * This contains two tests. One test where the frequency is below a threshold, where the cross section is equal to 1.0.
 * The other test has a cross-section < 1.
 *
 * ****************************************************************************************************************** */

void
test_compton_alpha (void)
{
  CU_ASSERT_DOUBLE_EQUAL_FATAL (compton_alpha (1e16), 1.0, EPSILON);
  CU_ASSERT_DOUBLE_EQUAL_FATAL (compton_alpha (1e17), 0.9964273280931072, EPSILON);
}

/** *******************************************************************************************************************
 *
 * @brief Test the Klein-Nishina cross-section formula
 *
 * @details
 *
 * Four tests are performed. Two tests for frequencies below the threshold, which should return the Thompson scattering
 * cross-section. The remaining two tests will be corrected by the KN formula.
 *
 * ****************************************************************************************************************** */

void
test_klein_nishina (void)
{
  /* These frequencies should be below the KM threshold and be equal to the
     Thompson cross-section */
  CU_ASSERT_DOUBLE_EQUAL_FATAL (klein_nishina (1e14), THOMPSON, EPSILON);
  CU_ASSERT_DOUBLE_EQUAL_FATAL (klein_nishina (1e15), THOMPSON, EPSILON);

  /* These frequency values will need the KN correction  --
     I think Epsilon needs to be a bit smaller in this case... */
  CU_ASSERT_DOUBLE_EQUAL_FATAL (klein_nishina (2e16), 6.650136664443343e-25, 1e-30);
  CU_ASSERT_DOUBLE_EQUAL_FATAL (klein_nishina (1e18), 6.546940139261951e-25, 1e-30);
}

/** *******************************************************************************************************************
 *
 * @brief Monte Carlo statistics from repeatedly Compton scattering a photon off thermal electrons
 *
 * @details
 *
 * A single-cell wind is set up with electron temperature t_e, and a photon of frequency freq (travelling along z,
 * in the local frame) is scattered nscat times with compton_scatter(). For each scatter we record
 *
 *   recoil = (f - 1) / f, where f is the electron rest-frame energy ratio returned by compton_scatter(), i.e. the
 *            probability of creating a k-packet in macro-atom mode
 *   gain   = w_after / w_before - 1, the fractional change in packet weight (energy)
 *   surv   = gain / f, the weight gain kept on average by packets that survive the k-packet choice in macro-atom
 *            mode (a packet survives with probability 1 / f), i.e. the energy the k-packet sink has to remove
 *
 * and return their means and standard errors. The globals that are modified are restored before returning.
 *
 * ****************************************************************************************************************** */

struct compton_mc_stats
{
  double recoil_mean, recoil_err;
  double gain_mean, gain_err;
  double surv_mean, surv_err;
};

static struct compton_mc_stats
run_compton_scatters (const double t_e, const double freq, const int rt_mode, const int nscat)
{
  struct compton_mc_stats stats;
  struct photon p;
  WindPtr wmain_save = wmain;
  PlasmaPtr plasmamain_save = plasmamain;
  const int rt_mode_save = geo.rt_mode;
  const int rel_mode_save = rel_mode;
  double f, recoil, gain;
  double sum_recoil = 0.0, sum_recoil2 = 0.0, sum_gain = 0.0, sum_gain2 = 0.0, sum_surv = 0.0, sum_surv2 = 0.0;
  int i;

  wmain = calloc (1, sizeof (wind_dummy));
  plasmamain = calloc (1, sizeof (plasma_dummy));
  wmain[0].nplasma = 0;
  plasmamain[0].nwind = 0;
  plasmamain[0].t_e = t_e;
  geo.rt_mode = rt_mode;
  rel_mode = REL_MODE_FULL;     /* the weight only changes under a Lorentz transform in full relativistic mode */

  for (i = 0; i < nscat; i++)
  {
    memset (&p, 0, sizeof (p));
    p.freq = freq;
    p.w = 1.0;
    p.lmn[2] = 1.0;
    p.grid = 0;
    p.frame = F_LOCAL;

    f = compton_scatter (&p);
    recoil = (f - 1.0) / f;
    gain = p.w - 1.0;

    sum_recoil += recoil;
    sum_recoil2 += recoil * recoil;
    sum_gain += gain;
    sum_gain2 += gain * gain;
    sum_surv += gain / f;
    sum_surv2 += (gain / f) * (gain / f);
  }

  stats.recoil_mean = sum_recoil / nscat;
  stats.recoil_err = sqrt (fmax (sum_recoil2 / nscat - stats.recoil_mean * stats.recoil_mean, 0.0) / nscat);
  stats.gain_mean = sum_gain / nscat;
  stats.gain_err = sqrt (fmax (sum_gain2 / nscat - stats.gain_mean * stats.gain_mean, 0.0) / nscat);
  stats.surv_mean = sum_surv / nscat;
  stats.surv_err = sqrt (fmax (sum_surv2 / nscat - stats.surv_mean * stats.surv_mean, 0.0) / nscat);

  free (wmain);
  free (plasmamain);
  wmain = wmain_save;
  plasmamain = plasmamain_save;
  geo.rt_mode = rt_mode_save;
  rel_mode = rel_mode_save;

  return stats;
}

/** *******************************************************************************************************************
 *
 * @brief The exact Klein-Nishina mean fractional energy loss per scatter, <1 - E'/E>, for an electron at rest
 *
 * @details
 *
 * Integrates over mu = cos(theta) with the KN differential cross section,
 * dsigma/dmu ~ r^2 (r + 1/r - 1 + mu^2), where r = E'/E = 1 / (1 + x (1 - mu)) and x = h nu / m_e c^2.
 * In the Thomson limit this tends to x (1 - 1.4 x).
 *
 * ****************************************************************************************************************** */

static double
kn_mean_fractional_loss (const double x)
{
  const int n = 20001;
  double mu, r, weight, factor;
  double sum_loss = 0.0, sum_weight = 0.0;
  int i;

  for (i = 0; i < n; i++)
  {
    mu = -1.0 + 2.0 * i / (n - 1);
    r = 1.0 / (1.0 + x * (1.0 - mu));
    weight = r * r * (r + 1.0 / r - 1.0 + mu * mu);
    factor = (i == 0 || i == n - 1) ? 0.5 : 1.0;        /* trapezium rule */
    sum_loss += factor * weight * (1.0 - r);
    sum_weight += factor * weight;
  }

  return sum_loss / sum_weight;
}

/** *******************************************************************************************************************
 *
 * @brief The expected mean fractional energy gain from thermal Doppler shifts, (4/3) <gamma^2 beta^2>
 *
 * @details
 *
 * The average is over the electron speed distribution sampled by compton_get_thermal_velocity(), i.e.
 * p(u) ~ u^2 exp(-u^2) for 0 < u < 5 with v = u sqrt(2 k T / m_e), capped at 0.5 c. To leading order this is
 * 4 theta + 20 theta^2, where theta = k T / m_e c^2. It applies in the Thomson limit (no recoil).
 *
 * ****************************************************************************************************************** */

static double
thermal_doppler_mean_gain (const double t_e)
{
  const int n = 20001;
  const double v_scale = sqrt (2.0 * BOLTZMANN * t_e / MELEC);
  double u, beta, weight, factor;
  double sum_gain = 0.0, sum_weight = 0.0;
  int i;

  for (i = 0; i < n; i++)
  {
    u = 5.0 * i / (n - 1);
    beta = fmin (u * v_scale / VLIGHT, 0.5);
    weight = u * u * exp (-u * u);
    factor = (i == 0 || i == n - 1) ? 0.5 : 1.0;
    sum_gain += factor * weight * (4.0 / 3.0) * beta * beta / (1.0 - beta * beta);
    sum_weight += factor * weight;
  }

  return sum_gain / sum_weight;
}

/** *******************************************************************************************************************
 *
 * @brief Test that the electron rest-frame recoil from compton_scatter() matches the Klein-Nishina mean energy loss
 *
 * @details
 *
 * The electrons are effectively at rest (t_e = 1 K), so (f - 1) / f is the fractional energy lost to recoil, and
 * <(f - 1) / f> should equal the exact KN mean fractional loss. This is the probability of creating a k-packet in
 * macro-atom mode, so it is also what sets the Compton heating delivered to k-packets. Frequencies cover the
 * dipole branch (x < 1e-4), either side of the switch to the KN branch, and the KN regime.
 *
 * ****************************************************************************************************************** */

void
test_compton_scatter_recoil (void)
{
  const double x_values[] = { 1e-5, 0.9e-4, 1.1e-4, 1e-3, 1e-2, 0.1, 0.5 };
  const int n_x = sizeof (x_values) / sizeof (x_values[0]);
  const int nscat = 50000;
  struct compton_mc_stats stats;
  double x, freq, expected;
  int i;

  for (i = 0; i < n_x; i++)
  {
    x = x_values[i];
    freq = x * MELEC * VLIGHT * VLIGHT / PLANCK;
    stats = run_compton_scatters (1.0, freq, RT_MODE_MACRO, nscat);
    expected = kn_mean_fractional_loss (x);
    printf ("\n  recoil: x = %8.2e  <(f-1)/f> = %10.4e +/- %8.2e  KN = %10.4e  ratio = %7.4f", x,
            stats.recoil_mean, stats.recoil_err, expected, stats.recoil_mean / expected);
    CU_ASSERT_DOUBLE_EQUAL (stats.recoil_mean, expected, 5.0 * stats.recoil_err + 0.005 * expected);
  }
  printf ("\n");
}

/** *******************************************************************************************************************
 *
 * @brief Test that in macro-atom mode the thermal Doppler shifts give a mean weight gain of (4/3) <gamma^2 beta^2>
 *
 * @details
 *
 * In macro-atom mode compton_dir() does not reduce the packet weight, so the only weight change comes from the two
 * Lorentz transforms into and out of the electron rest frame. The photon energy is low (x ~ 1e-8) so recoil is
 * negligible, and the mean weight gain should be (4/3) <gamma^2 beta^2> ~ 4 theta: the Compton cooling rate per
 * scatter that the k-packet sink (cooling_compton) is meant to balance.
 *
 * ****************************************************************************************************************** */

void
test_compton_scatter_doppler_macro (void)
{
  const double t_values[] = { 1e6, 1e7, 1e8 };
  const int n_t = sizeof (t_values) / sizeof (t_values[0]);
  const int nscat = 500000;
  const double freq = 1e12;
  struct compton_mc_stats stats;
  double expected, theta;
  int i;

  for (i = 0; i < n_t; i++)
  {
    stats = run_compton_scatters (t_values[i], freq, RT_MODE_MACRO, nscat);
    expected = thermal_doppler_mean_gain (t_values[i]);
    theta = BOLTZMANN * t_values[i] / (MELEC * VLIGHT * VLIGHT);
    printf ("\n  doppler (macro): T = %8.2e  <w'/w - 1> = %10.4e +/- %8.2e  (4/3)<g^2b^2> = %10.4e  "
            "4 theta = %10.4e  ratio = %7.4f", t_values[i], stats.gain_mean, stats.gain_err, expected, 4.0 * theta,
            stats.gain_mean / expected);
    CU_ASSERT_DOUBLE_EQUAL (stats.gain_mean, expected, 5.0 * stats.gain_err + 0.01 * expected);
  }
  printf ("\n");
}

/** *******************************************************************************************************************
 *
 * @brief Test that in two-level mode the mean fractional energy change per scatter is ~ 4 theta - x
 *
 * @details
 *
 * In RT_MODE_2LEVEL compton_dir() also divides the weight by f, so the packet loses the recoil energy and the mean
 * fractional change is the classic thermal Comptonisation result, (4/3) <gamma^2 beta^2> - <1 - E'/E>.
 * The cross terms are O(x theta), which the tolerance allows for.
 *
 * ****************************************************************************************************************** */

void
test_compton_scatter_2level (void)
{
  const double t_e = 1e7;
  const double x = 1e-3;
  const int nscat = 500000;
  struct compton_mc_stats stats;
  double expected, freq;

  freq = x * MELEC * VLIGHT * VLIGHT / PLANCK;
  stats = run_compton_scatters (t_e, freq, RT_MODE_2LEVEL, nscat);
  expected = thermal_doppler_mean_gain (t_e) - kn_mean_fractional_loss (x);
  printf ("\n  2level: T = %8.2e x = %8.2e  <w'/w - 1> = %10.4e +/- %8.2e  expected = %10.4e  ratio = %7.4f\n", t_e, x,
          stats.gain_mean, stats.gain_err, expected, stats.gain_mean / expected);
  CU_ASSERT_DOUBLE_EQUAL (stats.gain_mean, expected, 5.0 * stats.gain_err + 0.02 * expected);
}

/** *******************************************************************************************************************
 *
 * @brief Test that the thermal part of the net Compton energy exchange matches the Hazy/Guilbert alpha * beta
 *
 * @details
 *
 * Hazy (eq 6.5-6.7, fits to Guilbert 1986) gives the net photon energy change per unit path length as
 * sigma_T n_e (4 theta alpha beta - x alpha). We measure the temperature-dependent part of the net exchange from
 * two-level mode scatters off hot and cold electrons,
 *
 *   C_net = (sigma_KN / sigma_T) [<w'/w - 1>(T) - <w'/w - 1>(cold)] / 4 theta,
 *
 * and compare with alpha beta, multiplied by (4/3) <gamma^2 beta^2> / 4 theta to allow for the relativistic thermal
 * correction (about 4% at 5e7 K) that the non-relativistic Hazy expression leaves out. This is the coefficient that
 * total_comp() uses for the Compton cooling rate.
 *
 * ****************************************************************************************************************** */

void
test_compton_net_exchange_vs_hazy (void)
{
  const double t_e = 5e7;
  const double x_values[] = { 1e-2, 3e-2, 0.1 };
  const int n_x = sizeof (x_values) / sizeof (x_values[0]);
  const int nscat = 500000;
  struct compton_mc_stats cold, hot;
  double theta4, rel, x, freq, sk, c_net, c_net_err, expected;
  int i;

  theta4 = 4.0 * BOLTZMANN * t_e / (MELEC * VLIGHT * VLIGHT);
  rel = thermal_doppler_mean_gain (t_e) / theta4;

  for (i = 0; i < n_x; i++)
  {
    x = x_values[i];
    freq = x * MELEC * VLIGHT * VLIGHT / PLANCK;
    sk = klein_nishina (freq) / THOMPSON;
    cold = run_compton_scatters (1.0, freq, RT_MODE_2LEVEL, nscat);
    hot = run_compton_scatters (t_e, freq, RT_MODE_2LEVEL, nscat);
    c_net = sk * (hot.gain_mean - cold.gain_mean) / theta4;
    c_net_err = sk * sqrt (hot.gain_err * hot.gain_err + cold.gain_err * cold.gain_err) / theta4;
    expected = compton_alpha (freq) * compton_beta (freq) * rel;
    printf ("\n  net exchange: T = %8.2e x = %6.3f  C_net = %7.4f +/- %6.4f  alpha*beta*rel = %7.4f  ratio = %7.4f", t_e, x,
            c_net, c_net_err, expected, c_net / expected);
    CU_ASSERT_DOUBLE_EQUAL (c_net, expected, 5.0 * c_net_err + 0.02 * expected);
  }
  printf ("\n");
}

/** *******************************************************************************************************************
 *
 * @brief Test that the Compton cooling rate (k-packet sink) balances the Doppler gain kept by surviving packets
 *
 * @details
 *
 * In macro-atom mode a packet survives an electron scatter with probability 1 / f and keeps the Doppler weight gain,
 * so the energy created per unit path length is sigma_KN n_e <(w'/w - 1) / f>. This is balanced by the k-packet sink,
 * whose rate is the Compton cooling rate from total_comp(), sigma_T n_e alpha beta 4 theta. The surviving gain exceeds
 * the net thermal coefficient (see above) by the growth of the electron rest-frame recoil with temperature, which goes
 * into k-packets instead, so we expect
 *
 *   C_surv / (alpha beta rel) = 1 + a few percent,
 *
 * where C_surv = (sigma_KN / sigma_T) <(w'/w - 1) / f> / 4 theta. Before total_comp() included alpha this ratio was
 * about 0.75 at x = 0.1.
 *
 * ****************************************************************************************************************** */

void
test_compton_survivor_gain_vs_cooling (void)
{
  const double t_e = 5e7;
  const double x_values[] = { 1e-2, 3e-2, 0.1, 0.3 };
  const int n_x = sizeof (x_values) / sizeof (x_values[0]);
  const int nscat = 500000;
  struct compton_mc_stats stats;
  double theta4, rel, x, freq, sk, c_surv, c_surv_err, cooling, ratio, ratio_err;
  int i;

  theta4 = 4.0 * BOLTZMANN * t_e / (MELEC * VLIGHT * VLIGHT);
  rel = thermal_doppler_mean_gain (t_e) / theta4;

  for (i = 0; i < n_x; i++)
  {
    x = x_values[i];
    freq = x * MELEC * VLIGHT * VLIGHT / PLANCK;
    sk = klein_nishina (freq) / THOMPSON;
    stats = run_compton_scatters (t_e, freq, RT_MODE_MACRO, nscat);
    c_surv = sk * stats.surv_mean / theta4;
    c_surv_err = sk * stats.surv_err / theta4;
    cooling = compton_alpha (freq) * compton_beta (freq) * rel;
    ratio = c_surv / cooling;
    ratio_err = c_surv_err / cooling;
    printf ("\n  survivor gain: T = %8.2e x = %6.3f  C_surv = %7.4f +/- %6.4f  alpha*beta*rel = %7.4f  ratio = %7.4f", t_e, x,
            c_surv, c_surv_err, cooling, ratio);
    CU_ASSERT (ratio > 0.98 - 5.0 * ratio_err && ratio < 1.10 + 5.0 * ratio_err);
  }
  printf ("\n");
}

/** *******************************************************************************************************************
 *
 * @brief Create a CUnit test suite for Compton processes
 *
 * @details
 *
 * This function will create a test suite for functions related to Compton processes, e.g. scattering or thermal
 * velocity distribution.
 *
 * ****************************************************************************************************************** */

void
create_compton_test_suite (void)
{
  CU_pSuite suite = CU_add_suite ("Compton Processes", NULL, NULL);

  if (suite == NULL)
  {
    fprintf (stderr, "Failed to create `Compton Processes` suite\n");
    CU_cleanup_registry ();
    exit (CU_get_error ());
  }

  /* Add CPU tests to suite */
  if ((CU_add_test (suite, "Klein-Nisina Formula", test_klein_nishina) == NULL) ||
      (CU_add_test (suite, "Compton Alpha - heating cross section ", test_compton_alpha) == NULL) ||
      (CU_add_test (suite, "Compton Beta - cooling cross section", test_compton_beta) == NULL) ||
      (CU_add_test (suite, "Compton Formula", test_compton_func) == NULL) ||
      (CU_add_test (suite, "Compton Scatter - recoil vs Klein-Nishina", test_compton_scatter_recoil) == NULL) ||
      (CU_add_test (suite, "Compton Scatter - Doppler gain (macro)", test_compton_scatter_doppler_macro) == NULL) ||
      (CU_add_test (suite, "Compton Scatter - energy change (2level)", test_compton_scatter_2level) == NULL) ||
      (CU_add_test (suite, "Compton Scatter - net exchange vs Hazy alpha*beta", test_compton_net_exchange_vs_hazy) == NULL) ||
      (CU_add_test (suite, "Compton Scatter - survivor gain vs cooling rate", test_compton_survivor_gain_vs_cooling) == NULL))
  {
    fprintf (stderr, "Failed to add tests to `Compton Processes` suite\n");
    CU_cleanup_registry ();
    exit (CU_get_error ());
  }
}
