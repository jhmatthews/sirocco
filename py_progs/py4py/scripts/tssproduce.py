"""
Produces a set of CARAMEL and MEMECHO-format synthetic lightcurves from a set of models.
"""
# -*- coding: utf-8 -*-
import py4py.reverb.timeseries.output as tss_output
import py4py.reverb.timeseries.input as tss_import
import py4py.reverb.timeseries.process as tss_process

import astropy as ap
from astropy import units as u
# noinspection SpellCheckingInspection
from astropy.units import cds as ucds

import numpy as np

import pickle

import datetime

import sys

# ========== SETTINGS ==========
# -------- TF SETTINGS ---------
tf_lim_test = 9999999999
tf_delay_bins = 100
tf_wave = 6562.8
# ---- Seyfert Settings ----
spectrum_file_sey = "sey_100.spec"
suffix_sey = "sey"
bolometric_sey = 1.043e44
tf_line_sey = 28
databases_sey = {'min': {'path': "/home/swm1n12/bindata/sey_090", 'scale': 1.0/60, 'continuum': bolometric_sey*0.9},
                 'mid': {'path': "/home/swm1n12/bindata/sey_100", 'scale': 1.0/60, 'continuum': bolometric_sey},
                 'max': {'path': "/home/swm1n12/bindata/sey_110", 'scale': 1.0/60, 'continuum': bolometric_sey*1.1}}

# ---- QSO Settings ----
spectrum_file_qso = "qso_100.spec"
suffix_qso = "qso"
bolometric_qso = 1.043e46
tf_line_qso = 44
databases_qso = {'min': {'path': "/home/swm1n12/bindata/qso_090", 'scale': 1.0/50,  'continuum': bolometric_qso*0.9},
                 'mid': {'path': "/home/swm1n12/bindata/qso_100", 'scale': 1.0/100, 'continuum': bolometric_qso},
                 'max': {'path': "/home/swm1n12/bindata/qso_110", 'scale': 1.0/50,  'continuum': bolometric_qso*1.1}}

# ---- LIGHTCURVE SETTINGS -----
lightcurve_file = "light_1158.dat"
lightcurve_time_units = ucds.MJD
lightcurve_value_units = 1e-15 * u.erg / (u.s * u.angstrom * u.cm * u.cm)
# lightcurve_bolometric_correction = 9.0 * 5100.0 * u.angstrom
# We want a curve of the observed bolometric luminosity, so we need to rescale to 100 Pc
lightcurve_target_lum_qso = bolometric_qso * u.erg / (u.s * np.pi * np.power(u.parsec.to(u.cm)*100, 2) * u.cm * u.cm)
lightcurve_target_lum_sey = bolometric_sey * u.erg / (u.s * np.pi * np.power(u.parsec.to(u.cm)*100, 2) * u.cm * u.cm)


# ----- SPECTRUM SETTINGS ------
spectrum_bins_name = "Lambda"
spectrum_value_name = "A40P0.50"
# python spectra are per cm2 at 100 pc -> we need to divide ΔC by this value too
spectrum_value_units = u.angstrom * u.erg / (u.s * u.angstrom * (np.pi * np.power(u.parsec.to(u.cm)*100, 2) * u.cm * u.cm))
# spectrum_value_to_lightcurve_value = (np.pi * np.power(u.parsec.to(u.cm) * 100, 2)) * u.cm * u.cm / 1000
# spectrum_value_units = 1e-14 * u.erg / (u.s * u.angstrom * u.cm * u.cm)
spectrum_wave_units = u.angstrom

# ----- SPECTRUM (FULL) SETTINGS -----
# Used to create the non-continuum-subtracted spectra for MEMECHO
# As much continuum as possible! Can't go further due to the limits of Python.
spectrum_full_wave_range = [6200, 7000] * u.angstrom
spectrum_full_rebin_to = 100

# ----- SPECTRUM (SUBTRACTED) SETTINGS -----
# Used to create the continuum-subtracted spectra for CARAMEL
# As little continuum as possible! We want just the line to make it faster.
spectrum_line_wave_range = [6200, 7000] * u.angstrom
spectrum_line_subtract_range = [6300, 6850] * u.angstrom
spectrum_line_rebin_to = 70

# ------ SPECTRA TIMES ---------
spectra_times_file = "spectra_times.dat"
spectra_times_units = ucds.MJD
spectra_fudge_factor = 1

# ---- VISUALIZATION & OUTPUT SETTINGS ------
output_lightcurve_file = "out_lightcurve"
output_times_file = "out_times"
output_times_line_file = "out_times_line"
output_trailed_spec_file = "out_tspec"
output_spectra_file = "out_spectra"
output_animation_file = "out_anim"

is_reversed = False
visualise_outputs = True
visualise_animation = False
visualise_clean = False
visualise_rescaled_tfs = True
visualise_rescaled_tfs_max = 30
time_series_outputs = False

# --------- RESCALING ----------
delay_max = 90 * u.d
delta_continuum_range = 0.50

# --------- ERRORS ---------
error_ratio_to_variation = 0.020

# ------ PICKLES ------
pickle_tf_file = 'pickle_tf'
pickle_spectra_file = 'pickle_spectra'
pickle_times_file = 'pickle_times'
use_pickled_tf = True
use_pickled_times = True  # If set to true, won't change continuum range!
use_pickled_spectra = True

stop_at_min_max = False

# ===============
# Program begins!
# ===============

print("=== tssproduce started! ===")
print("Importing begins at: {}".format(datetime.datetime.now()))

# Import all data files
print("Importing spectra...")
spectrum_qso_line, continuum_fit_qso = tss_import.read_spectrum(spectrum_file_qso, spectrum_bins_name, spectrum_value_name,
                                                                frequency=False,
                                                                wave_units=spectrum_wave_units,
                                                                wave_name="Å",
                                                                value_units=spectrum_value_units,
                                                                value_name="erg s$^{-1}$ cm$^{-2}$ at 100 Pc",
                                                                limits=spectrum_line_wave_range,
                                                                subtract_continuum_with_mask=spectrum_line_subtract_range,
                                                                rebin_to=spectrum_line_rebin_to)
spectrum_qso_full = tss_import.read_spectrum(spectrum_file_qso, spectrum_bins_name, spectrum_value_name,
                                             frequency=False,
                                             wave_units=spectrum_wave_units,
                                             value_units=spectrum_value_units,
                                             value_name="erg s$^{-1}$ cm$^{-2}$ at 100 Pc",
                                             limits=spectrum_full_wave_range,
                                             rebin_to=spectrum_full_rebin_to)
spectrum_sey_line, continuum_fit_sey = tss_import.read_spectrum(spectrum_file_sey, spectrum_bins_name, spectrum_value_name,
                                                                frequency=False,
                                                                wave_units=spectrum_wave_units,
                                                                wave_name="Å",
                                                                value_units=spectrum_value_units,
                                                                value_name="erg s$^{-1}$ cm$^{-2}$ at 100 Pc",
                                                                limits=spectrum_line_wave_range,
                                                                subtract_continuum_with_mask=spectrum_line_subtract_range,
                                                                rebin_to=spectrum_line_rebin_to)
spectrum_sey_full = tss_import.read_spectrum(spectrum_file_sey, spectrum_bins_name, spectrum_value_name,
                                             frequency=False,
                                             wave_units=spectrum_wave_units,
                                             wave_name="Å",
                                             value_units=spectrum_value_units,
                                             value_name="erg s$^{-1}$ cm$^{-2}$ at 100 Pc",
                                             limits=spectrum_full_wave_range,
                                             rebin_to=spectrum_full_rebin_to)


print("Importing lightcurve file '{}'...".format(lightcurve_file))
lightcurve_qso = tss_import.read_lightcurve(lightcurve_file,
                                            time_units=lightcurve_time_units,
                                            value_units=lightcurve_value_units,
                                            time_name="MJD",
                                            value_name="erg s$^{-1}$",
                                            target_bolometric_luminosity=lightcurve_target_lum_qso,
                                            delta_continuum_range=delta_continuum_range)
lightcurve_sey = tss_import.read_lightcurve(lightcurve_file,
                                            time_units=lightcurve_time_units,
                                            value_units=lightcurve_value_units,
                                            time_name="MJD",
                                            value_name="erg s$^{-1}$",
                                            target_bolometric_luminosity=lightcurve_target_lum_sey,
                                            delta_continuum_range=delta_continuum_range)

print("Importing spectra timing file '{}'...".format(spectra_times_file))
spectra_times = tss_import.read_spectra_times(spectra_times_file,
                                              time_units=spectra_times_units,
                                              time_name="MJD")


# Produce a TF
if not use_pickled_tf:
    print("Generating Ψ begins at: {}".format(datetime.datetime.now()))
    tf_qso_full = tss_process.generate_tf(databases_qso, spectrum_qso_full, tf_delay_bins, tf_line_qso, tf_wave, 'qso_full', tf_lim_test)
    tf_qso_line = tss_process.generate_tf(databases_qso, spectrum_qso_line, tf_delay_bins, tf_line_qso, tf_wave, 'qso_line', tf_lim_test)
    tf_sey_full = tss_process.generate_tf(databases_sey, spectrum_sey_full, tf_delay_bins, tf_line_sey, tf_wave, 'sey_full', tf_lim_test)
    tf_sey_line = tss_process.generate_tf(databases_sey, spectrum_sey_line, tf_delay_bins, tf_line_sey, tf_wave, 'sey_line', tf_lim_test)
    print("Pickling Ψ for future use...")
    picklefile = open(pickle_tf_file+'_sey_full_tf.pickle', 'wb')
    pickle.dump(tf_sey_full, picklefile)
    picklefile.close()
    picklefile = open(pickle_tf_file+'_sey_line_tf.pickle', 'wb')
    pickle.dump(tf_sey_line, picklefile)
    picklefile.close()
    picklefile = open(pickle_tf_file+'_qso_full_tf.pickle', 'wb')
    pickle.dump(tf_qso_full, picklefile)
    picklefile.close()
    picklefile = open(pickle_tf_file+'_qso_line_tf.pickle', 'wb')
    pickle.dump(tf_qso_line, picklefile)
    picklefile.close()

else:
    print("Unpickling Ψ...")
    picklefile = open(pickle_tf_file+'_sey_full_tf.pickle', 'rb')
    tf_sey_full = pickle.load(picklefile)
    picklefile.close()
    picklefile = open(pickle_tf_file+'_sey_line_tf.pickle', 'rb')
    tf_sey_line = pickle.load(picklefile)
    picklefile.close()
    picklefile = open(pickle_tf_file+'_qso_full_tf.pickle', 'rb')
    tf_qso_full = pickle.load(picklefile)
    picklefile.close()
    picklefile = open(pickle_tf_file+'_qso_line_tf.pickle', 'rb')
    tf_qso_line = pickle.load(picklefile)
    picklefile.close()

# ========
# Convolve
# ========
spectra_qso_line = tss_process.generate_spectra_base(spectrum_qso_line, spectra_times)
spectra_qso_full = tss_process.generate_spectra_base(spectrum_qso_full, spectra_times)
spectra_sey_line = tss_process.generate_spectra_base(spectrum_sey_line, spectra_times)
spectra_sey_full = tss_process.generate_spectra_base(spectrum_sey_full, spectra_times)

if not use_pickled_times:
    print("Generating full lightcurve at: {}".format(datetime.datetime.now()))
    # Here we convert the simple lightcurve and series of sample times into a
    # high-resolution set of interpolated driving lightcurves and working times
    # that are at a resolution on or below that of the transfer function
    times_qso_full = tss_process.generate_times_and_delta_continuum(tf_qso_full, lightcurve_qso, delay_max)
    times_qso_line = tss_process.generate_times_and_delta_continuum(tf_qso_line, lightcurve_qso, delay_max)
    times_sey_full = tss_process.generate_times_and_delta_continuum(tf_sey_full, lightcurve_sey, delay_max)
    times_sey_line = tss_process.generate_times_and_delta_continuum(tf_sey_line, lightcurve_sey, delay_max)

    # Output the times to file
    times_sey_line.write(output_times_file+'_sey.dat', format='ascii', overwrite=True)
    times_qso_line.write(output_times_file+'_qso.dat', format='ascii', overwrite=True)

    # Then we pickle these times to disk for use later
    ap.io.misc.fnpickle(times_qso_full, pickle_times_file+'_qso_full_time.pickle')
    ap.io.misc.fnpickle(times_qso_line, pickle_times_file+'_qso_line_time.pickle')
    ap.io.misc.fnpickle(times_sey_full, pickle_times_file+'_sey_full_time.pickle')
    ap.io.misc.fnpickle(times_sey_line, pickle_times_file+'_sey_line_time.pickle')

else:
    print("Unpickling full lightcurve...")
    # Recover previously-generated times from disk for use again
    times_qso_full = ap.io.misc.fnunpickle(pickle_times_file+'_qso_full_time.pickle')
    times_qso_line = ap.io.misc.fnunpickle(pickle_times_file+'_qso_line_time.pickle')
    times_sey_full = ap.io.misc.fnunpickle(pickle_times_file+'_sey_full_time.pickle')
    times_sey_line = ap.io.misc.fnunpickle(pickle_times_file+'_sey_line_time.pickle')


# -------------
# Begin process
# -------------
if use_pickled_spectra:
    print("Unpickling time series...")
    # Recover previously-generated spectra from disk for use again
    spectra_qso_full = ap.io.misc.fnunpickle(pickle_spectra_file+'_qso_full.pickle')
    spectra_qso_line = ap.io.misc.fnunpickle(pickle_spectra_file+'_qso_line.pickle')
    spectra_sey_full = ap.io.misc.fnunpickle(pickle_spectra_file+'_sey_full.pickle')
    spectra_sey_line = ap.io.misc.fnunpickle(pickle_spectra_file+'_sey_line.pickle')
    spectra_qso_full_clean = ap.io.misc.fnunpickle(pickle_spectra_file+'_qso_full_clean.pickle')
    spectra_qso_line_clean = ap.io.misc.fnunpickle(pickle_spectra_file+'_qso_line_clean.pickle')
    spectra_sey_full_clean = ap.io.misc.fnunpickle(pickle_spectra_file+'_sey_full_clean.pickle')
    spectra_sey_line_clean = ap.io.misc.fnunpickle(pickle_spectra_file+'_sey_line_clean.pickle')
else:
    print("Generating time series begins at: {}".format(datetime.datetime.now()))
    # We generate the absolute minima and maxima possible for the spectrum given
    # the delta continuum range and the transfer functions.
    tss_process.generate_spectra_min_max(times_qso_full, tf_qso_full, spectra_qso_full, spectrum_qso_full, continuum_fit_qso)
    tss_process.generate_spectra_min_max(times_qso_line, tf_qso_line, spectra_qso_line, spectrum_qso_line)
    tss_process.generate_spectra_min_max(times_sey_full, tf_sey_full, spectra_sey_full, spectrum_sey_full, continuum_fit_sey)
    tss_process.generate_spectra_min_max(times_sey_line, tf_sey_line, spectra_sey_line, spectrum_sey_line)

    if stop_at_min_max:
        # If these are all we're interested in, just write them to disk for
        # looking at later
        spectra_qso_full.write(pickle_spectra_file+'_qso_full_range.dat', format='ascii', overwrite=True)
        spectra_qso_line.write(pickle_spectra_file+'_qso_line_range.dat', format='ascii', overwrite=True)
        spectra_sey_full.write(pickle_spectra_file+'_sey_full_range.dat', format='ascii', overwrite=True)
        spectra_sey_line.write(pickle_spectra_file+'_sey_line_range.dat', format='ascii', overwrite=True)
        sys.exit(1)

    # We now generate the spectra properly
    tss_process.generate_spectra_details(times_qso_line, tf_qso_line,
        spectra_qso_line, spectrum_qso_line, verbose=False)
    tss_process.generate_spectra_details(times_qso_full, tf_qso_full,
        spectra_qso_full, spectrum_qso_full, continuum_fit_qso, verbose=False)
    tss_process.generate_spectra_details(times_sey_line, tf_sey_line,
        spectra_sey_line, spectrum_sey_line, verbose=False)
    tss_process.generate_spectra_details(times_sey_full, tf_sey_full,
        spectra_sey_full, spectrum_sey_full, continuum_fit_sey, verbose=False)

    # Now we generate the errors! We want to set the per-pixel errors such that
    # the *total* error on the line is equal to error_ratio_to_variation
    # We calculate this for the QSO line case, as it's not polluted by the
    # continuum variation, then apply to the full QSO. Then we apply the errors
    # and keep clean copies of each.
    spectra_qso_line_clean = tss_process.generate_spectra_error(spectra=spectra_qso_line, error=error_ratio_to_variation)
    spectra_qso_full_clean = tss_process.copy_spectra_error(origin=spectra_qso_line, target=spectra_qso_full)
    spectra_sey_line_clean = tss_process.copy_spectra_error(origin=spectra_qso_line, target=spectra_sey_line, rescale=True)
    spectra_sey_full_clean = tss_process.copy_spectra_error(origin=spectra_sey_line, target=spectra_sey_full)

    # Pickle the spectra to be used later
    ap.io.misc.fnpickle(spectra_qso_full, pickle_spectra_file+'_qso_full.pickle')
    ap.io.misc.fnpickle(spectra_qso_line, pickle_spectra_file+'_qso_line.pickle')
    ap.io.misc.fnpickle(spectra_sey_full, pickle_spectra_file+'_sey_full.pickle')
    ap.io.misc.fnpickle(spectra_sey_line, pickle_spectra_file+'_sey_line.pickle')
    ap.io.misc.fnpickle(spectra_qso_full_clean, pickle_spectra_file+'_qso_full_clean.pickle')
    ap.io.misc.fnpickle(spectra_qso_line_clean, pickle_spectra_file+'_qso_line_clean.pickle')
    ap.io.misc.fnpickle(spectra_sey_full_clean, pickle_spectra_file+'_sey_full_clean.pickle')
    ap.io.misc.fnpickle(spectra_sey_line_clean, pickle_spectra_file+'_sey_line_clean.pickle')

    # Write the spectra out to an easily plotted file
    spectra_qso_full.write(output_spectra_file+'_qso_full.dat', format='ascii', overwrite=True)
    spectra_qso_line.write(output_spectra_file+'_qso_line.dat', format='ascii', overwrite=True)
    spectra_sey_full.write(output_spectra_file+'_sey_full.dat', format='ascii', overwrite=True)
    spectra_sey_line.write(output_spectra_file+'_sey_line.dat', format='ascii', overwrite=True)
    spectra_qso_full_clean.write(output_spectra_file+'_qso_full_clean.dat', format='ascii', overwrite=True)
    spectra_qso_line_clean.write(output_spectra_file+'_qso_line_clean.dat', format='ascii', overwrite=True)
    spectra_sey_full_clean.write(output_spectra_file+'_sey_full_clean.dat', format='ascii', overwrite=True)
    spectra_sey_line_clean.write(output_spectra_file+'_sey_line_clean.dat', format='ascii', overwrite=True)

    spectra_times_sey = tss_process.generate_times_line_emission(spectra_sey_line, spectra_times)
    spectra_times_qso = tss_process.generate_times_line_emission(spectra_qso_line, spectra_times)
    spectra_times_sey.write(output_times_line_file+'_sey.dat', format='ascii', overwrite=True)
    spectra_times_qso.write(output_times_line_file+'_qso.dat', format='ascii', overwrite=True)

if time_series_outputs:
    tss_output.write_caramel_data(lightcurve_qso, spectra_qso_line, spectra_times, suffix_qso)
    tss_output.write_memecho_data(lightcurve_qso, spectra_qso_full, spectra_times, suffix_qso)
    tss_output.write_caramel_data(lightcurve_sey, spectra_sey_line, spectra_times, suffix_sey)
    tss_output.write_memecho_data(lightcurve_sey, spectra_sey_full, spectra_times, suffix_sey)

# ==============================================================================
# OUTPUT VISUALIZATIONS
# ==============================================================================
# Generate trailed spectrogram
# ------------------------------------------------------------------------------
if visualise_outputs:
    print("Generating trailed spectrogram begins at: {}".format(datetime.datetime.now()))
    tss_output.trailed_spectrogram(spectra_qso_line, lightcurve_qso, spectra_times, output_trailed_spec_file+'_'+suffix_qso+'_line')
    tss_output.trailed_spectrogram(spectra_qso_full, lightcurve_qso, spectra_times, output_trailed_spec_file+'_'+suffix_qso+'_full')
    tss_output.trailed_spectrogram(spectra_sey_line, lightcurve_sey, spectra_times, output_trailed_spec_file+'_'+suffix_sey+'_line')
    tss_output.trailed_spectrogram(spectra_sey_full, lightcurve_sey, spectra_times, output_trailed_spec_file+'_'+suffix_sey+'_full')

    # tss_output.plot_spectra_rms([spectra_qso_line_clean, spectra_sey_line_clean],
    #                             ['out_specrms_qso', 'out_specrms_sey'])

    if visualise_clean:
        tss_output.trailed_spectrogram(spectra_qso_line_clean, lightcurve_qso, spectra_times, output_trailed_spec_file+'_'+suffix_qso+'_line_clean')
        tss_output.trailed_spectrogram(spectra_qso_full_clean, lightcurve_qso, spectra_times, output_trailed_spec_file+'_'+suffix_qso+'_full_clean')
        tss_output.trailed_spectrogram(spectra_sey_line_clean, lightcurve_sey, spectra_times, output_trailed_spec_file+'_'+suffix_sey+'_line_clean')
        tss_output.trailed_spectrogram(spectra_sey_full_clean, lightcurve_sey, spectra_times, output_trailed_spec_file+'_'+suffix_sey+'_full_clean')

if visualise_rescaled_tfs:
    # tf_qso_resc = tss_process.generate_tf(databases_qso, spectrum_qso_line, 25, tf_line_qso, tf_wave, 'out_qso_resc', 10000, dynamic_range=1.5)
    # tf_sey_resc = tss_process.generate_tf(databases_sey, spectrum_sey_line, 25, tf_line_sey, tf_wave, 'out_sey_resc', 10000, dynamic_range=1.5)

    tss_output.rescaled_rfs([tf_qso_line, tf_sey_line],
                            rescale_max_time=delay_max.value,
                            figure_max_time=delay_max.value/2.5,
                            keplerian={
                                'angle': 40, 'mass': 1.33e8, 'radius': [50, 100]
                            })

# ------------------------------------------------------------------------------
# Generate write_animation
# ------------------------------------------------------------------------------
if visualise_animation:
    print("Generating write_animation begins at: {}".format(datetime.datetime.now()))
    tss_output.write_animation(spectra_qso_full, lightcurve_qso, spectra_times, times_qso_full, output_animation_file + "_qso_full")
    tss_output.write_animation(spectra_qso_line, lightcurve_qso, spectra_times, times_qso_line, output_animation_file + "_qso_line")
    tss_output.write_animation(spectra_sey_full, lightcurve_sey, spectra_times, times_sey_full, output_animation_file + "_sey_full")
    tss_output.write_animation(spectra_sey_line, lightcurve_sey, spectra_times, times_sey_line, output_animation_file + "_sey_line")
    if visualise_clean:
        tss_output.write_animation(spectra_qso_full_clean, lightcurve_qso, spectra_times, times_qso_full, output_animation_file + "_qso_full_clean")
        tss_output.write_animation(spectra_qso_line_clean, lightcurve_qso, spectra_times, times_qso_line, output_animation_file + "_qso_line_clean")
        tss_output.write_animation(spectra_sey_full_clean, lightcurve_sey, spectra_times, times_sey_full, output_animation_file + "_sey_full_clean")
        tss_output.write_animation(spectra_sey_line_clean, lightcurve_sey, spectra_times, times_sey_line, output_animation_file + "_sey_line_clean")
