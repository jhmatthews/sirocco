"""
Timeseries internal processing module
"""
from py4py.reverb import TransferFunction, open_database
import numpy as np
from numpy.random import normal
from numpy import poly1d, floating
from numpy.typing import NDArray
import astropy as ap
import astropy.constants as apc
from astropy.table import Table
from astropy import units as u
from typing import Optional


def interpolation_across_range(
        x: NDArray[floating],
        y: NDArray[floating],
        x_int: float
) -> float:
    """
    Simple linear interpolation function

    Arguments:
        x: X values
        y: Y values
        x_int: X to find Y for
    Returns:
        Linear interpolation of Y for x_int
    """

    if x_int >= x[-1]:
        return y[-1]
    elif x_int <= x[0] == 0:
        return y[0]
    else:
        x_max_index = np.searchsorted(x, x_int)
        x_min_index = x_max_index - 1
        x_int_fraction = (x_int - x[x_min_index]) / (x[x_max_index] - x[x_min_index])
        y_diff = y[x_max_index] - y[x_min_index]
        return y[x_min_index] + (y_diff * x_int_fraction)


def generate_spectrum_bounds(spectrum: Table) -> NDArray[floating]:
    """
    Given a spectrum table with 'wave_min' and 'wave_max' columns,
    returns the full list of bins.

    Arguments:
        spectrum: The table in Python output format

    Returns:
        The bounds
    """
    # Take spectrum and produce full list of wavelength bins for it to feed to the TF
    bounds = list(spectrum["wave_min"])
    bounds.append(spectrum["wave_max"][-1])
    return np.array(bounds)


def generate_tf(
        databases: dict, spectrum: Table, delay_bins: int,
        line: int, wave: float, name: str, limit: int = 999999999,
        dynamic_range: float = 2, observer: Optional[int] = None,
        velocity: Optional[bool] = False
) -> TransferFunction:
    """
    Generates the response function for a system.

    Arguments:
        databases (dict):
            Dictionary of 'min', 'mid' and 'max' data, each containing a dictionary
            with 'path' (to the file), 'continuum' (the continuum used in creation)
            and 'scale' (number of spectral cycles used)
        spectrum (Table):
            Spectrum to template the wavelength bins off of
        delay_bins (int):
            Number of bins to bin delays by
        line (int):
            Python line number to select
        dynamic_range (float):
            Truncates the time axis to include 100 - 10^(-dynamic range)`% of the photons.
            E.g. `dynamic_range=2` includes 99% of photons. Used as some models have very long tails.
        wave (float):
            Frequency of the line selected (in A)
        name (str):
            Name of the output files.
        limit (int):
            Number of photons to limit the DB query to. Set low for testing.
        observer (int):
            Which observer (if several) to use.
        velocity (bool):
            Whether to plot in velocity space or not.

    Returns:
        TransferFunction: The response-mapped transfer function.

    Outputs:
        {name}_min.eps: Transfer function plot for minimum
        {name}.eps: Transfer function plot for midpoint
        {name}_max.eps: Transfer function plot for maximum
        {name}_resp.eps: Response function for minimum
    """
    db_mid = open_database(databases['mid']['path'], "root", "password")
    db_min = open_database(databases['min']['path'], "root", "password")
    db_max = open_database(databases['max']['path'], "root", "password")

    bounds = generate_spectrum_bounds(spectrum)

    tf_mid = TransferFunction(
        db_mid, name, continuum=databases['mid']['continuum'],
        wave_bins=(len(bounds)-1), delay_bins=delay_bins
    )
    if observer:
        tf_mid.spectrum(observer)

    tf_mid.line(line, wave).wavelength_bins(bounds).delay_dynamic_range(dynamic_range).run(
        scaling_factor=databases['mid']['scale'], limit=limit, verbose=True
    ).plot(velocity=velocity)

    tf_min = TransferFunction(
        db_min, name+'_min', continuum=databases['min']['continuum'], template=tf_mid
    ).run(
        scaling_factor=databases['min']['scale'], limit=limit
    ).plot(velocity=velocity)

    tf_max = TransferFunction(
        db_max, name+'_max', continuum=databases['max']['continuum'], template=tf_mid
    ).run(
        scaling_factor=databases['max']['scale'], limit=limit
    ).plot(velocity=velocity)

    tf_mid.response_map_by_tf(tf_min, tf_max).plot(
        response_map=True, name='resp', velocity=velocity
    )
    return tf_mid


def generate_spectra_base(spectrum: Table, spectra_times: Table) -> Table:
    """
    Generates the base spectra for each timestep.

    Arguments:
        spectrum (Table): The base, unmodified spectrum used for the output time series
        spectra_times(Table): Times to produce a spectrum for

    Returns:
        Table: With one spectrum per target spectra time, keyed by the times
    """
    spectra = ap.table.Table([spectrum['wave'], spectrum['value'], spectrum['error']])
    spectra.meta['bounds'] = generate_spectrum_bounds(spectrum)
    spectra['wave'].meta['name'] = spectrum['wave'].meta['name']
    spectra['value'].meta['name'] = spectrum['value'].meta['name']
    spectra['error'].meta['name'] = spectrum['error'].meta['name']
    spectra['value_min'] = spectrum["value"].copy(copy_data=True)
    spectra['value_max'] = spectrum["value"].copy(copy_data=True)

    for time in spectra_times['time']:
        # Add a base spectrum to the output spectra file
        spectra["value {}".format(time)] = spectra['value']
        spectra["value {}".format(time)].meta['time'] = time
    return spectra


def generate_times_and_delta_continuum(
        transfer_function: TransferFunction, lightcurve: Table, delay_max: Optional[float] = None
) -> Table:
    """
    Generates the timesteps to evaluate the TF at and the change in continuum at each

    Arguments:
        transfer_function (TransferFunction): The TF
        lightcurve(Table): The lightcurve base used for this.
        delay_max (float): The maximum delay to generate the delta continuum for

    Returns:
        times (np.array): Time domain broken down into small steps
    """
    # We need to evaluate at every bin-width's step
    delay_bins = (transfer_function.delay_bins() * u.s).to(lightcurve['time'].unit)

    if delay_max:
        # If we need to rescale this to a given maximum
        delay_bins *= (delay_max / delay_bins[-1])
    bin_width = (delay_bins[1] - delay_bins[0]).value

    # We need continuum in terms of delta from the mean
    times = ap.table.Table(
        [np.arange(lightcurve['time'][0], lightcurve['time'][-1] + bin_width, bin_width)], names=['time']
    )
    times['time'].unit = lightcurve['time'].unit
    times['C'] = np.zeros(len(times))
    times['C'].unit = lightcurve['value'].unit
    times['dC'] = np.zeros(len(times))
    times['dC'].unit = lightcurve['value'].unit
    times['dC%'] = np.zeros(len(times))
    times['time'].meta['name'] = lightcurve['time'].meta['name']
    times['C'].meta['name'] = lightcurve['value'].meta['name']
    times['dC'].meta['name'] = lightcurve['value'].meta['name']
    times.meta['delay_bins'] = delay_bins

    for step in range(0, len(times)):
        # Calculate the delta continuum from the 'current value and the starting value,
        # hence the pulse contribution to later timesteps
        times['C'][step] = interpolation_across_range(lightcurve['time'], lightcurve['value'], times['time'][step])
        times['dC'][step] = (times['C'][step] - lightcurve.meta['mean'].value)
        times['dC%'][step] = (times['dC'][step] / lightcurve.meta['mean'].value)
    return times


def generate_spectra_min_max(
    times: Table, transfer_function: TransferFunction, spectra: Table, spectrum: Table,
    continuum_fit: Optional[poly1d] = None
):
    """
    When passed a timeseries of spectra, finds the outermost extent of the flux envelope (minimum and maximum values),
    and modifies the timeseries of spectra to add that as columns.

    Arguments:
        times (Table): The list of times to take spectra and the associated delta continuum for that time
        transfer_function (TransferFunction): The response function of the system
        spectra (Table): The base spectrum in Python format
        spectrum (Table): The time series of spectra
        continuum_fit(Optional[poly1d]): A function that describes the background continuum as a function of wavelength

    Returns:
         None
    """
    delay_bins = transfer_function.delay_bins()
    dC_max = np.amax(times['dC'])
    dC_min = np.amin(times['dC'])

    pulses_added = []

    # First we generate 'extreme' spectra
    print("Generating 'extreme' spectra...".format(len(delay_bins)))
    for i in range(0, len(delay_bins)-1):
        response = transfer_function.response(delay_index=i)
        pulses_added.append((dC_max, i))

        for j in range(0, len(spectrum)):
            # Matching the format of spectra.c line:
            # x *= (freq * freq * 1e-8) / (dfreq * dd * C);
            dfreq = spectrum["freq_max"].quantity[j] - spectrum["freq_min"].quantity[j]
            invwave = (spectrum["freq"].quantity[j] / apc.c).to(1/spectrum["wave"].unit)
            spectra["value_min"][j] += (dC_min * response[j] * invwave * spectrum["freq"][j] / dfreq).value
            spectra["value_max"][j] += (dC_max * response[j] * invwave * spectrum["freq"][j] / dfreq).value

    np.savetxt("pulses_added_maximum.txt", pulses_added)

    if continuum_fit:
        dC_max = np.amax(times['dC%'])
        dC_min = np.amin(times['dC%'])
        for j in range(0, len(spectrum)):
            spectra["value_min"][j] += continuum_fit(spectrum["wave"][j]) * dC_min
            spectra["value_max"][j] += continuum_fit(spectrum["wave"][j]) * dC_max


def generate_spectra_details(
    times: Table, transfer_function: TransferFunction, spectra: Table, spectrum: Table,
    continuum_fit: Optional[poly1d] = None,
    calculate_error: Optional[bool] = False,
    error_over_variation: Optional[float] = 0.01,
    verbose: Optional[bool] = True
):
    """
    Generates synthetic spectra from a base spectrum, transfer function, and varying continuum.

    Arguments:
        times (Table):
            A table containing the timesteps the time-series should be evaluated over
        transfer_function (TransferFunction):
            The TF containing the response function to use
        spectra (Table):
            The output table for the finished timeseries
        spectrum (Table):
            A table containing the spectrum to use as the base for the time-series
        continuum_fit (poly1d):
            If this is a continuuum-subtracted timeseries, the function that describes the continuum
        calculate_error (Optional[bool]):
            Whether or not we should calculate the error on the timeseries steps
        error_over_variation (Optional[float]):
            If calculating error, what should be the target integrated line error over the line variation range
        verbose (Optional[bool]):
            Whether or not to print out info messages

    Returns:
        None.
    """
    delay_bins = times.meta['delay_bins']

    # Generate prefactors
    prefactor = np.zeros(len(spectrum))
    for j in range(0, len(spectrum)):
        dfreq = spectrum["freq_max"][j] - spectrum["freq_min"][j]
        invwave = (spectrum["freq"].quantity[j] / apc.c).to(1/spectrum["wave"].unit).value
        prefactor[j] = invwave * spectrum['freq'][j] / dfreq

    # For each timestep, we send out a 'pulse' of continuum
    print("Beginning {} time steps to generate {} spectra...".format(len(times), len(spectra.columns)-5))
    for step in range(0, len(times)-1):
        # For each time bin this pulse is smeared out over
        dC_abs = times['dC'][step]
        if verbose:
            print("Step {}: {:.1f}%".format(step+1, 100.0*(step+1)/len(times)))
        for i in range(0, len(delay_bins)-1):
            # Figure out what the bounds are for the region this pulse affects
            time_range = [times['time'][step] + delay_bins[i].value, times['time'][step] + delay_bins[i+1].value]
            response = transfer_function.response(delay_index=i)
            for column in spectra.colnames[5:]:
                # For each spectrum, if it occurs at a timestamp within this bin
                if time_range[0] <= spectra[column].meta['time'] < time_range[1]:
                    # Add this pulse's contribution to it
                    for j in range(0, len(spectrum)):
                        # Matching the format of spectra.c line:
                        # x *= (freq * freq * 1e-8) / (dfreq * dd * C);
                        spectra[column][j] += dC_abs * response[j] * prefactor[j]

    if calculate_error and error_over_variation:
        # If we're assigning errors
        L_t = np.zeros(len(spectra.colnames[5:]))
        for column, step in enumerate(spectra.colnames[5:]):
            L_t[step] = np.sum(column)
        dL = np.amax(L_t) - np.amin(L_t)
        error = dL / (error_over_variation * np.sqrt(len(spectra)))
        spectrum['error'] = error
        spectra['error'] = error

    if continuum_fit:
        # If we're adding a continuum change to this
        for column in spectra.colnames[5:]:
            # For each spectrum,
            dC_rel = interpolation_across_range(x=times['time'], y=times['dC%'],
                                                x_int=spectra[column].meta['time'])
            for j in range(0, len(spectra)):
                # Add the delta continuum to it
                spectra[column][j] += continuum_fit(spectrum["wave"][j]) * dC_rel


def generate_times_line_emission(spectra: Table, spectra_times: Table, verbose: bool = False):
    """
    Given a finished timeseries of spectra, produce the total line emission at each observation

    Arguments:
        spectra (Table): The time-series of spectra
        spectra_times (Table): The times at which the spectra are taken (i.e. the fake observations)
        verbose (bool): Whether or not to report on the total line variation

    Return:
        Table: A table with the line emission for each timestep
    """
    line_times = spectra_times.copy(copy_data=True)
    line_times['time'] = line_times['time'].quantity.to(u.s)
    line_times['line'] = np.zeros(len(line_times))

    # Integrated flux error dL = SQRT(dBin1^2 + dBin2^2 + ...)
    # dL = SQRT(N_Bins * dBin^2) = SQRT(N_Bins) * dBin, we divide by SQRT(N_Bins)
    line_times['line_error'] = np.zeros(len(line_times))
    line_times['line_error'] = np.sqrt(len(spectra)) * spectra['error'][0]

    # For each spectrum in the output, sum the total emission
    # This obviously only works for line spectra!
    for step in range(0, len(line_times)):
        line_times['line'][step] = np.sum(spectra[spectra.colnames[5+step]])

    if verbose:
        print("Variation is: {}".format(np.amax(line_times['line'])-np.amin(line_times['line'])))

    #
    # spec_max = np.amax(line_times['line'])
    # spec_min = np.amin(line_times['line'])

    # error = np.array(line_times['line_error']) / (spec_max - spec_min)
    # error = (error * 9)
    # value = (np.array(line_times['line']) - spec_min) / (spec_max - spec_min)
    # value = (value * 9) + 1
    #
    # line_times['line_error'] = error
    # line_times['line'] = value

    return line_times


def generate_spectra_error(spectra: Table,
                           error: float = 0.01,
                           fudge_factor: float = 1.0):
    """
    Given a timeseries of spectra, generates random errors of magnitude dependent on the luminosity variation within it.

    Arguments:
        spectra (Table):
            A timeseries of spectra.
        error (float):
            The ratio of error on the spectrum to the variation in luminosity within it.
        fudge_factor (float):
            A scaling factor.

    Returns:
        A copy of the spectra table with errors applied.
    """
    # We want integrated flux error / continuum variation <= error
    # err_L / var_l = error
    rhs = error
    print("Generating error, aiming for Line Error/Variation = {}".format(rhs))

    # So the error on the line should be err_L * var_l
    # Change in line is calculated from the spectra:
    line_array = np.zeros(len(spectra.colnames[5:]))
    for index, col_name in enumerate(spectra.colnames[5:]):
        line_array[index] = np.sum(spectra[col_name])

    var_l = np.amax(line_array) - np.amin(line_array)

    # Multiply by var_l to get error_L
    rhs *= var_l * fudge_factor

    # Now, error on the line = sqrt( error on bin 1^2, error on bin 2^2 ... )
    # Which is sqrt(number of bins) * error on bin
    # So divide through by sqrt(number of bins) to convert to single bin
    rhs /= np.sqrt(len(spectra))

    # We now have the actual error, stick on the spectrum
    spectra['error'] = rhs
    return apply_spectra_error(spectra)


def copy_spectra_error(origin: Table, target: Table, rescale: bool = False):
    """
    Spectra error are calculated from minimum to maximum line variation.
    This doesn't really work well in situations where the variation is low due to a mixed positive-negative
    response function. So we instead copy across the errors from another timeseries, rescaling if necessary.

    Arguments:
        origin (Table):
            The timeseries of spectra with errors.
        target (Table):
            The timeseries of spectra that we want to copy the errors to
        rescale (bool):
            Optionally, whether to rescale the errors to be proportionally the same, rather than absolutely the same.

    Return:
         Copy of the 'target' timeseries with errors applied
    """
    # Copy across errors and apply them
    if not rescale:
        target['error'] = origin['error'][0]

    else:
        max_origin = np.argmax(origin['value'])
        max_target = np.argmax(target['value'])
        rescaled_error = origin['error'][0] * target['value'][max_target] / origin['value'][max_origin]
        target['error'] = rescaled_error

    return apply_spectra_error(target)


def apply_spectra_error(spectra: Table):
    """
    Given a timeseries of spectra with an 'error' column, creates a copy and
    applies random normally-distributed errors to the values at each timestep.

    Arguments:
        spectra (Table): A timeseries of spectra, with timesteps as columns 5+.

    Returns:
        Table: A copy of the input spectra with the errors applied
    """
    # Now we have the final spectra, create clean copies then add the experimental errors
    clean_copy = spectra.copy(copy_data=True)

    for column in spectra.colnames[5:]:
        # For each spectrum
        if 'value' in column:
            for j in range(0, len(spectra)):
                # For each wavelength bin in each spectrum, add a random error
                spectra[column][j] += normal(scale=spectra['error'][j])

    return clean_copy
