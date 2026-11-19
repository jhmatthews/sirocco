"""
Contains the class used to create and manipulate reverberation maps from Python output files.

Example:

    For an existing delay output file called 'qso.delay_dump', to generate a TF plot for the C4
    line for a specific spectrum, with axis of velocity offset vs days, you would do::

        qso_conn = open_database('qso')
        tf_c4_1 = TransferFunction(
            qso_conn, continuum=1e43, wave_bins=100, delay_bins=100, filename='qso_c4_spectrum_1'
        )
        tf_c4_1.spectrum(1).line(443).run()
        tf_c4_1.plot(velocity=True, days=True)

    Given database queries can take a long time, it is advisable to pickle a TF that has been run
    so you can access it later on. Note, however: Once a TF has been restored from a pickle,
    you can no longer change the filters and re-run::

        with open('qso_c4_spectrum_1', 'wb') as file:
            pickle.dump(tf_c4_1, file)


"""
# -*- coding: utf-8 -*-
# pylint: disable=C0301
import numpy as np
import astropy.constants as apc  # pylint: disable=E1101
import time
import sys
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
import sqlalchemy
import sqlalchemy.ext.declarative
import sqlalchemy.exc
import sqlalchemy.orm
import sqlalchemy.orm.query
from sqlalchemy.engine import Engine
import matplotlib
from astropy import units as u
from typing import List, Optional, Union, Tuple

from py4py.array import calculate_fwhm, calculate_midpoints
from py4py.physics import keplerian_velocity, doppler_shift_wave, doppler_shift_vel


SECONDS_PER_DAY = 60 * 60 * 24
"""Constant used for rescaling data, that is probably superfluous and already present in Astropy"""


# ==============================================================================
# PHYSICS FUNCTIONS
# ==============================================================================
def calculate_delay(angle: float, phase: float, radius: float, days: bool = True) -> float:
    """
    Delay relative to continuum for emission from a point on the disk.

    Calculate delay for emission from a point on a keplerian disk, defined by
    its radius and disk angle, to an observer at a specified angle.

    Draw plane at r_rad_min out. Find x projection of disk position.
    Calculate distance travelled to that plane from the current disk position
    Delay relative to continuum is thus (distance from centre to plane)
    + distance from centre to point

    Args:
        angle (float):  Observer angle to disk normal, in radians
        phase (float):  Rotational angle of point on disk, in radians. 0 = in line to observer
        radius (float): Radius of the point on the disk, in m
        days (bool):    Whether the timescale should be seconds or days

    Returns:
        float:          Delay relative to continuum
    """
    vr_disk = np.array([radius*np.cos(phase), 0.0])
    vr_normal = np.array([np.sin(angle), np.cos(angle)])
    vr_plane = radius * vr_normal
    delay = np.dot((vr_plane - vr_disk), vr_normal) / apc.c.value

    if days:
        return delay / SECONDS_PER_DAY
    else:
        return delay


# ==============================================================================
# TRANSFER FUNCTION DEFINITION
# ==============================================================================
class TransferFunction:
    """
    Used to create, store and query emissivity and response functions
    """
    def __getstate__(self) -> dict:
        """
        Removes invalid data before saving to disk.

        Returns:
            dict: Updated internal dict, with references to external,
                  session-specific database things, removed.
        """
        state: dict = self.__dict__.copy()
        # Remove the unpicklable entries.
        del state['_session']
        del state['_query']
        del state['_database']
        return state

    def __setstate__(self, state: dict):
        """
        Restores the data from disk, and sets a flag to show this is a frozen TF.

        Args:
            state (dict): The unpickled object dict..
        """
        self.__dict__.update(state)
        self._unpickled = True


    def __init__(
                    self, database: sqlalchemy.engine.Connection, filename: str,
                    continuum: float, wave_bins: int = None, delay_bins: int = None,
                    template: 'TransferFunction' = None,
                    template_different_line: bool = False, template_different_spectrum: bool = False
            ):
        """
        Initialises the TF, optionally by templating off another TF.

        Sets up all the basic properties of the TF that are required to create
        it. It must be `.run()` to query the DB before it can itself be queried.
        If templating, it applies all the same filters that were applied to the
        template TF, unless explicitly told not to. Filters don't overwrite!
        They stack. So you can't simply call `.line()` to change the line the TF
        corresponds to if its template was a different line, unless you specify
        that the template was of a different line.

        Arguments:
            database (sqlalchemy.engine.Connection):
                The database to be queried for this TF.
            filename (string):
                The root filename for plots created for this TF.
            continuum (float):
                The continuum value associated with this TF. Central source + disk luminosity.
            wave_bins (int):
                Number of wavelength/velocity bins.
            delay_bins (int):
                Number of delay time bins.
            template (TransferFunction):
                Other TF to copy all filter settings from.
                Will match delay, wave and velocity bins exactly.
            template_different_line (bool):
                Is this TF going to share delay & velocity bins but have different wavelength bins?
            template_different_spectrum (bool):
                Is this TF going to share all specified bins but
                be taken on photons from a different observer.

        Todo:
            Consider making it impossible to apply filters after calling run().
        """
        assert (delay_bins is not None and wave_bins is not None) or template is not None,\
            "Must provide either resolutions or another TF to copy them from!"
        # self._query = database.query(Photon.Wavelength, Photon.Delay, Photon.Weight, Photon.X, Photon.Z)

        self._database = database
        session_maker = sqlalchemy.orm.sessionmaker(bind=self._database)
        self._session = session_maker()
        self._query = self._session.query(Photon.Wavelength, Photon.Delay, Photon.Weight)

        self._delay_dynamic_range = None
        self._velocity = None
        self._line_list = None
        self._line_wave = None
        self._line_num = None
        self._delay_range = None
        self._continuum = continuum
        self._filename = filename
        self._bins_wave_count = wave_bins
        self._bins_delay_count = delay_bins
        self._bins_vel = None
        self._bins_wave = None
        self._bins_delay = None
        self._emissivity = None
        self._response = None
        self._count = None
        self._wave_range = None
        self._spectrum = None
        self._unpickled = False

        if template is not None:
            # If we're templating off a pre-existing transfer function, copy over all the shared properties
            print("Templating '{}' off of '{}'...".format(self._filename, template._filename))
            # Regardless of what line we're templating off, we want to share the velocity and delay bins
            self._bins_wave_count = template._bins_wave_count
            self._bins_delay_count = template._bins_delay_count
            self._bins_vel = template._bins_vel
            self._bins_delay = template._bins_delay

            # Now we want to call all the same filter functions that've been applied to the template
            # (where appropriate)
            if template_different_line is False:
                # If we're templating off of the same line, we want the same wavelength bins
                self.wavelength_bins(template._bins_wave)
            if template._line_wave is not None and template_different_line is False:
                # If we're templating off the same line, record we're using that line
                self.line(template._line_num, template._line_wave)
            if template._velocity is not None:
                # If we're templating off a TF with velocity, record we're doing so
                self.velocities(template._velocity)
            if template._line_list is not None and template_different_line is False:
                # If we want the same bins for the same list of lines, record so
                self.lines(template._line_list)
            if template._spectrum is not None and template_different_spectrum is False:
                # If we want the same bins for the same spectrum, record so
                self.spectrum(template._spectrum)

    def spectrum(self, number: int) -> 'TransferFunction':
        """
        Constrain the TF to photons from a specific observer

        Args:
            number (int): Observer number from Python run

        Returns:
            TransferFunction: Self, so filters can be stacked
        """
        assert self._unpickled is False,\
            "TF restored from pickle! It cannot be rerun, filters cannot be applied."
        self._spectrum = number
        self._query = self._query.filter(Photon.Spectrum == number)
        return self

    def line(self, number: int, wavelength: float) -> 'TransferFunction':
        """
        Constrain the TF to only photons last interacting with a given line

        This includes being emitted in the specified line, or scattered off it

        Args:
            number (int): Python line number. Will vary based on data file!
            wavelength (float): Wavelength of the line in angstroms

        Returns:
            TransferFunction: Self, so filters can be stacked
        """
        assert self._unpickled is False,\
            "TF restored from pickle! It cannot be rerun, filters cannot be applied."
        self._line_wave = wavelength
        self._line_num = number
        self._query = self._query.filter(Photon.Resonance == number)
        return self

    def velocities(self, velocity: float) -> 'TransferFunction':
        """
        Constrain the TF to only photons with a range of Doppler shifts

        Args:
            velocity (float): Maximum doppler shift velocity in m/s. Applies
                              to both positive and negative Doppler shift

        Returns:
            TransferFunction: Self, so filters can be stacked
        """
        assert self._unpickled is False,\
            "TF restored from pickle! It cannot be rerun, filters cannot be applied."
        assert self._line_wave is not None,\
            "Cannot limit doppler shift around a line without specifying a line!"
        self._velocity = velocity
        self._query = self._query.filter(Photon.Wavelength >= doppler_shift_wave(self._line_wave, -velocity),
                                         Photon.Wavelength <= doppler_shift_wave(self._line_wave, velocity))
        return self

    def wavelengths(self, wave_min: float, wave_max: float) -> 'TransferFunction':
        """
        Constrain the TF to only photons with a range of wavelengths

        Args:
            wave_min (float): Minimum wavelength in angstroms
            wave_max (float): Maximum wavelength in angstroms

        Returns:
            TransferFunction: Self, so filters can be stacked
        """
        assert self._unpickled is False,\
            "TF restored from pickle! It cannot be rerun, filters cannot be applied."
        assert wave_min < wave_max,\
            "Minimum wavelength must be lower than maximum wavelength!"
        self._wave_range = [wave_min, wave_max]
        self._query = self._query.filter(Photon.Wavelength >= wave_min, Photon.Wavelength <= wave_max)
        return self

    def wavelength_bins(self, wave_range: np.ndarray) -> 'TransferFunction':
        """
        Constrain the TF to only photons with a range of wavelengths, and to a specific set of bins

        Args:
            wave_range (np.ndarray): Array of bins to use

        Returns:
            TransferFunction: Self, so filters can be stacked
        """
        assert self._unpickled is False,\
            "TF restored from pickle! It cannot be rerun, filters cannot be applied."
        assert len(wave_range) > 2,\
            "When providing an array, it must be of more than 2 entries! Use wavelength(min, max)."
        self._bins_wave = wave_range
        self._bins_wave_count = len(wave_range)-1
        self.wavelengths(self._bins_wave[0], self._bins_wave[-1])
        return self

    def lines(self, line_list: List[int]) -> 'TransferFunction':
        """
        Constrain the TF to only photons with a specific internal line number.
        This list number will be specific to the python atomic data file!

        Args:
            line_list (List[int]): List of lines

        Returns:
            TransferFunction: Self, so filters can be stacked
        """
        assert self._unpickled is False,\
            "TF restored from pickle! It cannot be rerun, filters cannot be applied."
        assert len(line_list) > 1,\
            "For a single line, use the 'line()' filter rather than 'lines()'!"
        self._line_list = line_list
        self._query = self._query.filter(Photon.Resonance.in_(line_list))
        return self

    def delays(self, delay_min: float, delay_max: float, days: bool = True) -> 'TransferFunction':
        """
        The delay range that should be considered when producing the TF.

        Args:
            delay_min (float): Minimum delay time (in seconds or days)
            delay_max (float): Maximum delay time (in seconds or days)
            days (bool): Whether or not the delay range has been provided in days

        Returns:
            TransferFunction: Self, so filters can be stacked
        """
        assert self._unpickled is False,\
            "TF restored from pickle! It cannot be rerun, filters cannot be applied."
        assert delay_min < delay_max,\
            "Minimum delay must be below maximum delay!"

        if days:
            self._delay_range = [delay_min * SECONDS_PER_DAY, delay_max * SECONDS_PER_DAY]
        else:
            self._delay_range = [delay_min, delay_max]
        self._query = self._query.filter(Photon.Delay > self._delay_range[0], Photon.Delay < self._delay_range[1])
        return self

    def delay_dynamic_range(self, delay_dynamic_range: float) -> 'TransferFunction':
        """
        If set, the TF will generate delay bins to cover this dynamic range of responses,
        i.e. (1 - 10^-ddr) of the delays.
        So a ddr of 1 will generate photons with delays up to 1 - (1/10) = the 90th percentile
        of delays. ddr=2 will give up to the 99th percentile, 3=99.9th percentile, etc.

        Arguably this is a bit of an ambiguous name

        Args:
            delay_dynamic_range (float): The dynamic range to be used when

        Returns:
            TransferFunction: Self, so filters can be stacked
        """
        assert self._unpickled is False,\
            "TF restored from pickle! It cannot be rerun, filters cannot be applied."
        assert delay_dynamic_range > 0,\
            "Cannot have a negative dynamic range!"
        self._delay_dynamic_range = delay_dynamic_range
        return self

    def cont_scatters(self, scat_min: int, scat_max: Optional[int] = None) -> 'TransferFunction':
        """
        Constrain the TF to only photons that have scattered min-max times via a
        continuum scattering process (e.g. electron scattering).

        Args:
            scat_min (int): Minimum number of continuum scatters
            scat_max (Optional[int]): Maximum number of continuum scatters, if desired

        Returns:
            TransferFunction: Self, so filters can be stacked
        """
        assert self._unpickled is False,\
            "TF restored from pickle! It cannot be rerun, filters cannot be applied."
        if scat_max is not None:
            assert scat_min < scat_max,\
                "Minimum continuum scatters must be below maximum scatters!"
        assert scat_min >= 0,\
            "Must select a positive number of continuum scatters"

        if scat_max is not None:
            self._query = self._query.filter(Photon.ContinuumScatters >= scat_min, Photon.ContinuumScatters <= scat_max)
        else:
            self._query = self._query.filter(Photon.ContinuumScatters == scat_min)
        return self

    def res_scatters(self, scat_min: int, scat_max: Optional[int] = None) -> 'TransferFunction':
        """
        Constrain the TF to only photons that have scattered min-max times via a
        resonant scattering process (e.g. line scattering).

        Args:
            scat_min (int): Minimum number of resonant scatters
            scat_max (Optional[int]): Maximum number of resonant scatters, if desired

        Returns:
            TransferFunction: Self, so filters can be stacked
        """
        assert self._unpickled is False,\
            "TF restored from pickle! It cannot be rerun, filters cannot be applied."
        if scat_max is not None:
            assert scat_min < scat_max,\
                "Minimum resonant scatters must be below maximum scatters!"
        assert scat_min >= 0,\
            "Must select a positive number of resonant scatters"

        if scat_max is not None:
            self._query = self._query.filter(Photon.ResonantScatters >= scat_min, Photon.ResonantScatters <= scat_max)
        else:
            self._query = self._query.filter(Photon.ResonantScatters == scat_min)
        return self

    def filter(self, *args) -> 'TransferFunction':
        """
        Apply a SQLalchemy filter directly to the content.

        Args:
            args: The list of filter arguments

        Returns:
            TransferFunction: Self, so filters can be stacked
        """
        assert self._unpickled is False,\
            "TF restored from pickle! It cannot be rerun, filters cannot be applied."
        self._query = self._query.filter(args)
        return self

    def response_map_by_tf(
        self, low_state: 'TransferFunction', high_state: 'TransferFunction',
        cf_low: float = 1.0, cf_high: float = 1.0
    ) -> 'TransferFunction':
        """
        Creates a response function for this transfer function by subtracting two transfer functions bracketing it.
        Requires two other completed transfer functions, bracketing this one in luminosity,
        all with matching wavelength/velocity and delay bins.

        Correction factors are there to account for things like runs that have been terminated early,
        e.g. if you request 100 spectrum cycles and stop (or Python dies) after 80,
        the total photon luminosity will only be 80/100. A correction factor allows you to bump this up.
        Arguably correction factors should be applied during the 'run()' method.

        Args:
            low_state (TransferFunction): A full, processed transfer function for a lower-luminosity system.
            high_state (TransferFunction): A full, processed transfer function for a higher-luminosity system.
            cf_low (float): Correction factor for low state. Multiplier to the whole transfer function.
            cf_high (float): Correction factor for high state. Multiplier to the whole transfer function.

        Returns:
            TransferFunction: Self, so plotting can be chained on.
        """
        assert self._emissivity is not None,\
            "You must run the TF query with '.run()' before response mapping it!"
        assert low_state._emissivity is not None and high_state._emissivity is not None,\
            "You must run the low and high state TF queries with '.run()' before response mapping using them!"
        assert np.array_equal(self._bins_wave, low_state._bins_wave) \
            and np.array_equal(self._bins_delay, low_state._bins_delay),\
            "Low state TF is binned differently to target TF! Cannot rescale using it."
        assert np.array_equal(self._bins_wave, high_state._bins_wave) \
            and np.array_equal(self._bins_delay, high_state._bins_delay),\
            "High state TF is binned differently to target TF! Cannot rescale using it."
        assert self._continuum is not None,\
            "TF missing continuum luminosity information!"
        assert low_state._continuum is not None,\
            "Low state TF missing continuum luminosity information!"
        assert high_state._continuum is not None,\
            "High state TF missing continuum luminosity information!"
        assert low_state._continuum <= self._continuum,\
            "Low state ionising luminosity greater than target TF ionising luminosity!"
        assert high_state._continuum >= self._continuum,\
            "High state ionising luminosity lower than target TF ionising luminosity!"

        continuum_difference = high_state._continuum - low_state._continuum

        # We divide the difference in TFs by the luminosity difference
        self._response = ((high_state._emissivity*high_state._continuum*cf_high) -
                          (low_state._emissivity*low_state._continuum*cf_low)) / continuum_difference
        return self

    def fwhm(self, response: bool = False, velocity: bool = True):
        """
        Calculates the full width half maximum of the delay-summed transfer function,
        roughly analogous to the line profile. Possibly meaningless for the response function?

        Args:
            response (bool): Whether to calculate the FWHM of the transfer or response function
            velocity (bool): Whether to return the FWHM in wavelength or velocity-space

        Returns:
            float: Full width at half maximum for the function.
                If the function is a doublet, this will not work properly.

        Todo:
            Catch doublets.
        """
        if response and self._response is None:
            raise RuntimeError("Must generate the response map using `response_map_from_tf()` before attempting to use it!")

        if velocity:
            midpoints = calculate_midpoints(self._bins_vel)
        else:
            midpoints = calculate_midpoints(self._bins_wave)

        if response:
            return calculate_fwhm(midpoints, np.sum(self._response, 0))

        return calculate_fwhm(midpoints, np.sum(self._emissivity, 0))

    def delay(
            self, response: bool = False, threshold: float = 0, bounds: float = None, days: bool = False
    ) -> Union[float, Tuple[float, float, float]]:
        """
        Calculates the centroid delay for the current data

        Args:
            response (bool): 
                Whether or not to calculate the delay from the response
            threshold (float): 
                Exclude all bins with value < the threshold fraction of the peak value.
                Standard value used in the reverb papers was 0.8.
            bounds (float): 
                Return the fractional bounds (i.e. `bounds=0.25`,
                the function will return `[0.5, 0.25, 0.75]`). Not implemented.
            days (bool): 
                Whether to return the delay in days or seconds

        Returns:
            Union[float, Tuple[float, float, float]]:
                Centroid delay, and lower and upper fractional bounds if bounds keyword provided

        Todo:
            Implement fractional bounds. Should just be able to call the centroid_delay function!
        """
        if response and self._response is None:
            raise RuntimeError("Must generate the response map using `response_map_from_tf()` before attempting to use it!")
        if threshold >= 1 or threshold < 0:
            raise ValueError("Threshold is a multiplier to the peak flux! It must be between 0 and 1")

        if response:
            data = np.sum(self._response, 1)
        else:
            data = np.sum(self._emissivity, 1)
        value_threshold = np.amax(data) * threshold
        delay_midpoints = calculate_midpoints(self._bins_delay)

        delay_weighted = 0
        value_total = 0
        for value, delay in zip(data, delay_midpoints):
            if value >= value_threshold:
                delay_weighted += value * delay
                value_total += value

        if days:
            return delay_weighted/(value_total * SECONDS_PER_DAY)

        return delay_weighted/value_total

    def delay_peak(self, response: bool = False, days: bool = False) -> float:
        """
        Calculates the peak delay for the transfer or response function,
        i.e. the delay at which the response is strongest.

        Args:
            response (bool): Whether or not to calculate the peak transfer or response function.
            days (bool): Whether to return the value in seconds or days.

        Returns:
            float: The peak delay.
        """
        if response and self._response is None:
            raise RuntimeError(
                "Must generate the response map using `response_map_from_tf()` before attempting to use it!"
            )

        data = self.transfer_function_1d(response=response, days=days)
        peak = data[np.argmax(data[:, 1]), 0]
        return peak

    def run(self, scaling_factor: float = 1.0, limit: int = None, verbose: bool = False) -> 'TransferFunction':
        """
        Performs a query on the photon DB and bins it.

        A TF must be run *after* all filters are applied and before any attempts
        to retrieve or process data from it. This can be a time-consuming call,
        on the order of 1 minute per GB of input file.

        Args:
            scaling_factor (float): 
                1/Number of cycles in the spectra file
            limit (int): 
                Number of photons to limit the TF to, for testing.
                Recommend testing filters on a small number of photons to begin with.
            verbose (bool): 
                Whether to output exactly what the query is.

        Returns:
            TransferFunction:   Self, for chaining commands
        """
        assert self._unpickled is False,\
            "TF restored from pickle! It cannot be rerun."
        assert self._emissivity is None,\
            "TF has already been run!"
        assert scaling_factor > 0,\
            "Negative scaling factors make no sense!"
        assert limit is None or limit > 0,\
            "Limit must either be zero or a positive number!"
        start = time.process_time()

        if verbose:
            if limit is not None:
                print("Limited to {} results...".format(limit))
            if self._velocity is not None:
                print("Limited to velocities -{} to +{}".format(self._velocity, self._velocity))
            if self._bins_wave is not None:
                print("Limited to preset wavelength bins from {} to {}".format(self._bins_wave[0], self._bins_wave[-1]))
            elif self._wave_range is not None:
                print("Limited to wavelengths {} to {}".format(self._wave_range[0], self._wave_range[1]))
            if self._line_num is not None:
                print("Limited to line {}, wavelength {}".format(self._line_num, self._line_wave))
            if self._spectrum is not None:
                print("Limited to spectrum {}".format(self._spectrum))
            if self._delay_range is not None:
                print("Limited to delays {} to {}".format(self._delay_range[0], self._delay_range[1]))

        if limit is None:
            data = np.asarray(self._query.all())
        else:
            data = np.asarray(self._query.limit(limit).all())

        assert len(data) > 0,\
            "No records found!"

        if verbose:
            print("Fetched {} records from '{}'...".format(len(data), self._filename))

        # Check if we've already got delay bins from another TF
        if self._bins_delay is None:
            # Data returned as Wavelength, Delay, Weight. Find min and max delays
            if self._delay_dynamic_range is not None:
                percentile = (1 - (10**(-self._delay_dynamic_range)))*100
                range_delay = [0, np.percentile(data[:, 1], percentile)]
                if verbose:
                    print(
                        "Delays up to the {} percentile value, {}d".format(
                            percentile, range_delay[1] / SECONDS_PER_DAY
                        )
                    )
            else:
                range_delay = [0, np.amax(data[:, 1])]

            self._bins_delay = np.linspace(
                range_delay[0], range_delay[1], self._bins_delay_count+1, endpoint=True, dtype=np.float64
            )

        # Check if we've already got wavelength bins from another TF
        if self._bins_wave is None:
            # If we have no velocity bins, this is a factory-fresh TF
            if self._bins_vel is None:
                # Data returned as Wavelength, Delay, Weight. Find min and max delays and wavelengths
                range_wave = [np.amin(data[:, 0]), np.amax(data[:, 0])]
            # If we do have velocity bins, this was templated off a different line
            # and we need to copy the velocities (but bins are in km! not m!)
            else:
                range_wave = [
                    doppler_shift_wave(self._line_wave, self._bins_vel[0] * 1000),
                    doppler_shift_wave(self._line_wave, self._bins_vel[-1] * 1000)
                ]
                print(
                    "Creating new wavelength bins from template, velocities from {:.2e}-{:.2e} to waves: {:.2f}-{:.2f}"
                    .format(self._bins_vel[0], self._bins_vel[-1], range_wave[0], range_wave[-1])
                )

            # Now create the bins for each dimension
            self._bins_wave = np.linspace(range_wave[0], range_wave[1],
                                          self._bins_wave_count+1, endpoint=True, dtype=np.float64)

        # Check if we've already got velocity bins from another TF and we have a line to center around
        if self._bins_vel is None and self._line_wave is not None:
            range_wave = [self._bins_wave[0], self._bins_wave[-1]]
            self._bins_vel = np.linspace(
                doppler_shift_vel(self._line_wave, range_wave[0]),
                doppler_shift_vel(self._line_wave, range_wave[-1]),
                self._bins_wave_count + 1, endpoint=True, dtype=np.float64
            )
            # Convert speed from m/s to km/s
            self._bins_vel = np.true_divide(self._bins_vel, 1000.0)

        # Now we bin the photons, weighting them by their photon weights for the luminosity
        self._emissivity, junk, junk = np.histogram2d(
            data[:, 1], data[:, 0], weights=data[:, 2], bins=[self._bins_delay, self._bins_wave]
        )
        # Keep an unweighted photon count for statistical error purposes
        self._count, junk, junk = np.histogram2d(
            data[:, 1], data[:, 0], bins=[self._bins_delay, self._bins_wave]
        )

        # Scaling factor! Each spectral cycle outputs L photons. If we do 50 cycles, we want a factor of 1/50
        self._emissivity *= scaling_factor
        # Scale to continuum luminosity
        self._emissivity /= self._continuum

        print("'{}' successfully run ({:.1f}s)".format(self._filename, time.process_time()-start))
        # Make absolutely sure this data is wiped as it's *HUGE*
        del data
        return self

    def _return_array(
        self, array: np.ndarray, delay: Optional[float] = None,
        wave: Optional[float] = None, delay_index: Optional[int] = None
    ) -> Union[int, float, np.ndarray]:
        """
        Internal function used by response(), emissivity() and count()

        Args:
            array (np.ndarray):    Array to return value from
            delay (Optional[float]):          Delay to return value for. Must provide this or delay_index.
            delay_index (Optional[int]):      Delay index to return value for. Must provide this or delay.
            wave (Optional[float]):           Wavelength to return value for

        Returns:
            Union[np.ndarray, float]: Either a subset of the array if only delay is provided,
                or the value of a single array element if delay and wavelength provided.

        Todo:
            Allow for only wavelength to be provided?
        """
        if delay is None and delay_index is None and wave is None:
            return array

        if delay is not None:
            if delay < self._bins_delay[0] or delay > self._bins_delay[-1]:
                if wave is None:
                    return np.zeros(self._bins_wave_count)
                else:
                    return 0
            delay_index = np.searchsorted(self._bins_delay, delay)

        elif delay_index is not None:
            if delay_index < 0 or delay_index > self._bins_delay_count:
                return 0

        if wave is None:
            return array[delay_index, :]

        return array[delay_index, np.searchsorted(self._bins_wave, wave)]

    def response_total(self) -> float:
        """
        Returns the total response.

        Returns:
            float: Total response.
        """
        if self._response is None:
            raise RuntimeError("Must generate the response map using `response_map_from_tf()` before attempting to use it!")
        return np.sum(self._response)

    def delay_bins(self) -> np.ndarray:
        """
        Returns the range of delays covered by this TF.

        Returns:
            np.ndarray: Array of the bin boundaries.
        """
        return self._bins_delay

    def response(
        self, delay: Optional[float] = None, wave: Optional[float] = None, delay_index: Optional[int] = None
    ) -> Union[float, np.ndarray]:
        """
        Returns the responsivity in either one specific wavelength/delay bin, or all wavelength bins
        for a given delay.

        Args:
            delay (Optional[float]):          Delay to return value for. Must provide this or delay_index.
            delay_index (Optional[int]):      Delay index to return value for. Must provide this or delay.
            wave (Optional[float]):           Wavelength to return value for.

        Returns:
            Union[int, np.ndarray]: Either the responsivity in one specific bin, or if wave is not specified
                the counts in each wavelength bin at this delay

        Todo:
            Allow for only wavelength to be provided?
        """
        if self._response is None:
            raise RuntimeError("Must generate the response map using `response_map_from_tf()` before attempting to use it!")
        return self._return_array(self._response, delay=delay, wave=wave, delay_index=delay_index)

    def emissivity(
        self, delay: Optional[float] = None, wave: Optional[float] = None, delay_index: Optional[int] = None
    ) -> Union[float, np.ndarray]:
        """
        Returns the emissivity in either one specific wavelength/delay bin, or all wavelength bins
        for a given delay.

        Args:
            delay (Optional[float]):          Delay to return value for. Must provide this or delay_index.
            delay_index (Optional[int]):      Delay index to return value for. Must provide this or delay.
            wave (Optional[float]):           Wavelength to return value for.

        Returns:
            Union[int, np.ndarray]: Either the emissivity in one specific bin, or if wave is not specified
                the counts in each wavelengthin bin at this delay

        Todo:
            Allow for only wavelength to be provided?
        """
        assert self._emissivity is not None,\
            "The TF has not been run! Use .run() to query the DB first."
        return self._return_array(self._emissivity, delay=delay, wave=wave, delay_index=delay_index)

    def count(
            self, delay: Optional[float] = None, wave: Optional[float] = None, delay_index: Optional[int] = None
    ) -> Union[int, np.ndarray]:
        """
        Returns the photon count in either one specific wavelength/delay bin, or all wavelength bins
        for a given delay.

        Args:
            delay (Optional[float]):          Delay to return value for. Must provide this or delay_index.
            delay_index (Optional[int]):      Delay index to return value for. Must provide this or delay.
            wave (Optional[float]):           Wavelength to return value for

        Returns:
            Union[int, np.ndarray]: Either the count in one specific bin, or if wave is not specified
                the counts in each wavelength bin at this delay

        Todo:
            Allow for only wavelength to be provided?
        """
        assert self._count is not None,\
            "The TF has not been run! Use .run() to query the DB first."
        assert delay_index is not None or delay is not None,\
            "You must provide a delay, or a delay index!"
        return self._return_array(self._count, delay=delay, wave=wave, delay_index=delay_index)

    def transfer_function_1d(self, response: bool = False, days: bool = True) -> np.ndarray:
        """
        Collapses the 2-d transfer/response function into a 1-d response function, and returns
        the bin midpoints and values in each bin for plotting.

        Args:
            response (bool): Whether or not to return the response function data
            days (bool): Whether the bin midpoints should be in seconds or days

        Returns:
            np.ndarray: A [bins, 2]-d array containing the midpoints of the delay bins,
                and the value of the 1-d transfer or response function in each bin.
        """
        if response and self._response is None:
            raise RuntimeError(
                "Must generate the response map using `response_map_from_tf()` before attempting to use it!"
            )

        if response:
            if days:
                return np.column_stack(
                    (calculate_midpoints(self._bins_delay / SECONDS_PER_DAY), np.sum(self._response, 1))
                )

            return np.column_stack((calculate_midpoints(self._bins_delay), np.sum(self._response, 1)))

        else:
            if days:
                return np.column_stack(
                    (calculate_midpoints(self._bins_delay / SECONDS_PER_DAY), np.sum(self._emissivity, 1))
                )

            return np.column_stack((calculate_midpoints(self._bins_delay), np.sum(self._emissivity, 1)))

    def plot(
            self, log: bool = False, normalised: bool = False, rescaled: bool = False, velocity: bool = False,
            name: str = None, days: bool = True, response_map=False, keplerian: dict = None,
            dynamic_range: int = None, rms: bool = False, show: bool = False,
            max_delay: Optional[float] = None,
            format: str = '.eps',
            return_figure: bool = False,
    ) -> Union['TransferFunction', Figure]:
        """
        Takes the data gathered by calling 'run' and outputs a plot

        Args:
            log (bool):
                Whether the plot should be linear or logarithmic.
            normalised (bool):
                Whether or not to rescale the plot such that the total emissivity = 1.
            rescaled (bool):
                Whether or not to rescale the plot such that the maximum emissivity = 1.
            velocity (bool):
                Whether the plot X-axis should be velocity (true) or wavelength (false).
            name (Optional[str]):
                The file will be output to 'tf_filename.eps'. May add the 'name' component
                to modify it to 'tf_filename_name.eps'. Useful for adding e.g. 'c4' or 'log'.
            days (bool):
                Whether the plot Y-axis should be in days (true) or seconds (false).
            response_map (bool):
                Whether to plot the transfer function map or the response function.
            keplerian (Optional[dict]):
                A dictionary describing the profile of a keplerian disk, the bounds of which will
                be overlaid on the plot. Arguments include
                angle (float) - Angle of disk to the observer,
                mass (float) - Mass of the central object in M_sol,
                radius (Tuple(float, float)) - Inner and outer disk radii, in $r_{ISCO}$.
                include_minimum_velocity - Whether or not to include the outer disk velocity profile
                (default no).
            dynamic_range (Optional[int]):
                If the plot is logarithmic,
                the dynamic range the colour bar should show. If not provided,
                will attempt to use the base dynamic range property, otherwise
                will default to showing 99.9% of all emissivity.
            max_delay (Optional[float]):
                The optional maximum delay to plot out to.
            rms (bool):
                Whether or not the line profile panel should show the root mean squared line profile.
            show (bool):
                Whether or not to display the plot to screen.
            format (str):
                The output file format. .eps by default.
            return_figure:
                If true, return the figure instead of platting it.

        Returns:
            TransferFunction: Self, for chaining outputs
        """
        assert response_map is False or self._response is not None,\
            "No data available for response map!"
        assert log is False or response_map is False,\
            "Cannot plot a logarithmic response map!"
        assert normalised is False or rescaled is False,\
            "Cannot be both normalised and rescaled!"
        assert self._bins_wave is not None,\
            "You must run the TF query with '.run()' before plotting it!"

        # matplotlib.rcParams["text.usetex"] = "True"
        matplotlib.rcParams.update({'font.size': 14})

        start = time.process_time()
        if name is not None:
            print("Plotting to file '"+self._filename+"_"+name+".eps'...")
        else:
            print("Plotting to file '"+self._filename+".eps'...")

        if dynamic_range is not None:
            log_range = dynamic_range
        elif self._delay_dynamic_range is not None:
            log_range = self._delay_dynamic_range
        else:
            log_range = 3

        # Set up the multiplot figure and axis
        fig, ((ax_spec, ax_none), (ax_tf, ax_resp)) = plt.subplots(
            2, 2, sharex='col', sharey='row',
            gridspec_kw={'width_ratios': [3, 1], 'height_ratios': [1, 3]}
        )
        ax_none.axis('off')
        ax_resp.invert_xaxis()
        fig.subplots_adjust(hspace=0, wspace=0)

        if response_map:
            ratio = np.sum(self._response)/np.sum(self._emissivity)
            ratio_exp = np.floor(np.log10(ratio))
            ratio_text = '\n'

            if ratio_exp < -1 or ratio_exp > 1:
                ratio_text_exp = r"{}{:.0f}{}".format("{", ratio_exp, "}")
                ratio_text += r"${:.2f}\times 10^{}$".format(ratio/(10**ratio_exp), ratio_text_exp)
            else:
                ratio_text += r"${:.3g}$".format(ratio)

            ax_tf.text(
                0.05, 0.95, r"$\frac{\Delta L}{L}/\frac{\Delta C}{C}=$"+ratio_text,
                transform=ax_tf.transAxes, fontsize=18, verticalalignment='top', horizontalalignment='left'
            )

        # Set the properties that depend on log and wave/velocity status
        cb_label = None
        cb_label_vars = r""
        cb_label_units = r""
        cb_label_scale = r""
        cb_map = "afmhot_r"

        # Copy the data for later modification.
        data_plot = None
        if response_map:
            data_plot = np.copy(self._response)
            print("Total response: {:.3e}".format(np.sum(data_plot)))
            psi_label = r"$\Psi_{R}$"
        else:
            data_plot = np.copy(self._emissivity)
            print("Total line: {:.3e}".format(np.sum(data_plot)))
            psi_label = r"$\Psi_{T}$"
        cb_label = psi_label

        # Set the xlabel and colour bar label - these differ if velocity or not
        x_bin_mult = 1
        if velocity:
            # We're rescaling the axis to e.g. 10^3 km/s but the colorbar is still in km/s
            # So when we scale by bin width, we need a multiplier on the bin widths
            oom = np.log10(np.amax(self._bins_vel))
            oom = oom - oom % 3
            bins_x = self._bins_vel/(10**oom)
            x_bin_mult = 10**oom
            ax_tf.set_xlabel(r'Velocity ($10^{:.0f}$ km s$^{}$)'.format(oom, '{-1}'))
            cb_label_vars = r"($v, \tau$)"
            cb_label_units = r"/ km s$^{-1}$"
        else:
            bins_x = self._bins_wave
            ax_tf.set_xlabel(r'Wavelength $\lambda$ ($\AA$)')
            cb_label_vars += r"($\lambda, \tau$)"
            cb_label_units = r"/ $\AA$"

        bins_x_midp = np.zeros(shape=self._bins_wave_count)
        for i in range(0, self._bins_wave_count):
            bins_x_midp[i] = (bins_x[i] + bins_x[i+1]) / 2

        # Set the ylabel and y bins for whether it's in days or seconds
        if days:
            bins_y = np.true_divide(self._bins_delay, float(SECONDS_PER_DAY))
            data_plot *= SECONDS_PER_DAY
            ax_tf.set_ylabel(r'Delay $\tau$ (days)')
            cb_label_units += r' d'
        else:
            bins_y = self._bins_delay
            ax_tf.set_ylabel(r'Delay $\tau$ (seconds)')
            cb_label_units += r' s'

        bins_y_midp = np.zeros(shape=self._bins_delay_count)
        for bin_y in range(0, self._bins_delay_count):
            bins_y_midp[bin_y] = (bins_y[bin_y] + bins_y[bin_y+1]) / 2

        # Rescale the values to be luminosity/km s^-1 d or /A d
        for bin_y in range(0, self._bins_delay_count):
            width_y = bins_y[bin_y+1] - bins_y[bin_y]
            for bin_x in range(0, self._bins_wave_count):
                width_x = bins_x[bin_x+1] - bins_x[bin_x]
                data_plot[bin_y][bin_x] /= (width_x * x_bin_mult * width_y)

        # Plot the spectrum and light curve, normalised
        data_plot_spec = np.sum(data_plot, 0)
        data_plot_resp = np.sum(data_plot, 1)
        exponent_spec = np.floor(np.log10(np.amax(data_plot_spec)))
        exponent_resp = np.floor(np.log10(np.amax(data_plot_resp)))
        exponent_resp_text = "{}{:.0f}{}".format("{", exponent_resp, "}")
        exponent_spec_text = "{}{:.0f}{}".format("{", exponent_spec, "}")

        ax_resp.plot(data_plot_resp/(10**exponent_resp), bins_y_midp, c='m')

        if velocity:
            ax_spec.set_ylabel(r'{}(v) $10^{}/$'.format(psi_label, exponent_spec_text)+r'km s$^{-1}$', fontsize=12)
        else:
            ax_spec.set_ylabel(r'{}($\lambda$) $10^{}/$'.format(psi_label, exponent_spec_text)+r'$\AA$', fontsize=12)

        if days:
            ax_resp.set_xlabel(r'{}($\tau$) $10^{}$ / d'.format(psi_label, exponent_resp_text))
        else:
            ax_resp.set_xlabel(r'{}($\tau$) $10^{}$ / s'.format(psi_label, exponent_resp_text))

        if response_map and rms:
            ax_spec.axhline(0, color='grey')
            ax_resp.axvline(0, color='grey')

            data_plot_rms = np.sqrt(np.sum(np.square(data_plot), 0) / self._bins_wave_count)
            exponent_rms = np.floor(np.log10(np.amax(data_plot_rms)))
            exponent_rms_text = "{}{:.0f}{}".format("{", exponent_rms, "}")
            maximum_spec = np.amax(data_plot_spec)/np.power(10, exponent_spec)
            maximum_rms = np.amax(data_plot_rms)/np.power(10, exponent_rms)
            data_plot_rms /= np.amax(data_plot_rms)
            data_plot_spec /= np.amax(data_plot_spec)

            ax_spec.plot(
                bins_x_midp, data_plot_rms, c='c',
                label=r'RMS {}(v)/{:.2f}$x10^{}$'.format(psi_label, maximum_rms, exponent_rms_text)
            )
            ax_spec.plot(
                bins_x_midp, data_plot_spec, c='m',
                label=r'{}(v)/{:.2f}$x10^{}$'.format(psi_label, maximum_spec, exponent_spec_text)
            )
            ax_spec.legend(loc='center left', bbox_to_anchor=(1, 0.5), frameon=False)

        elif response_map:
            ax_spec.axhline(0, color='grey')
            ax_resp.axvline(0, color='grey')
            ax_spec.plot(bins_x_midp, data_plot_spec/(10**exponent_spec), c='m')

        else:
            ax_spec.plot(bins_x_midp, data_plot_spec/(10**exponent_spec), c='m')

        # If this is a log plot, take log and correct label and limits
        if log:
            cb_max = np.log10(np.amax(data_plot))
            cb_min = cb_max-log_range
            cb_label = r"Log "+cb_label
            data_plot = np.ma.log10(data_plot)
        # Else just scale the data
        else:
            maxval = np.floor(np.log10(np.amax(data_plot)))
            data_plot /= np.power(10, maxval)
            cb_max = np.amax(data_plot)
            cb_min = np.amin(data_plot)
            dummy = "{}{:.0f}{}".format("{", maxval, "}")
            cb_label_scale = r" 10$^{}$".format(dummy)

        # If this is a response map, it may have a negative component and need a different plot
        if response_map:
            cb_max = np.amax([cb_max, np.abs(cb_min)])
            cb_min = -cb_max
            cb_map = 'RdBu_r'

        # Normalise or rescale the data. If doing neither, put units on cb.
        if normalised:
            data_plot /= np.sum(data_plot)
            cb_label_units = r""
            cb_label_scale = r""
        elif rescaled:
            data_plot /= np.amax(data_plot)
            cb_label_units = r""
            cb_label_scale = r""

        # Plot the main colourplot for the transfer function
        tf = ax_tf.pcolor(
            bins_x,
            bins_y,
            data_plot,
            vmin=cb_min, vmax=cb_max, cmap=cb_map
        )
        if not max_delay:
            ax_tf.set_ylim(bottom=bins_y[0], top=bins_y[-1])
        else:
            ax_tf.set_ylim(bottom=bins_y[0], top=max_delay)
        ax_tf.set_xlim(left=bins_x[0], right=bins_x[-1])
        ax_tf.set_aspect('auto')

        # Add lines for keplerian rotational outflows
        if keplerian is not None:
            resolution = 1000
            # Technically radians can output an array if k['angle'] is a list/tuple.
            # We want the exception to happen on this line if that's so for clarity.
            r_angle = float(np.radians(keplerian["angle"]))
            r_mass_bh = keplerian["mass"] * apc.M_sun.value
            r_rad_grav = (6 * apc.G.value * r_mass_bh / np.power(apc.c.value, 2))
            ar_wave = np.zeros(resolution)  # * u.angstrom
            ar_delay = np.zeros(resolution)  # * u.s
            ar_phase = np.linspace(0, np.pi*2, resolution)
            ar_rad = np.linspace(keplerian["radius"][0]*r_rad_grav, 20*keplerian["radius"][1]*r_rad_grav, resolution)
            ar_vel = np.zeros(resolution)
            r_rad_min = r_rad_grav * keplerian["radius"][0]
            r_rad_max = r_rad_grav * keplerian["radius"][1]
            r_vel_min = keplerian_velocity(r_mass_bh, r_rad_max)
            r_vel_max = keplerian_velocity(r_mass_bh, r_rad_min)
            include_minimum = keplerian.get('include_minimum_velocity', False)

            # ITERATE OVER INNER EDGE
            for r_phase, r_wave, r_delay, r_vel in np.nditer(
                    [ar_phase, ar_wave, ar_delay, ar_vel], op_flags=['readwrite']
            ):
                r_vel[...] = r_vel_max * np.sin(r_phase) * np.sin(r_angle) / (1e3 * x_bin_mult)
                r_wave[...] = doppler_shift_wave(self._line_wave, r_vel * 1e3 * x_bin_mult)
                r_delay[...] = calculate_delay(r_angle, r_phase, r_rad_min, u.day)
            if velocity:
                ax_tf.plot(ar_vel, ar_delay, '-', c='m')
            else:
                ax_tf.plot(ar_wave, ar_delay, '-', c='m')

            # # ITERATE OVER OUTER EDGE
            if include_minimum:
                for r_phase, r_wave, r_delay, r_vel in np.nditer(
                        [ar_phase, ar_wave, ar_delay, ar_vel], op_flags=['readwrite']
                ):
                    r_vel[...] = r_vel_min * np.sin(r_phase) * np.sin(r_angle) / (1e3 * x_bin_mult)
                    r_wave[...] = doppler_shift_wave(self._line_wave, r_vel * 1e3 * x_bin_mult)
                    r_delay[...] = calculate_delay(r_angle, r_phase, r_rad_max, u.day)
                if velocity:
                    ax_tf.plot(ar_vel, ar_delay, '-', c='m')
                else:
                    ax_tf.plot(ar_wave, ar_delay, '-', c='m')

            # ITERATE OVER BLUE BOUND
            for r_rad, r_wave, r_delay, r_vel in np.nditer(
                    [ar_rad, ar_wave, ar_delay, ar_vel], op_flags=['readwrite']
            ):
                r_rad = r_rad  # * u.m
                r_vel[...] = keplerian_velocity(r_mass_bh, r_rad) * np.sin(r_angle) / (1e3 * x_bin_mult)
                r_wave[...] = doppler_shift_wave(self._line_wave, r_vel * 1e3 * x_bin_mult)
                r_delay[...] = calculate_delay(r_angle, np.pi/2, r_rad, u.day)
            if velocity:
                ax_tf.plot(ar_vel, ar_delay, '-', c='m')
            else:
                ax_tf.plot(ar_wave, ar_delay, '-', c='m')

            # ITERATE OVER RED BOUND
            for r_rad, r_wave, r_delay, r_vel in np.nditer(
                    [ar_rad, ar_wave, ar_delay, ar_vel], op_flags=['readwrite']
            ):
                r_rad = r_rad  # * u.m
                r_vel[...] = -keplerian_velocity(r_mass_bh, r_rad) * np.sin(r_angle) / (1e3 * x_bin_mult)
                r_wave[...] = doppler_shift_wave(self._line_wave, r_vel * 1e3 * x_bin_mult)
                r_delay[...] = calculate_delay(r_angle, np.pi/2, r_rad, u.day)
            if velocity:
                ax_tf.plot(ar_vel, ar_delay, '-', c='m')
            else:
                ax_tf.plot(ar_wave, ar_delay, '-', c='m')

        cbar = plt.colorbar(tf, orientation="vertical")
        cbar.set_label(cb_label+cb_label_vars+cb_label_scale+cb_label_units)

        if name is None:
            plt.savefig("{}.{}".format(self._filename, format), bbox_inches='tight')
            print("Successfully plotted '{}.eps'({:.1f}s)".format(self._filename, time.process_time()-start))
        else:
            plt.savefig("{}_{}.{}".format(self._filename, name, format), bbox_inches='tight')
            print("Successfully plotted '{}_{}.eps'({:.1f}s)".format(self._filename, name, time.process_time()-start))

        if show:
            fig.show()

        plt.close(fig)

        if return_figure:
            return fig
        else:
            return self


# ==============================================================================
def open_database(file_root: str, user: str = None, password: str = None, batch_size: int = 25000) -> Engine:
    """
    Open or create a SQL database

    Will open a SQLite DB if one already exists, otherwise will create one from
    file. Note, though, that if the process is interrupted the code cannot
    intelligently resume - you must delete the half-written DB!

    Args:
        file_root (string): 
            Root of the filename (no '.db' or '.delay_dump')
        user (string):      
            Username. Here in case I change to PostgreSQL
        password (string):  
            Password. Here in case I change to PostgreSQL
        batch_size (int):   
            Number of photons to stage before committing. If
            too low, file creation is slow. If too high, get
            out-of-memory errors.

    Returns:
        Connection to the database opened
    """

    print("Opening database '{}'...".format(file_root))

    try:
        db_engine = sqlalchemy.create_engine(
            f"sqlite:///{file_root}.db"
        )
    except sqlalchemy.exc.SQLAlchemyError as e:
        print(e)
        sys.exit(1)

    # DOES IT ALREADY EXIST? ###
    session_maker = sqlalchemy.orm.sessionmaker(bind=db_engine)
    session = session_maker()

    start = time.process_time()

    try:
        session.query(Photon.Weight).first()
        # If so, we go with what we've found.
        print("Found existing filled photon database '{}'".format(file_root))

    except sqlalchemy.exc.SQLAlchemyError:
        # If not, we populate from the delay dump file. This bit is legacy!
        print("No existing filled photon database, reading from file '{}.delay_dump'".format(file_root))
        Base.metadata.create_all(db_engine)

        added = 0
        delay_dump = open("{}.delay_dump".format(file_root), 'r')
        for line in delay_dump:
            # For each line in this file, if it is not a comment
            if line.startswith('#'):
                continue
            try:
                # Try reading it in as a series of values.
                values = [float(i) for i in line.split()]
            except ValueError:
                print("Malformed line: '{}'".format(line))
                continue

            if len(values) != 14:
                # There should be 14 values per line in our base formatting!
                print("Incorrect number of values in line: '{}'".format(line))
                continue

            # Add the photo using the values. Some must be modified here; ideally, this would be done in Python.
            session.add(
                Photon(
                    Wavelength=values[2], Weight=values[3],
                    X=values[4], Y=values[5], Z=values[6],
                    ContinuumScatters=int(values[7]-values[8]), ResonantScatters=int(values[8]),
                    Delay=values[9],
                    Spectrum=int(values[10]), Origin=int(values[11] % 10),
                    Resonance=int(values[12]), LineResonance=int(values[13]),
                    Origin_matom=(values[11] > 9)
                )
            )
            added += 1
            if added > batch_size:
                # We commit in batches in order to avoid out-of-memory errors
                added = 0
                session.commit()

        session.commit()
        session.close()
        del session
        print("Successfully read in ({:.1f}s)".format(time.process_time()-start))

    return db_engine
# ==============================================================================


Base = sqlalchemy.ext.declarative.declarative_base()
"""Base class declared dynamically to bind to SQLalchemy"""

class Spectrum(Base):
    """
    The SQLalchemy table for the spectra. Unused.
    Could be removed but will break backward compatibility.
    Information required for this is not stored in the output files.

    # Todo: Implement or remove this table.
    """
    __tablename__ = "Spectra"
    id = sqlalchemy.Column(sqlalchemy.Integer, primary_key=True)
    angle = sqlalchemy.Column(sqlalchemy.Float)


class Origin(Base):
    """
    The SQLalchemy table for the photon origins. Unused.
    Could be removed but will break backward compatibility.
    Information required for this is not stored in the output files.

    # Todo: Implement or remove this table
    """
    __tablename__ = "Origin"
    id = sqlalchemy.Column(sqlalchemy.Integer, primary_key=True)
    name = sqlalchemy.Column(sqlalchemy.String)


class Photon(Base):
    """
    SQLalchemy class for a photon. Why are all the properties capitalised?
    Changing them to lowercase as would make sense breaks backwards compatibility.

    # ToDo: Change to lower case.
    """
    __tablename__ = "Photons"
    id = sqlalchemy.Column(sqlalchemy.Integer, primary_key=True, autoincrement=True)
    Wavelength = sqlalchemy.Column(sqlalchemy.Float)
    Weight = sqlalchemy.Column(sqlalchemy.Float)
    X = sqlalchemy.Column(sqlalchemy.Float)
    Y = sqlalchemy.Column(sqlalchemy.Float)
    Z = sqlalchemy.Column(sqlalchemy.Float)
    ContinuumScatters = sqlalchemy.Column(sqlalchemy.Integer)
    ResonantScatters = sqlalchemy.Column(sqlalchemy.Integer)
    Delay = sqlalchemy.Column(sqlalchemy.Float)
    Spectrum = sqlalchemy.Column(sqlalchemy.Integer)
    Origin = sqlalchemy.Column(sqlalchemy.Integer)
    Resonance = sqlalchemy.Column(sqlalchemy.Integer)
    LineResonance = sqlalchemy.Column(sqlalchemy.Integer)
    Origin_matom = sqlalchemy.Column(sqlalchemy.Boolean)


kep_sey = {"angle": 40, "mass": 1e7, "radius": [50, 2000]}
"""The default Keplerian outline settings for the Seyfert model"""

kep_qso = {"angle": 40, "mass": 1e9, "radius": [50, 20000]}
"""The default Keplerian outline settings for the QSO model"""
