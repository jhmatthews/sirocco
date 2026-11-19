"""
Array Functions

Standard functions used for manipulating arrays, e.g. to calculate full width half maxima, centroids or midpoints.
"""

from typing import Union, Tuple
from numpy.typing import NDArray
from numpy import floating
import numpy as np


def calculate_fwhm(
        midpoints: NDArray[floating],
        vals: NDArray[floating]
) -> float:
    """
    Calculate FWHM from arrays

    Taken from http://stackoverflow.com/questions/10582795/finding-the-full-width-half-maximum-of-a-peak
    I don't think this can cope with being passed a doublet or an array with no
    peak within it. Doublets will calculate FWHM from the HM of both!

    Args:
        midpoints: Array of bin midpoints
        vals: Array of bin values

    Returns:
        FWHM of the peak (should it exist!)
    """
    # Create 'difference' array by subtracting half maximum
    difference = vals - (np.amax(vals) / 2)
    # Find the points where the difference is positive
    indexes = np.where(difference > 0)[0]
    # The first and last positive points are the edges of the peak
    return abs(midpoints[indexes[-1]] - midpoints[indexes[0]])


def calculate_centroid(
    bins: np.ndarray, vals: np.ndarray, bounds: float = None
) -> Union[float, Tuple[float, float, float]]:
    """
    Returns the centroid position, with optional percentile bounds.

    Args:
        bins (np.ndarray):
            Array of bin bounds
        vals (np.ndarray):
            Array of bin values
        bounds (float):
            Fraction from 0-0.5. Percentile either side of the centroid to find (e.g. .2 -> 30%, 70%)

    Returns:
        Flux-weighted centroid, and if 'bounds' passed both lower and upper percentile bounds
    """
    centroid_total = np.sum(vals)
    centroid_position = np.sum(np.multiply(bins, vals))/centroid_total

    if bounds is not None:
        # If we're finding bounds
        bound_width = bounds/2
        bound_min = -1
        bound_max = -1
        # Find the upper bound
        value_total = 0
        for index, value in enumerate(vals):
            # Starting at 0, add the value in this bin to the running total
            value_total += value
            if value_total/centroid_total >= 0.5+bound_width:
                # If this total is > the bound we're looking for, record the bin and stop
                bound_max = bins[index]
                break
        # Find the lower bound
        value_total = centroid_total
        for index, value in enumerate(vals[::-1]):
            # Starting at the total value, subtract the value in this bin from the running total
            value_total -= value
            if value_total/centroid_total <= 0.5-bound_width:
                # If this total is < the bound we're looking for, record the bin and stop
                bound_min = bins[len(bins)-1-index]
                break
        # On reflection, they could both sum since I'm just iterating backwards.
        # Also, I could use zip() even though they're numpy arrays as zip works fine
        # if you don't want to modify the array entries.
        # Maybe go over this later, should be easy enough to test.

        # Return the centroid and the bins.
        # NOTE: If the value exceeds the bound range midway through a cell, it'll just return the min/max
        # for that cell as appropriate. This will overestimate the error on the centroid.
        return centroid_position, bound_min, bound_max
    else:
        return centroid_position


def calculate_midpoints(
        bins: NDArray[floating]
) -> NDArray[floating]:
    """
    Converts bin boundaries into midpoints

    Args:
        bins: Array of bin boundaries

    Returns:
        Array of bin midpoints (1 shorter!)
    """
    midpoints = np.zeros(shape=len(bins)-1)
    for i in range(0, len(bins)-1):
        midpoints[i] = (bins[i] + bins[i+1]) / 2
    return midpoints
