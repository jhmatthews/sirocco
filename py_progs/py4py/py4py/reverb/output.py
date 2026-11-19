"""
Reverberation mapping output module

This contains helper functions that bundle up the production of multiple transfer functions or the creation of
response functions.
"""
# -*- coding: utf-8 -*-
# pylint: disable=C0301
import numpy as np
from typing import List, Optional

from py4py.reverb import TransferFunction


def do_tf_plots(
    tf_list_inp: List[TransferFunction], dynamic_range: Optional[float] = None,
    keplerian: Optional[dict] = None, name: Optional[str] = None,
    file: Optional[str] = None, threshold: float = 0.8
):
    """
    Produces plots of the transfer functions for a list of provided TFs,
    with matching plotted dynamic ranges and name structures.
    Also optionally generates a file containing the centroid delays for each transfer function.

    Arguments:
        tf_list_inp (List[TransferFunction]): The transfer functions to plot.
        dynamic_range (Optional[float]): The dynamic range to plot them all across (see TransferFunction.plot).
        keplerian (Optional[dict]): The keplerian orbit parameters to overplot on all of them.
        name (Optional[str]): The name component to be added to the output plot filenames,
            e.g. [tf.name]_[name].eps. Useful for ending up with e.g. c4_with_keplerian.eps.
        file (Optional[str]): The output filename for the list of centroid delays.
        threshold (float): The peak flux threshold to use for calculating the centroid (see TransferFunction.delay).

    Outputs:
       ` {file}_tf_delay.txt` [Optional]
        `{tf.filename}_name.eps` [for each tf in tf_list_inp]
    """
    delays: List[float] = []
    for tf_inp in tf_list_inp:
        tf_inp.plot(velocity=True, keplerian=keplerian, log=False, name=name)
        tf_inp.plot(
            velocity=True, keplerian=keplerian, log=True,
            name=('log' if name is None else name+"_log"), dynamic_range=dynamic_range
        )
        delays.append(tf_inp.delay(threshold=threshold))

    if file is not None:
        print("Saving centroid transfer function delays to file: {}".format(file+"_tf_delay.txt"))
        np.savetxt(file+"_tf_delay.txt", np.array(delays, dtype='float'), header="Delay")


def do_rf_plots(
    tf_min: TransferFunction, tf_mid: TransferFunction, tf_max: TransferFunction,
    keplerian: Optional[dict] = None, name: Optional[str] = None, file: Optional[str] = None
):
    """
    Do response plot for a transfer function, optionally with keplerian disk lines on it.
    This will generate not just a response function from the two bracketing TFs, but one
    from the midpoint TF and the bracketing ones (e.g. Min-Max, Min-Mid, Mid-Max).
    Ideally, all three should be similar- if they are not, it suggests that the change
    in luminosity from min-max covers a point of inflection in d[Ionisation state]dL.

    Arguments:
        tf_min (TransferFunction): The low-state TF
        tf_mid (TransferFunction): The TF bracketed by low and high to produce the RF for.
        tf_max (TransferFunction): The high-state TF
        keplerian (Optional[dict]): The keplerian orbit parameters to overplot on all of them.
        name (Optional[str]): The name component to be added to the output plot filenames,
            e.g. [tf.name]_[name].eps. Useful for ending up with e.g. c4_with_keplerian.eps.
        file (Optional[str]): The output filename for the list of centroid delays.

    Outputs:
        `{tf_mid.name}_resp_mid.eps`
        `{tf_mid.name}_resp_low.eps`
        `{tf_mid.name}_resp_high.eps`
        `{file}_rf_delay.txt` [Optional]
    """
    if name is not None:
        name += '_'
    else:
        name = ''

    total_min: float = np.sum(tf_min._emissivity).item()
    total_mid: float = np.sum(tf_mid._emissivity).item()
    total_max: float = np.sum(tf_max._emissivity).item()

    calibration_factor: float = total_mid / ((total_min + total_max) / 2.0)

    tf_mid.response_map_by_tf(tf_min, tf_max, cf_low=1.0, cf_high=1.0).plot(
        velocity=True, response_map=True, keplerian=keplerian, name=name+"resp_mid"
    )
    rf_mid = tf_mid.delay(response=True, threshold=0.8)

    tf_mid.response_map_by_tf(tf_min, tf_mid, cf_low=calibration_factor, cf_high=1.0).plot(
        velocity=True, response_map=True, keplerian=keplerian, name=name+"resp_low"
    )
    rf_min = tf_mid.delay(response=True, threshold=0.8)
    tf_mid.response_map_by_tf(tf_mid, tf_max, cf_low=1.0, cf_high=calibration_factor).plot(
        velocity=True, response_map=True, keplerian=keplerian, name=name+"resp_high"
    )
    rf_max = tf_mid.delay(response=True, threshold=0.8)

    if file is not None:
        print("Saving RF plots to file: {}".format(file+"_rf_delay.txt"))
        np.savetxt(file+"_rf_delay.txt", np.array([rf_min, rf_mid, rf_max], dtype='float'), header="Delay")
