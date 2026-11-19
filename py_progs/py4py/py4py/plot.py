"""
Functions designed for plotting output files
"""
from typing import  List, Optional, Tuple

import numpy as np
from numpy import floating
from numpy.typing import NDArray
from matplotlib.figure import Figure
from matplotlib import pyplot as plt
from astropy.table import Table


def plot_spec(
        col: str,
        spectra: List[Table],
        names: List[str],
        log_x: bool = False,
        log_y: bool = False,
        scale_to: Optional[float] = None,
        lim_x: Optional[Tuple[float, float]] = None
) -> Figure:
    """
    Plots an array of python spectra imported as Tables.

    Example:

        To plot the 40* viewing angle of two spectra::

            plot_spec('A40P0.50', [spec_r01, spec_r10], ['1x', '10x'])

    Arguments:
        col: The column in the Astropy Table containing the spectra to plot.
        spectra: The list of Astropy Tables containing spectra to plot.
        names: The list of names of the spectra for the key.
        log_x: Whether the x-axis should be logarithmic or not.
        log_y: Whether the y-axis should be logarithmic or not.
        scale_to: Whether to rescale the plot to normalise to a wavelength.
        lim_x: If provided, the lower and upper bounds for the x-axis.

    Returns:
        The generated figure.
    """
    fig, ax = plt.subplots()
    ax.set_xlabel(r"Wavelength (\AA)")

    if scale_to:
        ax.set_ylabel(r"$F_{\lambda}/F_{"+str(scale_to)+r"\AA}$")
    else:
        ax.set_ylabel(r"Flux")
        scale_factor = 1

    max_y = -1
    min_y = 99999999

    for spectrum, name in zip(spectra, names):
        if scale_to:
            # Annoyingly, the file is sorted by frequency not wavelength
            # so we need to flip it to get searchsorted to run
            index_scale = np.searchsorted(spectrum['Lambda'][::-1], scale_to)
            scale_factor = spectrum[col][::-1][index_scale]

        if lim_x:
            min_y_curr = np.amin(spectrum[col][(spectrum['Lambda'] > lim_x[0]) & (spectrum['Lambda'] < lim_x[1])]) / scale_factor
            if min_y_curr < min_y:
                min_y = min_y_curr
            max_y_curr = np.amax(spectrum[col][(spectrum['Lambda'] > lim_x[0]) & (spectrum['Lambda'] < lim_x[1])]) / scale_factor
            if max_y_curr > max_y:
                max_y = max_y_curr

        if log_x and log_y:
            ax.loglog(spectrum['Lambda'], spectrum[col]/scale_factor, label=name)
        elif log_y:
            ax.semilogy(spectrum['Lambda'], spectrum[col]/scale_factor, label=name)
        elif log_x:
            ax.semilogx(spectrum['Lambda'], spectrum[col]/scale_factor, label=name)
        else:
            ax.plot(spectrum['Lambda'], spectrum[col]/scale_factor, label=name)

    if lim_x:
        ax.set_xlim(lim_x[0], lim_x[1])
        ax.set_ylim(min_y, max_y)

    ax.legend()
    return fig


def plot_dat(
        table: Table,
        grid_x: NDArray[floating],
        grid_z: NDArray[floating],
        title:str,
        label: str,
        volume: bool = True
):
    """
    Plots a given `py_wind` `.dat` file for a wind property.

    Example:

        To plot the H-I fraction of a model::

            plot_dat(
                table_h1_r01, x_r01, z_r01,
                'H-I, radius 1x', 'Log ion fraction', volume=False
            )

    Arguments:
        table: The Astropy Table containing the output of `py_wind` for a property.
        grid_x: The bounds of the X grid for the model.
        grid_z: The bounds of the Z grid for the model.
        title: The title of the plot.
        label: The colour bar label.
        volume: Whether to correct the value to be per-volume considering the model grid.

    Returns:
        The generated figure.
    """
    fig, ax = plt.subplots()
    fig.suptitle(title)
    ax.set_xlabel('Log X (cm)')
    ax.set_ylabel('Log Z (cm)')

    size = (len(grid_x)-1, len(grid_z)-1)
    data = np.reshape(table['var'], size)
    if volume:
        for xi in range(size[0]):
            for zi in range(size[1]):
                # We need to correct the per-area emission to a per-volume emission
                area = 2. * 2. * np.pi * (grid_x[xi+1]-grid_x[xi]) * (grid_z[zi+1]-grid_z[zi])
                area += 2. * np.pi * grid_x[xi+1] * (grid_z[zi+1]-grid_z[zi])
                area += 2. * np.pi * grid_x[xi]   * (grid_z[zi+1]-grid_z[zi])
                if area > 0:
                    data[xi, zi] = data[xi, zi] / area

    im = ax.pcolormesh(np.log10(grid_x), np.log10(grid_z), np.ma.log10(data.T))
    ax.set_xlim(14.5, 18)
    ax.set_ylim(13, 17)
    cbar = fig.colorbar(im, ax=ax).set_label(label)
    return fig


def plot_dat_many(
        tables: List[Table],
        grids_x: List[NDArray[floating]],
        grids_z: List[NDArray[floating]],
        x_lims: List[Tuple[float, float]],
        z_lims: List[Tuple[float, float]],
        titles: List[str],
        title: str,
        label: str,
        shared_y: bool = False,
        shared_cbar: bool = False,
        volume: bool = True,
        log: bool = True
    ):
    """
    Plots multiple `py_wind` `.dat` files of a wind property for multiple models.

    Example::

        To plot the H-I fraction of multiple models::

            plot_dat_many(
                [table_h1_r01, table_h1_r10, table_h1_r30],
                [x_r01, x_r10, x_r30], [z_r01, z_r10, z_r30],
                x_lims=[(14.5, 17.5), (15.5, 17.5), (16, 17.5)],
                z_lims=[(13, 17), (13, 17), (13, 17)],
                titles=['1x Radius', '10x Radius', '30x Radius'],
                title='H-I ion fraction', label='Log ion fraction',
                shared_y=True, volume=False, shared_cbar=True
            )

    Arguments:
        tables: The Astropy Tables containing the outputs of `py_wind` for a property.
        grids_x: The bounds of the X grid for each model.
        grids_z: The bounds of the Z grid for each model.
        x_lims: The lower and upper ranges of the x-axis for each subplot.
        z_lims: The lower and upper ranges of the z-axis for each subplot.
        titles: The titles of the subplots.
        title: The title of the plot.
        label: The colour bar label.
        shared_y: Whether the plots should share Y-axes.
        shared_cbar: Whether the plots should share the colour bar.
        volume: Whether to correct the value to be per-volume considering the model grid.
        log: Whether the colour scale should be logarithmic or not.

    Returns:
        The generated figure.
    """
    if shared_y:
        fig, axes = plt.subplots(1, len(tables), sharey='row')
    else:
        fig, axes = plt.subplots(1, len(tables))

    fig.suptitle(title, y=1.03)
    axes[0].set_ylabel('Log Z (cm)')
    images = []

    if shared_cbar:
        v_min = float('inf')
        v_max = -float('inf')

    for ax, table, grid_x, grid_z, xlim, zlim, title in zip(axes, tables, grids_x, grids_z, x_lims, z_lims, titles):
        ax.set_title(title)
        ax.set_xlabel('Log X (cm)')

        size = (len(grid_x)-1, len(grid_z)-1)
        data = np.reshape(table['var'], size)
        inwind = np.reshape(table['inwind'], size)

        if volume:
            for xi in range(size[0]):
                for zi in range(size[1]):
                    # We need to correct the per-area emission to a per-volume emission
                    area = 2. * 2. * np.pi * (grid_x[xi+1]-grid_x[xi]) * (grid_z[zi+1]-grid_z[zi])
                    area += 2. * np.pi * grid_x[xi+1] * (grid_z[zi+1]-grid_z[zi])
                    area += 2. * np.pi * grid_x[xi]   * (grid_z[zi+1]-grid_z[zi])
                    if area > 0:
                        data[xi, zi] = data[xi, zi] / area

        if shared_cbar:
            if data.max() > v_max:
                v_max = data.max()
            if data.min() < v_min:
                v_min = data.min()

        if log:
            image = ax.pcolormesh(np.log10(grid_x), np.log10(grid_z),
                                  np.ma.log10(data.T))
        else:
            image = ax.pcolormesh(np.log10(grid_x), np.log10(grid_z),
                                  np.ma.masked_where(inwind.T < 0, data.T))

        images.append(image)
        ax.set_xlim(xlim[0], xlim[1])
        ax.set_ylim(zlim[0], zlim[1])

    if not shared_cbar:
        for image, ax in zip(images, axes):
            cbar = fig.colorbar(image, ax=ax).set_label(label)
    else:
        for image, ax in zip(images, axes):
            image.clim = [v_min, v_max]

        cbar = fig.colorbar(images[-1], ax=axes[-1]).set_label(label)

    fig.tight_layout()
    if shared_y:
        fig.subplots_adjust(wspace=0)
    return fig
