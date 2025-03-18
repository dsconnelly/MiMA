import sys

import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

def make_plot(fname: str) -> None:
    """
    Plot the wind
    """

    names = ['ucomp', 'vcomp', 'gw_flux_x', 'gw_flux_y']
    labels = ['zonal wind', 'meridional wind', 'zonal flux', 'meridional flux']
    amaxes = [25, 25, 10, 10]
    factors = [1, 1, 1e3, 1e3]

    fig, axes = plt.subplots(nrows=2, ncols=2)
    fig.set_size_inches(9, 6)

    zipped = zip(names, labels, amaxes, factors, axes.flatten())
    with xr.open_dataset(fname, decode_times=False) as ds:
        ds = ds.mean(('time', 'lon'))
        lat = np.linspace(-90, 90, len(ds['lat']))
        p = ds['pfull'].values

        for name, label, amax, factor, ax in zipped:
            data = factor * ds[name].values
            ax.set_title(label)
            ax.invert_yaxis()

            ax.pcolormesh(
                lat, p, data,
                vmin=-amax, vmax=amax,
                shading='nearest',
                cmap='RdBu_r'
            )

    plt.tight_layout()
    plt.savefig('clim.png', dpi=400)

if __name__ == '__main__':
    fname = sys.argv[1]
    make_plot(fname)