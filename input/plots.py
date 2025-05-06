from typing import Optional

import sys

import cftime
import matplotlib.gridspec as gs
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

from matplotlib.colors import SymLogNorm

def plot_climatology(fname: str, season: Optional[str]=None) -> None:
    """Plot zonal-mean climatologies in relevant fields."""

    factors = [1, 1000, 86400] * 2
    amaxes = [75, 5, 20, 10, 0.5, 20]
    units = ['m / s', 'mPa', 'm / s / d'] * 2
    names = ['ucomp', 'gw_flux_x', 'gw_accel_x',
        'vcomp', 'gw_flux_y', 'gw_accel_y']
    
    widths = [4.5, 0.1] * 3
    fig = plt.figure(constrained_layout=True)
    spec = gs.GridSpec(2, 6, fig, width_ratios=widths)
    fig.set_size_inches(sum(widths), 2 * 3)
    axes, caxes = [], []

    for i in range(2):
        for j in range(3):
            axes.append(fig.add_subplot(spec[i, 2 * j]))
            caxes.append(fig.add_subplot(spec[i, 2 * j + 1]))

    with _open_dataset(fname) as ds:
        keep = ds['time.year'] > 1
        if season is not None:
            keep = keep & (ds['time.season'] == season)
            
        ds = ds.isel(time=keep).mean(('time', 'lon'))
        lat = np.linspace(-90, 90, len(ds['lat']))
        y = np.log10(ds['pfull'].values)
        ticks, labels = _get_yticks()

        zipped = zip(names, factors, amaxes, units, axes)
        for k, (name, factor, amax, unit, ax) in enumerate(zipped):
            if name.startswith('gw_accel'):
                norm = SymLogNorm(1e-3, vmin=-amax, vmax=amax)
                extras = {'norm' : norm}

            else:
                extras = {'vmin' : -amax, 'vmax' : amax}

            img = ax.pcolormesh(
                lat, y, factor * ds[name].values,
                shading='nearest', cmap='RdBu_r',
                **extras
            )

            ax.set_xlim(-90, 90)
            ax.set_xlabel('latitude')

            if k % 3 == 0:
                ax.set_yticks(ticks)
                ax.set_yticklabels(labels)
                ax.set_ylabel('pressure (hPa)')

            else:
                ax.set_yticks([])
                ax.set_yticklabels([])

            ax.invert_yaxis()
            plt.colorbar(img, cax=caxes[k])
            ax.set_title(f'{name} ({unit})')

    suffix = season if season else 'annual'
    oname = fname.split('.')[0] + f'-clim-{suffix}.png'
    plt.savefig(oname, dpi=400)

def plot_pruning(fname: str) -> None:
    """Plot pruning as a function of time."""

    n_days, dt = 360, 240
    dt_resample = 3 * 3600
    n_steps = 360 * int(n_days)

    ns = np.zeros((32, n_steps))
    ages = np.zeros((32, n_steps))
    j = np.zeros(32).astype(int)

    with open(fname) as f:
        for line in f:
            if not 'pruned' in line:
                continue

            pe, *parts = line.split()
            pe = int(pe)

            if 'nothing' in line:
                j[pe] = j[pe] + 1
                continue

            ns[pe, j[pe]] = float(parts[1])
            ages[pe, j[pe]] = float(parts[6])
            j[pe] = j[pe] + 1

    n_resample = dt_resample // dt
    ages = ages.reshape(32, -1, n_resample)
    ns = ns.reshape(32, -1, n_resample)

    ages = _weighted_mean(ages, ns, axis=2)
    ns = ns.sum(axis=2)

    ages = _weighted_mean(ages, ns, axis=0)
    ns = ns.mean(axis=0)

    fig, axes = plt.subplots(ncols=2)
    fig.set_size_inches(9, 3)

    days = np.linspace(0, n_days, len(ns))
    axes[1].plot(days, ages / 86400, color='k')
    axes[0].plot(days, ns, color='k')

    for ax in axes:
        ax.set_xlim(days.min(), days.max())
        ax.set_xlabel('integration day')

    axes[0].set_ylim(0, 20)
    axes[1].set_ylim(0, 3)

    axes[0].set_ylabel(f'prunes per {dt_resample // 3600} hours')
    axes[1].set_ylabel('age of pruned rays (days)')

    plt.tight_layout()
    plt.savefig('pruning.png', dpi=400)

def plot_qbo(fname: str) -> None:
    """Plot the QBO wind from a MiMA run."""

    widths = [4.5, 0.2]
    fig, (ax, cax) = plt.subplots(ncols=2, width_ratios=widths)
    fig.set_size_inches(sum(widths), 3)

    with _open_dataset(fname) as ds:
        spunup = ds['time.year'] > 2
        tropics = abs(ds['lat']) <= 15

        ds = ds.isel(
            time=spunup,
            lat=tropics
        ).mean(('lat', 'lon'))

        y = np.log10(ds['pfull'].values)
        ticks, labels = _get_yticks()

        time = cftime.date2num(
            ds['time'].values,
            'days since 0001-01-01'
        ) / 360

        img = ax.pcolormesh(
            time, y, ds['ucomp'].values.T,
            vmin=-50, vmax=50,
            shading='nearest',
            cmap='RdBu_r'
        )

        ax.set_xticks(np.linspace(time.min(), time.max(), 4))
        ax.set_xlim(time.min(), time.max())
        ax.set_xlabel('years')

        cbar = plt.colorbar(img, cax=cax)
        cbar.set_label('m / s')

        ax.set_yticks(ticks)
        ax.set_yticklabels(labels)
        ax.set_ylabel('pressure (hPa)')

        ax.set_title('tropical-mean $u$')
        ax.invert_yaxis()

    plt.tight_layout()
    oname = fname.split('.')[0] + '-qbo.png'
    plt.savefig(oname, dpi=400)

def _get_yticks() -> tuple[np.ndarray, list[float]]:
    """Get the log pressure coordinate and appropriate labels."""

    ticks = np.array([-1, 0, 1, 2, 3])
    labels = [float(f'{x:.2g}') for x in (10. ** ticks)]

    return ticks, labels
    
def _open_dataset(fname: str) -> xr.Dataset:
    """Utility for parsing MiMA time dimensions."""

    with xr.open_dataset(fname, decode_times=False) as ds:
        ds['time'] = cftime.num2date(
            ds['time'].values,
            units='days since 0001-01-01',
            calendar='360_day'
        )

        return ds
    
def _weighted_mean(a: np.ndarray, w: np.ndarray, axis: int) -> np.ndarray:
    """Take a weighted mean along a given axis."""

    w_sum = w.sum(axis=axis)
    
    return np.divide(
        (a * w).sum(axis=axis), w_sum,
        out = np.zeros_like(w_sum),
        where=(w_sum > 0)
    )

if __name__ == '__main__':
    task, *args = sys.argv[1:]
    func = globals()[f'plot_{task}']
    func(*args)