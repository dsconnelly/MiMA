from typing import Optional

import sys

import cftime
import numpy as np
import xarray as xr

from matplotlib import gridspec as gs, pyplot as plt
from matplotlib.colors import SymLogNorm

_DT = 240
_CALENDAR = '360_day'
_UNITS = 'days since 0001-01-01'

def plot_climatology(ds: xr.Dataset, season: Optional[str]=None) -> None:
    """Plot zonal-mean climatologies in relevant fields."""

    factors = [1, 1000, 86400] * 2
    amaxes = [75, 5, 20, 10, 0.5, 20]
    units = ['m / s', 'mPa', 'm / s / d'] * 2

    names = []
    for comp, coord in zip('uv', 'xy'):
        names.append(_get_uv_name(ds, comp))
        names.extend([f'gw_flux_{coord}', f'gw_accel_{coord}'])

    keep = [k for k, name in enumerate(names) if name in ds]
    names = list(map(names.__getitem__, keep))

    factors = list(map(factors.__getitem__, keep))
    amaxes = list(map(amaxes.__getitem__, keep))
    units = list(map(units.__getitem__, keep))
    
    n_cols = len(names) // 2
    widths = [4.5, 0.1] * n_cols
    fig = plt.figure(constrained_layout=True)
    spec = gs.GridSpec(2, n_cols * 2, fig, width_ratios=widths)
    fig.set_size_inches(sum(widths), 2 * 3)
    axes, caxes = [], []

    for i in range(2):
        for j in range(n_cols):
            axes.append(fig.add_subplot(spec[i, 2 * j]))
            caxes.append(fig.add_subplot(spec[i, 2 * j + 1]))

    lat = np.linspace(-90, 90, len(ds['lat']))
    xticks = np.linspace(-90, 90, 7)
    y = np.log10(ds['pfull'].values)
    yticks, ylabels = _get_yticks()

    keep = ds['time.year'] > 1
    if season is not None:
        keep = keep & (ds['time.season'] == season)

    ds = ds.isel(time=keep).mean(('time', 'lon'))
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
        ax.set_xticks(xticks)
        ax.set_xlabel('latitude')

        if k % 3 == 0:
            ax.set_yticks(yticks)
            ax.set_yticklabels(ylabels)
            ax.set_ylabel('pressure (hPa)')

        else:
            ax.set_yticks([])
            ax.set_yticklabels([])

        ax.invert_yaxis()
        plt.colorbar(img, cax=caxes[k])
        ax.set_title(f'{name} ({unit})')

    suffix = season if season else 'annual'
    plt.savefig(f'plots/climatology-{suffix}.png', dpi=400)

def plot_pruning(_, fname: str) -> None:
    """Plot pruning as a function of time."""

    n_days = 30
    dt_resample = 3 * 3600
    n_steps = (86400 // _DT) * n_days

    ns = np.zeros((32, n_steps))
    ages = np.zeros((32, n_steps))
    j = np.zeros(32).astype(int)

    with open(fname) as f:
        for line in f:
            if not 'pruned' in line:
                continue

            pe, *parts = line.split()
            pe = int(pe)

            if not ('nothing' in line):
                ns[pe, j[pe]] = float(parts[1])
                ages[pe, j[pe]] = float(parts[6])

            j[pe] = j[pe] + 1

    n_resample = dt_resample // _DT
    ages = ages.reshape(32, -1, n_resample)
    ns = ns.reshape(32, -1, n_resample)

    ages = _weighted_mean(ages, ns, 2)
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

    axes[0].set_ylim(0, 200)
    axes[1].set_ylim(0, 3)

    axes[0].set_ylabel(f'prunes per {dt_resample // 3600} hours')
    axes[1].set_ylabel('age of pruned rays (days)')

    plt.tight_layout()
    plt.savefig('plots/pruning.png', dpi=400)

def plot_qbo(ds: xr.Dataset) -> None:
    """Plot the QBO wind time series."""

    widths = [4.5, 0.2]
    fig, (ax, cax) = plt.subplots(ncols=2, width_ratios=widths)
    fig.set_size_inches(sum(widths), 3)

    tropics = abs(ds['lat']) <= 15
    ds = ds.isel(lat=tropics).mean(('lat', 'lon'))
    u = ds[_get_uv_name(ds)]

    years = cftime.date2num(ds['time'].values, _UNITS, calendar=_CALENDAR) / 360
    y = np.log10(ds['pfull'].values)
    yticks, ylabels = _get_yticks()
    
    img = ax.pcolormesh(
        years, y, u.values.T,
        vmin=-50, vmax=50,
        shading='nearest',
        cmap='RdBu_r'
    )

    signal = u.sel(pfull=10, method='nearest').values
    crossings = _get_zero_crossings(signal, years * 360) / 360
    ax.scatter(crossings, np.ones_like(crossings), color='k', marker='x')

    n_years = int(years.max()) + 1
    ax.set_xticks(np.linspace(0, n_years, n_years + 1))
    ax.set_xlim(0, n_years)
    ax.set_xlabel('years')

    ax.set_yticks(yticks)
    ax.set_yticklabels(ylabels)
    ax.set_ylabel('pressure (hPa)')

    cbar = plt.colorbar(img, cax=cax)
    cbar.set_label('m / s')

    ax.invert_yaxis()
    ax.set_title('tropical-mean $u$')

    plt.tight_layout()
    plt.savefig('plots/qbo.png', dpi=400)

def show_statistics(ds: xr.Dataset) -> None:
    """Print some statistics about stratospheric variability."""

    u = ds[_get_uv_name(ds)].sel(pfull=10, method='nearest')

    tropics = abs(ds['lat']) < 15
    u_qbo = u.isel(lat=tropics).mean(('lat', 'lon')).values
    days = cftime.date2num(u['time'], _UNITS, calendar=_CALENDAR)

    crossings = _get_zero_crossings(u_qbo, days)
    period = np.diff(crossings).mean() / 360
    print(f'    QBO period is {period:.3f} years')

    for season in ['DJF', 'JJA']:
        keep = (ds['time.season'] == season)
        keep = keep & (ds['time.year'] > 1)
        u_s = u.isel(time=keep)

        sign = 1 if season == 'DJF' else -1
        u_s = u_s.sel(lat=(sign * 60), method='nearest')
        u_s = u_s.mean(('lon', 'time')).item()

        name = 'boreal' if season == 'DJF' else 'austral'
        print(f'    {name} winter vortex speed is {u_s:.3f} m / s')

def _get_uv_name(ds: xr.Dataset, comp: str='u') -> str:
    """Get the zonal component of the wind from a file."""

    return [s for s in [f'{comp}comp', f'{comp}_gwf'] if s in ds][0]

def _get_yticks() -> tuple[np.ndarray, np.ndarray]:
    """Get log pressure coordinate ticks and labels for plots."""

    ticks = np.arange(-1, 4)
    labels = [float(f'{x:.2g}') for x in 10. ** ticks]

    return ticks, np.array(labels)

def _get_zero_crossings(u: np.ndarray, days: np.ndarray) -> np.ndarray:
    """Get appropriately smoothed zero crossings from the QBO wind."""

    u_hat = np.fft.rfft(u)
    freqs = np.fft.rfftfreq(len(u), days[1] - days[0])
    u_hat[freqs > 1 / 30] = 0
    u = np.fft.irfft(u_hat)

    return days[:-1][(u[:-1] > 0) & (u[1:] < 0)]

def _open_dataset(fname: str) -> xr.Dataset:
    """Open a MiMA output file and parse the time dimension correctly."""

    with xr.open_dataset(fname, decode_times=False) as ds:
        days = cftime.num2date(ds['time'].values, _UNITS, calendar=_CALENDAR)
        days = (days - days.min()).astype('timedelta64[D]').astype(int)
        ds['time'] = cftime.num2date(days, _UNITS, calendar=_CALENDAR)

    return ds

def _weighted_mean(a: np.ndarray, w: np.ndarray, axis: int) -> np.ndarray:
    """Take a weighted mean along a given axis."""

    w_sum = w.sum(axis=axis)

    return np.divide(
        (a * w).sum(axis=axis), w_sum,
        out=np.zeros_like(w_sum),
        where=(w_sum > 0)
    )

if __name__ == '__main__':
    fname, *tasks = sys.argv[1:]
    with _open_dataset(fname) as ds:
        for task in tasks:
            func_name, *args = task.split(':')
            globals()[func_name.replace('-', '_')](ds, *args)