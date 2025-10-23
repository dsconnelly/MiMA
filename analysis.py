import sys

from typing import Optional

import cftime
import matplotlib.gridspec as gs
import matplotlib.pyplot as plt
import numpy as np
import scipy.signal as signal
import xarray as xr

from scipy import stats

_AMAX = 70
_N_YEARS = 25
_SCALE = 1.2

_CALENDAR = '360_day'
_UNITS = 'days since 0001-01-01'

_JET_PRESSURE = 10
_QBO_PRESSURE = 20
_TROPICS_LAT = 10

_COMP_LEVEL = 10
_WINDOW = 30

def plot_composites(ds: xr.Dataset) -> None:
    """Plot composites of a particular QBO phase."""

    widths = [2, 4.5, 0.2]
    fig, axes = plt.subplots(3, len(widths), width_ratios=widths)
    fig.set_size_inches(sum(widths), 3 * 3)
    raxes, paxes, caxes = axes.T

    vnames = ['u', 'gw_flux_x', 'gw_accel_x']
    datas = {vname : _get_qbo(ds, vname) for vname in vnames}

    crossings = _get_zero_crossings(datas['u'], level=_COMP_LEVEL)
    days = cftime.date2num(ds['time'], _UNITS, calendar=_CALENDAR)
    keep = (_WINDOW < crossings) & (crossings < days.max() - _WINDOW)
    crossings = crossings[keep]

    print(f'Compositing {keep.sum()} events')

    yticks = np.arange(-1, 4)
    fmt = lambda v: str(v if v < 1 else int(v))
    ylabels = list(map(fmt, 10. ** yticks))
    y = np.log10(ds['pfull'].values)

    composites = [0] * 3
    for crossing in crossings:
        idx = (crossing - _WINDOW <= days) & (days <= crossing + _WINDOW)
        x = days[idx] - crossing

        for i, data in enumerate(datas.values()):
            composites[i] = composites[i] + data.isel(time=idx).values
            
    for i, factor in enumerate([1, 1000, 86400]):
        composites[i] = factor * composites[i] / len(crossings)

    mid = (len(x) - 1) // 2
    units = ['m / s', 'mPa', 'm / s / day']
    amaxes = [75, 2, 6]

    tasks = zip(raxes, paxes, caxes, composites, units, amaxes)
    for k, (rax, pax, cax, data, unit, amax) in enumerate(tasks):
        curve = data[mid]
        rax.plot(curve, y, color='k')

        img = pax.pcolormesh(
            x, y, data.T,
            vmin=-amax, vmax=amax,
            shading='nearest',
            cmap='RdBu_r'
        )

        cbar = plt.colorbar(img, cax=cax)
        cbar.set_ticks(np.linspace(-amax, amax, 5))
        cbar.set_label(unit)

        for ax in (rax, pax):
            ax.invert_yaxis()
            ax.set_yticks(yticks)
            ax.set_ylim(2, -1)

        if k == 2:
            rax.set_xlim(-0.3 * amax, 0.3 * amax)
        else:
            rax.set_xlim(-amax, amax)

        rax.set_xlabel(f'lag-0 mean ({unit})')
        rax.set_ylabel('pressure (hPa)')
        rax.set_yticklabels(ylabels)

        rax.tick_params('both', direction='in')
        rax.grid(color='lightgray')

        pax.axhline(np.log10(_COMP_LEVEL), color='k', ls='dashed')
        pax.set_xlabel('lag (days)')
        pax.set_yticklabels([])

    plt.tight_layout()
    plt.savefig('plots/composites.png', dpi=400, bbox_inches='tight')

def plot_qbo(ds: xr.Dataset) -> None:
    """Plot the QBO from an integration."""

    widths = [2, 4.5, 0.2]
    fig, axes = plt.subplots(3, len(widths), width_ratios=widths)
    fig.set_size_inches(sum(widths), 3 * 3)
    raxes, paxes, caxes = axes.T

    vnames = ['u', 'gw_flux_x', 'gw_accel_x']
    units = ['m / s', 'mPa', 'm / s / day']
    factors = [1, 1000, 86400]
    amaxes = [75, 2, 6]

    yticks = np.arange(-1, 4)
    fmt = lambda v: str(v if v < 1 else int(v))
    ylabels = list(map(fmt, 10. ** yticks))
    y = np.log10(ds['pfull'].values)

    years = cftime.date2num(ds['time'], _UNITS, calendar=_CALENDAR) / 360

    tasks = zip(raxes, paxes, caxes, vnames, units, factors, amaxes)
    for k, (rax, pax, cax, vname, unit, factor, amax) in enumerate(tasks):
        data = factor * _get_qbo(ds, vname=vname)
        rms = np.sqrt((data ** 2).mean('time'))

        rax.plot(rms, y, color='k')

        img = pax.pcolormesh(
            years, y, data.values.T,
            vmin=-amax, vmax=amax,
            shading='nearest',
            cmap='RdBu_r'
        )

        cbar = plt.colorbar(img, cax=cax)
        cbar.set_ticks(np.linspace(-amax, amax, 5))
        cbar.set_label(unit)

        rax.set_xlabel(f'RMS ({unit})')

        pax.set_xlabel('years')
        xticks = np.arange(0, _N_YEARS + 1, 5)
        pax.set_xlim(0, _N_YEARS)
        pax.set_xticks(xticks)

        for ax in (rax, pax):
            ax.invert_yaxis()
            ax.set_yticks(yticks)

        rax.set_yticklabels(ylabels)
        rax.set_ylabel('pressure (hPa)')
        pax.set_yticklabels([])

        rax.tick_params('both', direction='in')
        rax.grid(color='lightgray')
        rax.set_xlim(0, amax)

        pax.set_xlim(15, 20)
        rax.set_ylim(2, -1)
        pax.set_ylim(2, -1)

    plt.tight_layout()
    plt.savefig('plots/qbo.png', dpi=400, bbox_inches='tight')

def plot_source_levels(_) -> None:
    """Plot gravity wave sources as a function of latitude."""

    with open('INPUT/source-levels.in') as f:
        parse_line = lambda line: np.array(list(map(float, line.split())))
        lats, levels, *_ = map(parse_line, f)

    x = np.linspace(-90, 90, 180)
    y = np.interp(x, lats, levels)

    fig, ax = plt.subplots()
    fig.set_size_inches(4.5, 3)
    ax.plot(x, y / 1000, color='k')

    ax.set_xlim(-90, 90)
    ax.set_xlabel('latitude')
    ticks = np.linspace(-90, 90, 7)
    ax.set_xticks(ticks, labels=list(map(_format_lat, ticks)), rotation=30)

    ax.set_ylim(10, 20)
    ax.set_ylabel('source level (km)')

    ax.tick_params('both', direction='in')
    ax.grid(color='lightgray')

    plt.tight_layout()
    plt.savefig('plots/source-levels.png', dpi=400)

def plot_summary(ds: xr.Dataset, plot_func: str='pcolormesh') -> None:
    """Plot zonal wind climatologies and the QBO time series."""

    buffer = 0.4
    heights = [buffer, 2.5 - buffer, 0.68, 2.5 - buffer, buffer]
    widths = [3, 3, 0.15]

    fig = plt.figure(constrained_layout=False)
    fig.set_size_inches(sum(widths), sum(heights))

    spec = gs.GridSpec(
        len(heights), len(widths), fig,
        height_ratios=heights,
        width_ratios=widths,
        hspace=0, wspace=0.25
    )

    axes = [fig.add_subplot(spec[:2, j]) for j in (0, -2)]
    axes = axes + [fig.add_subplot(spec[3:, :-1])]
    cax = fig.add_subplot(spec[1:4, -1])

    lat = np.linspace(-90, 90, len(ds['lat']))
    y = np.log10(ds['pfull'].values)

    yticks = np.arange(-1, 4)
    fmt = lambda v: str(v if v < 1 else int(v))
    ylabels = list(map(fmt, 10. ** yticks))

    xticks = np.linspace(-90, 90, 5)
    xlabels = list(map(_format_lat, xticks))

    kwargs = dict(vmin=-_AMAX, vmax=_AMAX, cmap='RdBu_r')
    if plot_func == 'pcolormesh': kwargs['shading'] = 'nearest'

    for ax, season in zip(axes, ['DJF', 'JJA']):
        u = _get_season(ds, season).mean('time')
        getattr(ax, plot_func)(lat, y, u.values, **kwargs)

        ax.set_xlim(-90, 90)
        ax.set_xticks(xticks, labels=xlabels, rotation=10)

    u = _get_qbo(ds)
    years = cftime.date2num(ds['time'], _UNITS, calendar=_CALENDAR) / 360
    img = getattr(axes[-1], plot_func)(years, y, u.values.T, **kwargs)

    crossings = _get_zero_crossings(u) / 360
    y_c = np.log10(_QBO_PRESSURE) * np.ones_like(crossings)
    axes[-1].scatter(crossings, y_c, color='k', marker='x', s=(_SCALE * 50))

    axes[-1].set_xlabel('years')
    xticks = np.arange(0, _N_YEARS + 1, 5)
    axes[-1].set_xlim(0, _N_YEARS)
    axes[-1].set_xticks(xticks)

    cbar = plt.colorbar(img, cax=cax)
    cbar.set_ticks(np.linspace(-_AMAX, _AMAX, 5))
    cbar.set_label('m / s')

    for i, ax in enumerate(axes):
        ax.set_title(f'({chr(i + 97)})')
        ax.set_yticks(yticks)
        ax.invert_yaxis()

        if i != 1:
            ax.set_yticklabels(ylabels)
            ax.set_ylabel('pressure (hPa)')

        else:
            ax.set_yticklabels([])

    axes[-1].set_ylim(2, -1)
    plt.savefig(f'plots/summary-{plot_func}.png', dpi=400, bbox_inches='tight')

def show_statistics(ds: xr.Dataset) -> None:
    """Print some statistics to do with stratospheric variability."""

    crossings = _get_zero_crossings(_get_qbo(ds))
    mean, error = _get_statistics(np.diff(crossings) / 30)
    print(np.diff(crossings) / 30)
    print(f'    QBO period is {mean:.3f} +- {error:.3f} months')

    for season in ['DJF', 'JJA']:
        sign = 1 if season == 'DJF' else -1
        kwargs = dict(pfull=_JET_PRESSURE, lat=(sign * 60), method='nearest')
        name = 'boreal' if season == 'DJF' else 'austral'

        u = _get_season(ds, season).sel(**kwargs)
        g = u.groupby('time.year').mean().values
        speed, error = _get_statistics(g)

        print(
            f'    {name} winter vortex speed is ' +
            f'{speed:.3f} +- {error:.3f} m / s'
        )

def show_total_time(_) -> None:
    with open('slurm.out') as f:
        total, percent, count = 0, 0, 0

        for line in f:
            if 'MS-GWaM' not in line:
                continue

            if 'IO' in line:
                continue

            parts = line.strip().split()
            total = total + float(parts[3])
            percent = percent + float(parts[5])
            count = count + 1

    mean = total / count / 60
    percent = 100 * percent / count
    print(f'    MS-GWaM took {mean:.3f} min / yr ({percent:.3f}% of runtime)')

def _format_lat(v: float) -> str:
    """Format a latitude for display."""

    if v == 0:
        return '0'

    suffix = 'N' if v > 0 else 'S'
    return f'{int(abs(v))}{suffix}'

def _get_statistics(data: np.ndarray, conf: float=0.95) -> tuple[float, float]:
    """Get the mean an confidence interval width."""

    mean, error = data.mean(), stats.sem(data)
    h = error * stats.t.ppf((1 + conf) / 2, len(data) - 1)
    return mean, h

def _get_qbo(ds: xr.Dataset, vname: str='u') -> xr.DataArray:
    """Extract the QBO time series from a loaded dataset."""

    tropics = abs(ds['lat']) <= _TROPICS_LAT
    return ds[vname].isel(lat=tropics).mean(('lat', 'lon'))

def _get_season(ds: xr.Dataset, season: str) -> xr.DataArray:
    """Get the zonal wind climatology for a given season."""

    spunup = ds['time.year'] > 1
    keep = (ds['time.season'] == season) & spunup
    
    return ds['u'].isel(time=keep).mean('lon')

def _get_zero_crossings(
    u: xr.DataArray,
    level: Optional[float]=None,
    exclude: list[int]=[3]
) -> np.ndarray:
    """Smooth and extract zero crossings from the QBO wind."""

    if level is None:
        level = _QBO_PRESSURE

    days = cftime.date2num(u['time'], _UNITS, calendar=_CALENDAR)
    data = u.sel(pfull=level, method='nearest').values

    # fs = 1 / (days[1] - days[0])
    # sos = signal.butter(9, 1 / 120, output='sos', fs=fs)
    # data = signal.sosfilt(sos, data)

    idx = (data[:-1] > 0) & (data[1:] < 0)
    keep = [i not in exclude for i in range(idx.sum())]

    return days[:-1][idx][keep]

def _open_dataset(fname: str) -> xr.Dataset:
    """Open a MiMA output file and parse the time dimension correctly."""

    with xr.open_dataset(fname, decode_times=False) as ds:
        days = cftime.num2date(ds['time'].values, _UNITS, calendar=_CALENDAR)
        days = (days - days.min()).astype('timedelta64[D]').astype(int)
        ds['time'] = cftime.num2date(days, _UNITS, calendar=_CALENDAR)

    return ds

if __name__ == '__main__':
    fname, *tasks = sys.argv[1:]
    ds = _open_dataset(fname)

    for task in tasks:
        func_name, *args = task.split(':')
        func_name = func_name.replace('-', '_')
        globals()[func_name](ds, *args)
