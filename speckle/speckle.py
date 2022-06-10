import os
import sys
import numpy as np
import warnings
from astropy.io import fits
from astropy.utils.exceptions import AstropyWarning
import matplotlib.pyplot as pl
from matplotlib.colors import LogNorm
from seaborn.cm import mako
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import glob
from matplotlib import rcParams
rcParams["savefig.dpi"] = 150
rcParams["figure.dpi"] = 150
if sys.platform == 'linux':
    pl.rcParams['font.family'] = 'Liberation Sans'
elif sys.platform == 'darwin':
    pl.rcParams['font.family'] = 'Arial'

warnings.simplefilter('ignore', AstropyWarning)
warnings.simplefilter('ignore', UserWarning)


def normalize_01(im):

    return (im - im.min())/(im.max() - im.min())


class Speckle:

    def __init__(self, name, data_dir, inst="NESSI"):

        self.name = name
        self.data_dir = data_dir
        self.inst = inst
        if inst == "NESSI":
            self.blue_name = 'b'
            self.red_name = 'r'
            self.blue_wav = 562
            self.red_wav = 832
            self.skiprows = 19
        elif inst == "DSSI":
            self.blue_name = 'a'
            self.red_name = 'b'
            self.blue_wav = 692
            self.red_wav = 880
            self.skiprows = 29
        elif inst == "Zorro":
            self.blue_name = 'b'
            self.red_name = 'r'
            self.blue_wav = 562
            self.red_wav = 832
            self.skiprows = 29
        elif inst == "Alopeke":
            self.blue_name = '562'
            self.red_name = '832'
            self.blue_wav = 562
            self.red_wav = 832
            self.skiprows = 29

        self._load(name, data_dir)

    def _load(self, name, data_dir):

        print("Loading {} data".format(name))

        fp = glob.glob(os.path.join(data_dir, '{}*{}.fits'.format(name, self.blue_name)))[0]
        hl = fits.open(fp)
        self._fits_b = hl
        self._im_b = normalize_01(hl[0].data)
        self._hdr_b = hl[0].header

        fp = glob.glob(os.path.join(data_dir, '{}*{}.fits'.format(name, self.red_name)))[0]
        hl = fits.open(fp)
        self._fits_r = hl
        self._im_r = normalize_01(hl[0].data)
        self._hdr_r = hl[0].header

        fp = glob.glob(os.path.join(data_dir, '{}*{}.dat'.format(name, self.blue_name)))[0]
        self._cc_b = np.loadtxt(fp, skiprows=self.skiprows)

        fp = glob.glob(os.path.join(data_dir, '{}*{}.dat'.format(name, self.red_name)))[0]
        self._cc_r = np.loadtxt(fp, skiprows=self.skiprows)

        print("Data taken on {}".format(self.obs_date))


    @property
    def obs_date(self):
        return self._hdr_r['DATE-OBS']


    def plot(self, figsize=(5,3.5), title=None, fp=None, stretch=1, vrange=None, cmap=None, c1='#01cdfe', c2='#ff71ce'):

        # colors = ['navy', 'turquoise', 'darkorange']
        colors = [c1, c2]

        if vrange is None:
            vrange = 0.1, 99.9

        if cmap is None:
            cm = mako
        elif cmap == 'gray':
            cm = pl.cm.gray
        elif cmap == 'gray_r':
            cm = pl.cm.gray_r

        fontsize=30
        pl.style.use('seaborn-ticks')

        # rcParams['font.family'] = 'serif'
        rcParams['axes.facecolor'] = 'white'
        rcParams["xtick.direction"] = 'in'
        rcParams["ytick.direction"] = 'in'

        if title is None:
            title = self.name

        fig,ax = pl.subplots(1, 1, figsize=figsize)

        rho, theta = self._cc_b.T
        ax.plot(rho, theta, label=f'{self.blue_wav} nm', color=colors[0], lw=2)

        rho, theta = self._cc_r.T
        ax.plot(rho, theta, label=f'{self.red_wav} nm', color=colors[1], lw=2)
        ax.invert_yaxis()
        ax.xaxis.set_ticks_position('both')
        ax.yaxis.set_ticks_position('both')
        ax.minorticks_on()

        ax1 = inset_axes(ax, 4.4, 1.3, borderpad=2)
        vmin, vmax = np.percentile(self._im_b, vrange[0]), np.percentile(self._im_b, vrange[1])
        im = ax1.imshow(self._im_b, cmap=cm, norm=LogNorm(vmin=vmin, vmax=vmax), interpolation='none')
        ax1.xaxis.set_visible(False)
        ax1.yaxis.set_visible(False)
        ax1.set_title(f'{self.blue_wav} nm', fontsize=10)
#         pixscale = 0.0175649 # arcsec/pixel
        pixscale = self._hdr_b['PIXSCL']
        arcsec = int(round(1/pixscale))
        xl, yl = ax1.get_xlim(), ax1.get_ylim()
        xcoord = [xl[1]-1.3*arcsec, xl[1]-0.3*arcsec]
        ycoord = [yl[0]+0.07*np.diff(yl), yl[0]+0.07*np.diff(yl)]
        ax1.plot(xcoord, ycoord, color='white')
        if self.inst == "NESSI":
            ax1.text(xcoord[0]*0.98, ycoord[0]*0.95, '1 arcsec', color='white', fontsize=6)
        elif self.inst == "DSSI":
            ax1.text(xcoord[0]*1.08, ycoord[0]*0.95, '1 arcsec', color='white', fontsize=6)
        elif self.inst == "Zorro":
            ax1.text(xcoord[0]*1.12, ycoord[0]*0.95, '1 arcsec', color='white', fontsize=6)
        elif self.inst == "Alopeke":
            ax1.text(xcoord[0]*1.12, ycoord[0]*0.95, '1 arcsec', color='white', fontsize=6)

        ax2 = inset_axes(ax, 1.3, 1.3, borderpad=2)
        vmin, vmax = np.percentile(self._im_r, vrange[0]), np.percentile(self._im_r, vrange[1])
        im = ax2.imshow(self._im_r, cmap=cm, norm=LogNorm(vmin=vmin, vmax=vmax), interpolation='none')
        ax2.xaxis.set_visible(False)
        ax2.yaxis.set_visible(False)
        ax2.set_title(f'{self.red_wav} nm', fontsize=10)
#         pixscale = 0.0181887 # arcsec/pixel
        pixscale = self._hdr_r['PIXSCL']
        arcsec = int(round(1/pixscale))
        xl, yl = ax2.get_xlim(), ax2.get_ylim()
        xcoord = [xl[1]-1.3*arcsec, xl[1]-0.3*arcsec]
        ycoord = [yl[0]+0.07*np.diff(yl), yl[0]+0.07*np.diff(yl)]
        ax2.plot(xcoord, ycoord, color='white')
        if self.inst == "NESSI":
            ax2.text(xcoord[0]*0.98, ycoord[0]*0.95, '1 arcsec', color='white', fontsize=6)
        elif self.inst == "DSSI":
            ax2.text(xcoord[0]*1.08, ycoord[0]*0.95, '1 arcsec', color='white', fontsize=6)
        elif self.inst == "Zorro":
            ax2.text(xcoord[0]*1.12, ycoord[0]*0.95, '1 arcsec', color='white', fontsize=6)
        elif self.inst == "Alopeke":
            ax2.text(xcoord[0]*1.12, ycoord[0]*0.95, '1 arcsec', color='white', fontsize=6)

        yl = ax.get_ylim()
        ylim = (stretch * yl[0], yl[1])
        pl.setp(ax,
                title=title,
                xlabel='Separation [arcsec]',
                ylabel=r'$\Delta$mag',
                xlim=(rho.min(), rho.max()),
                ylim=ylim)

        ax.legend(loc='lower left', frameon=False)

        if fp is None:
            fp = '{}_{}.png'.format(self.name, self.obs_date)

        fig.tight_layout()
        fig.savefig(fp, dpi=400)
        pl.close()

        print("Wrote file: {}".format(fp))


def cli():

    import sys
    import time
    tick = time.time()

    import argparse

    parser = argparse.ArgumentParser(description="Plot speckle imaging data")
    parser.add_argument('-i', '--id', help='ID', type=str, default=None)
    parser.add_argument('-d', '--data_dir', help='Directory containing the data products',
        type=str, default='.')
    parser.add_argument('-s', '--stretch', help='y-axis stretch factor (default=1)', type=float, default=1)
    parser.add_argument('-n', '--name', help='Target name', type=str, default=None)
    parser.add_argument('-v', '--vrange', help='Log stretch vrange', type=str, default=None)
    parser.add_argument('-c', '--cmap', help='Color map name', type=str, default=None)
    parser.add_argument('-f', '--figsize', help='Figure size (comma-separated)', type=str, default='5,3.5')
    parser.add_argument('--inst', help='Instrument name', type=str, default='NESSI')
    args = parser.parse_args()

    if args.id is None:
        sys.exit('Must supply ID')

    id_ = args.id
    data_dir = args.data_dir
    stretch = args.stretch
    name = args.name
    vrange = args.vrange
    cmap = args.cmap
    inst = args.inst
    figsize = [float(i) for i in args.figsize.split(',')]

    if vrange is not None:
        vrange = list(map(float,vrange.split(',')))

    spkl = Speckle(id_, data_dir, inst=inst)
    spkl.plot(figsize=figsize, title=name, stretch=stretch, vrange=vrange, cmap=cmap)

    print("Script executed in {0:.1f} seconds\n".format(time.time() - tick))
