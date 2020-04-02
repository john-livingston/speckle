import os
import sys
import numpy as np
import warnings
from astropy.io import fits
from astropy.utils.exceptions import AstropyWarning
import matplotlib.pyplot as pl
from matplotlib.colors import LogNorm
from seaborn.cm import mako as cm
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

    def __init__(self, name, data_dir):

        self.name = name
        self.data_dir = data_dir

        self._load(name, data_dir)


    def _load(self, name, data_dir):

        print("Loading {} data".format(name))

        fp = glob.glob(os.path.join(data_dir, '{}*b.fits'.format(name)))[0]
        hl = fits.open(fp)
        self._fits_b = hl
        self._im_b = normalize_01(hl[0].data)
        self._hdr_b = hl[0].header

        fp = glob.glob(os.path.join(data_dir, '{}*r.fits'.format(name)))[0]
        hl = fits.open(fp)
        self._fits_r = hl
        self._im_r = normalize_01(hl[0].data)
        self._hdr_r = hl[0].header

        fp = glob.glob(os.path.join(data_dir, '{}*b.dat'.format(name)))[0]
        self._cc_b = np.loadtxt(fp, skiprows=19)

        fp = glob.glob(os.path.join(data_dir, '{}*r.dat'.format(name)))[0]
        self._cc_r = np.loadtxt(fp, skiprows=19)

        print("Data taken on {}".format(self.obs_date))


    @property
    def obs_date(self):
        return self._hdr_r['DATE-OBS']


    def plot(self, title=None, fp=None, stretch=1):

        # colors = ['navy', 'turquoise', 'darkorange']
        colors = ['#01cdfe', '#ff71ce']

        fontsize=30
        pl.style.use('seaborn-ticks')

        # rcParams['font.family'] = 'serif'
        rcParams['axes.facecolor'] = 'white'
        rcParams["xtick.direction"] = 'in'
        rcParams["ytick.direction"] = 'in'

        if title is None:
            title = self.name

        fig,ax = pl.subplots(1, 1, figsize=(5,3.5), sharex=True, sharey=True)

        rho, theta = self._cc_b.T
        ax.plot(rho, theta, label='562 nm', color=colors[0], lw=2)

        rho, theta = self._cc_r.T
        ax.plot(rho, theta, label='832 nm', color=colors[1], lw=2)
        ax.invert_yaxis()
        ax.xaxis.set_ticks_position('both')
        ax.yaxis.set_ticks_position('both')
        ax.minorticks_on()

        ax1 = inset_axes(ax, 4.4, 1.3, borderpad=2)
        vmin, vmax = np.percentile(self._im_b, 0.1), np.percentile(self._im_b, 99.9)
        im = ax1.imshow(self._im_b, cmap=cm, norm=LogNorm(vmin=vmin, vmax=vmax), interpolation='none')
        ax1.xaxis.set_visible(False)
        ax1.yaxis.set_visible(False)
        ax1.set_title('562 nm', fontsize=10)
        pixscale = 0.0175649 # arcsec/pixel
        arcsec = int(round(1/pixscale))
        xl, yl = ax1.get_xlim(), ax1.get_ylim()
        xcoord = [xl[1]-1.3*arcsec, xl[1]-0.3*arcsec]
        ycoord = [yl[0]+0.07*np.diff(yl), yl[0]+0.07*np.diff(yl)]
        ax1.plot(xcoord, ycoord, color='white')
        ax1.text(xcoord[0]*0.97, ycoord[0]*0.95, '1 arcsec', color='white', fontsize=6)

        ax2 = inset_axes(ax, 1.3, 1.3, borderpad=2)
        vmin, vmax = np.percentile(self._im_r, 0.1), np.percentile(self._im_r, 99.9)
        im = ax2.imshow(self._im_r, cmap=cm, norm=LogNorm(vmin=vmin, vmax=vmax), interpolation='none')
        ax2.xaxis.set_visible(False)
        ax2.yaxis.set_visible(False)
        ax2.set_title('832 nm', fontsize=10)
        pixscale = 0.0181887 # arcsec/pixel
        arcsec = int(round(1/pixscale))
        xl, yl = ax2.get_xlim(), ax2.get_ylim()
        xcoord = [xl[1]-1.3*arcsec, xl[1]-0.3*arcsec]
        ycoord = [yl[0]+0.07*np.diff(yl), yl[0]+0.07*np.diff(yl)]
        ax2.plot(xcoord, ycoord, color='white')
        ax2.text(xcoord[0]*0.97, ycoord[0]*0.95, '1 arcsec', color='white', fontsize=6)

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
