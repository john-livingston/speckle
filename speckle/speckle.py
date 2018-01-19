import os
import numpy as np
import warnings
from astropy.io import fits
from astropy.utils.exceptions import AstropyWarning
import matplotlib.pyplot as pl
from matplotlib.colors import LogNorm
from seaborn.cm import mako as cm
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import glob

warnings.simplefilter('ignore', AstropyWarning)
warnings.simplefilter('ignore', UserWarning)


def normalize_01(im):

    return (im - im.min())/(im.max() - im.min())


class Speckle:

    def __init__(self, epic, data_dir):

        self.epic = epic
        self.data_dir = data_dir

        self._load(epic, data_dir)


    def _load(self, epic, data_dir):

        print("Loading data for EPIC-{}".format(epic))

        fp = glob.glob(os.path.join(data_dir, '{}*b.fits'.format(epic)))[0]
        hl = fits.open(fp)
        self._im_b = normalize_01(hl[0].data)
        self._hdr_b = hl[0].header

        fp = glob.glob(os.path.join(data_dir, '{}*r.fits'.format(epic)))[0]
        hl = fits.open(fp)
        self._im_r = normalize_01(hl[0].data)
        self._hdr_r = hl[0].header

        fp = glob.glob(os.path.join(data_dir, '{}*b.dat'.format(epic)))[0]
        self._cc_b = np.loadtxt(fp, skiprows=19)

        fp = glob.glob(os.path.join(data_dir, '{}*r.dat'.format(epic)))[0]
        self._cc_r = np.loadtxt(fp, skiprows=19)

        print("Data taken on {}".format(self.obs_date))


    @property
    def obs_date(self):
        return self._hdr_r['DATE-OBS']


    def plot(self, fp=None):

        fig,ax = pl.subplots(1, 1, figsize=(5,3.5), sharex=True, sharey=True)

        rho, theta = self._cc_b.T
        ax.plot(rho, theta, label='562 nm')

        rho, theta = self._cc_r.T
        ax.plot(rho, theta, label='832 nm')
        ax.invert_yaxis()

        ax1 = inset_axes(ax, 4.4, 1.3, borderpad=2)
        im = ax1.imshow(self._im_b, cmap=cm, norm=LogNorm(), interpolation='none')
        ax1.xaxis.set_visible(False)
        ax1.yaxis.set_visible(False)
        ax1.set_title('562 nm', fontsize=10)

        ax2 = inset_axes(ax, 1.3, 1.3, borderpad=2)
        im = ax2.imshow(self._im_r, cmap=cm, norm=LogNorm(), interpolation='none')
        ax2.xaxis.set_visible(False)
        ax2.yaxis.set_visible(False)
        ax2.set_title('832 nm', fontsize=10)

        pl.setp(ax,
                title="EPIC-{}".format(self.epic),
                xlabel='Separation [arcsec]',
                ylabel=r'$\Delta$mag',
                xlim=(rho.min(), rho.max()))

        ax.legend(loc='lower left')

        if fp is None:
            fp = 'EPIC-{}_{}.png'.format(self.epic, self.obs_date)

        fig.tight_layout()
        fig.savefig(fp, dpi=400)
        pl.close()

        print("Wrote file: {}".format(fp))
