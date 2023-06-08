from pathlib import Path
import numpy as np

import scienceplots

from matplotlib import colors, pyplot as plt
import seaborn as sns
from mpl_toolkits.axes_grid1 import make_axes_locatable
import scienceplots
import matplotlib.colors as culs

plt.style.use(["science"])
plt.rcParams["font.size"] = "10.5"


def visualise_advection(x,y,time,fig_path):
    fig, ax = plt.subplots(figsize=(3,3))

    ax.scatter(x, y, color='k', s=0.05)
    ax.set_aspect(1)

    plt.savefig(f'{fig_path}/particles{time}.png', dpi=300)
    plt.close()


def plot_sketch(x,y,c, fig_string):
    fig, ax = plt.subplots(figsize=(3,3))
    cmap = sns.color_palette("Blues_r", as_cmap=True)

    lim = [1, 4]
    norm = culs.LogNorm(vmin=lim[0], vmax=lim[1])
    levels = np.linspace(lim[0], lim[1], 11)

    cs = ax.contourf(x, y, c,
            levels=levels,
            cmap=cmap,
            norm=norm,
            vmax=lim[1],
            vmin=lim[0],
            # extend='upper',
            )
    # divider = make_axes_locatable(ax)
    # cax = divider.append_axes('top', size="7%", pad=0.2)
    # fig.add_axes(cax)
    # cb = plt.colorbar(cs, cax=cax, orientation="horizontal")
    # cax.set_title(r"$c_{p}$")
    ax.set_aspect(1)

    plt.savefig(f'{fig_string}', dpi=300)
    plt.close()

def test(interpolator, fig_string='./figures/test.png'):
    x = np.linspace(-0.25, 1.75, 300)
    y = np.linspace(-0.25, 0.25, 400)
    xg, yg = np.meshgrid(x, y)

    newu = interpolator((1.25, xg, yg))

    fig, ax = plt.subplots(figsize=(3,3))
    cmap = sns.color_palette("seismic", as_cmap=True)

    cs = ax.contourf(xg, yg, newu,
            levels=np.linspace(-0.1,0.1,88),
            cmap=cmap,
            extend='both',
            )
    
    ax.set_aspect(1)

    plt.savefig(f'{fig_string}', dpi=300)
    plt.close()

