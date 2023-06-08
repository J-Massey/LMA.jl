import os
from pathlib import Path
from tkinter import Tcl
import numpy as np
from scipy.interpolate import RegularGridInterpolator
from tqdm import tqdm

import seaborn as sns
from matplotlib import pyplot as plt
from ftle_plots import visualise_advection, plot_sketch, test

plt.style.use(["science"])
plt.rcParams["font.size"] = "10.5"


def interpolator(x, y, t, dat):
    coords = (t, x, y)
    return RegularGridInterpolator(
        coords,
        (dat),
        method="linear",
        bounds_error=False,
        fill_value=None,
    )


def Jacobian(X, Y, dx, dy):
    # shapes
    nx, ny = X.shape

    # storage arrays
    J = np.empty([2, 2], float)
    FTLE = np.empty([nx, ny], float)

    for i in range(1, nx - 1):
        for j in range(1, ny - 1):
            J[0, 0] = (X[i + 1, j] - X[i - 1, j]) / (2 * dx)
            J[0, 1] = (X[i, j + 1] - X[i, j - 1]) / (2 * dy)
            J[1, 0] = (Y[i + 1, j] - Y[i - 1, j]) / (2 * dx)
            J[1, 1] = (Y[i, j + 1] - Y[i, j - 1]) / (2 * dy)

            # Green-Cauchy tensor
            C = np.dot(np.transpose(J), J)
            # Largest eigenvalue
            lamda = np.linalg.eigvals(C)
            FTLE[i, j] = max(lamda)
    return FTLE


def compute_ftle(interpolator, r, grid, t0, t1, dt):
    """
    Compute the finite time Lyoponov Exponent.
    Parameters
        interpolator : regular grid interpolator to return
                       value at any pos(t,X,Y)
        r            : 2D array of flattened X, Y coords.
        shape        : shape of particle seed
        grid         : [dx, dy] of particle mesh.
        t0           : start time
        t1           : end time
        dt           : time step for time integration
    """

    # FORWARD vs BACKWARD
    sign = 1
    if t1 < t0:
        sign = -1

    # integrate in time
    for time in np.arange(t0, t1, sign * dt):
        # visualise_advection(r[0][::8], r[1][::8], time, "./2-D/analysis/figures/advection")
        k1 = dt * interpolator((time, r[0], r[1]))
        k2 = dt * interpolator((time + 0.5 * dt, r[0] + 0.5 * k1, r[1] + 0.5 * k1))
        k3 = dt * interpolator((time + 0.5 * dt, r[0] + 0.5 * k2, r[1] + 0.5 * k2))
        k4 = dt * interpolator((time + dt, r[0] + k3, r[1] + k3))
        r += sign * (k1 + 2 * k2 + 2 * k3 + k4) / 6
        
    FTLE = Jacobian(r[0], r[0], grid[0], grid[1])

    return np.log(np.sqrt(FTLE[1:-1, 1:-1]))


def mwe():
    # mesh of prticles
    # desired resolution
    w, h = (1080, 920)
    # w, h = w // 4, h // 4
    x, dx = np.linspace(-0.2, 2.0, w, retstep=True)
    y, dy = np.linspace(-0.5, 0.5, h, retstep=True)
    X, Y = np.meshgrid(x, y)

    # build grid
    interp = interpolator(fs)

    # move particles once
    t0, t1, dt = 19, 18, 1

    # time increment
    n_snaps = 10
    tau = np.linspace(0, 19, n_snaps)
    for n in range(n_snaps):
        # integrate particle in time from various initial condition
        FTLE = compute_ftle(
            interp,
            np.array([X, Y]),
            [dx, dy],
            t0 - tau[n],
            t1 - tau[n],
            dt,
        )

        plot_sketch(x[1:-1], y[1:-1], np.clip(FTLE, 0, 4), f"./2-D/analysis/figures/ftle{n}.png")


def fns(data_dir):
    fnsv = [fn for fn in os.listdir(data_dir)
            if fn.startswith('smooth') and fn.endswith(f'.npy')]
    fnsv = Tcl().call('lsort', '-dict', fnsv)
    return fnsv


def collect_data(fns):
    resize_shape = np.load(f"{data_dir}/example_shape.npy")
    resize_shape = np.shape(resize_shape.squeeze())
    data = []
    for fn in tqdm(fns, desc="Loading data"):
        snap = np.load(f"{data_dir}/{fn}").squeeze()
        snap = np.resize(snap, resize_shape)
        data.append(snap)
    return np.array(data).squeeze()


def generate_grid(xlims, ylims):
    snap = np.load(f"{data_dir}/example_shape.npy")
    fsize = np.shape(snap.squeeze())
    x = np.linspace(*xlims, fsize[0])
    y = np.linspace(*ylims, fsize[1])
    # X, Y = np.meshgrid(x, y)
    t = np.linspace(0,4,len(fns(data_dir)))
    return x,y, t


if __name__ == "__main__":
    data_dir = "../vort-smooth-binary"
    fig_path = "./figures"

    xlims = (-0.2, 2.0)
    ylims = (-0.5, 0.5)
    x, y, t = generate_grid(xlims, ylims)
    dat = np.load(f"{data_dir}/flow-binary.npy")
    interp = interpolator(x,y,t,dat)
    print('success')

    x = np.linspace(-0.5, 2, 300)
    y = np.linspace(-0.5, 0.5, 400)
    xg, yg = np.meshgrid(x, y)
    newu = interp((5.25, xg, yg))

    fig, ax = plt.subplots(figsize=(3,3))
    cmap = sns.color_palette("seismic", as_cmap=True)

    cs = ax.contourf(xg, yg, newu,
            levels=np.linspace(-0.2,0.2,88),
            cmap=cmap,
            extend='both',
            )
    
    ax.set_aspect(1)

    plt.savefig(f'test.png', dpi=300)
    plt.close()

    # dat = collect_data(fns(data_dir))
    # print(dat.shape)
    # np.save(f"{data_dir}/flow-binary.npy", dat)

    # mwe()
