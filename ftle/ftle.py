from pathlib import Path
import numpy as np
from scipy.signal import savgol_filter
from skimage.measure import find_contours
from scipy.interpolate import RegularGridInterpolator

# from lotusvis.plot_flow import Plots
from lotusvis.flow_field import ReadIn
from lotusvis.assign_props import AssignProps
import scienceplots

import seaborn as sns
from matplotlib import pyplot as plt
from ftle_plots import visualise_advection, plot_sketch, test

plt.style.use(["science"])
plt.rcParams["font.size"] = "10.5"


def interpolator(fs):
    x, y, t = fs.x, fs.y, fs.t.mean(axis=3)
    coords = (t[:, 0][:, 0], x, y)
    return RegularGridInterpolator(
        coords,
        (np.einsum("ikj", fs.vorticity_z.mean(axis=3))),
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
        


def fluid_snap(sim_dir):
    fsim = ReadIn(sim_dir, "fluid", 1024, ext="vti")
    fs = fsim.snaps(part=False, save=True)
    fs = AssignProps(fs)
    return fs


if __name__ == "__main__":
    sim_dir = f"./2-D/2048"
    sim_dir = f"/home/masseyj/Workspace/lotus_projects/mode-one-roughness/data/outer-scaling/12000/12-2d"
    
    fig_path = f"/home/masseyj/Workspace/presentations/DisCoVor-2023/figures/12k-12-2d"
    fig_path = "./2-D/analysis"

    fs = fluid_snap(sim_dir)
    # test(fs, interpolator)
    mwe()
