import os
import time
import numpy as np
from tkinter import Tcl

from lotusvis.flow_field import ReadIn


def fluid_snap(sim_dir, fn, count):
    fsim = ReadIn(sim_dir, 'fluid', 4096, ext="vti")
    fsim.u_low_memory_saver(fn, count, "./data")
    fsim.v_low_memory_saver(fn, count, "./data")


def body_snap(sim_dir, fn, count):
    bsim = ReadIn(sim_dir, 'bodyF', 4096, ext="vti")
    bsim.save_sdf_low_memory(fn, count, "./data")


def bmask(count):
    fnsu, fnsv, fnsb = fns()
    for idx, (fnu, fnv, fnb) in enumerate(zip(fnsu, fnsv, fnsb)):
        u = np.load(os.path.join("./data", fnu))
        v = np.load(os.path.join("./data", fnv))
        b = np.load(os.path.join("./data", fnb))
        bmask = np.where(b <= 1, False, True)
        u = np.where(bmask, u, 0)
        print(u.shape)
        np.save(os.path.join("./data", f"u_{count}"), u)
        v = np.where(bmask, v, 0)
        np.save(os.path.join("./data", f"v_{count}"), v)
        count += 1
        # Now remove the files
        print(f"Removing {fnu}, {fnv}, {fnb}")
        try:
            os.remove(os.path.join("./data", fnu))
            os.remove(os.path.join("./data", fnv))
            os.remove(os.path.join("./data", fnb))
        except FileNotFoundError:
            pass
    return count


def fns():
    fnsu = [
        fn
        for fn in os.listdir("./data")
        if fn.startswith("fluid_u") and fn.endswith(f".npy")
    ]
    fnsu = Tcl().call("lsort", "-dict", fnsu)
    fnsv = [
        fn
        for fn in os.listdir("./data")
        if fn.startswith("fluid_v") and fn.endswith(f".npy")
    ]
    fnsv = Tcl().call("lsort", "-dict", fnsv)
    fnsb = [
        fn
        for fn in os.listdir("./data")
        if fn.startswith('bodyF') and fn.endswith(f".npy")
    ]
    fnsb = Tcl().call("lsort", "-dict", fnsb)
    return fnsu, fnsv, fnsb

# Specify the directory to monitor
directory_to_watch = "./lotus-data"


while True:
    count=0
    delete_count = 0
    for root, _, files in os.walk(directory_to_watch):
        # Process files
        for file in files:
            print("File created:", os.path.join(root, file))
            if root.endswith("datp"):
                # Buffer to allow finish writing
                time.sleep(0.1)
                # Sort the files
                dpdfs = [fp for fp in os.listdir(root)]
                dpdfs = Tcl().call("lsort", "-dict", dpdfs)
                for fn in dpdfs:
                    if fn.startswith("fluid") and fn.endswith(".pvti"):
                        path = os.path.join(root, fn)
                        fluid_snap(directory_to_watch, path, count)
                    if fn.startswith("bodyF") and fn.endswith(".pvti"):
                        path = os.path.join(root, fn)
                        body_snap(directory_to_watch, path, count)
                        count += 1
    for root, _, files in os.walk(directory_to_watch):
        for file in files:
            if (file.startswith("fluid") or file.startswith("bodyF")) and \
            (not file.endswith(".pvtr") and not file.endswith(".vtr") and not file.endswith("vtr.pvd")):
                file_path = os.path.join(root, file)
                os.remove(file_path)
                delete_count += 1
    bmask(count)

    print(f"Total files deleted: {delete_count}")

    # Sleep for a certain duration before checking again
    time.sleep(0.1)
