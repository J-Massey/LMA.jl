from pathlib import Path
import numpy as np
import os
from tkinter import Tcl

from lotusvis.flow_field import ReadIn
from lotusvis.assign_props import AssignProps

import os
import time
from watchdog.observers import Observer
from watchdog.events import FileSystemEventHandler


def fluid_snap(sim_dir, fn, count):
    fsim = ReadIn(sim_dir, fluid_ext, 4096, ext="vti")
    fsim.u_low_memory_saver(fn, count, save_path="./data")
    fsim.v_low_memory_saver(fn, count, save_path="./data")


def body_snap(sim_dir, fn, count):
    bsim = ReadIn(sim_dir, body_ext, 4096, ext="vti")
    bsim.save_sdf_low_memory(fn, count, save_path="./data")


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
        if fn.startswith(body_ext) and fn.endswith(f".npy")
    ]
    fnsb = Tcl().call("lsort", "-dict", fnsb)
    return fnsu, fnsv, fnsb


class FileHandler(FileSystemEventHandler):
    def __init__(self, parent_dir, file_extension):
        super(FileHandler, self).__init__()
        self.parent_dir = parent_dir
        self.count = 0
        self.delete_count = 0

    def on_created(self, event):
        if not event.is_directory:
            for dirpath, dirnames, filenames in os.walk(self.parent_dir):
                # for filename in filenames:
                if dirpath.endswith("datp"):
                    # Buffer to allow finish writing
                    time.sleep(0.1)
                    # Sort the files
                    dpdfs = [fp for fp in os.listdir(dirpath)]
                    dpdfs = Tcl().call("lsort", "-dict", dpdfs)
                    for fn in dpdfs:
                        if fn.startswith("fluid") and fn.endswith(".pvti"):
                            path = os.path.join(dirpath, fn)
                            fluid_snap(sim_dir, path, self.count)
                        if fn.startswith("bodyF") and fn.endswith(".pvti"):
                            path = os.path.join(dirpath, fn)
                            body_snap(sim_dir, path, self.count)
                            self.count += 1
            for dirpath, dirnames, filenames in os.walk(self.parent_dir):
                for filename in filenames:
                    if (filename.startswith("fluid") or filename.startswith("bodyF")) and \
                    (not filename.endswith(".pvtr") and not filename.endswith(".vtr") and not filename.endswith("vtr.pvd")):
                        file_path = os.path.join(dirpath, filename)
                        os.remove(file_path)
                        self.delete_count += 1
            bmask(self.count)

            print(f"Total files deleted: {self.delete_count}")


def watch_and_delete(parent_dir, file_extension):
    event_handler = FileHandler(parent_dir, file_extension)
    observer = Observer()
    observer.schedule(event_handler, parent_dir, recursive=False)
    observer.start()

    try:
        while True:
            time.sleep(0.1)
    except KeyboardInterrupt:
        observer.stop()

    observer.join()


if __name__ == "__main__":
    sim_dir = f"./lotus-data"
    fluid_ext = "fluid"
    body_ext = "bodyF"
    save_ext = "smooth"
    watch_and_delete(sim_dir, "bodyF.")
