
#!/usr/bin/env python3
# -*- coding: utf-8 -*-


def new_f90_res(n: int, re: float, k_lam: int, L: int):
    with open("lotus.f90","r") as fileSource:
        fileLines = fileSource.readlines()

    fileLines[11] = f"    real,parameter     :: Re = {re}\n"
    fileLines[13] = f"    real,parameter     :: L={L}, nu=L/Re\n"

    fileLines[19] = f"    real, parameter    :: A = 0.1*L, St_d = 0.3, k_x=0., k_z={float(k_lam)}, h_roughness=0.00\n"

    fileLines[23] = f"                          k_coeff = {1/0.9}, &\n"
    fileLines[27] = f"    integer            :: n(3), ndims={n}\n"

    # setup grid
    if k_lam % 2 ==0:
        fileLines[55] = f"      if(ndims==3) xg(3)%h = 4.\n"
    elif k_lam % 2 != 0:
        # If \k_lambda is odd we need to change the grid spacing so we don't have a half bump on the periodic boundary (L/8*delta_z = (k_lam/2+0.5)*L/k_lam)
        fileLines[55] = f"      if(ndims==3) xg(3)%h = {4*(1+1/k_lam)}\n"

    with open("lotus.f90","w") as fileOutput:
        fileOutput.writelines(fileLines)


