#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from lotus import run
from changef90 import new_f90_res
from pathlib import Path
import numpy as np


def run2d(L):
    new_f90_res(2, 10**5, 32, L)
    run(4, f'{cwd}/{L}')

if __name__ == "__main__":
    cwd = Path.cwd()
    Ls = 2**np.arange(10,14,1)
    # for L in Ls:
    #     run2d(L)
    run2d(512)
    

