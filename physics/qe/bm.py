#! /usr/bin/env python
"""
Get the VBM, CBM, direct and indirect band gap, similar to bm.m.

Usage: bm.py bandstr.pd num_val_bands

band data which must be in the following form:
kpath(1) band1(1) band2(1) band3(1) ... bandm(1)
kpath(2) band1(2) band2(2) band3(2) ... bandm(2)
......            ......
kpath(n) band1(n) band2(n) band3(n) ... bandm(n)
"""

import sys
import numpy as np


data_file = sys.argv[1]
num_val_bands = int(sys.argv[2])

# Load and split data
raw_data = np.loadtxt(data_file)
band_val = raw_data[:, 1:1+num_val_bands]
band_cond = raw_data[:, 1+num_val_bands:]

# Determine the indirect band gap
ik_vbm = np.argmax([max(row) for row in band_val])
ib_vbm = np.argmax(band_val[ik_vbm, :])
ik_cbm = np.argmin([min(row) for row in band_cond])
ib_cbm = np.argmin(band_cond[ik_cbm, :])
eg = band_cond[ik_cbm, ib_cbm] - band_val[ik_vbm, ib_vbm]
#print('VBM = %10.5f' % band_val[ik_vbm, ib_vbm])
#print('CBM = %10.5f' % band_cond[ik_cbm, ib_cbm])
print("IEg = %10.5f, (%4d, %4d) -> (%4d, %4d)" % 
      (eg, ib_vbm+1, ik_vbm+1, ib_cbm+1, ik_cbm+1))

# Determine the direct band gap
ik_vbm = np.argmin([min(band_cond[i, :]) - max(band_val[i, :])
                    for i in range(len(band_val))])
ib_vbm = np.argmax(band_val[ik_vbm, :])
ik_cbm = ik_vbm
ib_cbm = np.argmin(band_cond[ik_cbm, :])
eg = band_cond[ik_cbm, ib_cbm] - band_val[ik_vbm, ib_vbm]
print("DEg = %10.5f, (%4d, %4d) -> (%4d, %4d)" % 
      (eg, ib_vbm+1, ik_vbm+1, ib_cbm+1, ik_cbm+1))
