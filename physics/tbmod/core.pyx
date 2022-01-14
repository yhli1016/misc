# cython: language_level=3
import cython
import numpy as np
from libc.math cimport cos, sin


@cython.boundscheck(False)
@cython.wraparound(False)
def set_ham(double complex [:,::1] ham, int [:,::1] hop_ind,
        double complex [::1] hop_eng, double [::1] kpoint):
    cdef double TPI = 6.28318530717959
    cdef double phase
    cdef int ih, ii, jj

    for ih in range(hop_ind.shape[0]):
        phase = TPI * (kpoint[0] * hop_ind[ih, 0] + 
                       kpoint[1] * hop_ind[ih, 1] +
                       kpoint[2] * hop_ind[ih, 2])
        ii, jj = hop_ind[ih, 3], hop_ind[ih, 4]
        ham[ii, jj] = ham[ii, jj] + hop_eng[ih] * (cos(phase) + 1j * sin(phase))
