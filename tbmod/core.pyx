import cython
from libc.math cimport cos, sin


@cython.boundscheck(False)
@cython.wraparound(False)
def set_ham(complex [:, ::1] ham, int [:, ::1] ij, double [:, ::1] Rn,
            complex [:] tij, double [:] kpoint):
    cdef double PI2 = 6.28318530717959
    cdef double phase
    cdef int i, j, k
    for i in range(ham.shape[0]):
        for j in range(ham.shape[1]):
            ham[i][j] = 0.0
    for k in range(len(Rn)):
        i = ij[k, 0]
        j = ij[k, 1]
        phase =  PI2 * (kpoint[0] * Rn[k, 0] + kpoint[1] * Rn[k, 1] +
                        kpoint[2] * Rn[k, 2])
        ham[i, j] = ham[i, j] + tij[k] * (cos(phase) + 1j * sin(phase))
