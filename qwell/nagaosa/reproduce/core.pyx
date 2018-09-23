from scipy.special.cython_special import struve, y1
import cython


cdef double G(double a, double gamma):
    cdef double HALF_PI = 1.570796326794897
    cdef double abs_gamma = abs(gamma)
    cdef double x, H1, N1
    if abs_gamma < 1.0e-9:
        return 1.0
    else:
        x = 2 * abs_gamma / a
        H1 = struve(1, x)
        N1 = y1(x)
        return x * (HALF_PI * (H1 - N1) - 1)


@cython.boundscheck(False)
@cython.wraparound(False)
def integ_zh(double ze, double a, int num_series, double q, double l,
             double [::1] zgrid, double [::1] cos_zgrid):
    cdef int j, n
    cdef double zh, integral_zh, kernel, qn, zhn

    integral_zh = 0.0
    for j, zh in enumerate(zgrid):
        kernel = 0.0
        for n in range(-num_series, num_series+1):
            qn = q**abs(n)
            zhn = (-1)**n * zh + 2 * n * l
            kernel += qn * G(a, ze -zhn)
        integral_zh += cos_zgrid[j] * kernel
    return integral_zh
