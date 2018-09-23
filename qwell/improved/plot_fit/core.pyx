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
def integ_zh(double ze, double a, int num_series, double q, double k, double l,
             double [::1] zgrid, double [::1] sin_zgrid):
    cdef int j, n
    cdef double zh, integral_zh, kernel, qn, zpn, zmn

    integral_zh = 0.0
    for j, zh in enumerate(zgrid):
        kernel = G(a, ze - zh) + k * G(a, ze + zh)
        for n in range(1, num_series+1):
            qn = q**n
            zpn = -zh + 4 * n * l
            zmn =  zh + 4 * n * l
            kernel += qn * (k * G(a, ze - zpn) + G(a, ze + zpn) + 
                                G(a, ze - zmn) + G(a, ze + zmn) / k)
        integral_zh += sin_zgrid[j] * kernel
    return integral_zh
