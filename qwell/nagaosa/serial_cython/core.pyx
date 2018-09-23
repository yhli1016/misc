from scipy.special.cython_special import struve, y1
import cython

"""cdef declaration cannot be placed within an 'if' or 'for' statement."""

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
def calc_eb(double x,
            double mu, double epsilon_1, double epsilon_2,
            int num_grid, int num_series,
            double a_B, double R, double q, double l,
            double [::1] zgrid, double [::1] cos_zgrid):
    # Convert the unit of x from exciton Bohr radius to absolute Bohr
    cdef double a = x * a_B

    # Evaluation the kinetic term
    cdef double kinetic_energy = 1.0 / (2.0 * mu * a**2)

    # Evaluate the Coulomb attraction term
    cdef int i, j, n
    cdef double ze, zh, integral_ze, integral_zh, kernel
    cdef double qn, zhn

    integral_ze = 0.0
    for i, ze in enumerate(zgrid):
        integral_zh = 0.0
        for j, zh in enumerate(zgrid):
            kernel = 0.0
            for n in range(-num_series, num_series+1):
                qn = q**abs(n)
                zhn = zh * (-1)**n + 2 * n * l
                kernel += qn * G(a, ze -zhn)
            integral_zh += cos_zgrid[j] * kernel
        integral_ze += cos_zgrid[i] * integral_zh
    cdef double dS = (zgrid[1] - zgrid[0]) ** 2
    cdef double potential_energy = -2.0 / (epsilon_1 * l**2 * a) * integral_ze * dS

    # Evaluate the binding energy in effective Rydberg energy
    cdef double binding_energy = (kinetic_energy + potential_energy) / R
    return binding_energy
