import numpy as np
import scipy.linalg.lapack as spla
from scipy.special import jn


def _init(ham, nstep):
    # NOTE: it is found complex initial wave function is better for purely real
    # Hamiltonians, while real wave function works better with complex
    # Hamiltonian.
    if np.sum(np.abs(ham.imag)) > 1e-5:
        psi_0 = np.random.rand(ham.shape[0]) - 0.5
    else:
        psi_0 = np.exp(1j * 2 * np.pi * np.random.rand(ham.shape[0]))
    psi_0 /= np.linalg.norm(psi_0)
    rho = np.zeros(nstep, dtype=complex)
    rho[0] = np.vdot(psi_0, psi_0)
    return psi_0, rho


def diag(ham, dt, nstep):
    psi_0, rho = _init(ham, nstep)
    psi_t = psi_0.copy()

    e_mat, u_mat, info = spla.zheev(ham)
    e_mat = np.diag(np.exp(-1j * e_mat * dt))

    for it in range(1, nstep):
        psi_t = np.matmul(u_mat.T.conj(), psi_t)
        psi_t = np.matmul(e_mat, psi_t)
        psi_t = np.matmul(u_mat, psi_t)
        psi_t /= np.linalg.norm(psi_t)
        rho[it] = np.vdot(psi_0, psi_t)
    return rho


def cheb(ham, dt, nstep):
    psi_0, rho = _init(ham, nstep)
    psi_t = psi_0.copy()

    num_series = 10
    ck = np.array([2 * (-1j) ** k * jn(k, dt) for k in range(num_series)])
    ck[0] *= 0.5

    for it in range(1, nstep):
        phi_m2 = psi_t
        phi_m1 = np.matmul(ham, psi_t)
        psi = ck.item(0) * phi_m2 + ck.item(1) * phi_m1
        for k in range(2, num_series):
            phi = 2 * np.matmul(ham, phi_m1) - phi_m2
            psi += ck.item(k) * phi
            psi_t = psi
            phi_m2 = phi_m1
            phi_m1 = phi
        psi_t /= np.linalg.norm(psi_t)
        rho[it] = np.vdot(psi_0, psi_t)
    return rho


def cheb_sparse(ham, dt, nstep):
    psi_0, rho = _init(ham, nstep)
    psi_t = psi_0.copy()

    num_series = 10
    ck = np.array([2 * (-1j) ** k * jn(k, dt) for k in range(num_series)])
    ck[0] *= 0.5

    for it in range(1, nstep):
        phi_m2 = psi_t
        phi_m1 = ham * psi_t
        psi = ck.item(0) * phi_m2 + ck.item(1) * phi_m1
        for k in range(2, num_series):
            phi = 2 * ham * phi_m1 - phi_m2
            psi += ck.item(k) * phi
            psi_t = psi
            phi_m2 = phi_m1
            phi_m1 = phi
        psi_t /= np.linalg.norm(psi_t)
        rho[it] = np.vdot(psi_0, psi_t)
    return rho
