import numpy as np
import core
import propagate


class TBModel(object):
    """A simple yet efficient tight-binding model."""

    def __init__(self):
        self.num_orb = None
        self.hop_ind = None
        self.hop_eng = None
        self.ham = None

    def read_hr(self, hr_name, threshold=1.0e-4):
        """
        Read in seedname_hr.dat produced by Wannier90 and set up on-site
        energies and hopping terms.

        :param hr_name: filename of seedname_hr.dat
        :param threshold: float, threshold for the matrix elements
        :return: None.
        """
        # Parse seedname_hr.dat
        with open(hr_name, "r") as hr_file:
            hr_content = hr_file.readlines()
        hop_ind, hop_eng = [], []
        for line in hr_content:
            sline = line.split()
            Rij = (int(sline[0]), int(sline[1]), int(sline[2]),
                   int(sline[3])-1, int(sline[4])-1)
            tij = float(sline[5]) + 1j * float(sline[6])
            if abs(tij) >= threshold:
                hop_ind.append(Rij)
                hop_eng.append(tij)

        # Set up arrarys
        self.hop_ind = np.array(hop_ind, dtype=np.int32)
        self.hop_eng = np.array(hop_eng, dtype=np.complex128)
        self.num_orb = len(set(self.hop_ind[:, 3]))
        self.ham = np.zeros((self.num_orb, self.num_orb), dtype=np.complex128)


    def eval_energies(self, kpoints, orbital_list=None):
        """
        Get the energies for each kpoint.

        :param kpoints: (num_kpt, 3) float64 array
        :param orbital_list: list of integers
            indices of orbitals on which projections are evaluated
            COUNTED RROM 1
        :return: energies, (num_kpt, num_orb) float64 array
        :return: projection, (num_kpt, num_orn float64 array
        """
        # Allocate arrays
        num_kpt = len(kpoints)
        energies = np.zeros((num_kpt, self.num_orb), dtype=np.float64)
        projection = np.zeros((num_kpt, self.num_orb), dtype=np.float64)

        # Set the mask upon which to evaluate projection
        if orbital_list is not None:
            mask = np.array([True if i+1 in orbital_list else False
                             for i in range(self.num_orb)])
        else:
            mask = None

        # Loop over k-points to evaluate the energies and projections
        for ik, kpoint in enumerate(kpoints):
            self.ham *= 0.0
            core.set_ham(self.ham, self.hop_ind, self.hop_eng, kpoint)
            if orbital_list is not None:
                eigval, eigvec = np.linalg.eig(self.ham)
                energies[ik] = eigval.real
                for ib in range(self.num_orb):
                    projection[ik, ib] += np.sum(np.abs(eigvec[mask, ib])**2)
                    projection[ik, ib] /= np.sum(np.abs(eigvec[:, ib])**2)
            else:
                energies[ik] = np.sort(np.linalg.eigvals(self.ham)).real
        return energies, projection

    def propagate(self, kpoint, dt, nstep):
        """
        Evaluate eigenvalues at given kpoint using time-propagation method.

        This method is for prove-of-principle use.

        :param kpoint: 1*3 float64 array
            (ka, kb, kc), fractional coordinates of kpoint
        :param dt: float
            time step for propagation
            in Hartree atomic units (1 fs ~= 40 a.u.)
        :param step: integer
            total number of steps for propagation
        :return eng: 1*nstep float64 array
            x-axis for dos
        :return dos: 1*nstep float64 array
            density of states
        :return peaks: float64 array
            peaks in dos
        """
        EV2HAR = 1 / (2 * 13.60569253)

        # Set up the Hamiltonian.
        self.ham *= 0.0
        core.set_ham(self.ham, self.hop_ind, self.hop_eng, kpoint)
        self.ham *= EV2HAR

        # Propagate
        rho = propagate.cheb(self.ham, dt, nstep)

        # Add rho of -t
        rho_mt = np.conj(rho[1:][::-1])
        rho = np.hstack((rho_mt, rho))

        # Assemble dos
        dos = np.abs(np.fft.ifft(rho))
        nstep = len(dos)
        dos_1, dos_2 = dos[:nstep//2], dos[nstep//2:]
        dos = np.hstack((dos_2, dos_1))
        eng = 2 * np.pi / (dt * nstep) * \
            np.linspace(-nstep//2, nstep//2-1, nstep) / EV2HAR

        # Search for peaks in dos. 
        peaks = []
        for ie in range(1, dos.shape[0]-1):
            if dos.item(ie-1) < dos.item(ie) > dos.item(ie+1):
               peaks.append(eng.item(ie))
        peaks = np.array(peaks)
        return eng, dos, peaks
