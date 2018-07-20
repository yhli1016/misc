import numpy as np


def read_hr(file_name):
    global hr_terms
    with open(file_name, "r") as hr_file:
        hr_content = hr_file.readlines()
    for line in hr_content:
        line_split = line.split()
        Rn = np.array([float(line_split[0]), float(line_split[1]),
                       float(line_split[2])])
        ij = (int(line_split[3])-1, int(line_split[4])-1)
        tij = float(line_split[5]) + 1j * float(line_split[6])
        hr_terms.append([Rn, ij, tij])


def set_ham(k_point):
    global ham
    ham *= 0.0
    for hr_term in hr_terms:
        Rn = hr_term[0]
        ij = hr_term[1]
        tij = hr_term[2]
        hij = tij * np.exp(1j * 2 * np.pi * np.dot(k_point, Rn))
        ham[ij] += hij


def eval_energies(k_points):
    global ham, hr_terms, num_orbital
    num_k_point = len(k_points)
    energies = np.zeros((num_k_point, num_orbital), dtype="float")
    for ik, k_point in enumerate(k_points):
        print("ik = %d" % ik)
        set_ham(k_point)
        energies[ik] = np.sort(np.linalg.eigvals(ham)).real
    return energies


def main():
    # Parameters
    kx_min = 0.35
    kx_max = 0.65
    ky_min = 0.0
    ky_max = 0.3
    nx = 20
    ny = 20

    # Create k-mesh
    kx_list = np.linspace(kx_min, kx_max, nx)
    ky_list = np.linspace(ky_min, ky_max, ny)
    k_mesh = [np.array((kx, ky, 0.0)) for kx in kx_list for ky in ky_list]

    # Evaluate energies
    full_energies = eval_energies(k_mesh)

    # Output
    with open("energies.dat", "w") as output:
        for ik, k_point in enumerate(k_mesh):
            kx = k_point[0]
            ky = k_point[1]
            energy_1 = full_energies[ik, 17]
            energy_2 = full_energies[ik, 18]
            output.write("%12.5f%12.5f%12.5f%12.5f\n" %
                         (kx, ky, energy_1, energy_2))


# Load hr.dat
hr_terms = []
read_hr("hop.dat")

# Create the Hamiltonian
# We define it as global variable to avoid recreating it every time when calling
# set_ham. Hopefully it may save some time.
num_orbital = 30
ham = np.zeros((num_orbital, num_orbital), dtype="complex")

# Run the code
main()
