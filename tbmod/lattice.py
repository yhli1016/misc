"""Functions for lattice operations."""

from math import sin, cos, sqrt, pi
import numpy as np


def gen_lattice_vectors(lattice_parameters):
    """
    Generate lattice vectors from given lattice parameters.

    :param lattice_parameters: list containing the six lattice parameters
            [a, b, c, alpha, beta, gamma]
    :return:
        lattice_vectors: a 3*3 array containing the lattice vectors
    """
    assert len(lattice_parameters) == 6
    a = lattice_parameters[0]
    b = lattice_parameters[1]
    c = lattice_parameters[2]
    alpha = lattice_parameters[3] * 1.0 / 180 * pi
    beta = lattice_parameters[4] * 1.0 / 180 * pi
    gamma = lattice_parameters[5] * 1.0 / 180 * pi
    lattice_vectors = np.zeros((3, 3))
    lattice_vectors[0, :] = [a, 0, 0]
    lattice_vectors[1, :] = [b*cos(gamma), b*sin(gamma), 0]
    lattice_vectors[2, 0] = c * cos(beta)
    lattice_vectors[2, 1] = c * (cos(alpha) - cos(beta)*cos(gamma)) / sin(gamma)
    lattice_vectors[2, 2] = c * sqrt(1 + 2*cos(alpha)*cos(beta)*cos(gamma)
                                     - cos(alpha)**2 - cos(beta)**2
                                     - cos(gamma)**2) / sin(gamma)
    return lattice_vectors


def gen_reciprocal_vectors(lattice_vectors):
    """
    Generate reciprocal lattice vectors from real-space lattice vectors.

    :param lattice_vectors: 3*3 array
    :return: reciprocal_vectors: 3*3 array
    """
    reciprocal_vectors = np.zeros((3, 3))

    # NOTE: the volume here should have a sign, i.e. DO NOT add any abs()
    # here. Otherwise the results will be wrong when volume < 0.
    a0 = lattice_vectors[0]
    a1 = lattice_vectors[1]
    a2 = lattice_vectors[2]
    volume = np.dot(np.cross(a0, a1), a2)
    reciprocal_vectors[0] = np.cross(a1, a2)
    reciprocal_vectors[1] = np.cross(a2, a0)
    reciprocal_vectors[2] = np.cross(a0, a1)
    reciprocal_vectors *= (2 * pi / volume)

    # Alternatively, you can use the algorithm below which strictly follows the
    # definition of reciprocal lattice vectors.
    #product = 2 * pi * np.eye(3)
    #for i in range(3):
    #    reciprocal_vectors[i] = np.linalg.solve(lattice_vectors, product[i]);

    return reciprocal_vectors


def cart2frac(lattice_vectors, cartesian_coordinates):
    """
    Convert Cartesian coordinates to fractional coordinates.

    :param lattice_vectors: 3*3 array containing the coordinates of lattice
            vectors
    :param cartesian_coordinates: N*3 array containing the Cartesian coordinates
            of atoms
    :return: fractional_coordinates: N*3 array
    """
    fractional_coordinates = np.zeros(cartesian_coordinates.shape)
    conversion_matrix = np.linalg.inv(lattice_vectors.T)
    for i, row in enumerate(cartesian_coordinates):
        fractional_coordinates[i] = np.matmul(conversion_matrix, row.T)
    return fractional_coordinates


def frac2cart(lattice_vectors, fractional_coordinates):
    """
    Convert fractional coordinates to Cartesian coordinates.

    :param lattice_vectors: 3*3 array
    :param fractional_coordinates: N*3 array
    :return: cartesian coordinates: N*3 array
    """
    cartesian_coordinates = np.zeros(fractional_coordinates.shape)
    conversion_matrix = lattice_vectors.T
    for i, row in enumerate(fractional_coordinates):
        cartesian_coordinates[i] = np.matmul(conversion_matrix, row.T)
    return cartesian_coordinates
