import numpy as np
from lattice import frac2cart, cart2frac, get_lattice_volume, gen_reciprocal_vectors


lat_vec = np.array([[1.0, 1.0, -1.0],
                    [2.0, -1.5, 1.0],
                    [0.0, 0.1, 2.1]])
frac_coord = np.array([[0.5, 0.0, 0.0, 0.7, 6.0],
                       [1.0, 1.0, 0.0, 0.1, -0.1],
                       [0.0, 1.0, 1.0, -1.2, 0.7]])

print("lat_vec")
print(lat_vec)
print("frac_coord")
print(frac_coord)

cart_coord = frac2cart(lat_vec.T, frac_coord.T).T
print("cart_coord")
print(cart_coord)

delta = cart2frac(lat_vec.T, cart_coord.T).T - frac_coord
print("abs(delta)")
print(np.abs(delta))

volume = get_lattice_volume(lat_vec.T)
print("volume")
print(volume)

recip_vec = gen_reciprocal_vectors(lat_vec.T).T
print("recip_vec")
print(recip_vec)
