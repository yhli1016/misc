import numpy as np
from math import sqrt


def f(x, y):
    assert x.shape[0] == y.shape[0]
    n = x.shape[0]
    sum_x = np.sum(x)
    sum_y = np.sum(y)
    sum_xy = np.sum(x * y)
    return n * sum_xy - sum_x * sum_y


# Prepare data
data = np.array([
    [5.916, 2.432],
    [7.823, 2.398],
    [10.316, 2.384],
    [14.274, 2.366],
    [19.846, 2.357],
    [26.151, 2.348]])

x = 1 / data[:, 0]**2
y = data[:, 1]

# Evaluate the corrcoeff
corrcoef = f(x, y) / sqrt(f(x, x) * f(y, y))
print(corrcoef)
print(np.corrcoef(x, y)[0, 1])
