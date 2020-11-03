import numpy as np
np.set_printoptions(precision=20)

mu = 5.
sigma = 3.
a = 1.
x_vec = np.arange(1., 10., 1)
y_vec = a * np.exp(-(x_vec - mu)**2 / (2. * sigma**2))
print(x_vec)
print(y_vec)
