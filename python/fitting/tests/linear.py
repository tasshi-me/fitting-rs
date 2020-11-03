import numpy as np
np.set_printoptions(precision=20)

# positive slope
a = 1.2
b = 3.
x_vec = np.arange(1., 10., 1)
y_vec = a * x_vec + b
print(x_vec)
print(y_vec)

# negative slope
a = -1.2
b = 3.
x_vec = np.arange(1., 10., 1)
y_vec = a * x_vec + b
print(x_vec)
print(y_vec)

# zero slope
a = 0.
b = 3.
x_vec = np.arange(1., 10., 1)
y_vec = a * x_vec + b
print(x_vec)
print(y_vec)
