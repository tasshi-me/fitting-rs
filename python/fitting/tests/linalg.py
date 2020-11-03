import numpy as np
np.set_printoptions(precision=20)

# has one solution 1
a = [[2., 1., -3., -2.],
     [2., -1., -1., 3.],
     [1., -1., -2., 2.],
     [-1., 1., 3., -2.]]

b = [-4., 1., -3., 5.]

x = np.linalg.solve(a, b)

assert all(x == [1., 2., 2., 1.]), "Invalid"
print(x)

# has one solution 2
a = [[2., 1., -3.], [2., -1., -1.], [1., -1., -2.], [-1., 1., 3.]]

b = [-2., -2., -5., 7.]

x = np.linalg.solve(a, b)

assert all(x == [1., 2., 2.]), "Invalid"
print(x)

# has inf solutions
a = [[2., 1., -3., -2.], [2., -1., -1., 3.], [1., -1., -2., 2.]]

b = [4., 1., -3.]

x = np.linalg.solve(a, b)

assert all(x == [1., 2., 2.]), "Invalid"
print(x)
