# pylint: disable=invalid-name

"""
gaussian
"""
import numpy as np

# Gaussian distribution


def val(x, mu, sigma, a):
    """
    Gaussian distribution
    """
    y = a * np.exp(-(x - mu)**2 / (2. * sigma**2))
    return y

# Gaussian fitting with Caruana's algorithm


def gaussian_fit_caruanas(x, y):
    """
    Gaussian fitting with Caruana's algorithm
    """
    sum_x = sum(x)
    sum_x_pow2 = sum(x**2)
    sum_x_pow3 = sum(x**3)
    sum_x_pow4 = sum(x**4)
    a = [[len(x), sum_x, sum_x_pow2], [sum_x, sum_x_pow2, sum_x_pow3], [
        sum_x_pow2, sum_x_pow3, sum_x_pow4]]

    sum_log_y = sum(np.log(y))
    sum_x_log_y = sum(np.log(y)*x)
    sum_x_pow2_log_y = sum(np.log(y)*(x**2))
    b = [sum_log_y, sum_x_log_y, sum_x_pow2_log_y]

    ans = np.linalg.solve(a, b)

    mu = -ans[1] / (2. * ans[2])
    sigma = np.sqrt(-1. / (2. * ans[2]))
    a = np.exp(ans[0] - (ans[1]**2 / (4. * ans[2])))

    return [mu, sigma, a]

# Gaussian fitting with Guo's algorithm


def gaussian_fit_guos(X, Y):
    """
    Gaussian fitting with Guo's algorithm
    """
    sum_y_pow2 = sum(Y**2)
    sum_x_y_pow2 = sum(X*(Y**2))
    sum_x_pow2_y_pow2 = sum((X**2) * (Y**2))
    sum_x_pow3_y_pow2 = sum((X**3) * (Y**2))
    sum_x_pow4_y_pow2 = sum((X**4) * (Y**2))
    a = [
        [sum_y_pow2, sum_x_y_pow2, sum_x_pow2_y_pow2],
        [sum_x_y_pow2, sum_x_pow2_y_pow2, sum_x_pow3_y_pow2],
        [sum_x_pow2_y_pow2, sum_x_pow3_y_pow2, sum_x_pow4_y_pow2]
    ]

    sum_y_pow2_log_y = sum((Y**2) * np.log(Y))
    sum_x_y_pow2_log_y = sum(X * (Y**2) * np.log(Y))
    sum_x_pow2_y_pow2_log_y = sum((X**2) * (Y**2) * np.log(Y))
    b = [sum_y_pow2_log_y, sum_x_y_pow2_log_y, sum_x_pow2_y_pow2_log_y]

    ans = np.linalg.solve(a, b)

    mu = -ans[1] / (2. * ans[2])
    sigma = np.sqrt(-1. / (2. * ans[2]))
    a = np.exp(ans[0] - (ans[1]**2 / (4. * ans[2])))

    return [mu, sigma, a]

# Gaussian fitting with Guo's algorithm with iterative procedure


def gaussian_fit_guos_iterative_procedure(X, Y, iterations):
    """
    Gaussian fitting with Guo's algorithm with iterative procedure
    """
    for _ in range(iterations):
        sum_y_pow2 = sum(Y**2)
        sum_x_y_pow2 = sum(X*(Y**2))
        sum_x_pow2_y_pow2 = sum((X**2) * (Y**2))
        sum_x_pow3_y_pow2 = sum((X**3) * (Y**2))
        sum_x_pow4_y_pow2 = sum((X**4) * (Y**2))
        a = [
            [sum_y_pow2, sum_x_y_pow2, sum_x_pow2_y_pow2],
            [sum_x_y_pow2, sum_x_pow2_y_pow2, sum_x_pow3_y_pow2],
            [sum_x_pow2_y_pow2, sum_x_pow3_y_pow2, sum_x_pow4_y_pow2]
        ]

        sum_y_pow2_log_y = sum((Y**2) * np.log(Y))
        sum_x_y_pow2_log_y = sum(X * (Y**2) * np.log(Y))
        sum_x_pow2_y_pow2_log_y = sum((X**2) * (Y**2) * np.log(Y))
        b = [sum_y_pow2_log_y, sum_x_y_pow2_log_y, sum_x_pow2_y_pow2_log_y]

        ans = np.linalg.solve(a, b)
        # prev_Y = Y
        Y = np.exp(ans[0]+(ans[1]*X)+(ans[2]*(X**2)))
        # print(Y - prev_Y)

    mu = -ans[1] / (2. * ans[2])
    sigma = np.sqrt(-1. / (2. * ans[2]))
    a = np.exp(ans[0] - (ans[1]**2 / (4. * ans[2])))

    return [mu, sigma, a]
