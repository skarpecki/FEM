from math import sqrt
from math import pow
from pprint import pprint as pp

def func_test(x,y):
    return pow(x,2)*y

def func_1(x,y):
    return -2*(pow(x,2))*y + 2*x*y + 4

def func_2(x,y):
    return -5 *pow(x,2)*y + 2*x*pow(y,2) + 10

def _fill_psx(val, no_of_y_points, no_of_x_points, x_factor, y_factor):
    psx = []
    for y in range(no_of_y_points):
        for x in range(no_of_x_points):
            psx.append((-val + x_factor * x * val, -val + y_factor * y * val))
            #psx.append((-val + 2 * x * val, -val + 2 * y * val)) - for 4 points
            #psx.append((-val + x * val, -val + y * val)) - for 9 points
    return psx



def gaussian_quadrature_2d_4pt(func):
    psx = _fill_psx(1/sqrt(3), 2, 2, 2, 2)
    result = 0
    for x, y in psx:
        result += func(x,y)
    return result


def gaussian_quadrature_2d_9pt(func):
    psx = _fill_psx(sqrt(3/5), 3, 3, 1, 1)
    result = 0
    for x, y in psx:
        wx = 5/9
        wy = 5/9
        if x == 0:
            wx = 8/9
        if y == 0:
            wy = 8/9
        result += func(x, y)*wx*wy
    return result

if __name__ == '__main__':
    print(gaussian_quadrature_2d_9pt(func_2))
    print(gaussian_quadrature_2d_4pt(func_1))
