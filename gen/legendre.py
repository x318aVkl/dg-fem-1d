import sympy as sp
import numpy as np
from math import factorial

MAX_DEGREE = 12


def gen_degree(degree):
    t = sp.symbols("t")

    N = degree + 1

    poly = []
    dpoly = []

    for i in range(N):
        p = 1.0
        p = 1.0 / (2.0**i * factorial(i)) * sp.diff((t**2 - 1)**i, t, i)

        poly.append(p)
        dpoly.append(sp.diff(p, t))
    
    return poly, dpoly, t


def gen_roots(degree):
    poly, _, t = gen_degree(degree)

    q = sp.roots(poly[-1], t)
    return np.array([float(qi) for qi in q])



def format_func(file, d, p):
    for i in range(d+1):
        file.write(f"if (i == {i}) {{\n")

        s = sp.ccode(sp.expand(p[i]))
        file.write(f"    return {s};\n")
        file.write(f"}}\n")



with open("include/basis/legendre_shape.h", "wt") as file:

    d = MAX_DEGREE

    p, dp, t = gen_degree(d)
    format_func(file, d, p)


with open("include/basis/legendre_grad.h", "wt") as file:

    d = MAX_DEGREE

    p, dp, t = gen_degree(d)
    format_func(file, d, dp)

