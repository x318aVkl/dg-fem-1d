import sympy as sp
import numpy as np

from legendre import gen_roots


MAX_DEGREE = 12



def gen_degree(degree):
    t = sp.symbols("t")

    N = degree + 1

    #ts = np.linspace(-1, 1, N)

    ts = np.sort(gen_roots(N))

    poly = []
    dpoly = []

    for i in range(N):
        p = 1.0
        for j in range(N):
            if i != j:
                p = p * (t - ts[j]) / (ts[i] - ts[j])

        poly.append(p)
        dpoly.append(sp.diff(p, t))
    
    return poly, dpoly


def format_func(file, d, p):
    if d == 0:
        file.write("if (d == 0) {\n")
    else:
        file.write(f"}} else if (d == {d}) {{\n")

    for i in range(d+1):
        file.write(f"    if (i == {i}) {{\n")
        s = sp.ccode(sp.expand(p[i]))
        file.write(f"        return {s};\n")
        file.write(f"    }}\n")



with open("include/basis/lagrange_shape.h", "wt") as file:

    for d in range(MAX_DEGREE):
        p, dp = gen_degree(d)
        format_func(file, d, p)

    file.write("}\n")


with open("include/basis/lagrange_grad.h", "wt") as file:

    for d in range(MAX_DEGREE):
        p, dp = gen_degree(d)
        format_func(file, d, dp)

    file.write("}\n")

