
from gen.legendre import gen_degree

import sympy as sp
import numpy as np


N = 12

polys, dpolys, t = gen_degree(N)



with open("include/basis/quadrature_point.h", "wt") as file:
    with open("include/basis/quadrature_weight.h", "wt") as file_w:
        for d in range(1, N):

            poly = polys[d]
            dpoly = dpolys[d]

            q = sp.roots(poly, t)
            q = np.array([float(qi) for qi in q])

            w = q * 0.0

            for i in range(len(w)):
                xi = q[i]
                w[i] = 2.0 / ((1.0 - xi*xi) * (float(dpoly.subs(t, xi))**2))

            print(d, q)

            file.write(f"if (d == {d-1}) {{\n")
            file_w.write(f"if (d == {d-1}) {{\n")

            for i in range(d):
                file.write(f"    if (q == {i}) {{\n")
                file.write(f"        return {q[i]};\n")
                file.write("    }\n")

                file_w.write(f"    if (q == {i}) {{\n")
                file_w.write(f"        return {w[i]};\n")
                file_w.write("    }\n")

            file.write("}\n")
            file_w.write("}\n")


