import sympy as sp
import numpy as np
import dolfin as df


class ManufacturedSolution:
    def __init__(self, geo_map, f):
        self.t = geo_map.r_ref["t"]
        self.s = geo_map.r_ref["s"]
        self.map = dict()
        self.geo_map = geo_map
        self.geodict = geo_map.map
        self.f = f
        self.initialize()

    def initialize(self):
        # Manufactured solution:
        print("Computing Manufactured Solution.")
        self.map["psi"] = self.f
        self.map["nuhat"] = (
            self.geodict["K^tt"] * sp.diff(self.f, self.t, self.t)
            + 2 * self.geodict["K^ts"] * sp.diff(self.f,
                                                 self.s, self.t)
            + self.geodict["K^ss"]*sp.diff(self.f, self.s, self.t)
            - self.geodict["K^tt"]
            * self.geodict["G^t_tt"] * sp.diff(self.f, self.t)
            - 2 * self.geodict["K^st"]
            * self.geodict["G^t_st"] * sp.diff(self.f, self.t)
            - self.geodict["K^ss"]
            * self.geodict["G^t_ss"] * sp.diff(self.f, self.t)
            - self.geodict["K^tt"]
            * self.geodict["G^s_tt"] * sp.diff(self.f, self.s)
            - 2 * self.geodict["K^st"]
            * self.geodict["G^s_st"] * sp.diff(self.f, self.s)
            - self.geodict["K^ss"]
            * self.geodict["G^s_ss"] * sp.diff(self.f, self.s)
        )
        self.map["nu"] = (
            self.geodict["g^tt"] * sp.diff(self.f, self.t, self.t)
            + 2 * self.geodict["g^st"] * sp.diff(self.f,
                                                 self.s, self.t)
            + self.geodict["g^ss"] * sp.diff(self.f, self.s, self.s)
            + (1/self.geodict["sqrt_g"]) * (
                self.geodict["g^tt"]
                * sp.diff(self.geodict["sqrt_g"], self.t)
                * sp.diff(self.f, self.t)
                + self.geodict["g^st"]
                * sp.diff(self.geodict["sqrt_g"], self.t)
                * sp.diff(self.f, self.s)
                + self.geodict["g^st"]
                * sp.diff(self.geodict["sqrt_g"], self.s)
                * sp.diff(self.f, self.t)
                # A problem occurs in the following THREE lines of code (when run for EllipsoidMap)!
                + self.geodict["g^ss"]
                * sp.diff(self.geodict["sqrt_g"], self.s)  # This line is the problem
                * sp.diff(self.f, self.s)
                )
            )

        self.evalf = dict()
        for key in self.map.keys():
            try:
                self.evalf[key] = sp.lambdify([self.t, self.s], self.map[key], "numpy")
            except:
                simp = self.map[key].doit()
                print("Lamdify failed for key: ", key)
                print("With expression:", simp)
                exit()

    def eval(self, key):
        v = self.evalf[key](self.geo_map.r_ref_vals["t"],
                            self.geo_map.r_ref_vals["s"])
        if isinstance(v, int) or isinstance(v, float):
            return v*np.ones_like(self.t_vals)
        else:
            return v

    def get_function(self, key):
        f = df.Function(self.geo_map.S_ref)
        f.rename("{}MMS".format(key), "tmp")
        F = self.eval(key)
        f.vector()[:] = F
        return f

    def psi(self):
        return self.get_function("psi")

    def laplacian(self):
        return self.get_function("nu")

    def curvplacian(self):
        return self.get_function("nuhat")
