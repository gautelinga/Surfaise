import dolfin as df
import numpy as np
import random


# Class representing the intial conditions.
class GenericIC(df.UserExpression):
    # Base class
    def __init__(self, u_, **kwargs):
        self.size = len(u_)
        super().__init__(**kwargs)

    def eval(self, values, x):
        for i in range(self.size):
            values[i] = 0.0

    def value_shape(self):
        return (self.size,)


class RandomIC(GenericIC):
    def __init__(self, u_, amplitude=1.0, dims=1, **kwargs):
        self.amplitude = amplitude
        self.dims = dims
        super().__init__(u_, **kwargs)

    def eval(self, values, x):
        super().eval(values, x)
        for d in range(self.dims):
            values[d] = self.amplitude*(2*random.random()-1)

class NematicPerturbationIC(GenericIC):
    def __init__(self, u_, x0=0, y0=0, r=1, theta0=0, amplitude=1.0, dims=1, **kwargs):
        self.amplitude = amplitude
        self.dims = dims
        self.x0 = x0
        self.y0 = y0
        self.r = r
        self.t0 = theta0
        super().__init__(u_, **kwargs)

    def eval(self, values, x):
        super().eval(values, x)
        dtheta = self.amplitude*(2*random.random()-1) * np.heaviside(self.r**2-(x[0]-self.x0)**2-(x[1]-self.y0)**2,0)
        for d in range(self.dims):
            values[d] = 0
            if d == 0:
                values[d] = 2 * np.cos(self.t0+dtheta)*np.sin(self.t0+dtheta)
            if d == 2:
                values[d] = 2 * (np.cos(self.t0+dtheta)**2 - 0.5)

class StripedIC(GenericIC):
    def __init__(self, u_, alpha=0.0, q0=1./np.sqrt(2),
                 amplitude=0.5, **kwargs):
        self.alpha = alpha
        self.q0 = q0
        self.amplitude = amplitude
        super().__init__(u_, **kwargs)

    def eval(self, values, x):
        super().eval(values, x)
        x_n = x[0]*np.cos(self.alpha) + x[1]*np.sin(self.alpha)
        values[0] = self.amplitude*np.sin(self.q0*x_n)


# Stripes around
class AroundStripedIC(StripedIC):
    def __init__(self, u_, q0=1./np.sqrt(2),
                 amplitude=0.5, **kwargs):
        super().__init__(alpha=np.pi/2, q0=q0, amplitude=amplitude, **kwargs)


# Stripes along
class AlongStripedIC(StripedIC):
    def __init__(self, u_, q0=1./np.sqrt(2),
                 amplitude=0.5, **kwargs):
        super().__init__(alpha=0.0, q0=q0, amplitude=amplitude, **kwargs)


# Initial conditions for manufactured solution
class MMSIC(GenericIC):
    def __init__(self, u_, geo_map, **kwargs):
        self.map = geo_map
        super().__init__(u_, **kwargs)

    def eval(self, values, x):
        super().eval(values, x)
        # values[0] = 0.5*np.sin(x[0]/np.sqrt(2))
        # values[0] = self.map.psiMMS
        values[0] = (np.sin(x[1]/np.sqrt(2)))**2 + (np.sin(x[0]/np.sqrt(2)))**2


# Circular stripes
class CircularStripedIC(StripedIC):
    def eval(self, values, x):
        # angle alpha is not used for now
        super().eval(values, x)
        x_n = np.sqrt(x[0]**2+x[1]**2)
        values[0] = self.amplitude*np.sin(self.q0*x_n)
        # values[0] = np.sin(1.15*np.sqrt(x[0]**2+x[1]**2)/np.sqrt(2))
