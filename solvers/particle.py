#Lamb-Oseen vortex particle

import numpy as np 
from scipy.special import expn

class LambOseenVortexParticle(object):

    def __init__(self, Gamma, t0, nu, x0, y0, uInf, vInf):
        self.Gamma = Gamma
        self.t0 = t0
        self.t = 0
        self.nu = nu
        self.x0 = x0
        self.y0 = y0
        self.x = x0
        self.y = y0
        self.uInf = uInf
        self.vInf = vInf


    def vorticity(self, xy):
        rx = xy[:,0] - self.x
        ry = xy[:,1] - self.y
        e = 1.0e-16
        rSq = rx*rx + ry*ry + e
        a = 4.0*self.nu*(self.t + self.t0)
        c = self.Gamma/(4.0*(np.pi)*self.nu*(self.t + self.t0))
        g = np.exp(-rSq/a)
        omega = c*g

        return omega

    def vorticity_initial_blobs(self, xy, sigma):
        rx = xy[:,0] - self.x
        ry = xy[:,1] - self.y
        e = 1.0e-16
        rSq = rx*rx + ry*ry + e
        a = 4.0*self.nu*(self.t + self.t0 - (sigma**2/(2.0*self.nu)))
        c = self.Gamma/(4.0*(np.pi)*self.nu*(self.t + self.t0 - (sigma**2/(2*self.nu))))
        g = np.exp(-rSq/a)
        omega = c*g

        return omega

    def velocity(self, xy):
        rx = xy[:,0] - self.x
        ry = xy[:,1] - self.y
        e = 1.0e-16
        rSq = rx*rx + ry*ry + e
        a = 4.0*self.nu*(self.t + self.t0)
        c = self.Gamma/(2.0*(np.pi)*rSq)*self.exponent(rSq, a)

        u = -c*ry + self.uInf
        v = c*rx + self.vInf

        return u, v

    def velocity_gradient(self, xy):
        rx = xy[:,0] - self.x
        ry = xy[:,1] - self.y
        e = 1.0e-16
        rSq = rx*rx + ry*ry + e
        a = 4.0*self.nu*(self.t + self.t0)       
        c = self.Gamma/(np.pi*rSq*rSq)

        u, v = self.velocity(xy)

        dudx = c*rx*ry
        dudy = c*(ry*ry - 0.5*rSq)
        dvdx = c*(0.5*rSq - rx*rx)
        dvdy = -c*rx*ry

        g = self.exponent(rSq, a)
        dgdx, dgdy = self.exponent_gradient(rx, ry, rSq, a)

        dudx = dudx*g + (u - self.uInf)*dgdx
        dudy = dudy*g + (u - self.uInf)*dgdy
        dvdx = dvdx*g + (v - self.vInf)*dgdx
        dvdy = dvdy*g + (v - self.vInf)*dgdy

        return dudx, dudy, dvdx, dvdy

    # def pressure(self, u ,v):
    #     pUnsteady = -((u - self.uInf)*self.uInf + (v - self.vInf)*self.vInf)
    #     pDyn = 0.5*(u*u + v*v)
    #     pDynInf = 0.5*(self.uInf*self.uInf + self.vInf*self.vInf)
    #     return pDynInf - pUnsteady - pDyn

    def pressure(self, xy):
        rx = xy[:, 0] - self.x
        ry = xy[:, 1] - self.y
        rSq = rx*rx + ry*ry + 1e-16
        a = 4.0*self.nu*(self.t + self.t0)
        rhoSq = rSq/a

        p = -0.125*(self.Gamma/np.pi)**2*(expn(2, 2.0*rhoSq)/rSq - 2.0*expn(2, rhoSq)/rSq + 1.0/rSq)

        return p   

    def pressure_gradient(self, xy, u, v):
        dudx, dudy, dvdx, dvdy = self.velocity_gradient(xy)

        dpdx = -(dudx*(u - self.uInf) + dudy*(v - self.vInf))
        dpdy = -(dvdx*(u - self.uInf) + dvdy*(v - self.vInf))

        return dpdx, dpdy

    def exponent(self, rSq, a):
        return 1 - np.exp(-rSq/a)

    def exponent_gradient(self, rx, ry, rSq, a):
        c = 2.0*np.exp(-rSq/a)/a
        gx = c*rx
        gy = c*ry

        return gx, gy

    def evolve(self, Dt):
        self.x = self.x + self.uInf*Dt
        self.y = self.y + self.vInf*Dt

        self.t = self.t + Dt



