import numpy as np
import math as math
import matplotlib.pyplot as plt

class Cross_section:
    """
    Atributes:
        name: [str]
        coordinates: [lst]
    """
    def __init__(self, name, coordinates, n=200):
        self.name = name
        self.coord = np.asarray(coordinates)
        self.n = n


    def shape(self):
        if np.array_equal(self.coord[0], self.coord[-1]):
            return self.coord
        else:
            return np.append(self.coord, [self.coord[0]], axis=0)
        return

    @property
    def height(self):
        # x, y = self.shape.T
        # h = y.max() - y.min()
        return self.shape()[:,1].max()-self.shape()[:,1].min()

    @property
    def dz(self):
        return self.height/self.n

    def z_arr(self):
        arr = np.arange(self.shape()[:, 1].min() + self.dz / 2, self.height, self.dz)
        return arr

    def coor(self):
        """
        calculates x- and z-coordinates of each location
        """
        zd_ = np.zeros(1)
        xd_ = np.zeros(1)
        for i in range(self.shape().shape[0]-1):
            if self.shape()[i, 1] < self.shape()[i+1, 1]:
                zd_sel = self.z_arr()[np.nonzero(np.where(((self.z_arr() >= self.shape()[i, 1]) & (self.z_arr() <= self.shape()[i + 1, 1])), self.z_arr(), 0))]
                xd_sel = np.interp(zd_sel, (self.shape()[i, 1], self.shape()[i + 1, 1]), (self.shape()[i, 0], self.shape()[i + 1, 0]))
                zd_ = np.append(zd_, zd_sel)
                xd_ = np.append(xd_, xd_sel)
            elif self.shape()[i, 1] > self.shape()[i+1, 1]:
                zd_sel = self.z_arr()[np.nonzero(np.where(((self.z_arr() <= self.shape()[i, 1]) & (self.z_arr() >= self.shape()[i+1, 1])), self.z_arr(), 0))]
                xd_sel = np.interp(zd_sel, (self.shape()[i+1, 1], self.shape()[i, 1]), (self.shape()[i+1, 0], self.shape()[i, 0]))
                zd_ = np.append(zd_, zd_sel)
                xd_ = np.append(xd_, xd_sel)
            else:
                pass
        zd_ = zd_[1:]
        xd_ = xd_[1:]
        ind = np.lexsort((xd_, zd_))
        zd_ = zd_[ind]
        return np.vstack((xd_[ind], zd_[::-1])).T

    def width_at_z(self):
        """" werkt alleen voor figuren die niet meer dan 2 x-waarden hebben voor één z-waarde"""
        w = np.copy(self.coor())
        w[1::2, 0] -= w[0::2, 0]
        return w[1::2, :]

    def dA(self):
        dA = np.copy(self.width_at_z())
        dA[:, 0] *= self.dz
        return dA

    @property
    def area(self):
        return sum(self.dA()[:, 0])

    @property
    def center_of_gravity(self):
        return np.sum(self.dA()[:, 0]*self.dA()[:,1])/np.sum(self.dA(), axis=0)[0]

    def static_moment(self):
        smd = np.copy(self.dA())
        smd[:, 1] -= self.center_of_gravity
        smd[:, 0] *= smd[:, 1]
        sm = np.copy(smd)
        for x in range(np.size(sm, axis=0)):
            sm[x, 0] = np.sum(smd[:x, 0])
        return sm

    @property
    def moment_of_inertia(self):
        """geldt alleen voor symetrische doorsneden"""
        day = np.copy(self.dA())
        day[:,1] -= self.center_of_gravity
        day[:,1] *= day[:,1]
        inertia = day[:,0]*day[:,1]
        return np.sum(inertia)

    def sig_bend_unit(self):
        sig_unit = np.vstack(((self.z_arr()-self.center_of_gravity)/self.moment_of_inertia, self.z_arr()-self.center_of_gravity)).T
        return sig_unit

    def sig_norm_unit(self):
        sig_unit = np.vstack((np.full_like(self.z_arr() , 1/np.sum(self.dA()[:, 0])), self.z_arr()-self.center_of_gravity)).T
        return sig_unit

    def tau_unit(self):
        tau_unit = np.copy(self.static_moment())
        tau_unit[:,0] = tau_unit[:,0]/(self.width_at_z()[:,0]*self.moment_of_inertia)
        return tau_unit


def principle_stress(sig_x, tau, sig_z=0):
    sig1 = (sig_x+sig_z)/2+(((sig_x-sig_z)/2)**2+tau**2)**0.5
    sig2 = (sig_x+sig_z)/2-(((sig_x-sig_z)/2)**2+tau**2)**0.5
    theta = math.degrees(math.atan((2*tau/(sig_x-sig_z))/2))
    return sig1, sig2, theta


def principle_stress_section(sig_x, tau, sig_z=None):
    """" insert only stresses not z_location"""
    if sig_z:
        sig1 = (sig_x + sig_z) / 2 + (((sig_x - sig_z) / 2) ** 2 + tau ** 2) ** 0.5
        sig2 = (sig_x + sig_z) / 2 - (((sig_x - sig_z) / 2) ** 2 + tau ** 2) ** 0.5
        theta = np.degrees(np.arctan((2 * tau / (sig_x - sig_z)) / 2))
        return sig1, sig2, theta
    sig_z = np.zeros(sig_x.shape)
    sig1 = (sig_x + sig_z) / 2 + (((sig_x - sig_z) / 2) ** 2 + tau ** 2) ** 0.5
    sig2 = (sig_x + sig_z) / 2 - (((sig_x - sig_z) / 2) ** 2 + tau ** 2) ** 0.5
    theta = np.degrees(np.arctan((2 * tau / (sig_x - sig_z)) / 2))
    return sig1, sig2, theta
