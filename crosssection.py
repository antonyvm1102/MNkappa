import numpy as np
import matplotlib.pyplot as plt

# class Cross_section():
#     def __init__(self,x_coor,y_coor):
#         x.coor = self.x_coor

def rectangle(b,h):
    return np.array([[0,0], [0,h], [b,h], [b,0]])

# def t_beam():
#     pass
# def i_beam():
#     pass

#solid

def y_arr(h):
    dy = h/200
    return np.arange(h - dy / 2, dy / 2, dy)

def dA_arr(y,h):
    dy = h / 200
    return np.full(y.size, dy * b)

def calc_h(cross_section):
    return y.max()

def eps_c(eps_top, x, y):
    return (y-x)/x * eps_top

def stress_c(eps,material):

def force_c(stress,da):
    return stress*da

#reinforcement
def eps_r(eps_top, x, d):
    return

def sig_r(eps, material):

    return sig_r

def force_r(sigma, As):
    return sigma * As

# general

def moment(force_c,y_c,force_r,d_r):
    """ reinforcemen: [[As,d],[As,d]]"""
    return np.sum(force_c * y_c) + sum(force_r*d_r)


b, h = 500, 1000

cs = rectangle(b, h)

dA =
eps_top = 5

np.zeros(y.size)



plt.fill(x, y, color=grey)
plt.show()

print(dy)