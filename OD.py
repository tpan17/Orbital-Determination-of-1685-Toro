#Tiffany Pan
#OD
#7/22/2015

from __future__ import division
from math import *
import visual
from visual import *

#


#angle correcting for quadrant
def rad_ang(sin_ang, cos_ang):
    inv_sine = arcsin(sin_ang)
    inv_cos = arccos(cos_ang)
    if sin_ang > 0 and sin_ang < 1:
        if cos_ang > 0 and cos_ang < 1:
            return inv_sine #first quadrant
    if sin_ang > 0 and sin_ang < 1:
        if cos_ang < 0 and cos_ang > -1:
            return inv_cos #second quadrant
    if sin_ang < 0 and sin_ang > -1:
         if cos_ang < 0 and cos_ang > -1:
            value = inv_cos - 2 * inv_sine
            return value #third quadrant
    if sin_ang < 0 and sin_ang > -1:
        if cos_ang > 0 and cos_ang < 1:
            ang = 2 * pi + inv_sine
            return ang #fourth quadrant
    if sin_ang == 0:
        if cos_ang == 1 or cos_ang == -1:
            return inv_cos #x axis
    if cos_ang == 0:
        if sin_ang == 1 or sin_ang == -1:
            return inv_sine #y axis
#^finding angle value with sin and cos of angle (for future use)


#CALCULATE a
def calc_a(r, r_dot):
    value1 = 2 / mag(r)
    value2 = dot(r_dot, r_dot)
    a = 1/ (value1 - value2)
    return a
a = calc_a(r, r_dot)
print "a = ",
print a



#CALCULATE e
def calc_e(r, r_dot, a):
    cross_product = cross(r, r_dot)
    value1 = (mag(cross_product)) ** 2
    value2 = 1 - (value1 / a)
    e = sqrt(value2)
    return e
e = calc_e(r, r_dot, a)
print "e = ",
print e



#CALCULATE I
def calc_I(r, r_dot):
    cross_prod = cross(r, r_dot)
    z_comp = cross_prod.z
    value1 = mag(cross_prod)
    value2 = z_comp / value1
    I = degrees(arccos(value2))
    return I
I = calc_I(r, r_dot)
print "I = ",
print I



#CALCULATE O
def calc_O(r, r_dot, I):
    cross_prod = cross(r, r_dot)
    x_comp = cross_prod.x
    y_comp = cross_prod.y
    z_comp = cross_prod.z
    value1 = mag(cross_prod)
    if z_comp > 0:
        sin_O1 = x_comp / (value1 * sin(radians(I)))
        cos_O2 = y_comp / (-value1 * sin(radians(I)))
        return degrees(rad_ang(sin_O1, cos_O2))
    if z_comp < 0:
        sin_O3 = x_comp / (-value1 * sin(radians(I)))
        cos_O4 = y_comp / (value1 * sin(radians(I)))
        return degrees(rad_ang(sin_O3, cos_O4))
O = calc_O(r, r_dot, I)
print "O = ",
print O



#CALCULATE w
def calc_f(a, e, r, r_dot):
    cross_prod = cross(r, r_dot)
    magn = mag(cross_prod)
    value1 = (a * (1 - e**2) / mag(r) - 1)
    cos_f = (1 / e) * value1
    value2 = (a * (1 - e**2)) / magn
    value3 = dot(r, r_dot) / (e * mag(r))
    sin_f = value2 * value3
    f = degrees(rad_ang(sin_f, cos_f))
    return f
f = calc_f(a, e, r, r_dot)
#^calculate f to use in calculating w


def calc_w(f, r, I, O):
    X = r.x
    Z = r.z
    sin_fPlusw = Z / (mag(r) * sin(radians(I)))
    value1 = X / mag(r)
    value2 = cos(radians(I)) * sin_fPlusw * sin(radians(O))
    value3 = value1 + value2
    cos_fPlusw = value3 / cos(radians(O))
    fPlusw = degrees(rad_ang(sin_fPlusw, cos_fPlusw))
    w = fPlusw - f
    while w < 0:
        w += 360
    while w > 360:
        w -= 360
    return w
w = calc_w(f, r, I, O)
print "w = ",
print w



#CALCULATE M0
def calc_M(e, r, a):
    value1 = 1 / e
    value2 = 1 - (mag(r) / a)
    cos_E = value1 * value2
    E = degrees(arccos(cos_E))
    if f >= 180 and f <= 360:
        E = 360 - E
    M = radians(E) - e * sin(radians(E))
    return M
M = calc_M(e, r, a)
print M


def calc_Mfinal(a, M):
    k = .01720209895
    n = k / (a ** (3/2))
    M0 = degrees(M + n * (2451545.0 - 2456842.5))
    while M0 < 0:
        M0 += 360
    while M0 > 360:
        M0 -= 360
    return M0
#^accounting for change in time
M0 = calc_Mfinal(a, M)
print "M0 = ",
print M0


