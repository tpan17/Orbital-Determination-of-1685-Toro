#Tiffany Pan
#OD
#7/22/2015

from __future__ import division
from math import *
import visual
from visual import *
import numpy as np

#INPUT FILE
OD_input = open("odinput1.txt") #raw_input("\nEnter file name(ODinput.txt): "))
OD_input = OD_input.readlines()

o_elements = np.arange(len(OD_input) * 8)
o_elements.shape = (len(OD_input), 8)
print o_elements

#RA & DEC
RAs = []
Decs = []
times = []
solar_vecs = []
earth_vel_vecs = []
#^designates empty lists for RAs, Decs, observation times, and vectors to be appended to
for line in OD_input:
    line = line.strip() #takes away line breaks
    line = line.split(';') #removes semicolons and separates each item
    line[0] = float(line[0][0:2]) + float(line[0][3:5])/60 + float(line[0][6:])/3600
    RAs.append(line[0])
    if line[1][1] == "-":
        line[1] = float(line[1][1:4]) - float(line[1][5:7])/60 - float(line[1][8:])/3600
    else:
        line[1] = float(line[1][1:3]) + float(line[1][4:6])/60 + float(line[1][7:])/3600
    Decs.append(line[1])
    #^converts RAs to decimal hours and Decs to decimal degrees & appends to corresponding lists
    times.append(line[2])
    solar_vecs.append(line[3])
    earth_vel_vecs.append(line[4])
    #^appends different values to separate lists
ra = []
for i in range(len(RAs)):
    RA = RAs[i]
    ra.append(RA)

dec = []
for i in range(len(Decs)):
    Dec = Decs[i]
    dec.append(Dec)



#JULIAN DATE
def calc_JD(times):
    JDs = []
    for i in range(len(times)):
        Y = float(times[i][1:5])
        M = float(times[i][6:8])
        D = float(times[i][9:11])
        UT = float(times[i][12:14]) + float(times[i][15:17])/60 + float(times[i][18:])/3600
        term1 = 367 * Y
        value1 = (7 * (Y + int((M + 9) / 12)) / 4 )
        term2 = int(value1)
        term3 = int(275 * M / 9)
        J0 = term1 - term2 + term3 + D + 1721013.5
        JD = J0 + UT/24
        JDs.append(JD)
    return JDs
JDs = calc_JD(times)
#^converts observation times to Julian Date and returns in a list
time = []
for i in range(len(JDs)):
    t_ = JDs[i]
    time.append(t_)
    tneg1 = time[0]
    t1 = time[len(time)-1]



#SOLAR VECTOR
def solar_vector(solar_vecs):
    v = []
    for vector in solar_vecs:
        vector = vector.split(',')
        vectors = []
        for i in range(len(vector)):
            i = float(vector[i])
            vectors.append(i)
        v.append(vectors)
    return v
solar_vector = solar_vector(solar_vecs)
#^function separates solar vectors and puts them into a list
R = []
for i in range(len(solar_vector)):
    s_vector = vector(solar_vector[i])
    R.append(s_vector)
    Rneg1 = R[0]
    R1 = R[len(R)-1]



#EARTH VELOCITY VECTOR
def earth_vel_vector(earth_vel_vecs):
    v = []
    for vector in earth_vel_vecs:
        vector = vector.split(',')
        vectors = []
        for i in range(len(vector)):
            i = float(vector[i])
            vectors.append(i)
        v.append(vectors)
    return v
earth_vel_v = earth_vel_vector(earth_vel_vecs)
#^function separates earth velocity vectors and puts them into a list
E_vel = []
for i in range(len(earth_vel_v)):
    e_vec = vector(earth_vel_v[i])
    E_vel.append(e_vec)
    E_velneg1 = E_vel[0]
    E_vel1 = E_vel[len(E_vel)-1]



#P-HAT
def p_hat(ra, dec):
    p_h = []
    for i in range(len(ra)):
        rad_alpha = radians(ra[i]*15)
        rad_dec = radians(dec[i])
        i_hat = cos(rad_alpha) * cos(rad_dec)
        j_hat = sin(rad_alpha) * cos(rad_dec)
        k_hat = sin(rad_dec)
        p_hat = vector(i_hat, j_hat, k_hat)
        p_h.append(p_hat)
    return p_h
p_hat_ = p_hat(ra, dec)
p_hat_b = []
for i in range(len(p_hat_)):
    value = p_hat_[i]
    p_hat_b.append(value)
    p_hat_bneg1 = p_hat_b[0]
    p_hat_b1 = p_hat_b[len(p_hat_b)-1]


#STELLAR ABERRATION
c = 173.144633
def p_hat_corrected(p_hat_b, E_vel):
    p = []
    for i in range(len(p_hat_b)):
        numer = c * p_hat_b[i] - E_vel[i]
        denom = mag(c * p_hat_b[i] - E_vel[i])
        p_hat_new = numer / denom
        p.append(p_hat_new)
    return p
p_ = p_hat_corrected(p_hat_b, E_vel)
p_hat = []
for i in range(len(p_)):
    value = p_[i]
    p_hat.append(value)
    p_hatneg1 = p_hat[0]
    p_hat1 = p_hat[len(p_hat)-1]
print p_hat
print p_hat1



#k
k = .01720209895




#PREITERATION

#magnitudes of p's
def calc_m_pneg1(R0, p_hat0, a1, a3):
    pneg1_num = a1 * dot(cross(Rneg1, p_hat0), p_hat1) - dot(cross(R0, p_hat0), p_hat1) + a3 * dot(cross(R1, p_hat0), p_hat1)
    pneg1_denom = a1 * dot(cross(p_hatneg1, p_hat0), p_hat1)
    pneg1 = pneg1_num / pneg1_denom
    return pneg1
def calc_m_p0(R0, p_hat0, a1, a3):
    p0_num = a1 * dot(cross(p_hatneg1, Rneg1), p_hat1) - dot(cross(p_hatneg1, R0), p_hat1) + a3 * dot(cross(p_hatneg1, R1), p_hat1)
    p0_denom = -dot(cross(p_hatneg1, p_hat0), p_hat1)
    p0 = p0_num / p0_denom
    return p0 
def calc_m_p1(R0, p_hat0, a1, a3):
    p1_num = a1 * dot(cross(p_hat0, Rneg1), p_hatneg1) - dot(cross(p_hat0, R0), p_hatneg1) + a3 * dot(cross(p_hat0, R1), p_hatneg1)
    p1_denom = a3 * dot(cross(p_hat0, p_hat1), p_hatneg1)
    p1 = p1_num / p1_denom
    return p1
#r's (vectors)
def calc_rneg1(m_pneg1, p_hatneg1, Rneg1):
    rneg1 = m_pneg1 * p_hatneg1 - Rneg1
    return rneg1
def calc_r0(m_p0, p_hat0, R0):
    r0 = m_p0 * p_hat0 - R0
    return r0
def calc_r1(m_p1, p_hat1, R1):
    r1 = m_p1 * p_hat1 - R1
    return r1
#r dot (vector)
def calc_r0dot(rneg1, r0, r1):
    value1 = (r0 - rneg1)/(-Tneg1)
    value2 = (r1 - r0)/(T1)
    r0_dot = 1/2 * (value1 + value2)
    return r0_dot
#f(T)
def funct_f(T, r0, r0_dot):
    m_r0 = mag(r0)
    term2 = (T**2) / (2 * m_r0**3) 
    term3 = (T**3 * dot(r0, r0_dot)) / (2 * m_r0**5)
    f = 1 - term2 + term3
    return f
#g{T}
def funct_g(T, r0, r0_dot):
    m_r0 = mag(r0)
    term2 = (T**3) / (6 * m_r0**3)
    g = T - term2
    return g
#r(T)
def funct_r(f, g, r0, r0_dot):
    r_T = f * r0 + g * r0_dot
    return r_T
    
#for loop that iterates over combinations of observations
for t in range(1, len(JDs)-1):
    #T
    Tneg1 = k * (tneg1 - JDs[t])
    T0 = k * (t1 - tneg1)
    T1 = k * (t1 - JDs[t])
    #a1 & a3
    a1 = T1 / T0
    a3 = - Tneg1 / T0
    #magnitude of p
    m_pneg1 = calc_m_pneg1(R[t], p_hat[t], a1, a3)
    m_p0 = calc_m_p0(R[t], p_hat[t], a1, a3)
    m_p1 = calc_m_p1(R[t], p_hat[t], a1, a3)
    #r's (vector)
    rneg1 = calc_rneg1(m_pneg1, p_hatneg1, Rneg1)
    r0 = calc_r0(m_p0, p_hat[t], R[t])
    r1 = calc_r1(m_p1, p_hat1, R1)
    #r dot (vector)
    r0_dot = calc_r0dot(rneg1, r0, r1)
    #f(T)
    fneg1 = funct_f(Tneg1, r0, r0_dot)
    f0 = funct_f(T0, r0, r0_dot)
    f1 = funct_f(T1, r0, r0_dot)
    #g(T)
    gneg1 = funct_g(Tneg1, r0, r0_dot)
    g0 = funct_g(T0, r0, r0_dot)
    g1 = funct_g(T1, r0, r0_dot)
    #r(T) vector
    rT_neg1 = funct_r(fneg1, gneg1, r0, r0_dot)
    rT_0 = funct_r(f0, g0, r0, r0_dot)
    rT_1 = funct_r(f1, g1, r0, r0_dot)


#vector p's
rneg1_new = rT_neg1
r1_new = rT_1
pneg1 = rneg1_new + Rneg1
p1 = r1_new + R1

def calc_p_hatneg1(pneg1):
    p_hatneg1 = pneg1 / mag(pneg1)
    return p_hatneg1
p_hatneg1_new = calc_p_hatneg1(pneg1)

def calc_p_hat1(p1):
    p_hat1 = p1 / mag(p1)
    return p_hat1
p_hat1_new = calc_p_hat1(p1)


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




#ORBITAL ELEMENTS
###CALCULATE a
##def calc_a(r, r_dot):
##    value1 = 2 / mag(r)
##    value2 = dot(r_dot, r_dot)
##    a = 1/ (value1 - value2)
##    return a
####a = calc_a(r, r_dot)
####print "\na = ",
####print a

#CALCULATE e
##def calc_e(r, r_dot, a):
##    cross_product = cross(r, r_dot)
##    value1 = (mag(cross_product)) ** 2
##    value2 = 1 - (value1 / a)
##    e = sqrt(value2)
##    return e
####e = calc_e(r, r_dot, a)
####print "e = ",
####print e

###CALCULATE I
##def calc_I(r, r_dot):
##    cross_prod = cross(r, r_dot)
##    z_comp = cross_prod.z
##    value1 = mag(cross_prod)
##    value2 = z_comp / value1
##    I = degrees(arccos(value2))
##    return I
####I = calc_I(r, r_dot)
####print "I = ",
####print I

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
##O = calc_O(r, r_dot, I)
##print "O = ",
##print O

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
##f = calc_f(a, e, r, r_dot)
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
##w = calc_w(f, r, I, O)
##print "w = ",
##print w

#CALCULATE M0
def calc_M(e, r, a, f):
    value1 = 1 / e
    value2 = 1 - (mag(r) / a)
    cos_E = value1 * value2
    E = degrees(arccos(cos_E))
    if f >= 180 and f <= 360:
        E = 360 - E
    M = radians(E) - e * sin(radians(E))
    return M
##M = calc_M(e, r, a)
##print "M(t) = ",
##print degrees(M)

te = 2457000.0
for t in range(1, len(JDs)-1):
    t0 = time[t]
def calc_Mfinal(a, M):
    k = .01720209895
    n = k / (a ** (3/2))
    M0 = degrees(M + n * (te - t0))
    while M0 < 0:
        M0 += 360
    while M0 > 360:
        M0 -= 360
    return M0
#^accounting for change in time
##M0 = calc_Mfinal(a, M)
##print "M(te) = ",
##print M0

def OD(r, r_dot):
    value1 = 2 / mag(r)
    value2 = dot(r_dot, r_dot)
    a = 1/ (value1 - value2)
    cross_product = cross(r, r_dot)
    value1 = (mag(cross_product)) ** 2
    value2 = 1 - (value1 / a)
    e = sqrt(value2)
    cross_prod = cross(r, r_dot)
    z_comp = cross_prod.z
    value1 = mag(cross_prod)
    value2 = z_comp / value1
    I = degrees(arccos(value2))
    return I
    return a, e, I, O, f, w, M, M0



#ITERATION
#formulas to find new r0, r0 dot, a1, & a3
def new_r0(rneg1_new, r1_new, fneg1, f1, gneg1, g1):
    numer = g1 * rneg1_new - gneg1 * r1_new
    denom = fneg1 * g1 - f1 * gneg1
    r0 = numer / denom
    return r0
def new_r0_dot(rneg1_new, r1_new, fneg1, f1, gneg1, g1):
    numer = f1 * rneg1_new - fneg1 * r1_new
    denom = f1 * gneg1 - fneg1 * g1
    r0_dot = numer / denom
    return r0_dot
def new_a1(fneg1, f1, gneg1, g1):
    a1 = g1 / (g1 * fneg1 - gneg1 * f1)
    return a1
def new_a3(fneg1, f1, gneg1, g1):
    a3 = -gneg1 / (g1 * fneg1 - gneg1 * f1)
    return a3

e = radians(23.5)
#iteration
total = 0
prev_a1 = 100000
prev_a3 = 100000
for t in range(1, len(JDs)-1):
    #magnitude of p
    m_pneg1 = calc_m_pneg1(R[t], p_hat[t], a1, a3)
    m_p0 = calc_m_p0(R[t], p_hat[t], a1, a3)
    m_p1 = calc_m_p1(R[t], p_hat[t], a1, a3)
    #r's (vector)
    rneg1 = calc_rneg1(m_pneg1, p_hatneg1, Rneg1)
    r0 = calc_r0(m_p0, p_hat[t], R[t])
    r1 = calc_r1(m_p1, p_hat1, R1)
    #r dot (vector)
    r0_dot = calc_r0dot(rneg1, r0, r1)
    #f(T)
    fneg1 = funct_f(Tneg1, r0, r0_dot)
    f0 = funct_f(T0, r0, r0_dot)
    f1 = funct_f(T1, r0, r0_dot)
    #g(T)
    gneg1 = funct_g(Tneg1, r0, r0_dot)
    g0 = funct_g(T0, r0, r0_dot)
    g1 = funct_g(T1, r0, r0_dot)
    #r(T) vector
    rT_neg1 = funct_r(fneg1, gneg1, r0, r0_dot)
    rT_0 = funct_r(f0, g0, r0, r0_dot)
    rT_1 = funct_r(f1, g1, r0, r0_dot)
    while abs(float(a1)-float(prev_a1)) >= .00000001:
        while abs(float(a3)-float(prev_a3)) >= .00000001:
            prev_a1 = a1
            prev_a3 = a3
            r0 = new_r0(rT_neg1, rT_1, fneg1, f1, gneg1, g1)
            r0_dot = new_r0_dot(rT_neg1, rT_1, fneg1, f1, gneg1, g1)
            a1 = new_a1(fneg1, f1, gneg1, g1)
            a3 = new_a3(fneg1, f1, gneg1, g1)
            r = r0.x*vector(1, 0, 0) + (r0.y*cos(e) + r0.z*sin(e))*vector(0, 1, 0) + (-r0.y*sin(e) + r0.z*cos(e))*vector(0, 0, 1)
            r_dot = r0_dot.x*vector(1, 0, 0) + (r0_dot.y*cos(e) + r0_dot.z*sin(e))*vector(0, 1, 0) + (-r0_dot.y*sin(e) + r0_dot.z*cos(e))*vector(0, 0, 1)
            o_elements[t-1,0:] = OD(r, r_dot)
            total += 1
print total
print o_elements

###r and r-dot to ecliptic coordinates
##def r_(r0):
##    e = radians(23.5)
##    r = r0.x*vector(1, 0, 0) + (r0.y*cos(e) + r0.z*sin(e))*vector(0, 1, 0) + (-r0.y*sin(e) + r0.z*cos(e))*vector(0, 0, 1)
##    return r
##
##
##def r_dot_(r0_dot):
##    e = radians(23.5)
##    r_dot = r0_dot.x*vector(1, 0, 0) + (r0_dot.y*cos(e) + r0_dot.z*sin(e))*vector(0, 1, 0) + (-r0_dot.y*sin(e) + r0_dot.z*cos(e))*vector(0, 0, 1)
##    return r_dot






