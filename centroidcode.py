# Problem 2: Semi-Automated Centroid Calculation
from math import *
import numpy as np

# part b
import pyfits
import matplotlib.pyplot as plt
import matplotlib
import warnings

warnings.filterwarnings("ignore")

print "To calculate the centroid of an object, identify its coordinates, close the image display,\
and input the x and y coordinates upon request. The program will additionally ask how many points \
should be used in its centroid calculation. \n"


imgData = pyfits.getdata("1948_OA.00000012.15H27M52.8S_20D07M48SS.REDUCED.FIT")
matplotlib.rcParams['image.origin'] = 'lower'
plt.imshow(imgData, vmin = imgData.mean(), vmax = 2 * imgData.mean())
plt.gray()
#plt.show()

n = 5 #sqrt(int(raw_input("How many total points do you want to include in the centroid calculation (9, 25, 49, etc.)? ")))

x = int(float(raw_input("Enter x-coordinate: ")))
y = int(float(raw_input("Enter y-coordinate: ")))

imgD = imgData[y-(int(n/2)):y+int(n/2)+1:1,x-int(n/2):x+int(n/2)+1:1]
Int = imgD.sum()
print imgD

# calculating the x coord of centroid
def xIntensity(imgD):
    inten1 = imgD.copy()
    for i in range(0,len(imgD)/2):
        inten1[::,i:i+1:] *= x - abs((int(n/2) - i))
    inten1[::,n/2:n/2+1:] *= x
    for i in range(len(imgD)/2 + 1,len(imgD)):
        inten1[::,i:i+1:] *= x + abs((int(n/2) - i))
    return inten1 
xInt = xIntensity(imgD)
cent_x = xInt.sum() / Int
print cent_x # x coord of centroid


# calculating the y coord of centroid
def yIntensity(imgD):
    inten2 = imgD.copy()
    for i in range(0,len(imgD)/2):
        inten2[i:i+1:,::] *= y - abs((int(n/2) - i))
    inten2[n/2:n/2+1:,::] *= y
    for i in range(len(imgD)/2 + 1,len(imgD)):
        inten2[i:i+1:,::] *= y + abs((int(n/2) - i))
    return inten2
yInt = yIntensity(imgD)
cent_y = yInt.sum() / Int
print cent_y # y coord of centroid


# x uncertainty
def get_uncertaintyX(cent_x, Int, imgD):
    intenX = imgD.sum(axis = 0)
    numer = 0
    for i in range(0,len(imgD)):
        distance = cent_x - (x - (int(n/2) - i))
        distance = distance ** 2
        numer += intenX[i] * distance
    uncertainty = sqrt(numer / (Int * (Int-1)))
    return uncertainty
x_uncert = get_uncertaintyX(cent_x, Int, imgD)
print x_uncert

#y uncertainty
def get_uncertaintyY(cent_y, Int, imgD):
    intenY = imgD.sum(axis = 1)
    numer = 0
    for i in range(0,len(imgD)):
        distance = cent_y - (y - (int(n/2) - i))
        distance = distance ** 2
        numer += intenY[i] * distance
    uncertainty = sqrt(numer / (Int * (Int-1)))
    return uncertainty
y_uncert = get_uncertaintyY(cent_y, Int, imgD)
print y_uncert

#total uncertainty
def total_uncert(x_uncert, y_uncert):
    uncert = sqrt(x_uncert**2 + y_uncert**2)
    return uncert
print total_uncert(x_uncert, y_uncert)



