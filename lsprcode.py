from math import *
from decimal import Decimal 
import numpy as np
import pyfits
import matplotlib.pyplot as plt
import matplotlib
import warnings

warnings.filterwarnings("ignore")

print "To find the right ascension and declination of an object, identify its (x,y) coordinates, \
close the image display, and input the x and y coordinates upon request. \
The program will additionally ask for a file of reference stars, the center of the field of view of the image,\
and if you want to use the flattened RA and Dec or not. \n"

imgData = pyfits.getdata("1948_OA.00000012.15H27M52.8S_20D07M48SS.REDUCED.FIT")
matplotlib.rcParams['image.origin'] = 'lower'
plt.imshow(imgData, vmin = imgData.mean(), vmax = 2 * imgData.mean())
plt.gray()
plt.show()
#^imports and displays the FITS image

x = float(raw_input("Enter x-coordinate: "))
y = float(raw_input("Enter y-coordinate: "))
#^asks for x and y coordinates of an object in the FITS file

reference_stars = open(raw_input("\nEnter file name(referencestars.txt): ")) #tested with "referencestars.txt"
reference_stars = reference_stars.readlines()
N = len(reference_stars) #sets number of reference stars to N

center_A = radians(15.0 * float(raw_input("\nEnter the RA of center of field of view (in hours, ex 15.72386): "))) #tested with 15.72386
center_D = radians(float(raw_input("Enter the Dec of center of field of view (in degrees, ex -23.97258): "))) #tested with -23.97258
#^coordinates of center of field (in radians)

L = 3910 * .6 / .0068
#^effective focal length with reducer expressed in # of pixels
x_total = 0
y_total = 0
xsqr_total = 0
ysqr_total = 0
xy_total = 0
ra_total = 0
rax_total = 0
ray_total = 0
dec_total = 0
decx_total = 0
decy_total = 0
#^original totals for each variable set to 0 to be added to in for loop

answer = raw_input("\nDo you want to use a flattened RA and Dec? ").lower()
#^allows user to choose if they want to use flat alpha/delta or not

sum_sigma_a = 0
sum_sigma_d = 0

for line in reference_stars:
    line = line.strip("\n") #takes away line breaks 
    line = line.split(',') #removes commas and separates each item
    for i in range(len(line)):
        line[i] = float(line[i]) #turns each item to number
    for i in range(2, 3):
        line[i] = radians(float(line[i])*15)
    for i in range(3, len(line)):
        line[i] = radians(float(line[i])) #turns RA and Dec into radians to be used in trig calculations
    if answer == "yes":
        numer_1 = cos(line[3])*sin(line[2]-center_A)
        H = sin(line[3])*sin(center_D)+cos(line[3])*cos(center_D)*cos(line[2]-center_A)
        numer_2 = sin(line[3])*cos(center_D)-cos(line[3])*sin(center_D)*cos(line[2]-center_A)
        flat_a = numer_1 / H - line[0] / L 
        flat_d = numer_2 / H - line[1] / L #need to convert ra and dec into flat alpha and delta
        x_total += line[0] #sums all x values
        y_total += line[1] #sums all y values
        xy_total += (line[0] * line[1]) #sums all x's times y's
        ra_total += flat_a #sums all flat RAs
        dec_total += flat_d #sums all flat Decs
        xsqr_total += (line[0] ** 2) #sums all x squares
        ysqr_total += (line[1] ** 2) #sums all y squares
        rax_total += (flat_a * line[0]) #sums all flat RAs times corresponding x's
        ray_total += (flat_a * line[1]) #sums all flat RAs times corresponding y's
        decx_total += (flat_d * line[0]) #sums all flat Decs times corresponding x's
        decy_total += (flat_d * line[1]) #sums all flat Decs times corresponding y's
    #^flattens alpha and delta
    if answer == "no":
        x_total += line[0] #sums all x values
        y_total += line[1] #sums all y values
        xy_total += (line[0] * line[1]) #sums all x's times y's
        ra_total += line[2] #sums all flat RAs
        dec_total += line[3] #sums all flat Decs
        xsqr_total += (line[0] ** 2) #sums all x squares
        ysqr_total += (line[1] ** 2) #sums all y squares
        rax_total += (line[2] * line[0]) #sums all flat RAs times corresponding x's
        ray_total += (line[2] * line[1]) #sums all flat RAs times corresponding y's
        decx_total += (line[3] * line[0]) #sums all flat Decs times corresponding x's
        decy_total += (line[3] * line[1]) #sums all flat Decs times corresponding y's
    #^does not flatten alpha or delta
    if answer != "yes" and answer != "no":
        print "None" #if user inputs an answer that isn't no or yes

#alpha equation
alpha = np.array([[ra_total], [rax_total], [ray_total]])

matrix_alpha = np.array([[N, x_total, y_total],
                         [x_total, xsqr_total, xy_total],
                         [y_total, xy_total, ysqr_total]])
#^3 x 3 matrix with sums
inverse_matrix_alpha = np.linalg.inv(matrix_alpha)

print "\nThe following are the values of the plate constants b1, a11, and a12 in a 3 x 1 matrix:"
alpha_elements = np.dot(inverse_matrix_alpha, alpha)
print alpha_elements #multiplied alpha by inverse of matrix_alpha to obtain [b1, a11, and a12] in 3 x 1 matrix

#delta equation
delta = np.array([[dec_total], [decx_total], [decy_total]])

matrix_delta = np.array([[N, x_total, y_total],
                         [x_total, xsqr_total, xy_total],
                         [y_total, xy_total, ysqr_total]])
#^3 x 3 matrix with sums
inverse_matrix_delta = np.linalg.inv(matrix_delta)

print "\nThe following are the values of the plate constants b2, a21, and a22 in a 3 x 1 matrix:"
delta_elements = np.dot(inverse_matrix_delta, delta)
print delta_elements #multiplied delta by inverse of matrix_delta to obtain [b2, a21, and a22] in 3 x 1 matrix

sum_sigma_a = 0
sum_sigma_d = 0
#unflattening flat alpha and delta to calculate uncertainties
for line in reference_stars:
    line = line.strip("\n") #takes away line breaks 
    line = line.split(',') #removes commas and separates each item
    for i in range(len(line)):
        line[i] = float(line[i]) #turns each item to number
    for i in range(2, 3):
        line[i] = radians(float(line[i])*15)
    for i in range(3, len(line)):
        line[i] = radians(float(line[i])) #turns RA and Dec into radians to be used in trig calculations
    if answer == "yes":
        flat_alpha = alpha_elements[0] + alpha_elements[1]*x + alpha_elements[2]*y #calculating flat alpha
        flat_delta = delta_elements[0] + delta_elements[1]*x + delta_elements[2]*y #calculating flat delta
        d = cos(center_D) - flat_delta * sin(center_D) #delta used to unflatten reference stars/asteroid
        gamma = sqrt((flat_alpha**2) + d**2) #gamma used to unflatten reference stars/asteroid
        unflat_alpha = degrees(center_A + atan(-flat_alpha/d)) / 15 #unflat alpha in decimal hours
        unflat_delta = degrees(atan((sin(center_D) + flat_delta * cos(center_D))/gamma)) #unflat delta in decimal degrees
        sum_sigma_a += (degrees(line[2])/15-unflat_alpha)**2 #sum for alpha uncertainty with flattening
        sum_sigma_d += (degrees(line[3])-unflat_delta)**2 #sum for delta uncertainty with flattening
    if answer == "no":
        fit_alpha = degrees(alpha_elements[0] + alpha_elements[1]*x + alpha_elements[2]*y) / 15 #calculating fit alpha in decimal hours
        fit_delta = degrees(delta_elements[0] + delta_elements[1]*x + delta_elements[2]*y) #calculating fit delta in decimal degrees
        sum_sigma_a += (degrees(line[2])/15-fit_alpha)**2 #sum for alpha uncertainty without flattening
        sum_sigma_d += (degrees(line[3])-fit_delta)**2 #sum for delta uncertainty without flattening
if answer == "yes":
    print "\nRight ascension of the object: ",
    print unflat_alpha
    print "Declination of the object: ",
    print unflat_delta
if answer == "no":
    print "\nRight ascension of the object: ", 
    print fit_alpha
    print "Declination of the object: ", 
    print fit_delta
sigma_alpha = sqrt(sum_sigma_a/(N-3)) #alpha uncertainty
sigma_delta = sqrt(sum_sigma_d/(N-3)) #delta uncertainty
print "\nThe following are the RA and Dec uncertainties, respectively: "
print sigma_alpha
print sigma_delta
#^prints uncertainties

