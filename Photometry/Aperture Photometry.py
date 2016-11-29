#Tiffany Pan
#Aperture Photometry Redux
#7/13/15

from math import *
import pyfits
import matplotlib.pyplot as plt
import matplotlib
import warnings

warnings.filterwarnings("ignore")

print "Identify the coordinates of 5 reference stars and your asteroid, close the image display,\
and input the x and y coordinates upon request along with each star's apparent magnitude.\
The program will additionally ask how many points should be used in its centroid calculation. \n"
print "NOTES: The chosen star and their values are comments in the program\n\
The best radii (aperture, inner annulus, outer annulus): 1.75, 4, 8\n\
The asteroid magnitude takes a while to load :-)"

#IMAGE
imgData = pyfits.getdata("1948_OA.00000012.15H27M52.8S_20D07M48SS.REDUCED.FIT")
matplotlib.rcParams['image.origin'] = 'lower'
plt.imshow(imgData, vmin = imgData.mean(), vmax = 2 * imgData.mean())
plt.gray()
plt.show()
#^displays FITS image



#USER INPUT
#xy_mag1 = [float(raw_input("\nEnter x coordinate of 1st reference star: ")), float(raw_input("Enter y coordinate of 1st reference star: ")), float(raw_input("Enter apparent magnitude of 1st reference star: "))]
#xy_mag2 = [float(raw_input("\nEnter x coordinate of 2nd reference star: ")), float(raw_input("Enter y coordinate of 2nd reference star: ")), float(raw_input("Enter apparent magnitude of 2nd reference star: "))]
#xy_mag3 = [float(raw_input("\nEnter x coordinate of 3rd reference star: ")), float(raw_input("Enter y coordinate of 3rd reference star: ")), float(raw_input("Enter apparent magnitude of 3rd reference star: "))]
#xy_mag4 = [float(raw_input("\nEnter x coordinate of 4th reference star: ")), float(raw_input("Enter y coordinate of 4th reference star: ")), float(raw_input("Enter apparent magnitude of 4th reference star: "))]
#xy_mag5 = [float(raw_input("\nEnter x coordinate of 5th reference star: ")), float(raw_input("Enter y coordinate of 5th reference star: ")), float(raw_input("Enter apparent magnitude of 5th reference star: "))]
#xy_magA = [float(raw_input("\nEnter x coordinate of your asteroid: ")), float(raw_input("Enter y coordinate of your asteroid: "))]
#^user specifies (x,y) & magnitude of 5 chosen reference stars and the asteroid (placed in separate lists)

xy_mag1 = [661.02, 231.998, 15.97]
xy_mag2 = [661.02, 192.2, 15.61]
xy_mag3 = [664.811, 385.504, 15.93]
xy_mag4 = [728.298, 329.597, 16.04]
xy_mag5 = [842.953, 303.065, 15.58]
xy_magA = [631.938, 256.001, 16.16]
#^values for chosen reference stars and asteroid


#CENTROID CODE
n = sqrt(int(raw_input("\nHow many total points do you want to include in the centroid calculation (9, 25, 49, etc.)? ")))
#^user specifies number of points to use in centroid; n is set to be the sqrt of the total number inputted

def pixel_box(xy_mag):
    imgD = imgData[xy_mag[1]-(int(n/2)):xy_mag[1]+int(n/2)+1:1, xy_mag[0]-int(n/2):xy_mag[0]+int(n/2)+1:1]
    return imgD
#^function that extracts a box of pixel intensity values of specified size around the inputted x and y values
imgD1 = pixel_box(xy_mag1) #box of pixel intensities for reference star 1
Int1 = imgD1.sum() #sum of all intensity values in box for reference star 1
imgD2 = pixel_box(xy_mag2)
Int2 = imgD2.sum()
imgD3 = pixel_box(xy_mag3)
Int3 = imgD3.sum()
imgD4 = pixel_box(xy_mag4)
Int4 = imgD4.sum()
imgD5 = pixel_box(xy_mag5)
Int5 = imgD5.sum()
imgDA = pixel_box(xy_magA)
IntA = imgDA.sum()


#calculating the x coord of centroid
def xIntensity(imgD, xy_mag):
    inten1 = imgD.copy()
    for i in range(0,len(imgD)/2):
        inten1[::,i:i+1:] *= xy_mag[0] - abs((int(n/2) - i))
    inten1[::,n/2:n/2+1:] *= xy_mag[0]
    for i in range(len(imgD)/2 + 1,len(imgD)):
        inten1[::,i:i+1:] *= xy_mag[0] + abs((int(n/2) - i))
    return inten1 
xInt1 = xIntensity(imgD1, xy_mag1)
xInt2 = xIntensity(imgD2, xy_mag2)
xInt3 = xIntensity(imgD3, xy_mag3)
xInt4 = xIntensity(imgD4, xy_mag4)
xInt5 = xIntensity(imgD5, xy_mag5)
xIntA = xIntensity(imgDA, xy_magA)
#^array of intensities weighted with x value for ref stars and asteroid
cent_x1 = xInt1.sum() / Int1
cent_x2 = xInt2.sum() / Int2
cent_x3 = xInt3.sum() / Int3
cent_x4 = xInt4.sum() / Int4
cent_x5 = xInt5.sum() / Int5
cent_xA = xIntA.sum() / IntA
#^x coord centroid values for ref stars and asteroid
print cent_x1

#calculating the y coord of centroid
def yIntensity(imgD, xy_mag):
    inten2 = imgD.copy()
    for i in range(0,len(imgD)/2):
        inten2[i:i+1:,::] *= xy_mag[1] - abs((int(n/2) - i))
    inten2[n/2:n/2+1:,::] *= xy_mag[1]
    for i in range(len(imgD)/2 + 1,len(imgD)):
        inten2[i:i+1:,::] *= xy_mag[1] + abs((int(n/2) - i))
    return inten2
yInt1 = yIntensity(imgD1, xy_mag1)
yInt2 = yIntensity(imgD2, xy_mag2)
yInt3 = yIntensity(imgD3, xy_mag3)
yInt4 = yIntensity(imgD4, xy_mag4)
yInt5 = yIntensity(imgD5, xy_mag5)
yIntA = yIntensity(imgDA, xy_magA)
#^array of intensities weighted with y value for ref stars and asteroid
cent_y1 = yInt1.sum() / Int1
cent_y2 = yInt2.sum() / Int2
cent_y3 = yInt3.sum() / Int3
cent_y4 = yInt4.sum() / Int4
cent_y5 = yInt5.sum() / Int5
cent_yA = yIntA.sum() / IntA
#^y coord centroid values for ref stars and asteroid
print cent_y1


#APERTURE AREA
r_star = float(raw_input("\nEnter the value of the aperture radius: "))
#^user-chosen radius of aperture (tested with 1.75)
def pixel_boxap(cent_x, cent_y, r_star):
    imgAR = imgData[cent_y-(r_star + 5):cent_y+(r_star + 5)+1:1, cent_x-(r_star + 5):cent_x+(r_star + 5)+1:1]
    return imgAR
def apt_area(cent_x, cent_y):
    total = 0
    number_aperture = 0
    pixel_boxAP = pixel_boxap(cent_x, cent_y, r_star)
    for i in range(len(pixel_boxAP)):
        diff_y = i-float(len(pixel_boxAP)/2)
        for j in range(len(pixel_boxAP)):
            diff_x = j-float(len(pixel_boxAP)/2)
            distance = sqrt((diff_y)**2 + (diff_x)**2)
            value = pixel_boxAP[i, j]
            if distance <= r_star:
                total += value
                number_aperture += 1
    return [total, number_aperture]
#^function loops through all pixels & returns (in a list) the summed pix count & number of pixels in circular aperture of each star
#apt_area(cent_x1, cent_y1)[0] = summed pix count in apt for ref star 1
#apt_area(cent_x1, cent_y1)[1] = total number of pixels in aperture for ref star 1
#same goes for the rest of ref stars and asteroid



#ANNULUS AREA
r_Iannulus = float(raw_input("\nEnter the value of the inner annulus radius: "))
r_Oannulus = float(raw_input("\nEnter the value of the outer annulus radius: "))
def pixel_boxann(cent_x, cent_y, r_Oannulus):
    imgAN = imgData[cent_y-(r_Oannulus + 5):cent_y+(r_Oannulus + 5)+1:1, cent_x-(r_Oannulus + 5):cent_x+(r_Oannulus + 5)+1:1]
    return imgAN
#^user-chosen radius of inner annulus and outer annulus (tested with 3.5 and 8, respectively)
def annulus_area(cent_x, cent_y):
    total = 0
    number_annulus = 0
    pixel_boxAN = pixel_boxann(cent_x, cent_y, r_Oannulus)
    for i in range(len(pixel_boxAN)):
        diff_y = i-float(len(pixel_boxAN)/2)
        for j in range(len(pixel_boxAN)):
            diff_x = j-float(len(pixel_boxAN)/2)
            distance = sqrt((diff_y)**2 + (diff_x)**2)
            value = pixel_boxAN[i, j]
            if distance >= r_Iannulus:
                if distance <= r_Oannulus:
                    total += value
                    number_annulus += 1
    return [total, number_annulus]
#^function loops through all pixels & returns (in a list) the summed pix count & number of pixels in circular annulus of each star
#annulus_area(cent_x1, cent_y1)[0] = summed pix count in annulus for ref star 1
#annulus_area(cent_x1, cent_y1)[1] = total number of pixels in annulus for ref star 1
#same goes for the rest of ref stars and asteroid



#BACKGROUND VALUES
avg_backperpix1 = annulus_area(cent_x1, cent_y1)[0] / annulus_area(cent_x1, cent_y1)[1]
avg_backperpix2 = annulus_area(cent_x2, cent_y2)[0] / annulus_area(cent_x2, cent_y2)[1]
avg_backperpix3 = annulus_area(cent_x3, cent_y3)[0] / annulus_area(cent_x3, cent_y3)[1]
avg_backperpix4 = annulus_area(cent_x4, cent_y4)[0] / annulus_area(cent_x4, cent_y4)[1]
avg_backperpix5 = annulus_area(cent_x5, cent_y5)[0] / annulus_area(cent_x5, cent_y5)[1]
avg_backperpixA = annulus_area(cent_xA, cent_yA)[0] / annulus_area(cent_xA, cent_yA)[1]
#^sets the average background-per-pixel value for each star & asteroid by dividing pix count by num of pixels



#STAR & ASTEROID VALUES
star1 = apt_area(cent_x1, cent_y1)[0] - avg_backperpix1
star2 = apt_area(cent_x2, cent_y2)[0] - avg_backperpix2
star3 = apt_area(cent_x3, cent_y3)[0] - avg_backperpix3
star4 = apt_area(cent_x4, cent_y4)[0] - avg_backperpix4
star5 = apt_area(cent_x5, cent_y5)[0] - avg_backperpix5
Asteroid = apt_area(cent_xA, cent_yA)[0] - avg_backperpixA
#^sets the pixel count for each star & asteroid by subtracting background from star
mag1 = xy_mag1[2]
mag2 = xy_mag2[2]
mag3 = xy_mag3[2]
mag4 = xy_mag4[2]
mag5 = xy_mag5[2]
#^known magnitudes for each star as variables (from earlier lists)



#CALCULATING CONSTANT
def const_using_mag(mag, star):
    const = mag + 2.5 * log10(star)
    return const
#^function returns calculated constant for each star using magnitude and pixel count

def constant(mag1, mag2, mag3, mag4, mag5):
    avg_const = (const_using_mag(mag1, star1)+const_using_mag(mag2, star2)+
                 const_using_mag(mag3, star3)+const_using_mag(mag4, star4)+
                 const_using_mag(mag5, star5)) / 5
    return avg_const
#^function averages the five constants calculated for each star to obtain one average value
avg_constant = constant(mag1, mag2, mag3, mag4, mag5)
#^avg_constant value



#CALCULATING MAGNITUDE
def magnitude(obj, avg_constant):
    magnitude = -2.5 * log10(obj) + avg_constant
    return magnitude
#^function calculates magnitude of asteroid using pixel count and average constant
print "\nThe magnitude of the asteroid: ",
print magnitude(Asteroid, avg_constant)


#CALCULATING BACKGROUND MAG
def mag_arcsec(backperpix):
    avg_backperarcsec2 = backperpix * (.59787 ** 2)
    mags_arcsec2 = magnitude(avg_backperarcsec2, avg_constant)
    return mags_arcsec2
#^returns mag of background sky by converting background per pix value to background per arcsec^2 and then finding the magnitude
print "\nThe magnitude of the background sky (in mags-per-square-arcsecond): ",
print mag_arcsec(avg_backperpixA)
