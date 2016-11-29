#Calculating Jackknife Uncertainties
import numpy as np
from math import *

teamData = np.array([[1.38185773447,0.450131251259,9.26858905551,274.095364618,125.605511774,144.287025793],
          [1.36274838244,0.437235514259,9.31467701989,274.027361239,127.63030769,136.835532108],
          [1.34641348317,0.418046571591,9.20339738992,274.197027382,128.641885809,130.319248216],
          [1.31159809825,0.401313596534,9.08612645232,274.36763645,131.582549192,118.176785391],
          [1.30305051024,0.393222941684,8.89168489932,274.665398878,131.595829935,115.647880947]])
teamAverage = np.array([0.0,0.0,0.0,0.0,0.0,0.0])
result = np.array([0.0,0.0,0.0,0.0,0.0,0.0])
def column(matrix, i):
    return [row[i] for row in matrix]
print column(teamData, 0)
for i in range(0, len(teamData)+1):
    teamAverage += column(teamData, i)[i-1]
    teamAverage = teamAverage/5
for i in range(0, len(teamData)+1):
    result += (column(teamData, i)[i-1]-teamAverage)**2
    result = result*4/5
for i in range(0, len(teamData)+1):
    result[i] = sqrt(result[i])
#Uncertainty in terms of standard deviation
print "Uncertainty of Semimajor Axis: " +str(result[0])
print "Uncertainty of Eccentricity: " +str(result[1])
print "Uncertainty of Inclination: " +str(result[2])
print "Uncertainty of Longitude of Ascending Node: "+str(result[3])
print "Uncertainty of Argument of Perihelion: " +str(result[4])
print "Uncertainty of Mean Anomaly: " +str(result[5])
