"""
.. versionadded:: 1.1.0
   This demo depends on new features added to contourf3d.
"""

from subprocess import Popen, PIPE, check_call, call
import shlex

from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
from matplotlib import cm
from scipy.interpolate import griddata
import numpy as np

import sys

ALPHA = 0
BETA = 1
ENERGY = 2
KINETIC = 3
POTENTIAL = 4
VARIANCE = 5
ERROR = 6


filename = str(sys.argv[1])
minEnergy = 0
maxEnergy = 210000000
#varianceCutoff = float(sys.argv[2])
#errorCutoff = float(sys.argv[2])


#print filename

import csv
with open(filename) as file:
	reader = csv.reader(file,delimiter=' ')
	parsedFile = []
	alpha1 = []
	beta1 = []
	energy1 = []
	kinetic1 = []
	potential1 = []
	variance1 = []
	error1 = []
	virialRatio1 = []
	for row in reader:
		parsedRow = []
		for cell in row:
			if (cell != ''):
				parsedRow.append(float(cell))
		parsedFile.append(parsedRow)
		
		alpha1.append(parsedRow[ALPHA])
		beta1.append(parsedRow[BETA])

		kinetic1.append(parsedRow[KINETIC])
		potential1.append(parsedRow[POTENTIAL])
		virialRatio1.append(parsedRow[KINETIC] / parsedRow[POTENTIAL])

		energy1.append(parsedRow[ENERGY])
		variance1.append(parsedRow[VARIANCE])
		error1.append(parsedRow[ERROR])


variance1.sort()
error1.sort()
virialRatio1.sort()

maxVariance = 1000000000000000
maxError    = 1000000000000000

#print "Max variance and error:"
#print int(len(variance1)*varianceCutoff)
#print len(variance1)
#maxVariance = variance1[int(len(variance1)*varianceCutoff)]
#print maxVariance
#maxError = error1[int(len(error1)*errorCutoff)]
#print maxError

import random

with open(filename) as file:
	reader = csv.reader(file,delimiter=' ')
	parsedFile = []
	alpha = []
	beta = []
	energy = []
	kinetic = []
	potential = []
	variance = []
	error = []
	virialRatio = []
	for row in reader:
		parsedRow = []
		for cell in row:
			if (cell != ''):
				parsedRow.append(float(cell))
#		parsedFile.append(parsedRow)
		ratio = parsedRow[KINETIC] / parsedRow[POTENTIAL]
#		if (parsedRow[ENERGY] < 20):
#			print parsedRow

		if (abs(parsedRow[ERROR]) < maxError and abs(parsedRow[VARIANCE]) < maxVariance and parsedRow[ENERGY] > minEnergy and parsedRow[ENERGY] < maxEnergy and parsedRow[KINETIC] > 0 ):
			alpha.append(parsedRow[ALPHA])
			beta.append(parsedRow[BETA])
			energy.append(parsedRow[ENERGY])
			kinetic.append(parsedRow[KINETIC])
			potential.append(parsedRow[POTENTIAL])
			virialRatio.append(ratio)
			variance.append(parsedRow[VARIANCE])
			error.append(parsedRow[ERROR])
virialRatio.sort()
#print " "
#print " "
#print "Virial Ratio: "
#print virialRatio
m = min(energy)
#for i in range(len(virialRatio)):
#	if (energy[i] == m):
#		print " "
#	print i, "  ", kinetic[i], " & ", potential[i], " & ", virialRatio[i], " :: ", energy[i]
#print min(virialRatio)
#print energy
#k = 39
#print kinetic[k]
#print potential[k]
#print virialRatio[k]
print energy
print m



"""
data = [(60, 5, '121'), (61, 5, '103'), (62, 5, '14.8'), (63, 5, '48.5'), (64, 5, '57.5'), (65, 5, '75.7'), (66, 5, '89.6'), (67, 5, '55.3'), (68, 5, '63.3'), (69, 5, '118'), (70, 5, '128'), (71, 5, '105'), (72, 5, '115'), (73, 5, '104'), (74, 5, '134'), (75, 5, '123'), (76, 5, '66.3'), (77, 5, '132'), (78, 5, '145'), (79, 5, '115'), (80, 5, '38.2'), (81, 5, '10.4'), (82, 5, '18.4'), (83, 5, '87'), (84, 5, '86.7'), (85, 5, '78.9'), (86, 5, '89.9'), (87, 5, '108'), (88, 5, '57.1'), (89, 5, '51.1'), (90, 5, '69.1'), (91, 5, '59.8'), (60, 6, '48.9'), (61, 6, '33.3'), (62, 6, '-19.2'), (63, 6, '-17.5'), (64, 6, '-6.5'), (65, 6, '75.7'), (66, 6, '89.6'), (67, 6, '55.3'), (68, 6, '99.8'), (69, 6, '156'), (70, 6, '141'), (71, 6, '54.1'), (72, 6, '66.1'), (73, 6, '98.9'), (74, 6, '155'), (75, 6, '146'), (76, 6, '111'), (77, 6, '132'), (78, 6, '145'), (79, 6, '97.3'), (80, 6, '101'), (81, 6, '59.4'), (82, 6, '70.4'), (83, 6, '142'), (84, 6, '145'), (85, 6, '140'), (86, 6, '56.9'), (87, 6, '77.8'), (88, 6, '21.1'), (89, 6, '27.1'), (90, 6, '48.1'), (91, 6, '41.8')]
x, y, z = zip(*data)
z = map(float, z)
grid_x, grid_y = np.mgrid[min(x):max(x):100j, min(y):max(y):100j]
grid_z = griddata((x, y), z, (grid_x, grid_y), method='cubic')

print z
print " "
print " "
rpint grid_z
"""




"""
print alpha
print beta
grid_alpha, grid_beta = np.mgrid[min(alpha):max(alpha):30j  , min(beta):max(beta):30j]
grid_energy = griddata((alpha,beta), energy, (grid_alpha, grid_beta), method='cubic')

fig = plt.figure()
ax = fig.gca(projection='3d')
#fig,ax = plt.subplots()
#cset = ax.contourf(grid_alpha, grid_beta, grid_energy, 100, zdir='z', offset=min(energy) - 0.1*(max(energy)-min(energy)), cmap=cm.coolwarm  )

ax.plot_surface(grid_alpha, grid_beta, grid_energy, rstride=1, cstride=1, linewidth=0.1, alpha=0.7, cmap=cm.coolwarm)

A = []
for i in range(len(alpha)):
	for j in range(len(alpha)):
		A.append(alpha[i])
#ax.plot(alpha,energy)
#print min(energy)
ax.set_title("6N, 2D, repulsion, Jastrow, omega=0.1\n300k cycles, min(energy)="+str(min(energy)))
#plt.savefig(filename+".png")

#cset = ax.contourf(grid_alpha, grid_beta, grid_energy, zdir='x', offset=max(alpha) + 0.1*(max(alpha)-min(alpha)), cmap=cm.coolwarm  )
#cset = ax.contourf(grid_alpha, grid_beta, grid_energy, zdir='y', offset=max(beta), cmap=cm.coolwarm  )

ax.set_xlabel('Alpha')
ax.set_ylabel('Beta')
ax.set_zlabel('Energy')

if (int(sys.argv[2]) == 1):
	for ii in xrange(0,360,1):
		print ii
		ax.view_init(elev=30., azim=ii)
		plt.savefig("imgs/"+filename+"_"+str(ii)+".png")
command = "ffmpeg -i imgs/"+str(filename)+"_%d.png "+str(filename)+".mp4"
args = shlex.split(command)
print args
Popen(args)


plt.show()
"""






"""
fig = plt.figure()
ax = fig.gca(projection='3d')
X, Y, Z = axes3d.get_test_data(0.05)
ax.plot_surface(X, Y, Z, rstride=8, cstride=8, alpha=0.3)
cset = ax.contourf(X, Y, Z, zdir='z', offset=-100, cmap=cm.coolwarm)
cset = ax.contourf(X, Y, Z, zdir='x', offset=-40, cmap=cm.coolwarm)
cset = ax.contourf(X, Y, Z, zdir='y', offset=40, cmap=cm.coolwarm)

ax.set_xlabel('X')
ax.set_xlim(-40, 40)
ax.set_ylabel('Y')
ax.set_ylim(-40, 40)
ax.set_zlabel('Z')
ax.set_zlim(-100, 100)

plt.show()
"""