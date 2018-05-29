from __future__ import division
from functions import bessel,E_r,E_phi,E_z,energy_density
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
from pylab import *
from matplotlib.widgets import Slider
import argparse


''' Everything is in Gaussian Units, i.e. c = 1'''
lDefault = 2
omegaDefault = 1
c=1
krat = [1/np.sqrt(2),1/np.sqrt(2)]

''' The k vectors determine the propagation in radial and z directions.
    Their (omega/c) = k_f = sqrt(k_r^2 + k_z^2)
    I've chosen a 1/sqrt(2) term so that the propagate equally in both 
    directions '''

ds = 0.1
x = np.arange(-20,20,ds)
y = np.arange(-20,20,ds)
mx,my = np.meshgrid(x,y)

def besselPlot():
    ''' Shows a 2D plot of the Bessel Function as a Function of x'''
    plt.figure()
    plt.plot(x,bessel(x,0,k_r,l))
    plt.xlabel('x (m)')
    plt.ylabel('Bessel Intesity')
    plt.title(r'Bessel Function $J_%d(r)$'%l)
    plt.show()
    
def showEz():
    ''' Shows the Real Part of the Z component of the Electric Field,
		to change this just change the value for z = np.real(E...) in the drawE function.
		
		I don't really like it all that much. '''
    ''' Contour Plot Styles'''
    contourfillflag = True
    contourfillonlyflag = False
    contourlevels = 20
    contourlabelsflag = False
    contourlevelsbarflag = True
    gridflag = True
    
    fig,ax = plt.subplots()
    def drawE_z(i):
        z = np.real(E_z(mx,my,0,k_r,k_z,omega,l,i,alpha))
        extent = (-20, 20, -20, 20)
        ax.imshow(z, extent=extent)
        if contourfillflag:
            if contourlevels == True:
                if not contourfillonlyflag:
                    c1 = ax.contour(x,y,z,colors='k')
                c2 = ax.contourf(x,y,z)
            else:
                if not contourfillonlyflag:
                    c1 = ax.contour(x,y,z,contourlevels,colors='k')
                c2 = ax.contourf(x,y,z,contourlevels)
            if contourlevelsbarflag and i == 0:
                cb2 = fig.colorbar(c2)
        else:
            if contourlevels == True:
                c1 = ax.contour(x,y,z)
            else:
                c1 = ax.contour(x,y,z,contourlevels)
            if contourlevelsbarflag and i == 0:
                fig.colorbar(c1)
        if contourlabelsflag and ((not contourfillflag) or (not contourfillonlyflag)) and i == 0:
            l1 = ax.clabel(c1)

    def animateE_z(i):
        ax.clear()
        drawE_z(i)
        return ax,
        
    drawE_z(0)
    ani = animation.FuncAnimation(fig, animateE_z, np.arange(1, 200), interval=1, blit=True)
    plt.show()


def slider(fieldname,l):
	fig, ax = plt.subplots()
	dt = 1e-1
	time = np.arange(0.0, 10.0, dt)
	ptype = fieldname.lower()
	
	if ptype == 'ez':
		title = (r'$E_z$')
		f = np.zeros([len(time),len(mx),len(my)])
		for t in range(len(time)):
			f[t,:,:] = np.real(E_z(mx,my,0,k_r,k_z,omega,l,t*dt,alpha))
	# In case I want to add in more sliders

	cs=ax.contourf(mx,my,f[0,:,:])
	contour_axis = plt.gca()
	contour_axis.set_title(title)
	cbar = fig.colorbar(cs)
	contour_axis.set_aspect('equal')

	axmax = plt.axes([0.2, 0.01, 0.65, 0.03])  #slider location and size
	stime = Slider(axmax, 'Time',0, time.max(),0,'%1.2f',)  #slider properties

	def update(val):
		contour_axis.clear()
		contour_axis.contourf(mx,my,f[int(stime.val/dt),:,:])
		contour_axis.set_title(r'$E_z$')
		plt.draw()                   
	stime.on_changed(update)

	plt.show()

def avg_energy_density():
	u = energy_density(mx,my,k_r,k_z,omega,l,alpha)
	fig,ax = plt.subplots()
	cax = ax.contourf(mx,my,u)
	ax.set_aspect('equal')
	ax.set_title(r'Time Averaged Energy Density $J_%d$'%(l)+str('\n')+\
		str(r'$\omega: %.1f, c:%.1f$'%(omega,c)))
	cbar = fig.colorbar(cax)
	ax.grid()
	plt.show()
	
parser = argparse.ArgumentParser(description='Which plots to make.')
parser.add_argument('-p','--plot', metavar='PLOT', action ="store",
                    default='', help='the plot to make')
parser.add_argument('-l','--lterm', metavar = 'L',type = int, nargs = "+",
                    default = lDefault, help = "What l value do you want to test.")
parser.add_argument('-w','--omega', metavar = 'OMEGA',type = float, nargs = "+",
					default = omegaDefault, help = "Angular Frequency, should be less than 1.")              

args = parser.parse_args()
plot = args.plot
l = args.lterm
omega = args.omega

if omega == omegaDefault:
	k_r = (omega/c)*krat[0]
	k_z = (omega/c)*krat[1]
else:
	k_r = (omega[0]/c)*krat[0]
	k_z = (omega[0]/c)*krat[1]
	omega = omega[0]
if l != lDefault:
	l = l[0]
alpha = np.arctan(k_r/k_z)

print("\nPlot: %s"%args.plot)
if np.sqrt(krat[0]*krat[0]+krat[1]*krat[1]) > 1:
	print('Fix your kratio bro')

if args.plot.lower() == '':
	print('This is my help menu :)')
	print('If you do not know how to run this script you have come to the right place.\n')
	print('The plot types are\n1.\tbessel\n2.\tintensity\n3.\tslider\n')
	print('Just type -p and the name of the plot you want to make.\n')
	print('You can also change the l term with the -l [int] command.\n')
	print('You can also change the angular frequency with the -w [float] command.\n')
	print('Ex. python runfunctions.py - p slider -l 6 -w 0.23 ')
elif args.plot.lower() == 'bessel':
	besselPlot()
	print('This is just the Bessel Function')
elif args.plot.lower() == 'intensity':
	avg_energy_density()
elif args.plot.lower() == 'slider':
	if l == lDefault:
		slider('ez',l)
	else:
		slider('ez',l[0])
elif args.plot.lower() == 'showez':
	showEz()
print('\n\nAngular Momentum/l term:%d'%l)
print('Angular Frequency:%.3f'%omega)
print('k_r: %.3f, k_z: %.3f'%(krat[0],krat[1]))

