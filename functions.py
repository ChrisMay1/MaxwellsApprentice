from __future__ import division
import numpy as np
import math

def bessel(x,y,kr,p):
    ''' Solves for the Bessel function you desire '''
    r = np.sqrt(x*x + y*y)
    J = 0
    mmax = 25
    for m in range(mmax):
        J += ((kr*r/2)**(2*m+p))*((-1)**(m))/(math.factorial(m)*math.gamma(m+p+1))
    return J

def E_r(x,y,z,k_r,k_z,omega,l,t,alpha):
	phi = np.arctan(y/x)
	exp = np.cos(alpha)*np.exp(1j*(k_z*z-omega*t-l*phi))
	if l == 0:
		Jsum = -bessel(x,y,k_r,1)
	else:
		Jsum = (bessel(x,y,k_r,l+1) - bessel(x,y,k_r,l-1))/2
	return Jsum*exp
    
def E_phi(x,y,z,k_r,k_z,omega,l,t,alpha):
	phi = np.arctan(y/x)
	r = np.sqrt(x*x + y*y)
	exp = np.cos(alpha)*np.exp(1j*(k_z*z-omega*t-l*phi))
	J = bessel(x,y,k_r,l)
	return 1j*l*J*exp/(r*k_r)
    
def E_z(x,y,z,k_r,k_z,omega,l,t,alpha):
	phi = np.arctan(y/x)
	return 1j*np.sin(alpha)*bessel(x,y,k_r,l)*np.exp(1j*(k_z*z-omega*t-l*phi))

def B_r(x,y,z,k_r,k_z,omega,l,t,alpha):
	phi = np.arctan(y/x)
	r = np.sqrt(x*x+y*y)
	exp = np.exp(1j*(k_z*z-omega*t-l*phi))
	return bessel(x,y,k_r,l)*exp
	
def B_phi(x,y,z,k_r,k_z,omega,l,t,alpha):
	phi = np.arctan(y/x)
	exp = np.exp(1j*(k_z*z-omega*t-l*phi))
	if l == 0:
		Jsum = -bessel(x,y,k_r,1)
	else:
		Jsum = (bessel(x,y,k_r,l+1) - bessel(x,y,k_r,l-1))/2
	return Jsum*exp
	
def energy_density(x,y,k_r,k_z,omega,l,alpha):
	term1 = np.sin(alpha)*bessel(x,y,k_r,l)**2
	if l == 0:
		term2 = (1 + np.cos(alpha)*np.cos(alpha))*(2*bessel(x,y,k_r,1)**2)
	else:
		term2 = (1 + np.cos(alpha)*np.cos(alpha))*(bessel(x,y,k_r,l+1)**2 +bessel(x,y,k_r,l-1)**2)
	return term1*term1/(4*np.pi) + term2/(16*np.pi)
	
