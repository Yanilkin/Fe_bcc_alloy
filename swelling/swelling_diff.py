import matplotlib.pyplot as plt
import numpy as np
import math
from scipy.integrate import odeint

import sys
sys.path.append('/home/yanilkin/spparks_code_generator/code')
import defect_properties as defpr

volume = 1.000000e-18 # in Litres
volume = volume/1000
a = 3.14e-10 # lattice constant in m
b = math.sqrt(3)/2*a # Burger's vector sqrt(3)/2*a
Va = (a*a*a)/2.0 # Atomic volume
na = volume/Va # Number of atoms
pi = 3.14 # Pui number

#Parameters
G = 5e-4
Ks = 1.2e13*0 # Surface sink
Kd = 5e13 # Dislocation network sink
Kiv = 4*3.14*(3.1e-10+2.9e-10)/Va # Recombination rate constant
zid = 3.0 # Bias factor
zvd = 1.0 # Bias factor
#Cp = 6e23 # From Cockerman for T=300C
Cp = 1e21*1e0 # From Igata for T=1150K 	
Ro = 0.2 # Initial radius in nanometers


Ediss = 2.21 #v5c2
N = 4
v1C = 1e-4

def Difi(T):
	return 1e-8*np.exp(-0.4*11600/T)

def Difv(T):
	E_vac_mig = 1.3	#1.1	#eV
	d_vac = a*math.sqrt(3.)/2
	nu = 1.5e13#*55	#3.6e12
	return 8./6*d_vac**2*nu*np.exp(-E_vac_mig*11600/T)

def zib(R):
#	return 1+4*a/R/1e-9
	R = R - Ro
#	return 2.5-1.0*R*R/(1000+R*R)
	return 1.0	

def cv(T,R):
	Ki = Ks+ zid*Kd + zib(R)*4*pi*(R*1e-9)*Cp # Radius in nanometers
	Kv = Ks+ zvd*Kd + 4*pi*(R*1e-9)*Cp # Radius in nanometers
	Di = Difi(T)
	Dv = Difv(T)
	return np.sqrt(np.power(Ki*Di/(2*Kiv*(Di+Dv)),2)+G*Ki*Di/(Kiv*(Di+Dv)*Kv*Dv))-Ki*Di/(2*Kiv*(Di+Dv))

def ci(T,R):
	Ki = Ks+ zid*Kd + zib(R)*4*pi*(R*1e-9)*Cp # Radius in nanometers
	Kv = Ks+ zvd*Kd + 4*pi*(R*1e-9)*Cp # Radius in nanometers
	Di = Difi(T)
	Dv = Difv(T)
	return np.sqrt(np.power(Kv*Dv/(2*Kiv*(Di+Dv)),2)+G*Kv*Dv/(Kiv*(Di+Dv)*Ki*Di))-Kv*Dv/(2*Kiv*(Di+Dv))

def cveq(T,R):
	Ef = 3.0 #eV
	return 0*np.exp(-Ef*11600/T)

#T = 2000
#R = 1
#print 1/R*1e9*(Difv(T)*(cv(T,R)-0*cveq(T))-Difi(T)*ci(T,R))

def f(R,t):
	T = 900
	"""Grow rate equation for single pore"""
	dR = 1.0/R*1e18*(Difv(T)*(cv(T,R)-cveq(T,R))-zib(R)*Difi(T)*ci(T,R)) # Radius in nanometers
	return dR

if __name__ == "__main__":
	# Time integration of grow equation
	t = np.linspace(0,100/G,10000)
	results = odeint(f,Ro,t)
	R = results[:,0]
	#R = np.linspace(1,10,100)
	#T = 1150
	plt.plot(t*G,100*4.0*pi/3*np.power(R*1e-9,3)*Cp,color='b') #Swelling in % Initial radius 1.0 nm
	#plt.plot(t*G,100*1.0/3*np.power(R*1e-9,3)*Cp,color='b') #Swelling in % Initial radius 1.0 nm
	#plt.plot(t*G,Difv(T)*cv(T,R)-Difi(T)*ci(T,R),color='m',label='Dvcv-Dici')
	#plt.plot(R,(Difv(T)*cv(T,R)-Difi(T)*ci(T,R)),color='m',label='Dvcv-Dici')

	print t[1.0/100*10000]*G, 100*4.0*pi/3*np.power(R[1.0/100*10000]*1e-9,3)*Cp
	print t[10.0/100*10000]*G, 100*4.0*pi/3*np.power(R[10.0/100*10000]*1e-9,3)*Cp
	print t[99.0/100*10000]*G, 100*4.0*pi/3*np.power(R[99.0/100*10000]*1e-9,3)*Cp
	#plt.text(1e-3,1e-6,'Gt',fontsize=14)
	#plt.yscale('log')
	plt.xlabel('Dose, dpa')
	plt.ylabel('Swelling, %')
	plt.savefig('Swelling.png')
	plt.show()
	plt.close()

	plt.plot(t*G,R,color='b') #Dependence of radius 
	plt.ylabel('Radius, nm')
	plt.xlabel('Dose, dpa')
	plt.savefig('Radius.png')
	plt.show()
	plt.close()

	T = np.linspace(500,2000,100)
	plt.plot(T,Difi(T)*ci(T,1)*zib(1),color='b',label='Dici*zib')
	plt.plot(T,Difi(T)*ci(T,1),color='g',label='Dici')
	plt.plot(T,Difv(T)*cv(T,1),color='r',label='Dvcv')
#	plt.plot(T,Difv(T)*cveq(T,1),color='g',label='Dvcv')
	plt.plot(T,Difv(T)*(cv(T,1)-cveq(T,1))-Difi(T)*ci(T,1)*zib(1),color='m',label='Dv(cv-cveq)-Dici')
	#plt.text(1e-3,1e-6,'Gt',fontsize=14)
	plt.axis([500,2000,1e-22,1e-17])
	plt.yscale('log')
	plt.legend(loc='best')
	plt.xlabel('Temperature, K')
	plt.ylabel('Point defect flux')
	#plt.title('kd = 5e13')
	#plt.savefig('eq_flux_dn.png')
	#plt.show()
