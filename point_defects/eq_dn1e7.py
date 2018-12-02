import matplotlib.pyplot as plt
import numpy as np
import math

volume = 1.000000e-18 # in Litres
volume = volume/1000
a = 2.86e-10 # lattice constant in m
b = math.sqrt(3)/2*a # Burger's vector sqrt(3)/2*a
Va = (a*a*a)/2.0 # Atomic volume
na = volume/Va # Number of atoms
pi = 3.14

#Parameters
G = 1e-7
Ks = 1.2e13*0
Kd =1e14
Kiv = 4*3.14*(3.1e-10+2.9e-10)/Va
zi = 1.1
zv = 1.0
R = 10.0
Cp = 1e21*1e0*0
Ki = Ks + zi*Kd + 4*pi*(R*1e-9)*Cp
Kv = Ks + zv*Kd + 4*pi*(R*1e-9)*Cp

def Difi(T):
	return 1.5e12/6.0*1.43**2*1e-20*np.exp(-0.9*11600/T)

def Difv(T):
	E_vac_mig = 1.15	#1.1	#eV
	d_vac = a*math.sqrt(3.)/2
	nu = 1.5e13#*55	#3.6e12
#	return 8./6*d_vac**2*nu*np.exp(-E_vac_mig*11600/T)
	return 1.5e12/6.0*1.43**2*1e-20*np.exp(-E_vac_mig*11600/T)

def cv(T):
	Di = Difi(T)
	Dv = Difv(T)
	return np.sqrt(np.power(Ki*Di/(2*Kiv*(Di+Dv)),2)+G*Ki*Di/(Kiv*(Di+Dv)*Kv*Dv))-Ki*Di/(2*Kiv*(Di+Dv))

def ci(T):
	Di = Difi(T)
	Dv = Difv(T)
	return np.sqrt(np.power(Kv*Dv/(2*Kiv*(Di+Dv)),2)+G*Kv*Dv/(Kiv*(Di+Dv)*Ki*Di))-Kv*Dv/(2*Kiv*(Di+Dv))

def cveq(T):
	Ef = 1.8 #eV
	return np.exp(-Ef*11600/T)


T = np.linspace(500,2500,100)
plt.plot(T,ci(T),color='b',label='ci')
plt.plot(T,cv(T),color='r',label='cv')
plt.plot(T,G/Ki/Difi(T),ls='--',color='b',label='ci*')
plt.plot(T,G/Kv/Difv(T),ls='--',color='r',label='cv*')
plt.plot(T,cveq(T),color='g',label='cv_eq')

#plt.text(1e-3,1e-6,'Gt',fontsize=14)
plt.yscale('log')
#plt.axis([1e-5,1e3,1e-8,1e-2])
plt.axis([500,2000,1e-13,1e-2])
plt.legend(loc='best')
plt.xlabel('Temperature, K')
plt.ylabel('Atomic concentration, 1/at')
#plt.title('kd = 5e13')
plt.savefig('eq_concentration_dn_G1e7.png')
plt.show()
plt.close()

#plt.plot(1000.0/T,Difi(T)*ci(T),color='b',label='Dici')
#plt.plot(1000.0/T,Difv(T)*cv(T),color='r',label='Dvcv')
#plt.plot(1000.0/T,np.sqrt(Difv(T)*G/Kiv*Ki/Kv),ls='--',color='r')
#plt.plot(1000.0/T,G/Kv+T*0,ls='--',color='r')

plt.plot(T,Difi(T)*ci(T),color='b',label='Dici')
plt.plot(T,Difv(T)*cv(T),color='r',label='Dvcv')
plt.plot(T,Difv(T)*cveq(T),color='g',label='Dvcveq')
plt.plot(T,np.sqrt(Difv(T)*G/Kiv*Ki/Kv),ls='--',color='r')
plt.plot(T,G/Kv+T*0,ls='--',color='r')

#plt.plot(T,Difv(T)*cveq(T),color='g',label='Dvcv')
#plt.plot(1000.0/T,Cp*R*(Difv(T)*(cv(T)-0*cveq(T))-Difi(T)*ci(T))/1e23,color='m',label='Dv(cv-cveq)-Dici')
#plt.text(1e-3,1e-6,'Gt',fontsize=14)
plt.axis([500,2000,1e-24,1e-20])
plt.yscale('log')
plt.legend(loc='best')
plt.xlabel('Temperature, K')
plt.ylabel('Point defect flux, m2/s')
#plt.title('kd = 5e13')
plt.savefig('eq_flux_dn_G1e7.png')
plt.show()
