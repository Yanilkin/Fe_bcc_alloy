import matplotlib.pyplot as plt
import numpy as np

a = 3.14e-10 # lattice constant in m
Va = (a*a*a)/2.0 # Atomic volume
pi = 3.14

G = 5e-4
Kd = 5e13 # Dislocation network sink
Kiv = 4*3.14*(3.1e-10+2.9e-10)/Va # Recombination rate constant
zid = 3.0 # Bias factor
zvd = 1.0 # Bias factor
Cp = 1e21
R =120*1e-9


def Difv(T):
	E_vac_mig = 1.3	#1.1	#eV
	d_vac = a*np.sqrt(3.)/2
	nu = 1.5e13#*55	#3.6e12
	return 8./6*d_vac**2*nu*np.exp(-E_vac_mig*11600/T)

T = np.linspace(500,1400,100)

plt.plot(1000.0/T,100*np.sqrt(Difv(T)/Kiv/G)*Kd*(zid-zvd)*1.0, ls='--', color='b')
plt.plot(1000.0/T,100*np.sqrt(Difv(T)/Kiv/G)*Kd*(zid-zvd)*10.0, ls='--', color='r')
plt.plot(1000.0/T,100*np.sqrt(Difv(T)/Kiv/G)*Kd*(zid-zvd)*100.0, ls='--', color='g')

Qi = zid*Kd/(4*pi*Cp*R)
Qv = zvd*Kd/(4*pi*Cp*R)
print Qi

plt.plot(1000.0/T,100*(zid-zvd)/zid*Qi/(1+Qi)/(1+Qv)*100.0+T*0, ls='--', color='g')


dpa1 = np.loadtxt('Swelling_1dpa.dat',skiprows=1)
dpa10 = np.loadtxt('Swelling_10dpa.dat',skiprows=1)
dpa100 = np.loadtxt('Swelling_100dpa.dat',skiprows=1)

plt.plot(1000.0/dpa1[:,0],dpa1[:,1],marker='o', color='b',label='1dpa')
plt.plot(1000.0/dpa10[:,0],dpa10[:,1],marker='o', color='r',label='10dpa')
plt.plot(1000.0/dpa100[:,0],dpa100[:,1],marker='o', color='g',label='100dpa')
plt.xlabel('1/Temperature, 1000/K')
plt.ylabel('Swelling, %')
plt.yscale('log')
plt.legend(loc='best')
plt.savefig('Swelling_T.png')
plt.show()