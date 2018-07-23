import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate

data = np.genfromtxt('snls_3rdyear_lcparams.txt', names=True)



total_z = (1+data['zcmb'])*(1+data['zhel'])-1

a = 1/(1+total_z)

# Cosmological paramters (these we could sample??)
Omega_m0 = 0.3089
H_0 = 67.74 #km/s/Mpc
c = 299792 #km/s

# Calculating the DL from the mu, We should probably use a more precise
# calculation of MU. meaning a more accurate Absolute magnitude

mu = data ['mb'] + 19.5
DL_mu = 10**((mu-25)/5)

D_H = 4514.94


# This is their method, assuming everything they did is correct
Y = a* DL_mu/D_H



# Now to calculate the error
total_z_err = total_z * data['dz']*(1/(1+data['zcmb'])**2 + 1/(1+data['zhel'])**2)**(1/2)

a_err = 1/(1+total_z)**2 *total_z_err

DL_mu_err = 2**(mu/5-5) * 5**(mu/5 -6) * np.log(10) * data['dmb']

Y_err = Y * ((a_err/a)**2 + (DL_mu_err/DL_mu)**2)**(1/2)

LT_RM_err = np.zeros(Y.shape[0])
sum_err2 = 0

LT_RM = np.zeros(Y.shape[0])
sum = 0 
i = 1
while i < Y.shape[0]:
	sum = sum + a[i]*(Y[i] -Y[i-1])
	sum_err2 = sum_err2 + (a[i]*(Y[i]-Y[i-1]))**2 * ((a_err[i]/a[i])**2 + (Y_err[i]**2 + Y_err[i-1]**2)/(Y[i]-Y[i-1])**2)
	LT_RM[i] = 1-sum
	LT_RM_err[i] = sum_err2 **(1/2)
	i = i+1
LT_RM[0] = 1 - a[0]*Y[0]


# Now if we use a classical version of the luminosity distance see fore
# example 1601.01451

def E(x):
	return (Omega_m0*(1+x)**3 + (1-Omega_m0))**(-1/2)

integral = np.zeros(data['zcmb'].shape[0])
i = 0
while i< data['zcmb'].shape[0]:
	integra = integrate.quad(E, 0, data['zcmb'][i])
	integral[i] = integra[0]
	i= i+1

DL_z = (1+data['zhel'])*c/H_0 * integral

Y_z = a* DL_z/D_H

LT_RMz = np.zeros(Y.shape[0])
sum = 0 
i = 1
while i < Y_z.shape[0]:
	sum = sum + a[i]*(Y_z[i] -Y_z[i-1])
	LT_RMz[i] = 1-sum
	i = i+1
LT_RMz[0] = 1

fig1 = plt.figure()
#plt.plot(LT_RM, a, 'r.', label='Their method')
plt.errorbar(LT_RM, a, xerr = LT_RM_err, fmt= '-r', label='Their method')
plt.plot(LT_RMz, a, 'b.', label='Their method, dl calculated from z')
plt.title('Their lookback time method')
plt.legend()
plt.xlabel('t')
plt.ylabel('a(t)')
fig1.savefig('theirmethod.png')
plt.show()
# Now lets calculate the classical lookback time!

def E2(x):
	return 1/((1+x)*(Omega_m0*(1+x)**3 + (1-Omega_m0))**(1/2))

integral = np.zeros(data['zcmb'].shape[0])
i = 0
while i< data['zcmb'].shape[0]:
	integrab = integrate.quad(E2, 0, data['zcmb'][i])
	integral[i] = integrab[0]
	i = i +1

t_l = 1/H_0 * integral /(1+data['zhel']) #Not sure if the last part here is needed. 

# Now if we look at their derivation, we should note, that they are
# not actually calculating the lookback time. so to get the lookback time
# we can compare to the classical value we do the following:

t_lRM = 1/(H_0) *(1-LT_RM)
t_lRMz = 1/(H_0) *(1-LT_RMz)

fig2 = plt.figure()
plt.plot(t_l, a, 'g.', label = 'traditional lookback time')
plt.plot(t_lRM, a, 'r.', label = 'RM lookback time corrected')
plt.plot(t_lRMz, a, 'b.', label = 'RM lookback time corrected, dl from z')

plt.title('alternative lookback time methods')
plt.legend()
plt.xlabel('t')
plt.ylabel('a(t)')
fig2.savefig('correctedlookbacktime.png')


#What if we are calculating t_l wrong?
integral = np.zeros(total_z.shape[0])
i = 0
while i<total_z.shape[0]:
	integrac = integrate.quad(E2, 0 ,total_z[i])
	integral[i] = integrac[0]
	i = i+1
	
t_l2 = 1/H_0 * integral

fig3 = plt.figure()
plt.plot(t_l, a, 'g.', label = 'lookback time calculated with zcmb in integral')
plt.plot(t_l2, a, 'k.', label = 'lookback time calculate with total z in integral')
plt.plot(t_lRMz, a, 'b.', label = 'Their lookback time, corrected, and with trad. Dl')

plt.title('which Z to integrate')
plt.legend()
plt.xlabel('t')
plt.ylabel('a(t)')
fig3.savefig('whichzvalue.png')
#plt.show()
