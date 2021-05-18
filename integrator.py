import scipy.integrate as integrate
import numpy as np
import matplotlib.pyplot as plt

#-------------- Legend ---------------

sigma_xs = np.arange(1,10,0.5)
sigma_zs = np.arange(1,10,0.5)
sigma_x = 0.
sigma_z = 0.
phi = np.pi/2.
x_e = 0.584
r_a = 0.003

prefactor = 0.

#-------------------------------------

def P(r,theta):
  
  return r * np.exp(-0.5*r**2 * ((np.sin(theta) / sigma_x)**2 + (np.cos(theta) / sigma_z)**2 ) + r*x_e * ( (np.sin(phi) * np.sin(theta)/ sigma_x**2) + (np.cos(phi) * np.cos(theta) / sigma_z**2) ) )
  
r_lower = 0
r_upper = r_a
theta_lower = 0
theta_upper = 2*np.pi

SX, SZ = np.meshgrid(sigma_xs, sigma_zs)

results = np.zeros_like(SX, dtype=float)

for i in range(len(sigma_xs)):
  for j in range(len(sigma_zs)):
    sigma_x = sigma_xs[i]
    sigma_z = sigma_zs[j]
    prefactor = np.exp(-0.5 * x_e**2 * ( (np.sin(phi) / sigma_x)**2 + (np.cos(phi) / sigma_z)**2 ) )/(2*np.pi * sigma_x * sigma_z)
    results[i][j] = prefactor * integrate.dblquad(P, r_lower, r_upper, theta_lower, theta_upper)[0]
    

print(results)

plt.pcolormesh(SX, SZ, results)
plt.colorbar()
plt.ylabel("$\\sigma_{z'} [km]$")
plt.xlabel("$\\sigma_{x'} [km]$")
plt.title("Probability of collision")
#plt.yticks(sigma_zs, sigma_zs)
#plt.xticks(sigma_xs, sigma_xs)
#plt.ylim([sigma_zs[0]-1, sigma_zs[-1]])
#plt.xlim([sigma_xs[0]-1, sigma_xs[-1]])
plt.show()

#result = integrate.dblquad(P, r_lower, r_upper, theta_lower, theta_upper)
#print("Result : ", prefactor * result[0])
#print(" Numerical Error : ", result[1])

  



