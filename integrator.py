import scipy.integrate as integrate
import numpy as np

#-------------- Legend ---------------

sigma_x = 1.
sigma_z = 1.
phi = 1.
x_e = 1.
r_a = 1.

prefactor = np.exp(-0.5 * x_e**2 * ( (np.sin(phi) / sigma_x)**2 + (np.cos(phi) / sigma_z)**2 ) )/(2*np.pi * sigma_x * sigma_z)

#-------------------------------------

def P(r,theta):
  
  return r * np.exp(-0.5*r**2 * ((np.sin(theta) / sigma_x)**2 + (np.cos(theta) / sigma_z)**2 ) + r*x_e * ( (np.sin(phi) * np.sin(theta)/ sigma_x**2) + (np.cos(phi) * np.cos(theta) / sigma_z**2) ) )
  
r_lower = 0
r_upper = r_a
theta_lower = 0
theta_upper = 2*np.pi
result = integrate.dblquad(P, r_lower, r_upper, theta_lower, theta_upper)

print(prefactor * result[0])
print(result[1])

  



