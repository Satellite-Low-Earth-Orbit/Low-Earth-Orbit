import scipy.integrate as integrate
import numpy as np
import matplotlib.pyplot as plt

#-------------- Legend ---------------

sigma_xs = np.arange(1,10,0.5)
sigma_zs = np.arange(1,10,0.5)
sigma_x = 3.
sigma_z = 3.
phi = np.pi/2.
#x_e = 0.584
x_es = np.linspace(0, 21.5, 500)
#r_a = 0.01
r_as = np.linspace(0.01, 0.2, 6)

prefactor = 0.

earth_r = 6371
orbit_height = 500
orbit_r = (earth_r + orbit_height)*10**3
G = 6.67408 * 10**-11  # Gravitational constant
M = 5.972 * 10**24     # Mass of the Earth
T = np.sqrt(4*np.pi**2 * orbit_r**3 / (G*M)) # orbit period
N = 1000 # number of satellites in the cluster
S = 2*np.pi*orbit_r / N  # spacing between the satellites in the cluster

#-------------------------------------

print(T)

def P(r,theta):
  
  return r * np.exp(-0.5*r**2 * ((np.sin(theta) / sigma_x)**2 + (np.cos(theta) / sigma_z)**2 ) + r*x_e * ( (np.sin(phi) * np.sin(theta)/ sigma_x**2) + (np.cos(phi) * np.cos(theta) / sigma_z**2) ) )
  
def prob_collision():
  
  # -------- Monte Carlo sampling -----------
  
  n_samples = 100000
  x_e_samples = np.random.uniform(n_samples) * S/2.0
  temp = 0.
  
  # ---------- Monte Carlo Averaging ---------
  
  for i in range(n_samples):
  
    x_e = x_e_samples[i]
    prefactor = np.exp(-0.5 * x_e**2 * ( (np.sin(phi) / sigma_x)**2 + (np.cos(phi) / sigma_z)**2 ) )/(2*np.pi * sigma_x * sigma_z)
    temp += prefactor * integrate.dblquad(P, r_lower, ra, theta_lower, theta_upper)[0]
  
  return temp/n_samples

def p_survival_per_day(prob_collision):

  n_encounters_per_day = 2./(T / (60. * 60. * 24.) )
  
  p_survival_per_day = (1 - prob_collision) ** n_encounters_per_day

  return p_survival_per_day
  
  
def p_survival_per_year(p_survival_per_day):

  return p_survival_per_day ** 365.

r_lower = 0
#r_upper = r_a
theta_lower = 0
theta_upper = 2*np.pi

results = np.zeros((len(r_as), len(x_es)))

for j, ra in enumerate(r_as):
  for i, elem in enumerate(x_es):
    x_e = elem
    prefactor = np.exp(-0.5 * x_e**2 * ( (np.sin(phi) / sigma_x)**2 + (np.cos(phi) / sigma_z)**2 ) )/(2*np.pi * sigma_x * sigma_z)
    prob = integrate.dblquad(P, r_lower, ra, theta_lower, theta_upper)[0]
    results[j,i] = prefactor * prob
  
for j in range(len(r_as)):
  plt.plot(x_es, results[j], label="$r_a$ = "+str(np.round(r_as[j],2)))
plt.xlabel("$x_e$")
plt.ylabel("P")
plt.legend()
plt.show()

"""SX, SZ = np.meshgrid(sigma_xs, sigma_zs)

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
plt.show()

#result = integrate.dblquad(P, r_lower, r_upper, theta_lower, theta_upper)
#print("Result : ", prefactor * result[0])
#print(" Numerical Error : ", result[1])

"""



