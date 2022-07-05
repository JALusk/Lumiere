import matplotlib.pyplot as plt
import math
import numpy as np
import scipy.optimize

#Get data from file
with open('luminosity_2018hna_nouv.txt', 'r') as f:
    lines = f.readlines()
    luminosity = []
    time = []
    for i in lines:
        time.append(float(i.split(',')[1]))
        luminosity.append(float(i.split(',')[3]))
    f.close()

#Find log of the luminosity
log_luminosity = [math.log10(x) for x in luminosity]

#Change days to seconds
time_s = [x*86400 for x in time]

#Linear fit to post-plateau decline
decline_luminosity = []
decline_time_s = []
decline_indices = [ind for ind, x in enumerate(luminosity) if ind > 35 if x < (1*10**41.2)]
for i in decline_indices:
    decline_luminosity.append(luminosity[i])
    decline_time_s.append(time_s[i])

decline_luminosity = np.asarray(decline_luminosity)
decline_time_s = np.asarray(decline_time_s)

#print("Decline lum: ", decline_luminosity)
#print("Decline time: ", decline_time_s)
#print("Luminosity: ", luminosity)
#print("Time in sec: ", time_s)

#Calculate nickel mass
g1 = -1.32*10**-6
g2 = -1.02*10**-7
ni_mass_sol = 0.1
sol_to_g = 1.989*10**33

def co_decay_erg_per_s(t, ni_mass_sol):
       s = ((3.90*10**10)*(math.e**(g1*t)) + (6.78*10**9)*(math.e**(g2*t) - math.e**(g1*t)))
       ni_lum_erg_per_s = s*ni_mass_sol*sol_to_g
       return ni_lum_erg_per_s

popt, pocov = scipy.optimize.curve_fit(co_decay_erg_per_s, decline_time_s, decline_luminosity, p0=ni_mass_sol*sol_to_g)

ni_mass_sol = popt[0]
print("Nickel mass in solar masses: ", ni_mass_sol)
print("Nickel mass in grams: ", ni_mass_sol*sol_to_g)

ni_lum_erg_per_s = co_decay_erg_per_s(decline_time_s, ni_mass_sol)
log_lum_fit = []
for i in ni_lum_erg_per_s:
    log_lum_fit.append(math.log10(i))

#Plot lightcurve
plt.scatter(time_s, log_luminosity)
plt.plot(decline_time_s, log_lum_fit, label = 'Linear fit for nickel-56 mass')
plt.title("Time since explosion (seconds) vs. log base 10 of luminosity")
plt.xlabel('Time since explosion (seconds)')
plt.ylabel('Log(luminosity)')
plt.savefig('../Desktop/2018hna_lightcurve_nouv.png')
plt.show()
