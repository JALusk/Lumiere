import matplotlib.pyplot as plt
import math
import numpy as np

#Get data from file
with open('luminosity_2018hna.txt', 'r') as f:
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
decline_luminosity = luminosity[35:]
decline_time = time_s[35:]
fit = np.polyfit(decline_time, decline_luminosity, 1)
line = np.poly1d(fit)
print("Linear fit to post-plateau tail: ", line)

#Calculate nickel mass
g1 = 1.32*10**-6
g2 = 1.02*10**-7
s = ((3.90*10**10)*(e**(g2*t)) + (6.78*10**9)(e**(-g2*t) - e**(g1*t)))

#Plot lightcurve
plt.scatter(time_s, log_luminosity)
plt.plot(decline_time, line(decline_time))
plt.title("Time since explosion (seconds) vs. log base 10 of luminosity")
plt.xlabel('Time since explosion (seconds)')
plt.ylabel('Log(luminosity)')
plt.savefig('../Desktop/2018hna_lightcurve.png')
