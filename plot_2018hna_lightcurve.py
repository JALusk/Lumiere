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
print("Time: ", time)
print("Luminosity: ", luminosity)
#Find log of the luminosity
log_luminosity = [math.log10(x) for x in luminosity]
print("Log of luminosity: ", log_luminosity)

#Linear fit to post-plateau decline
decline_luminosity = log_luminosity[35:]
decline_time = time[35:]
fit = np.polyfit(decline_time, decline_luminosity, 1)
line = np.poly1d(fit)
print("Line: ", line)
print("Decline time: ", decline_time)
print("Fit: ", fit)
print("Decline luminosity: ", decline_luminosity)

#Plot lightcurve
plt.scatter(time, log_luminosity)
plt.plot(decline_time, line(decline_time))
plt.title("Time since explosion (days) vs. log base 10 of luminosity")
plt.xlabel('Time since explosion (days)')
plt.ylabel('Log(luminosity)')
plt.savefig('../Desktop/2018hna_lightcurve.png')
