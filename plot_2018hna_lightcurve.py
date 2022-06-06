import matplotlib.pyplot as plt
import math

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

#Plot lightcurve
plt.scatter(time, log_luminosity)
plt.title("Time since explosion (days) vs. log base 10 of luminosity")
plt.xlabel('Time since explosion (days)')
plt.ylabel('Log(luminosity)')
plt.savefig('../Desktop/2018hna_lightcurve.png')
