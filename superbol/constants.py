"""Constants for use by the bolometric correction routine
"""

# Coefficients for polynomial fit to bolometric correction - color relation
coeff_BminusV = [-0.823, 5.027, -13.409, 20.133, -18.096, 9.084, -1.950]
coeff_VminusI = [-1.355, 6.262, -2.676, -22.973, 35.524, -15.340]
coeff_BminusI = [-1.096, 3.038, -2.246, -0.497, 0.7078, 0.576, -0.713,
                 0.239, -0.027]

# Ranges of validity for polynomial fits
min_BminusV = -0.2
max_BminusV = 1.65
min_VminusI = -0.1
max_VminusI = 1.0
min_BminusI = -0.4
max_BminusI = 3.0

# RMS errors in polynomial fits
rms_err_BminusV = 0.113
rms_err_VminusI = 0.109
rms_err_BminusI = 0.091

# Zeropoint for use in the calculation of bolometric magnitude
mbol_zeropoint = 11.64
