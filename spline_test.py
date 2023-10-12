import math
import numpy as np
from scipy import interpolate

from superbol import fqbol
from superbol import mag2flux

#low
flux01 = mag2flux.MonochromaticFlux(98, 2, 1, 0)
flux02 = mag2flux.MonochromaticFlux(198, 2, 2, 0)
flux03 = mag2flux.MonochromaticFlux(148, 2, 3, 0)
#high
flux04 = mag2flux.MonochromaticFlux(102, 2, 1, 0)
flux05 = mag2flux.MonochromaticFlux(202, 2, 2, 0)
flux06 = mag2flux.MonochromaticFlux(152, 2, 3, 0)
sed_low = [flux01, flux02, flux03]
sed_high = [flux04, flux05, flux06]

low_result = fqbol.SplineIntegralCalculator().calculate(sed_low)
high_result = fqbol.SplineIntegralCalculator().calculate(sed_high)

print("Low bound to spline integrated flux: ", low_result)
print("High bound to spline integrated flux: ", high_result)

