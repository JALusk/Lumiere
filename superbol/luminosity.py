import math

from superbol.bc_polynomial import calc_bolometric_correction as bc
from superbol.bc_polynomial import quadrature_sum
from superbol.constants import mbol_zeropoint


def calc_Fbol(color_value, color_err, color_type, v_magnitude,
              v_magnitude_err):
    """Calculates the bolometric flux of a Type II-P supernova.

    Args:
        color_value (float): B-V, V-I, or B-I color of the supernova in
            magnitudes (corrected for reddening and extinction from
            the host and MWG.)
        color_err (float): Uncertainty in the photometric color.
        color_type (str): String signifying which color color_value
            represents. Valid values are "BminusV" for B-V, "VminusI"
            for V-I, and "BminusI" for B-I.
        v_magnitude (float): Photometric magnitude in the V band, corrected for
            host + MWG extinction.
        v_magnitude_err (float): Uncertainty in the V band magnitude after
            correction for host + MWG extinction.

    Returns:
        tuple: A tuple containing the bolometric flux used when calculating the
        bolometric luminosity, and the uncertainty in that number.

        (Fbol, uncertainty)

        (-999, -999) if the bolometric correction calculated from the
        color_value and color_type is -999 (which means the observed
        color is outside the range of validity of the polynomial fit.)
    """
    bolometric_correction, bc_err = bc(color_value, color_err, color_type)

    if bolometric_correction == -999:
        Fbol = -999
        Fbol_uncertainty = -999
    else:
        Fbol = 10**(-0.4 *
                    (bolometric_correction + v_magnitude + mbol_zeropoint))
        Fbol_uncertainty = (0.4 * math.log(10) * Fbol *
                            math.sqrt(bc_err**2 + v_magnitude_err**2))

    return Fbol, Fbol_uncertainty


def calc_4piDsquared(distance, distance_err):
    """Calculates :math:`4\\pi D^2`, to convert flux to luminosity.

    Args:
        distance (float): The distance to the supernova in centimeters.
        distance_err (float): The uncertainty in the distance to the supernova.

    Returns:
        tuple: A tuple containing the :math:`4 \\pi D^2`, and the uncertainty of this number.

        (4piDsquared, uncertainty)
    """
    fourPiDsquared = 4.0 * math.pi * distance**2.0
    fourPiDsquared_uncertainty = 8.0 * math.pi * distance * distance_err

    return fourPiDsquared, fourPiDsquared_uncertainty


def calc_Lbol(color_value, color_err, color_type, v_magnitude, v_magnitude_err,
              distance, distance_err):
    """Calculates the bolometric luminosity of a Type II-P Supernova.

    Args:
        color_value (float): B-V, V-I, or B-I color of the supernova in
            magnitudes (corrected for reddening and extinction from
            the host and MWG.)
        color_err (float): Uncertainty in the photometric color.
        color_type (str): String signifying which color color_value
            represents. Valid values are "BminusV" for B-V, "VminusI"
            for V-I, and "BminusI" for B-I.
        v_magnitude (float): Photometric magnitude in the V band, corrected for
            host + MWG extinction.
        v_magnitude_err (float): Uncertainty in the V band magnitude after
            correction for host + MWG extinction.
        distance (float): The distance to the supernova in centimeters.
        distance_err (float): The uncertainty in the distance to the supernova.

    Returns:
        tuple: A tuple containing the bolometric luminosity in ergs per second,
        and the uncertainty in that value.

        (Lbol, uncertainty)

        (-999, -999) if the bolometric correction is -999 (which means
        the observed color value is outside the range of vaidity of the
        polynomial fit used to determine the bolometric correction.)
    """
    Fbol, Fbol_err = calc_Fbol(color_value, color_err, color_type, v_magnitude,
                               v_magnitude_err)
    fourPiDsquared, fourPiDsquared_err = calc_4piDsquared(distance,
                                                          distance_err)

    if Fbol == -999:
        Lbol = -999
        Lbol_uncertainty = -999
    else:
        Lbol = Fbol * fourPiDsquared
        Lbol_uncertainty = quadrature_sum(fourPiDsquared * Fbol_err,
                                          Fbol * fourPiDsquared_err)

    return Lbol, Lbol_uncertainty
