import math

from itertools import groupby

def group_magnitudes(magnitudes, keyfunc=math.floor):
    """Group magnitudes by applying keyfunc() to each magnitude.time.

    Magnitudes with the same return value for keyfunc(magnitude.time) will be 
    part of the same group."""
    grouped_magnitudes = [list(it) for k, it in groupby(magnitudes, lambda x: keyfunc(x.time))]
    return grouped_magnitudes
