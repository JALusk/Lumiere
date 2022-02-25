import matplotlib.pyplot as plt

from superbol import extinction, read_osc

"""
TODO Need to create a class that contains the fluxes, magnitudes and etc. that are calculated from 
an osc collection
This class should also contain a function that opens a graph to view the flux and its details
"""
class FluxCurve:
    """constructor for file path,
    constructor for photometry given
    constructor for url given
    """
    photometries = []
    luminosities = []
    times = []

    def __init__(self, photometries):
        self.photometries = photometries
        self.luminosities
        # fluxes = []
        # for photometry_dict in photometries:
        #     try:
        #         observed_magnitude = read_osc.get_observed_magnitude(photometry_dict)
        #         extinction_value = extinction.get_extinction_by_name(extinction_table, observed_magnitude.band)
        #         observed_magnitude.magnitude = extinction.correct_observed_magnitude(observed_magnitude, extinction_value)
        #         fluxes.append(observed_magnitude.convert_to_flux())
        #     except:
        #         pass
        
        # distance = lum.Distance(3.0E7 * 3.086E18, 7.0E6 * 3.086E18)
        # self.luminosities = lightcurve.calculate_lightcurve(fluxes, distance, lqbol.calculate_qbol_flux)


    def draw_curve(self):
        x=[1,3,5,7]
        y=[2,4,6,1]
        print('drawing curve')
        plt.plot(x,y)
        plt.show()
    
    def sum_fn(self):
        print('yeah im working. what ab it')