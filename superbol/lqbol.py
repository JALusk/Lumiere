class SpectralEnergyDistribution(object):

    def sort_by_wavelength(self, fluxes):
        return sorted(fluxes, key=lambda flux: flux.wavelength)
