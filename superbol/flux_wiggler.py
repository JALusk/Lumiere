import random

from superbol import lightcurve

#Make a copy of the existing list of fluxes
def copy_flux_list(fluxes, distance, flux_calculator):
    flux_list_copy = lightcurve.calculate_lightcurve(fluxes, distance, flux_calculator)
    #Should this call lightcurve.calculate_lightcurve or another function?
    #Is there a more useful + accurate name than flux_list_copy? 
    print("Copied flux list: ", flux_list_copy)
    return flux_list_copy

#Wiggle them
def wiggle_fluxes(flux_list_copy):
    for flux in flux_list_copy:
        top = flux.flux + flux.flux_uncertainty
        bottom = flux.flux - flux.flux_uncertainty
        #Where does df come from? lightcurve.calculate_lightcurve doesn't retunr an error bar...
        flux.flux = random.uniform(bottom, top)
    print("Wiggled fluxes: ", flux_list_copy)
    return flux_list_copy

