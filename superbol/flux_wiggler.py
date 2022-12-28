import random

from superbol import lightcurve

#Make a copy of the existing list of fluxes
def copy_flux_list(sed):
    flux_list_copy = [None] * len(sed)
    for i in range(len(sed)):
        flux_list_copy[i] = sed[i]
    for flux in flux_list_copy:
        print("Copied flux list: ", flux)
    return flux_list_copy

#Wiggle them
def wiggle_fluxes(flux_list_copy):
    for flux in flux_list_copy:
        top = flux.flux + flux.flux_uncertainty
        print("Top: ", top)
        bottom = flux.flux - flux.flux_uncertainty
        print("Bottom; ", bottom)
        #Where does df come from? lightcurve.calculate_lightcurve doesn't retunr an error bar...
        flux.flux = random.uniform(bottom, top)
    #print("Wiggled fluxes: ", flux_list_copy)
    return flux_list_copy

