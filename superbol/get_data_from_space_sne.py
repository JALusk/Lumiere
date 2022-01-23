import requests
# from superbol.FlxuCurve import FluxCurve

def get_supernova_photometry(name):
    return requests.get(f'https://api.astrocats.space/{name}/photometry/').json()[name]["photometry"]

SUPERNOVA_NAME = 'SN2014J'
supernova_photometry = get_supernova_photometry(SUPERNOVA_NAME)
print(supernova_photometry)
# fluxCurve = FluxCurve(photometries=supernova_photometry)
# fluxCurve.draw_curve()
# # print(photometry)
