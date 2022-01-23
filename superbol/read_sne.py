import unittest
import requests

from superbol import lightcurve, lqbol, lum, read_osc

def get_supernova_photometry(name):
    return requests.get(f'https://api.astrocats.space/{name}/photometry/').json()[name]["photometry"]
