import requests

def get_supernova_lum_dist(name):
    return requests.get(f'https://api.astrocats.space/{name}/lumdist').json()[name]["lumdist"]

def get_supernova_photometry(name):
    return requests.get(f'https://api.astrocats.space/{name}/photometry/').json()[name]["photometry"]
