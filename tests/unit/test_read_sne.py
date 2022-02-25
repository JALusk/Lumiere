import json
import os
import unittest
import numpy as np
from .context import superbol
import requests
from superbol import read_sne
from httmock import urlmatch, all_requests, HTTMock, response

SUPERNOVA_NAME = 'SN2000cb'
DIRNAME = os.path.dirname(__file__)
SUPERNOVA_FILE_PATH = os.path.join(
    DIRNAME, f'../../data/{SUPERNOVA_NAME}.json')


@all_requests
def google_mock(url, request):
	content = json.load(open(SUPERNOVA_FILE_PATH))
	return response(200, content, None, None, 5, request)

class TestReadSne(unittest.TestCase):
    def test_get_photometry(self):
        with open(SUPERNOVA_FILE_PATH) as f, HTTMock(google_mock):
            req_array = np.array(read_sne.get_supernova_photometry(SUPERNOVA_NAME))
            file_array = np.array(json.load(f)[SUPERNOVA_NAME]['photometry'])
            assert np.array_equal(req_array, file_array)

    def test_get_photometry(self):
        with open(SUPERNOVA_FILE_PATH) as f, HTTMock(google_mock):
            req_array = np.array(read_sne.get_supernova_lum_dist(SUPERNOVA_NAME))
            file_array = np.array(json.load(f)[SUPERNOVA_NAME]['lumdist'])
            assert np.array_equal(req_array, file_array)
