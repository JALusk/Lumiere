import unittest
from .context import superbol
from superbol import mag2flux

class TestObservation(unittest.TestCase):

    def test_observation(self):
        my_obs = mag2flux.Observation()
