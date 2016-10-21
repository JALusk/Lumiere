import unittest
import os
from mock import Mock
from .context import superbol
import superbol.sn as sn
import tables as tb
from pkg_resources import resource_filename

class TestSNInitialization(unittest.TestCase):
   
    def setUp(self):
        self.sn_name = "SN1987A"
        self.my_sn = sn.SN(self.sn_name)

    def test_sn_name_initialization(self):
        expected = self.sn_name
        result = self.my_sn.name
        self.assertEqual(expected, result)

    def test_sn_source_initialization_when_no_source_given(self):
        expected = None
        result = self.my_sn.source
        self.assertEqual(expected, result)

    def test_sn_source_initialization_when_source_given(self):
        source = '/path/to/source'
        my_sn = sn.SN(self.sn_name, source)
        expected = source
        result = my_sn.source
        self.assertEqual(expected, result)

class TestSNDefaultSourceOpening(unittest.TestCase):

    def setUp(self):
        self.sn_name = "sn1998a"
        self.default_source = None
        self.my_sn = sn.SN(self.sn_name, self.default_source)

    def test_sn_open_default_source_file_returns_pytables_file_object(self):
        h5file = self.my_sn.open_source_h5file()
        result = isinstance(h5file, tb.file.File)
        h5file.close()
        self.assertTrue(result)

    def test_sn_open_default_source_file_returns_default_sn_data_file(self):
        h5file = self.my_sn.open_source_h5file()
        expected = resource_filename('superbol', 'data/sn_data.h5')
        result = h5file.filename
        h5file.close()
        self.assertEqual(expected, result)

class TestSNCustomSourceOpening(unittest.TestCase):

    def setUp(self):
        self.sn_name = "sn1998a"
        self.custom_source = "tests/test_data.h5"
        self.my_sn = sn.SN(self.sn_name, self.custom_source)

    def test_sn_open_custom_source_file_returns_pytables_file_object(self):
        h5file = self.my_sn.open_source_h5file()
        result = isinstance(h5file, tb.file.File)
        h5file.close()
        self.assertTrue(result)

    def test_sn_open_custom_source_file_returns_correct_filename(self):
        h5file = self.my_sn.open_source_h5file()
        expected = self.custom_source
        result = h5file.filename
        h5file.close()
        self.assertEqual(expected, result)
