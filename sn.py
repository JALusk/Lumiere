import tables as tb
hdf5_filename = './hdf5/sn_data.h5'

class SN(object):
    """A supernova is the explosion that ends the life of a star

    The SN needs to be conatained within the HDF5 database before it is used
    by SNoBoL. Once there, simply create a supernova by calling the constructor
    with the name of the SN as a string of the form "sn[YEAR][Letter(s)]"

    For example:
    sn1987A = SN('sn1987a')
    sn1999em = SN('sn1999em')

    Attributes
    ----------
    name : Name of the supernova, "sn" followed by the year of first observation
           along with a letter designating the order of observation in that
           year. "sn1987a" was the first SN observed in 1987. "sn2000cb" was the
           eightieth SN observed in 2000.
    """

    def __init__(self, name):
        """Initializes the SN with supplied value for [name]"""
        self.name = name
        
        # Load SN data from HDF5 file
        self.read_hdf5()

    def read_hdf5(self):
        h5file = tb.open_file(hdf5_filename, 'r')
       
        # Get the desired node from the HDF5 file
        sn_node = h5file.get_node('/sn', self.name)
