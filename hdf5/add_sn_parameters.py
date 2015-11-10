from tables import *

h5file = open_file('sn_data.h5', 'a')

sn_name = 'sn1998a'

class Parameters(IsDescription):
    heliocentric_z = Float64Col()
    heliocentric_z_err = Float64Col()
    heliocentric_z_ref = StringCol(19)
    distance_Mpc = Float64Col()
    distance_Mpc_err = Float64Col()
    distance_Mpc_ref = StringCol(19)
    Av_gal = Float64Col()
    Av_gal_ref = StringCol(19)
    Av_host = Float64Col()
    Av_host_ref = StringCol(19)
    explosion_JD = Float64Col()
    explosion_JD_err = Float64Col()
    explosion_JD_ref = StringCol(19)

table = h5file.create_table('/sn/' + sn_name, 'parameters', Parameters)

parameters = table.row

parameters['heliocentric_z'] = 
parameters['heliocentric_z_err'] =
parameters['heliocentric_z_ref'] =
parameters['distance_Mpc'] =
parameters['distance_Mpc_err'] =
parameters['distance_Mpc_ref'] =
parameters['Av_gal'] =
parameters['Av_gal_ref'] =
parameters['Av_host'] =
parameters['Av_host_ref'] =
parameters['explosion_JD'] =
parameters['explosion_JD_err'] =
parameters['explosion_JD_ref'] =

parameters.append()

table.flush()
