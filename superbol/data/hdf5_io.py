import tables as tb

h5file = tb.open_file('sn_data.h5', 'a')

def get_phot_table_description():
    # The data fields for this will be the same as those for the other SNe
    phot_description = h5file.root.sn.sn1998a.phot.description
    return phot_description

def get_parameters_table_description():
    # The data fields for this will be the same as those for the other SNe
    parameters_description = h5file.root.sn.sn1998a.parameters.description
    return parameters_description

def make_new_sn_group(sn_name):
    # Make a new group to hold the SN data, if it doesn't already exist
    if h5file.__contains__("/sn/"+sn_name):
        print("Found existing HDF5 group /sn/"+sn_name)
        print("Skipping group creation")
    else:
        print("Creating HDF5 group /sn/"+sn_name)
        group = h5file.create_group(h5file.root.sn, sn_name)

def make_new_sn_phot_table(sn_name):
    phot_description = get_phot_table_description()
    # Make a table to hold the photometry, if it doesn't aleady exist
    if h5file.__contains__("/sn/"+sn_name+"/phot"):
        print("Found existing HDF5 table /sn/"+sn_name + "/phot")
        print("Skipping table creation")
        phot_table = h5file.get_node("/sn/"+sn_name+"/phot")
        return phot_table
    else:
        print("Creating HDF5 table /sn/"+sn_name+"/phot")
        phot_table = h5file.create_table("/sn/"+sn_name, "phot", phot_description)
        return phot_table

def make_new_sn_parameters_table(sn_name):
    parameters_description = get_parameters_table_description()
    # Make a table to hold the photometry, if it doesn't aleady exist
    if h5file.__contains__("/sn/"+sn_name+"/parameters"):
        print("Found existing HDF5 table /sn/"+sn_name + "/parameters")
        print("Skipping table creation")
        parameters_table = h5file.get_node("/sn/"+sn_name+"/parameters")
        return parameters_table
    else:
        print("Creating HDF5 table /sn/"+sn_name+"/parameters")
        parameters_table = h5file.create_table("/sn/"+sn_name, "parameters", parameters_description)
        return parameters_table

def make_new_filter_entry(filter_name, filter_eff_wl, filter_flux_zeropoint,
                          note, ref):
    # If the user provided eff_wl and flux_zeropoint, make a new filter entry
    new_filter = filter_table.row
    
    largest_id = set_new_filter_id()
    new_filter['filter_id'] = largest_id
    new_filter['name'] = filter_name
    new_filter['eff_wl'] = filter_eff_wl
    new_filter['flux_zeropoint'] = filter_flux_zeropoint
    new_filter['note'] = note
    new_filter['ref'] = ref
    
    return new_filter, largest_id

def append_new_filter(new_filter):
    new_filter.append()

def write_filter_table():
    filter_table.flush()
    
def set_new_filter_id():
    # Give new filter a unique ID larger than the largest ID in the HDF5 file
    # NOTE: JLusk <jeremy.lusk@ou.edu> This will cause problems if two
    # users try to merge HDF5 changes back in to the github repo. Find a
    # better way to do this in the future!
    largest_id = max([x['filter_id'] for x in filter_table.iterrows()]) + 1
    print("Adding new filter with id", largest_id)
    return largest_id

def get_old_filter_id(filter_name):
    # If no filter is specified, use the parameters already stored in SuperBoL
    filter_table = h5file.root.filters
    filter_id = min([x['filter_id'] for x in filter_table.iterrows() if x['name'].decode('ascii') == filter_name])
    return filter_id

def make_new_observation_entry(filter_id, jd, magnitude, uncertainty, reference, note, phot_table):
    new_observation = phot_table.row

    new_observation['filter_id'] = filter_id
    new_observation['jd'] = jd
    new_observation['magnitude'] = magnitude
    new_observation['uncertainty'] = uncertainty
    new_observation['reference'] = reference
    new_observation['note'] = note
    new_observation['telescope_id'] = 0

    return new_observation

def append_new_observation(new_observation):
    new_observation.append()

def write_photometry_table():
    phot_table.flush()
    
def add_sn_filter_photometry(sn_name, filter_name, filter_eff_wl,
                             filter_flux_zeropoint, ref, note, jd, mag, err):

    make_new_sn_group(sn_name)
    phot_table = make_new_sn_phot_table(sn_name)

    filter_table = h5file.root.filters

    if filter_eff_wl is not None and filter_flux_zeropoint is not None:
        new_filter, filter_id = make_new_filter_entry(filter_name,
                                                      filter_eff_wl,
                                                      filter_flux_zeropoint,
                                                      note, ref)
        append_new_filter(new_filter)
        write_filter_table()        
    else:
        filter_id = get_old_filter_id(filter_name)

    for i in range(len(jd)):
        new_observation = phot_table.row
        new_observation['filter_id'] = filter_id
        new_observation['jd'] = jd[i]
        new_observation['magnitude'] = mag[i]
        new_observation['uncertainty'] = err[i]
        new_observation['reference'] = ref
        new_observation['note'] = note
        new_observation['telescope_id'] = 0

        new_observation.append()

    phot_table.flush()

    h5file.close()

def add_sn_parameters(sn_name, Av_gal, Av_gal_ref, Av_host, Av_host_ref,
                      distance_Mpc, distance_Mpc_err, distance_Mpc_ref,
                      explosion_JD, explosion_JD_err, explosion_JD_ref,
                      heliocentric_v_kms, heliocentric_v_kms_err,
                      heliocentric_v_kms_ref):

    make_new_sn_group(sn_name)
    parameters_table = make_new_sn_parameters_table(sn_name)

    new_parameters = parameters_table.row
    new_parameters['Av_gal'] = Av_gal
    new_parameters['Av_gal_ref'] = Av_gal_ref
    new_parameters['Av_host'] = Av_host
    new_parameters['Av_host_ref'] = Av_host_ref
    new_parameters['distance_Mpc'] = distance_Mpc
    new_parameters['distance_Mpc_err'] = distance_Mpc_err
    new_parameters['distance_Mpc_ref'] = distance_Mpc_ref
    new_parameters['explosion_JD'] = explosion_JD
    new_parameters['explosion_JD_err'] = explosion_JD_err
    new_parameters['explosion_JD_ref'] = explosion_JD_ref
    new_parameters['heliocentric_v_kms'] = heliocentric_v_kms
    new_parameters['heliocentric_v_kms_err'] = heliocentric_v_kms_err
    new_parameters['heliocentric_v_kms_ref'] = heliocentric_v_kms_ref
    
    new_parameters.append()

    parameters_table.flush()

    h5file.close()
