import argparse
import sys
import hdf5_io

parser = argparse.ArgumentParser(description='Add SN photometry to SuperBoL')
parser.add_argument('input_files', metavar='filename', type=str, nargs='+',
                    help='One or more edited copies of sn_filter.dat')
args = parser.parse_args()

for filename in args.input_files:
    with open(filename, 'r') as f:
        read_phot_data = False
        jd = []
        mag = []
        err = []
        # Grab the photometry from the data file.
        for line in f:
            # Grab the supernova name
            if line.startswith('sn_name'):
                sn_name = line.strip().split()[2]
                if sn_name[0] == '\'' and sn_name[-1] == '\'':
                    sn_name = sn_name[1:-1]
                    
            # Grab the filter name
            if line.startswith('filter_name'):
                filter_name = line.strip().split()[2]
                if filter_name[0] == '\'' and filter_name[-1] == '\'':
                    filter_name = filter_name[1:-1]
                    
            # Grab the effective wavelength of the filter, or set to None
            if line.startswith('filter_eff_wl'):
                filter_eff_wl = line.strip().split()[2]
                if filter_eff_wl == 'None':
                    filter_eff_wl = None
                else:
                    filter_eff_wl = float(filter_eff_wl)
                    
            # Grab the flux zeropoint of the filter, or set to None
            if line.startswith('filter_flux_zeropoint'):
                filter_flux_zeropoint = line.strip().split()[2]
                if filter_flux_zeropoint == 'None':
                    filter_flux_zeropoint = None
                else:
                    filter_flux_zeropoint = float(filter_flux_zeropoint)
                    
            # Grab the reference
            if line.startswith('ref'):
                ref = line.strip().split()[2]
                if ref[0] == '\'' and ref[-1] == '\'':
                    ref = ref[1:-1]
                    
            # Grab the note
            if line.startswith('note'):
                note = line.strip().split()[2:]
                note = " ".join(note).lstrip('\'').rstrip('\'')

            # Finally, grab the photometric data
            if line.strip() == "# START_DATA":
                read_phot_data = True
            elif line.strip() == "# END_DATA":
                read_phot_data = False
            elif read_phot_data:
                obs = [float(x) for x in line.strip().split()]
                jd.append(obs[0])
                mag.append(obs[1])
                err.append(obs[2])

        # Have user check over data
        print "Please examine the input data for any errors:\n"
        print "sn_name = ", sn_name
        print "filter_name = ", filter_name
        print "filter_eff_wl = ", filter_eff_wl
        print "filter_flux_zeropoint = ", filter_flux_zeropoint
        print "reference = ", ref
        print "note = ", note
        print "\n"
        print "Data:"
        for i in range(len(jd)):
            print jd[i], mag[i], err[i]

        # Have user affirm adding data to HDF5 file
        print "\n"
        print "Enter above information into HDF5 file? (y/n)"
        s = raw_input('--> ')
        if s == 'y':
            hdf5_io.add_sn_filter_photometry(sn_name, filter_name,
                                             filter_eff_wl,
                                             filter_flux_zeropoint, ref, note,
                                             jd, mag, err)
        elif s == 'n':
            "Make changes to .dat file and re-run script. Exiting..."
            sys.exit()
        else:
            sys.stdout.write('Please respond with \'y\' or \'n\'.\n')
