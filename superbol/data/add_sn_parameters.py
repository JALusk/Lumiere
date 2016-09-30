import argparse
import sys
import hdf5_io

parser = argparse.ArgumentParser(description='Add SN parameters to SuperBoL')
parser.add_argument('input_files', metavar='filename', type=str, nargs='+',
                    help='One or more edited copies of sn_parameters.dat')
args = parser.parse_args()

for filename in args.input_files:
    with open(filename, 'r') as f:
        for line in f:
            # Grab the supernova name
            if line.startswith('sn_name'):
                sn_name = line.strip().split()[2]
                if sn_name[0] == '\'' and sn_name[-1] == '\'':
                    sn_name = sn_name[1:-1]
                    
            # Grab the total MWG extinction in the V-band
            if line.startswith('Av_gal '):
                Av_gal = line.strip().split()[2]
                    
            # Grab the reference for the MWG extinction in the V-band
            if line.startswith('Av_gal_ref'):
                Av_gal_ref = line.strip().split()[2:]
                Av_gal_ref = " ".join(Av_gal_ref).lstrip('\'').rstrip('\'')

            # Grab the total host extinction in the V-band
            if line.startswith('Av_host '):
                Av_host = line.strip().split()[2]
                    
            # Grab the reference for the host extinction in the V-band
            if line.startswith('Av_host_ref'):
                Av_host_ref = line.strip().split()[2:]
                Av_host_ref = " ".join(Av_host_ref).lstrip('\'').rstrip('\'')

            # Grab the distance to the SN in megaparsecs
            if line.startswith('distance_Mpc '):
                distance_Mpc = line.strip().split()[2]
                    
            # Grab the uncertainty in the distance to the SN in megaparsecs
            if line.startswith('distance_Mpc_err'):
                distance_Mpc_err = line.strip().split()[2]
            
            # Grab the reference for the host extinction in the V-band
            if line.startswith('distance_Mpc_ref'):
                distance_Mpc_ref = line.strip().split()[2:]
                distance_Mpc_ref = " ".join(distance_Mpc_ref).lstrip('\'').rstrip('\'')

            # Grab the eplosion JD
            if line.startswith('explosion_JD '):
                explosion_JD = line.strip().split()[2]

            # Grab the uncertainty in the explosion JD
            if line.startswith('explosion_JD_err'):
                explosion_JD_err = line.strip().split()[2]

            # Grab the reference for the explosion JD
            if line.startswith('explosion_JD_ref'):
                explosion_JD_ref = line.strip().split()[2:]
                explosion_JD_ref = " ".join(explosion_JD_ref).lstrip('\'').rstrip('\'')

            # Grab the heliocentric velocity in km/s 
            if line.startswith('heliocentric_v_kms '):
                heliocentric_v_kms = line.strip().split()[2]

            # Grab the heliocentric velocity in km/s
            if line.startswith('heliocentric_v_kms_err'):
                heliocentric_v_kms_err = line.strip().split()[2]

            # Grab the reference for the heliocentric velocity
            if line.startswith('heliocentric_v_kms_ref'):
                heliocentric_v_kms_ref = line.strip().split()[2:]
                heliocentric_v_kms_ref = " ".join(heliocentric_v_kms_ref).lstrip('\'').rstrip('\'')

        # Have user check over data
        print("Please examine the input data for any errors:\n")
        print("sn_name = ", sn_name)
        print("distance_Mpc = ", distance_Mpc)
        print("distance_Mpc_err = ", distance_Mpc_err)
        print("distance_Mpc_ref = ", distance_Mpc_ref)
        print("Av_gal = ", Av_gal)
        print("Av_gal_ref = ", Av_gal_ref)
        print("Av_host = ", Av_host)
        print("Av_host_ref = ", Av_host_ref)
        print("explosion_JD = ", explosion_JD)
        print("explosion_JD_err = ", explosion_JD_err)
        print("explosion_JD_ref = ", explosion_JD_ref)
        print("heliocentric_v_kms = ", heliocentric_v_kms)
        print("heliocentric_V_kms_err = ", heliocentric_v_kms_err)
        print("heliocentric_v_kms_ref = ", heliocentric_v_kms_ref)

        # Have user affirm adding data to HDF5 file
        print("\n")
        print("Enter above information into HDF5 file? (y/n)")
        s = input('--> ')
        if s == 'y':
            hdf5_io.add_sn_parameters(sn_name, Av_gal, Av_gal_ref, Av_host,
                                      Av_host_ref, distance_Mpc,
                                      distance_Mpc_err, distance_Mpc_ref,
                                      explosion_JD, explosion_JD_err,
                                      explosion_JD_ref, heliocentric_v_kms,
                                      heliocentric_v_kms_err,
                                      heliocentric_v_kms_ref)
        elif s == 'n':
            "Make changes to .dat file and re-run script. Exiting..."
            sys.exit()
        else:
            sys.stdout.write('Please respond with \'y\' or \'n\'.\n')
