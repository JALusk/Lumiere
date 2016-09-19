# Formatting your data for SuperBoL

In order to make the process of adding your own photometric data to SuperBoL
as simple as possible, please follow these steps:

1. Make a copy of `sn_filter.dat` and rename it
    - for example, `sn2016a_U.dat`
2. Repeat for each filter used to observe that supernova
3. Edit your new .dat files, filling in the variables and adding your photometry
4. Make a copy of `sn_parameters.dat` and rename it
    - for example, `sn2016a_parameters.dat`
5. Edit this file, filling in the variables

# Adding your data to the HDF5 file

Once you have edited copies of the `sn_filter.dat` and `sn_parameters.dat`
files, you need to run the included scripts to actually include that data
in the HDF5 file.

1. To add the parameters to the HDF5 file, run

`python add_sn_parameters.py your_sn_parameters_file.dat`

where `your_sn_parameters_file.dat` is the name of your edited copy of
`sn_parameters.dat`

This will prompt you to look over the data and confirm it is as you expect.
Once you verify the data and type `y`, your data will be added to the HDF5
file.

2. To add your photometry to the HDF5 file, run

`python add_sn_filter.py your_sn_filter1.dat your_sn_filter2.dat ...`

where the filenames that follow the script will be parsed and added in order.
For each set of data, you will be asked to verify the data and confirm it is
as you expect.
Once you verify the data and type `y`, your data will be added to the HDF5
file.
