import numpy as np

def zip_photometry(photometry, dt=0.3):
    sorted_phot = sort_photometry(photometry)
    binned_phot = bin_JDs(sorted_phot, dt)
    zipped_phot = combine_repeated_detections(binned_phot)

    return zipped_phot


def sort_photometry(photometry):
    return np.sort(photometry, order='jd')

def bin_JDs(sorted_photometry, dt=0.3):
    """Combine groups of observations separated by more than dt days"""
    
    break_index = 0
    deltas = np.diff(sorted_photometry['jd'])
    num_detections = len(sorted_photometry['jd'])

    for i, delta in enumerate(deltas):
        if delta > dt:
            # Found a break from one night to the next
            # Group the previous night's observations together
            group = sorted_photometry['jd'][break_index:i+1]
            mean_jd = np.mean(group)
            # Replace the JD's of the group elements with the mean of the group
            for j in range(len(group)):
                sorted_photometry['jd'][break_index + j] = mean_jd

            break_index = i + 1

    last_group = sorted_photometry['jd'][break_index - num_detections:num_detections]
    mean_last_group = np.mean(last_group)

    for k in range(len(last_group)):
        sorted_photometry['jd'][break_index + k] = mean_last_group

    return sorted_photometry

def combine_repeated_detections(binned_phot):
    """Average together observations taken with the same filter on the same
    JD (after binning)"""
    
    JD, jdindex = np.unique(binned_phot['jd'], return_index = True)
    zipped_photometry = np.empty(np.shape(binned_phot[0]), dtype = binned_phot.dtype)
    # Split into sub-arrays for each days' detections
    for subarray in np.split(binned_phot, jdindex[1:]):
        jd = subarray['jd'][0]
        name, nameindex, count = np.unique(subarray['name'], return_index = True, return_counts = True)
        # Make an array of filter names which appear more than once
        repeated_filters = name[count > 1]
        subarray_no_repeats = np.delete(subarray, [n for n in range(len(subarray['name'])) if subarray['name'][n] in repeated_filters])
        # Combine the repeated detections together by averaging magnitudes
        for repeated_filter in repeated_filters:
            average_mag = np.mean(subarray['magnitude'][np.where(subarray['name'] == repeated_filter)])
            err = np.sqrt(np.sum(subarray['uncertainty'][np.where(subarray['name'] == repeated_filter)]**2))
            averaged_detection = np.array([(subarray['jd'][0], repeated_filter, average_mag, err)], dtype=binned_phot.dtype)
                                                            
            subarray_no_repeats = np.append(subarray_no_repeats, averaged_detection)
        zipped_photometry = np.append(zipped_photometry, subarray_no_repeats)
    zipped_photometry = np.delete(zipped_photometry, 0)

    return zipped_photometry
