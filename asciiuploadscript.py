# Author: Igor Andreoni andreoni@caltech.edu, modified by Kaustav Das and Tawny Sit

import requests
import glob
import re
import numpy as np
from astropy.io import ascii, fits
from astropy.table import Table
from astropy.time import Time
import json

# Copy the API token from your Fritz account
token = '8fea7814-af0d-42e8-ae9e-deff36e86bc4'

spec_files = glob.glob('*.spec')
for file in spec_files:
    with open(file, 'r') as f:
        lines = f.readlines()
    cleaned_lines = [line for line in lines if 'NaN' not in line]
    with open(file, 'w') as f:
        f.writelines(cleaned_lines)
print('are you done?')

"""
Default observers and reducers (these are NOT necessary to be set,
since they will be used only when the -d option is called).
User IDs can be found at: https://fritz.science/api/user
"""
default_observer = [1296,1268] # change
default_reducer = [155] # this is my ID


def api(method, endpoint, data=None):
    """API request"""

    headers = {'Authorization': f'token {token}'}
    response = requests.request(method, endpoint, json=data, headers=headers)
    return response


def str2bool(v):
    if v.lower() in ('yes', 'true', 't', 'y', '1', 'Yes', 'True'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0', 'No', 'False'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')


def check_source_exists(source):
    """check if the source exists"""

    response = api('HEAD', f'https://fritz.science/api/sources/{source}')

    if response.status_code == 200:
        print(f"Source {source} was found on Fritz")
    else:
        print(f"Source {source} does not exist on Fritz!!")

    return response


def get_groups(source):
    """Get the groups a source belongs to"""

    response = api('GET',
                   f'https://fritz.science/api/sources/{source}'
                   )

    if response.status_code == 200:
        groups = response.json()['data']['groups']
    else:
        print(f'HTTP code: {response.status_code}, {response.reason}')

    return groups


def upload_spectrum(spec, observers, reducers, group_ids=[], date=None,
                    inst_id=3, ztfid=None, meta=None):
    """
    Upload a spectrum to the Fritz marshal
    ----
    Parameters
    spec astropy table object
        table with wavelength, flux, fluxerr
    observers list of int
        list of integers corresponding to observer usernames
    reducers list of int
        list of integers corresponding to reducer usernames
    group_ids list of int
        list of IDs of the groups to share the spectum with
    date str
         date (UT) of the spectrum
    inst_id int
        ID of the instrument (default DBSP)
    ztfid str
        ID of the ZTF source, e.g. ZTF20aclnxgz
    """

    # Remove meta from spec, they are already stored
    #print(meta)
    #spec._meta = meta # do not remove meta

    #spec._meta = None
    #Get rid of NaNs in lpipe output since fritz does not like NaNs and inf
    good_rows = ~np.isnan(spec['flux'])
    usespec = spec[good_rows]
    good_rows = ~np.isinf(spec['flux'])
    usespec = spec[good_rows]
    good_rows = ~np.isnan(spec['fluxerr'])
    usespec = spec[good_rows]
    good_rows = ~np.isinf(spec['fluxerr'])
    usespec = spec[good_rows]

    print(type(usespec))

    filename_ = ztfid+'_'+str(date)[:10]+'.ascii'


    from astropy.io import ascii
    from astropy.table import Table


    spec = usespec

    
    meta_str = '\n'.join([f"# {key}: {value}" for key, value in spec.meta.items()])

    with open(f'./{filename_}', 'w') as f:
        # Write the metadata
        f.write(meta_str + '\n')
        # Write the table data
        ascii.write(spec, f, format='no_header', overwrite=True)
        #ascii.write(usespec, f'./{filename_}', format='no_header', overwrite=True)

    print(ztfid, type(ztfid))
    print(filename_, type(filename_))

    

    with open(f'./{filename_}') as f:
                        data={
                        'wave_column': 0,
                        'flux_column': 1,
                        'fluxerr_column': 2,
                        'ascii': f.read(),
                        'group_ids': [1],
                        "observed_by": observers,
                        #"altdata": meta,
                        "observed_at": str(date),
                        "instrument_id": inst_id,
                        "reduced_by": reducers,
                        #"origin": "",
                        "obj_id": str(ztfid),
                        "filename": str(filename_)
    
                        # here add all the fields you already have here that you might need for bookkeeping:
                        # so like the observers, reducers, ... and other things.
                        }
                


    print(data['filename'])
    response = api('POST',
                   'https://fritz.science/api/spectrum/ascii',
                   data = data)

    print(f'HTTP code: {response.status_code}, {response.reason}')
    if response.status_code == 400:
        print(f'JSON response: {response.json()}')
           
    #data = {
    #        "observed_by": observers,
    #        #"group_ids": group_ids,
    #        "group_ids": [1],
    #        "altdata": meta,
    #        "observed_at": str(date),
    #       "fluxes": list(usespec['flux']),
    #        "errors": list(usespec['fluxerr']),
    #        "wavelengths": list(usespec['wavelength']),
    #        "instrument_id": inst_id,
    #        "reduced_by": reducers,
    #        "origin": "",
    #        "obj_id": ztfid
    #        }
    #response = api('POST',
    #               'https://fritz.science/api/spectra',
    #               data=data)

 

def crop_spectrum(spectrum):
    """
    Crop the spectra
    ----
    Parameters
        spec astropy.table object
            table with the spectrum
    ----
    Returns
        cropped spectrum
    """

    # Define the ranges to crop the spectra
    default_crop = str2bool(input("Are you happy with the default crop \
between 3,200A and 10,000A? \n"))
    if default_crop is True:
        crop_range = (3200, 10000)
    else:
        happy = False
        while happy is False:
            raw_range = input(f"Input the desired range (two numbers, \
comma-separated): \n").split(",")
            crop_range = tuple([int(s) for s in raw_range])
            print(f"Your range for cropping will be {crop_range}")
            happy = str2bool(input(f"Please confirm that this \
is the range you want [y/n]:\n"))
    spectrum = spectrum[spectrum['wavelength'] < np.max(crop_range)]
    spectrum = spectrum[spectrum['wavelength'] > np.min(crop_range)]

    return spectrum


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Upload spectra to Fritz.')
    parser.add_argument('ID', nargs='+', help='Spectra filenames')
    parser.add_argument('--date', dest='date', type=str,
                        required=True, help='Date of the observations (UT), \
for example 2020-11-10T00:00:00 or 2021-01-12')
    parser.add_argument('--inst', dest='inst_id', type=int,
                        required=True, help='Instrument ID, \
for example. inst_id = 3 for DBSP, inst_id = 7 for LRIS. \
Instrument IDs can be found at: \
https://fritz.science/api/instrument')
    parser.add_argument('-d', action='store_true',
                        help='Use default reducer ID and observer; \
please customize the code if you want to use this (unrequired) option')
    args = parser.parse_args()

    # Observers, reducers
    if args.d is True:
        #default observer
        observers = default_observer
        #default reducer
        reducers = default_reducer
    else:
        print("users IDs can be found at: https://fritz.science/api/user")
        observers = input("Enter comma-separated IDs (int) \
of the observers:\n")
        observers = observers.split(",")
        reducers = input("Enter comma-separated IDs (int) of the reducers:\n")
        reducers = reducers.split(",")

    # For each ID, check which files are available
    for source_filename in args.ID:
        files = glob.glob(f"{source_filename}")
        if len(files) == 0:
            print(f"No files named {source_filename} found")
            filename = input(f"Enter the correct spectrum file name \
or enter 'c' to skip this source and continue: \n")
        else:
            filename = files[0]
        if filename == 'c':
            continue
        elif filename[-4:] == 'fits':
            hdul = fits.open(filename)
            spec = hdul[-1].data
            
            # Fix the table format
            spec = Table(hdul[-1].data)
            spec.rename_column("wave", "wavelength")
            spec.rename_column("sigma", "fluxerr")
            
            # build header
            for i, hud in enumerate(hdul): # saves info from all headers - might be too much/redundant, but better safe than sorry
                header = hdul[i].header
                if i==0: # contains coadd file paths and software versions (possibly useful for debugging)
                    fileext = 'VERSION'
                    spec.meta = {'header': {fileext: dict(header)}}
                else: # index headers by file extension name (eg. RED00XX/BLUE00XX = header from raw exposures, RED/BLUE/SPLICED = reduced spectrum)
                    fileext = header['EXTNAME']
                    spec.meta['header'][str(fileext)] = dict(header)
                    try:
                        del spec.meta['header'][str(fileext)]['COMMENT'] # this is just info about the FITS file format, not needed or relevant
                    except:
                        pass

        else:
            # Read the file as ascii
            spec = ascii.read(filename)
            header = None
            spec.rename_column("col1", "wavelength")
            spec.rename_column("col2", "flux")
            # Uncertainty for DBSP
            if len(spec.colnames) > 2 and args.inst_id == 3:
                spec.rename_column("col3", "fluxerr")
            elif len(spec.colnames) > 2 and args.inst_id == 7:
                spec.rename_column("col4", "fluxerr")
            elif len(spec.colnames) == 2:
                spec["fluxerr"] = np.zeros(len(spec))

        # Crop DBSP spectra
        if args.inst_id == 3:
            print(f"The current wavelength range is between \
{np.round(np.min(spec['wavelength']))} and \
{np.round(np.max(spec['wavelength']))}.")
            do_crop = str2bool(input("Do you want to crop the spectrum? \n"))
            if do_crop is True:
                spec = crop_spectrum(spec)

        # Metadata
        # For LRIS:
        if args.inst_id == 7:
            meta = spec._meta['comments']
        # other instruments:
        else:
            try:
                meta = spec.meta['header']
            except KeyError:
                print("WARNING! Unable to find header metadata")
                print("No metadata will be uploaded")
                meta = None
        # Extract the source filename
        if not "ZTF" in source_filename:
            source = input(f"No 'ZTF' found in the file name, please enter \
the name of the source for the spectrum {source_filename}:\n")
        else:
            span = re.search("ZTF", source_filename).span()
            source = source_filename[span[0]:span[0]+12]
        response = check_source_exists(source)
        if response.status_code != 200:
            print(f"Skipping {source}...")
            continue

        # Groups
        groups = get_groups(source)
        group_ids = []
        if len(groups) == 0:
            print(f"{source} is not saved in any group")
        else:
            print(f"{source} was saved in the following groups:")
        for g in groups:
            print(f"{g['id']}, {g['name']}")
            group_ids.append(g['id'])
        print(f"The spectrum will be sent to groups {group_ids}")
        ok_groups = str2bool(input("Are you happy with that?\n") or 'y')
        if ok_groups is False:
            group_ids = input("Please enter comma-separated group IDs (int)\n")
            group_ids = group_ids.split(",")
            group_ids = [int(g) for g in group_ids]

        # Upload
        print(f"Uploading spectrum {source_filename} for source {source} \
to group IDs {group_ids}")
        #new_fname = ztfid+'_'+str(date)[:10]+'.ascii'
        #spec.meta['comments']  = list([key + '=' +  str(header[key]) for key in list(header.keys())])
        #spec.write(new_fname, format='ascii.no_header', overwrite=True)
        #print(meta)
        upload_spectrum(spec, observers, reducers, group_ids=group_ids,
                        date=Time(args.date).iso, inst_id=args.inst_id,
                        ztfid=source, meta=None)

