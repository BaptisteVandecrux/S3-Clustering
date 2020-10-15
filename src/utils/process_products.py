#!/bin/python3
"""Download and process product from wishlist."""

from netCDF4 import Dataset
from shapely import geometry as geo
import numpy as np
import pandas
import json
import os
from os import path
import argparse
from concurrent.futures import ProcessPoolExecutor
import urllib.parse
import urllib.request
import csv
import getpass
import sys

def load_area_of_interest(geojson_path):
    """Load a geojson for the area of interest."""
    with open(geojson_path) as f:
        area_of_interest = geo.shape(json.load(f)["geometry"])
    return area_of_interest

def check_masks(rootgrp, key, downsample_mask=None):
    """Check if masked_array contains data."""
    if np.all(rootgrp.variables[key].mask):
        print("Warning: No %s data on track" % key)
    elif downsample_mask and np.all(rootgrp.variables[key][downsample_mask].mask):
        print("Warning: No %s data on mask" % key)

def convert_eastward_lon(lon):
    """Convert eastward 0:360 degrees longitude to -180:180 eastward longitude."""
    return lon - 360 * (lon // 180)

def downsample(rootgrp, indexes, field="sig0_ice_sheet_ku", index="index_1hz_meas_20_ku"):
    """Downsample the 20Hz data into 1Hz"""
    n_short = indexes.size
    downsampled_field = np.ma.masked_all(n_short)
    for i, j in enumerate(indexes):
        idx = np.where(rootgrp.variables[index][:] == j) # perf: 38.4%
        d = rootgrp.variables[field][idx].compressed() # perf: 58.1%
        if d.size != 0:
            downsampled_field[i] = np.mean(d)
            downsampled_field.mask[i] = False
    return downsampled_field

def get_coord_mask(lon, lat, area_of_interest):
    """Return the mask indicating if the coordinates in arrays lon and lat are in greenland.

        :type lon: Array[float]
        :param lon: Array of longitudes. Should be in the [-180, 180] format, positive eastward.

        :type lat: Array[float]
        :param lat: Array of latitudes

        :type area_of_interest: shapely.geometry.shape Object
        :param area_of_interest: Geographical Area to filter

        :rtype: Array[Bool]
    """
    return np.array([False if area_of_interest.contains(geo.Point(convert_eastward_lon(lon[i]), lat[i])) else True for i in range(lon.size)])

def check_size_mismatch(rootgrp, key, ref):
    """Check if array has the reference size."""
    if rootgrp.variables[key][:].shape[0] != ref:
        print("Warning: %s size mismatch: %d instead of %d" % (key, rootgrp.variables[key][:].shape[0], ref))

def fetch_local(sat, cycle_number,  data_folder):
    """Download the data for given satellite and cycle number on ESA's OData API."""

    if not path.isdir(data_folder):
        os.mkdir(data_folder)
    if not path.isdir(path.join(data_folder, "%scycle%03d/" % (sat, cycle_number))):
        os.mkdir(path.join(data_folder, "%scycle%03d/" % (sat, cycle_number)))

    new_dict = []
    knames = ["Id", "Name", "Downloaded"]
    i = 0
    with open(path.join(WISHLIST_FOLDER, "%scycle%03d.csv" % (sat, cycle_number)), "r") as csvfile:
        reader = csv.DictReader(csvfile, fieldnames=knames)
        rows = list(reader)
        N = len(rows)
        for row in rows[1:]:
            print("\rGetting %s %d n°%d/%d" % (sat, cycle_number, i, N), end="")
            year = row["Name"][16:20]
            month = row["Name"][20:22]
            day = row["Name"][22:24]

            path_org = "/eodata/Sentinel-3/SRAL/SR_2_LAN/"+year+'/'+month+'/'+ day+'/'+row["Name"]+'.SEN3/enhanced_measurement.nc'
            path_save = DATA_FOLDER+"/%scycle%03d/%s_enhanced.nc" % (sat, 
                                                                 cycle_number, 
                                                                 row["Name"])
            try: 
                os.symlink(path_org,path_save)
            except:
                print('(!)')
                print('Error with '+row["Name"])
                print(sys.exc_info()[0])
                
            new_dict.append(row)
            i += 1

    with open(path.join(WISHLIST_FOLDER, "%scycle%03d.csv" % (sat, cycle_number)), "w") as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=knames)
        writer.writeheader()
        for row in new_dict:
            writer.writerow(row)
    print(" Done!")
    
def fetch_all(data_folder, sat=None):
    """Download all the data in the wishlists."""
    for f in os.listdir(WISHLIST_FOLDER):
        sat_wish = f[:3]
        try: 
            cycle = int(f[8:11])

            if sat is None or sat == sat_wish:
                fetch_local(sat_wish, cycle, data_folder)
        except:
            print('Error with: ' +f)
            
def download(sat, cycle_number, username, password, data_folder):
    """Download the data for given satellite and cycle number on ESA's OData API."""

    if not path.isdir(data_folder):
        os.mkdir(data_folder)
    if not path.isdir(path.join(data_folder, "%scycle%03d/" % (sat, cycle_number))):
        os.mkdir(path.join(data_folder, "%scycle%03d/" % (sat, cycle_number)))

    password_mgr = urllib.request.HTTPPasswordMgrWithDefaultRealm()
    password_mgr.add_password(None, "https://scihub.copernicus.eu/dhus/odata/", username, password)
    handler = urllib.request.HTTPBasicAuthHandler(password_mgr)

    opener = urllib.request.build_opener(handler)
    # ...and install it globally so it can be used with urlopen.
    urllib.request.install_opener(opener)
    new_dict = []
    knames = ["Id", "Name", "Downloaded"]
    i = 0
    with open(path.join(WISHLIST_FOLDER, "%scycle%03d.csv" % (sat, cycle_number)), "r") as csvfile:
        reader = csv.DictReader(csvfile, fieldnames=knames)
        rows = list(reader)
        N = len(rows)
        for row in rows:
            print("\rGetting %s %d n°%d/%d" % (sat, cycle_number, i, N), end="")
            retry = 10
            while retry > 0 and row["Id"] != "Id" and not path.isfile(path.join(data_folder, "%scycle%03d/%s_enhanced.nc" % (sat, cycle_number, row["Name"]))):
                try:
                    url = "https://scihub.copernicus.eu/dhus/odata/v1/Products('%s')/Nodes('%s.SEN3')/Nodes('enhanced_measurement.nc')/$value" % (row["Id"], row["Name"])
                    page = urllib.request.urlopen(url)
                    content = page.read()
                    f = open(path.join(DATA_FOLDER, "%scycle%03d/%s_enhanced.nc" % (sat, cycle_number, row["Name"])), "wb")
                    f.write(content)
                    f.close()
                except Exception as e:
                    retry -= 1
                    print(e)
            if retry != 0 and row["Id"] != "Id":
                row["Downloaded"] = True
            new_dict.append(row)
            i += 1

    with open(path.join(WISHLIST_FOLDER, "%scycle%03d.csv" % (sat, cycle_number)), "w") as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=knames)
        writer.writeheader()
        for row in new_dict:
            writer.writerow(row)
    print(" Done!")

def download_all(username, password, data_folder, sat=None):
    """Download all the data in the wishlists."""
    for f in os.listdir(WISHLIST_FOLDER):
        sat_wish = f[:3]
        cycle = int(f[8:11])
        if sat is None or sat == sat_wish:
            download(sat_wish, cycle, username, password, data_folder)
      
def preprocess_file(filename, cycle_number, filenumber, number_rows, filetotal, area_of_interest, ice_filter=True, waveform=False):
    """Preprocess data for a track file.

        :type filename: string
        :param filename: path to file to process

        :type cycle_number: Int
        :param cycle_number: Cycle number of the data being processed

        :type filenumber: Int
        :param filenumber: Number of the file being processed

        :type number_rows: Int
        :param number_rows: Current number of rows in dataframe

        :type filetotal: Int
        :param filetotal: Total number of files

        :type area_of_interest: shapely.geometry.shape Object
        :param area_of_interest: Geographical Area to filter

        :type ice_filter: Bool
        :param ice_filter: If set, points not flagged as ice by ice_sheet_snow_facies_flag_01_ku will be filtered out

        :rtype: Dataframe
    """
    # Open enhanced measurements file
    rootgrp = Dataset(filename, "r", format="NETCDF4") # perf: 15.2%

    # Check for size mismatch
    n_01 = rootgrp.variables["time_01"][:].size
    check_size_mismatch(rootgrp, "lat_01", n_01)
    check_size_mismatch(rootgrp, "lon_01", n_01)
    check_size_mismatch(rootgrp, "tb_238_01", n_01)
    check_size_mismatch(rootgrp, "tb_365_01", n_01)
    check_size_mismatch(rootgrp, "tb_238_quality_flag_01", n_01)
    check_size_mismatch(rootgrp, "tb_365_quality_flag_01", n_01)
    check_size_mismatch(rootgrp, "ice_sheet_snow_facies_flag_01_ku", n_01)

    n_20_ku = rootgrp.variables["time_20_ku"][:].size
    check_size_mismatch(rootgrp, "sig0_ice_sheet_20_ku", n_20_ku)

    # check_size_mismatch(rootgrp, "sig0_leading_edge_ice_20_ku", n_20_ku)
    # check_size_mismatch(rootgrp, "width_leading_edge_ice_20_ku", n_20_ku)
    # check_size_mismatch(rootgrp, "slope_first_trailing_edge_ice_20_ku", n_20_ku)
    # check_size_mismatch(rootgrp, "slope_second_trailing_edge_ice_20_ku", n_20_ku)
    # check_size_mismatch(rootgrp, "peakiness_1_20_ku", n_20_ku)
    # check_size_mismatch(rootgrp, "peakiness_2_20_ku", n_20_ku)
    check_size_mismatch(rootgrp, "waveform_20_ku", n_20_ku)
    check_size_mismatch(rootgrp, "waveform_qual_ice_20_ku", n_20_ku)

    n_20_c = rootgrp.variables["time_20_c"][:].size
    check_size_mismatch(rootgrp, "sig0_leading_edge_ice_20_plrm_ku", n_20_c)
    check_size_mismatch(rootgrp, "width_leading_edge_ice_20_plrm_ku", n_20_c)
    check_size_mismatch(rootgrp, "slope_first_trailing_edge_ice_20_plrm_ku", n_20_c)
    check_size_mismatch(rootgrp, "slope_second_trailing_edge_ice_20_plrm_ku", n_20_c)
    check_size_mismatch(rootgrp, "peakiness_1_20_plrm_ku", n_20_c)
    check_size_mismatch(rootgrp, "waveform_20_plrm_ku", n_20_c)

    check_size_mismatch(rootgrp, "sig0_ice_sheet_20_c", n_20_c)

    # check_size_mismatch(rootgrp, "sig0_leading_edge_ice_20_c", n_20_c)
    # check_size_mismatch(rootgrp, "width_leading_edge_ice_20_c", n_20_c)
    # check_size_mismatch(rootgrp, "slope_first_trailing_edge_ice_20_c", n_20_c)
    # check_size_mismatch(rootgrp, "slope_second_trailing_edge_ice_20_c", n_20_c)
    # check_size_mismatch(rootgrp, "peakiness_1_20_c", n_20_c)
    # check_size_mismatch(rootgrp, "peakiness_2_20_c", n_20_c)
    check_size_mismatch(rootgrp, "waveform_20_c", n_20_c)

    # Check quality flags
    mask = np.zeros(n_01, dtype=np.bool)
    mask = mask | (rootgrp.variables["tb_238_quality_flag_01"][:] == 1)
    mask = mask | (rootgrp.variables["tb_365_quality_flag_01"][:] == 1)

    # Filter geographic data
    geo_mask = get_coord_mask(rootgrp.variables["lon_01"][:], rootgrp.variables["lat_01"][:], area_of_interest) # perf: 12.6%
    mask = mask | geo_mask

    # Add data masks
    mask = mask | rootgrp.variables["time_01"][:].mask
    mask = mask | rootgrp.variables["lat_01"][:].mask
    mask = mask | rootgrp.variables["lon_01"][:].mask
    mask = mask | rootgrp.variables["tb_238_01"][:].mask
    mask = mask | rootgrp.variables["tb_365_01"][:].mask
    mask = mask | rootgrp.variables["ice_sheet_snow_facies_flag_01_ku"][:].mask

    # Filter if the ice_sheet_snow_facies is not set
    if ice_filter:
        mask = mask | (rootgrp.variables["ice_sheet_snow_facies_flag_01_ku"][:] == 0)

    # Reduced index
    indexes_01 = np.where(~mask)[0]

    if indexes_01.size == 0:
        return pandas.DataFrame()

    # Reducing arrays before downsampling for performance
    time_01 = rootgrp.variables["time_01"][indexes_01]
    lat_01 = rootgrp.variables["lat_01"][indexes_01]
    lon_01 = rootgrp.variables["lon_01"][indexes_01]
    tb_238_01 = rootgrp.variables["tb_238_01"][indexes_01]
    tb_365_01 = rootgrp.variables["tb_365_01"][indexes_01]
    ice_sheet_snow_facies_flag_01_ku = rootgrp.variables["ice_sheet_snow_facies_flag_01_ku"][indexes_01]

    if waveform:
        waveform_20_ku = rootgrp.variables["waveform_20_ku"][rootgrp.variables["index_first_20hz_meas_01_ku"][indexes_01], :]
        waveform_20_c = rootgrp.variables["waveform_20_c"][rootgrp.variables["index_first_20hz_meas_01_c"][indexes_01], :]
        # waveform_20_plrm_ku = rootgrp.variables["waveform_20_plrm_ku"][rootgrp.variables["index_first_20hz_meas_01_ku"][indexes_01], :]

    sig0_leading_edge_ice_20_ku = rootgrp.variables["sig0_leading_edge_ice_20_ku"][rootgrp.variables["index_first_20hz_meas_01_ku"][indexes_01]]
    width_leading_edge_ice_20_ku = rootgrp.variables["width_leading_edge_ice_20_ku"][rootgrp.variables["index_first_20hz_meas_01_ku"][indexes_01]]
    slope_first_trailing_edge_ice_20_ku = rootgrp.variables["slope_first_trailing_edge_ice_20_ku"][rootgrp.variables["index_first_20hz_meas_01_ku"][indexes_01]]
    slope_second_trailing_edge_ice_20_ku = rootgrp.variables["slope_second_trailing_edge_ice_20_ku"][rootgrp.variables["index_first_20hz_meas_01_ku"][indexes_01]]
    peakiness_1_20_ku = rootgrp.variables["peakiness_1_20_ku"][rootgrp.variables["index_first_20hz_meas_01_ku"][indexes_01]]
    peakiness_2_20_ku = rootgrp.variables["peakiness_2_20_ku"][rootgrp.variables["index_first_20hz_meas_01_ku"][indexes_01]]
    #
    sig0_leading_edge_ice_20_c = rootgrp.variables["sig0_leading_edge_ice_20_c"][rootgrp.variables["index_first_20hz_meas_01_c"][indexes_01]]
    width_leading_edge_ice_20_c = rootgrp.variables["width_leading_edge_ice_20_c"][rootgrp.variables["index_first_20hz_meas_01_c"][indexes_01]]
    slope_first_trailing_edge_ice_20_c = rootgrp.variables["slope_first_trailing_edge_ice_20_c"][rootgrp.variables["index_first_20hz_meas_01_c"][indexes_01]]
    slope_second_trailing_edge_ice_20_c = rootgrp.variables["slope_second_trailing_edge_ice_20_c"][rootgrp.variables["index_first_20hz_meas_01_c"][indexes_01]]
    peakiness_1_20_c = rootgrp.variables["peakiness_1_20_c"][rootgrp.variables["index_first_20hz_meas_01_c"][indexes_01]]
    peakiness_2_20_c = rootgrp.variables["peakiness_2_20_c"][rootgrp.variables["index_first_20hz_meas_01_c"][indexes_01]]
    waveform_qual_ice_20_ku = rootgrp.variables["waveform_qual_ice_20_ku"][rootgrp.variables["index_first_20hz_meas_01_c"][indexes_01]]

    sig0_leading_edge_ice_20_plrm_ku = rootgrp.variables["sig0_leading_edge_ice_20_plrm_ku"][rootgrp.variables["index_first_20hz_meas_01_c"][indexes_01]]
    width_leading_edge_ice_20_plrm_ku = rootgrp.variables["width_leading_edge_ice_20_plrm_ku"][rootgrp.variables["index_first_20hz_meas_01_c"][indexes_01]]
    slope_first_trailing_edge_ice_20_plrm_ku = rootgrp.variables["slope_first_trailing_edge_ice_20_plrm_ku"][rootgrp.variables["index_first_20hz_meas_01_c"][indexes_01]]
    slope_second_trailing_edge_ice_20_plrm_ku = rootgrp.variables["slope_second_trailing_edge_ice_20_plrm_ku"][rootgrp.variables["index_first_20hz_meas_01_c"][indexes_01]]
    peakiness_1_20_plrm_ku = rootgrp.variables["peakiness_1_20_plrm_ku"][rootgrp.variables["index_first_20hz_meas_01_c"][indexes_01]]



    # Downsample 20Hz data
    sig0_ice_sheet_01_ku = downsample(rootgrp, indexes_01, "sig0_ice_sheet_20_ku", "index_1hz_meas_20_ku") # perf: 14.1%
    sig0_ice_sheet_01_c = downsample(rootgrp, indexes_01, "sig0_ice_sheet_20_c", "index_1hz_meas_20_c") # perf: 13.8%

    rootgrp.close()
    # Create new mask for downsampled data
    mask = sig0_ice_sheet_01_c.mask | sig0_ice_sheet_01_ku.mask

    # Reducing index
    indexes_01 = np.where(~mask)[0]

    if indexes_01.size == 0:
        return pandas.DataFrame()

    ## Return all fields of interest
    df = pandas.DataFrame({
        "time_01": time_01[indexes_01],
        "lat_01": lat_01[indexes_01],
        "lon_01": convert_eastward_lon(lon_01[indexes_01]),
        "tb_238_01": tb_238_01[indexes_01],
        "tb_365_01": tb_365_01[indexes_01],
        "ice_sheet_snow_facies_flag_01_ku": ice_sheet_snow_facies_flag_01_ku[indexes_01],

        "sig0_ice_sheet_01_ku": sig0_ice_sheet_01_ku[indexes_01],
        "sig0_ice_sheet_01_c": sig0_ice_sheet_01_c[indexes_01],

        "sig0_leading_edge_ice_20_ku": sig0_leading_edge_ice_20_ku[indexes_01],
        "width_leading_edge_ice_20_ku": width_leading_edge_ice_20_ku[indexes_01],
        "slope_first_trailing_edge_ice_20_ku": slope_first_trailing_edge_ice_20_ku[indexes_01],
        "slope_second_trailing_edge_ice_20_ku": slope_second_trailing_edge_ice_20_ku[indexes_01],
        "peakiness_1_20_ku": peakiness_1_20_ku[indexes_01],
        "peakiness_2_20_ku": peakiness_2_20_ku[indexes_01],
        "waveform_qual_ice_20_ku": waveform_qual_ice_20_ku[indexes_01],

        "sig0_leading_edge_ice_20_c": sig0_leading_edge_ice_20_c[indexes_01],
        "width_leading_edge_ice_20_c": width_leading_edge_ice_20_c[indexes_01],
        "slope_first_trailing_edge_ice_20_c": slope_first_trailing_edge_ice_20_c[indexes_01],
        "slope_second_trailing_edge_ice_20_c": slope_second_trailing_edge_ice_20_c[indexes_01],
        "peakiness_1_20_c": peakiness_1_20_c[indexes_01],
        "peakiness_2_20_c": peakiness_2_20_c[indexes_01],

        "sig0_leading_edge_ice_20_plrm_ku": sig0_leading_edge_ice_20_plrm_ku[indexes_01],
        "width_leading_edge_ice_20_plrm_ku": width_leading_edge_ice_20_plrm_ku[indexes_01],
        "slope_first_trailing_edge_ice_20_plrm_ku": slope_first_trailing_edge_ice_20_plrm_ku[indexes_01],
        "slope_second_trailing_edge_ice_20_plrm_ku": slope_second_trailing_edge_ice_20_plrm_ku[indexes_01],
        "peakiness_1_20_plrm_ku": peakiness_1_20_plrm_ku[indexes_01]
    })

    # Save waveform data
    if waveform:
        df = df.join(pandas.DataFrame(waveform_20_ku[indexes_01, :]).add_prefix("waveform_20_ku_"))
        df = df.join(pandas.DataFrame(waveform_20_c[indexes_01, :]).add_prefix("waveform_20_c_"))
        # df = df.join(pandas.DataFrame(waveform_20_plrm_ku[0, indexes_01, :]).add_prefix("waveform_20_plrm_ku_"))

    print("\rCycle %d: File %d/%d: %d rows, %d added" % (cycle_number, filenumber, filetotal, number_rows, df.size), end="")

    return df

def preprocess_file_map(args):
    """Like preprocess_file but takes a array of parameters instead of multiple parameters"""
    return preprocess_file(*args)


def preprocess_all(data_folder, preproc_folder, area_of_interest, n_threads=1, ice_filter=True, waveform=False, sat=None):
    """Preprocess all cycle. A satellite can be defined with sat."""
    cycles = os.listdir(data_folder)
    cycles.sort()
    for cycle in cycles:
        cyclesat = cycle[:3]
        cyclenum = int(cycle[8:])
        if sat is not None:
            if cyclesat == sat:
                preprocess_cycle(cyclesat, cyclenum, data_folder, preproc_folder, area_of_interest, n_threads, ice_filter, waveform)
        else:
            preprocess_cycle(cyclesat, cyclenum, data_folder, preproc_folder, area_of_interest, n_threads, ice_filter, waveform)


def preprocess_cycle(satellite, cycle_number, data_folder, preproc_folder, area_of_interest, n_threads=1, ice_filter=True, waveform=False):
    """Preprocess data for a given cycle and satellite.

        :type satellite: str
        :param satellite: Codename of the satellite to processed

        :type cycle_number: int
        :param cycle_number: Cycle number to process

        :type data_folder: str
        :param data_folder: Path to the folder containing the satellite data

        :type preproc_folder: str
        :param preproc_folder: Path to the output folder to save the preprocessed data into.

        :type area_of_interest: shapely.geometry.shape Object
        :param area_of_interest: Geographical Area to filter

        :type n_threads: int
        :param n_threads: Number of threads to use. If below 2, will not use multithreading.

        :rtype: None
    """
    if not path.isdir(preproc_folder):
        os.mkdir(preproc_folder)

    # Get all files
    files = os.listdir(os.path.join(data_folder, "%scycle%03d" % (satellite, cycle_number)))

    # Iterate over files
    n = len(files)
    if n_threads < 2:
        i = 1
        df = pandas.DataFrame()
        for f in files:
            filedf = preprocess_file(os.path.join(data_folder, "%scycle%03d" % (satellite, cycle_number), f), cycle_number, i, df.size, n, area_of_interest, ice_filter, waveform)
            df = pandas.concat([df, filedf], ignore_index=True)
            i += 1
    else:
        with ProcessPoolExecutor(max_workers=n_threads) as executor:
            dfs = executor.map(preprocess_file_map, [[os.path.join(data_folder, "%scycle%03d" % (satellite, cycle_number), f), cycle_number, i, 0, n, area_of_interest, ice_filter, waveform] for i, f in enumerate(files)])
        df = pandas.concat(dfs, ignore_index=True)
        df = df.sort_values(by="time_01", ignore_index=True)

    # Save data
    df.to_hdf(os.path.join(preproc_folder, "%scycle%03d.hdf" % (satellite, cycle_number)), key="%scycle%03d" % (satellite, cycle_number))



CURRENT_FOLDER = path.dirname(os.path.abspath(__file__))
PREPROC_FOLDER = path.join(CURRENT_FOLDER, "../../data/PREPROC")
DATA_FOLDER = path.join(CURRENT_FOLDER, "../../data/PRODUCT")
WISHLIST_FOLDER = path.join(CURRENT_FOLDER, "../../data/wishlists")
GEOJSON_PATH = path.join(CURRENT_FOLDER, "../assets/maritime_grl.geo.json")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Preprocess Sentinel-3 data to filter a geographic area and invalid data and keep only the fields needed.')
    parser.add_argument('--platform', type=str, help='Codename for the platform used. Can be either S3A or S3B.', default=None)
    parser.add_argument('--cycle', type=int, help='Satellite cycle number.', default=None)
    parser.add_argument('--geojson', type=str, help='Path to the shapely geojson describing the region to study.', default=GEOJSON_PATH)
    parser.add_argument('--outfolder', type=str, help='Folder where to store the results.', default=PREPROC_FOLDER)
    parser.add_argument('--infolder', type=str, help='Folder where the data product are stored.', default=DATA_FOLDER)
    parser.add_argument('--threads', type=int, help='Number of threads to run.', default=1)
    parser.add_argument('--no_waveform', action='store_false', help='Ignore waveforms.')
    parser.add_argument('--no_ice_filter', action='store_false', help='Deactivate the filtering based on snow facies classification.')
    parser.add_argument('--download', action='store_true', help='Download the selected cycle.')
    parser.add_argument('--fetch_local', action='store_true', help='Fetch files on eodata folder for the selected cycle.')


    args = parser.parse_args()

    if args.platform is None and args.cycle is not None:
        raise Exception("Cannot define cycle without platform")

    if args.download:
        username = input("Username: ")
        password = getpass.getpass("Password: ")
        if args.cycle is not None:
            download(args.platform, args.cycle, username, password, args.infolder)
        else:
            download_all(username, password, args.infolder, sat=args.platform)
        exit()
        
    if args.fetch_local:
        if args.cycle is not None:
            fetch_local(args.platform, args.cycle, args.infolder)
        else:
            fetch_all(args.infolder, sat=args.platform)
        exit()

    if args.cycle is not None:
        preprocess_cycle(args.platform, args.cycle, args.infolder, args.outfolder, load_area_of_interest(args.geojson), args.threads, args.no_ice_filter, args.no_waveform)
    else:
        preprocess_all(args.infolder, args.outfolder, load_area_of_interest(args.geojson), args.threads, args.no_ice_filter, args.no_waveform, sat=args.platform)

