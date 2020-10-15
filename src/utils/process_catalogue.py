#!/bin/python3
"""Extract the products needed from the catalogue."""
import os
from os import path
import csv
import pandas as pd

CURRENT_FOLDER = path.dirname(os.path.abspath(__file__))
CATALOG_PATH = path.join(CURRENT_FOLDER, "../../data/catalogueview")
WISHLIST_PATH = path.join(CURRENT_FOLDER, "../../data/wishlists")

def process_catalog(sat):
    """Process the catalog to only keep the requiered products."""
    knames = []
    with open(path.join(CATALOG_PATH, "%s_catalogue_filtered.csv" % sat), 'w') as csvfile:
        knames = ['Id', 'Name', 'ContentLength', 'IngestionDate',
                  'ContentDate:Start', 'ContentDate:End', 'Checksum:Algorithm',
                  'Checksum:Value', 'cycle_number', 'orbit_number', 'length'
                 ]
        writer = csv.DictWriter(csvfile, fieldnames=knames)
        writer.writeheader()

        # Iterate over all folders
        years = os.listdir(path.join(CATALOG_PATH, sat))
        years.sort()
        for year in years:
            months = os.listdir(path.join(CATALOG_PATH, sat, year))
            months.sort()
            for month in months:
                print("\rProcessing %s %s %s." % (sat, year, month), end="")
                for filename in os.listdir(path.join(CATALOG_PATH, sat, year, month)):
                    with open(path.join(CATALOG_PATH, sat, year, month, filename)) as csvfile:
                        reader = csv.DictReader(csvfile)
                        for row in reader:
                            # Filter product for SR_2_LAN_ products and NT_003 processor
                            if row["Name"][:13] == "%s_SR_2_LAN_" % sat and row["Name"][91:94] == "003":
                                row["length"] = row["Name"][64:68]
                                row["cycle_number"] = row["Name"][69:72]
                                row["orbit_number"] = row["Name"][73:76]
                                writer.writerow(row)
        print(" Catalogue for %s processed!" % sat)

# The list below are the relative number of the orbits of S3A and S3B passing over Greenland.
# They were manually selected by making a request on the dhus platform.
# https://scihub.copernicus.eu/dhus/#/home
ASCENDING_ORBITS = [
    1, 11, 12, 13, 14, 15, 26, 27, 28, 29, 40, 41, 42, 43, 44, 54, 55, 56, 57,
    58, 68, 69, 70, 71, 72, 83, 84, 85, 86, 97, 98, 99, 100, 101, 111, 112, 113,
    114, 115, 125, 126, 127, 128, 129, 140, 141, 142, 143, 154, 155, 156, 157,
    158, 168, 169, 170, 171, 172, 182, 183, 184, 185, 186, 197, 198, 199, 200,
    211, 212, 213, 214, 215, 225, 226, 227, 228, 229, 239, 240, 241, 242, 243,
    254, 255, 256, 257, 268, 269, 270, 271, 272, 282, 283, 284, 285, 286, 297,
    298, 299, 300, 311, 312, 313, 314, 315, 325, 326, 327, 328, 329, 339, 340,
    341, 342, 343, 354, 355, 356, 357, 368, 369, 370, 371, 372, 382, 383, 384,
    385
]

DESCENDING_ORBITS = [
    10, 11, 12, 13, 14, 24, 25, 26, 27, 28, 39, 40, 41, 42, 53, 54, 55, 56, 67,
    68, 69, 70, 71, 81, 82, 83, 84, 85, 96, 97, 98, 99, 110, 111, 112, 113, 124,
    125, 126, 127, 128, 138, 139, 140, 141, 142, 153, 154, 155, 156, 167, 168,
    169, 170, 171, 181, 182, 183, 184, 185, 195, 196, 197, 198, 199, 210, 211,
    212, 213, 224, 225, 226, 227, 228, 238, 239, 240, 241, 242, 252, 253, 254,
    255, 256, 267, 268, 269, 270, 281, 282, 283, 284, 285, 295, 296, 297, 298,
    299, 309, 310, 311, 312, 313, 324, 325, 326, 327, 338, 339, 340, 341, 342,
    352, 353, 354, 355, 356, 367, 368, 369, 370, 381, 382, 383, 384
]

def make_wishlist(sat):
    """Look at the filtered catalogue csv file and find the tracks passing over Greenland."""
    if not path.isdir(WISHLIST_PATH):
        os.mkdir(WISHLIST_PATH)

    df = pd.read_csv(path.join(CATALOG_PATH, "%s_catalogue_filtered.csv" % sat))
    cycle_grouped = df.groupby("cycle_number")
    for cycle, df_cycle in cycle_grouped:
        print("\rProcessing %s cycle %d." % (sat, cycle), end="")
        with open(path.join(WISHLIST_PATH, "%scycle%03d.csv" % (sat, cycle)), "w") as csvfile:
            knames = ["Id", "Name"]
            writer = csv.DictWriter(csvfile, fieldnames=knames)
            writer.writeheader()
            orbit_grouped = df_cycle.groupby("orbit_number")
            for orbit_n, df_orbit in orbit_grouped:
                df_orbit = df_orbit.sort_values("ContentDate:Start")
                # If there is 2 products, then the first is descending and the second is ascending track
                if len(df_orbit.index) == 2:
                    if orbit_n in DESCENDING_ORBITS:
                        writer.writerow({"Id": df_orbit["Id"].values[0], "Name": df_orbit["Name"].values[0]})
                    if orbit_n in ASCENDING_ORBITS:
                        writer.writerow({"Id": df_orbit["Id"].values[1], "Name": df_orbit["Name"].values[1]})
                # If there is another number of products, we don't know which track it is, so we take all
                elif orbit_n in DESCENDING_ORBITS or orbit_n in ASCENDING_ORBITS:
                    for i in range(len(df_orbit.index)):
                        writer.writerow({"Id": df_orbit["Id"].values[i], "Name": df_orbit["Name"].values[i]})
    print(" Wishlists for %s processed!" % sat)

# S3A
if not path.isfile(path.join(CATALOG_PATH, "S3A_catalogue_filtered.csv")):
    process_catalog("S3A")

# S3B
if not path.isfile(path.join(CATALOG_PATH, "S3B_catalogue_filtered.csv")):
    process_catalog("S3B")

make_wishlist("S3A")
make_wishlist("S3B")
