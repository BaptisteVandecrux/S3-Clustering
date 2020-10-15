"""Includes small tools used in the project."""
from datetime import datetime
from mpl_toolkits.basemap import Basemap
from matplotlib import pyplot as plt
import pandas
import os
from os import path
import numpy as np


CURRENT_FOLDER = path.dirname(os.path.abspath(__file__))
PREPROC_FOLDER = path.join(CURRENT_FOLDER, "../../data/PREPROC_no_waveform")

def load_data():
    dfs = []

    for f in os.listdir(PREPROC_FOLDER):
        if f[2] == "A":
            print("Loading %s" % f)
            dfs.append(pandas.read_hdf(path.join(PREPROC_FOLDER, f)))
            dfs[-1]["cycle"] = int(f[8:11])

    df = pandas.concat(dfs)
    df["tb_ratio"] = (df["tb_238_01"] - df["tb_365_01"]) / (df["tb_238_01"] + df["tb_365_01"])
    df["tb_avg"] = (df["tb_238_01"] + df["tb_365_01"]) / 2
    df["ku_c_diff"] = df["sig0_ice_sheet_01_ku"] - df["sig0_ice_sheet_01_c"]
    df.reset_index(drop=True, inplace=True)

    return df

def convert_eastward_lon(lon):
    """Convert eastward 0:360 degrees longitude to -180:180 eastward longitude."""
    return lon - 360 * (lon // 180)


TIMESTAMP_OFFSET = datetime.strptime("2000-01-01 00:00:00", "%Y-%m-%d %H:%M:%S").timestamp()

def ESA_time_to_str(ts):
    """Convert ESA timestamp to string."""
    return datetime.fromtimestamp(ts + TIMESTAMP_OFFSET).ctime()


def ESA_time_to_datetime(ts):
    """Convert ESA timestamp to datetime."""
    return datetime.fromtimestamp(ts + TIMESTAMP_OFFSET)

def str_to_ESA_time(datestring):
    """Convert date string to ESA timestamp."""
    return datetime.strptime(datestring, "%Y-%m-%d %H:%M:%S").timestamp() - TIMESTAMP_OFFSET

def drawGreenland(data, lon_field, lat_field, field, field_label):
    """Draws a map of Greenland with the given data or field."""
    if isinstance(field, str):
        vmin = data[field].quantile(q=0.05)
        vmax = data[field].quantile(q=0.95)
        field_val = data[field].to_numpy()
        if field_label=='':
            field_label = field
    else:
        field_val = field
        vmin = np.quantile(field_val, q=0.05)
        vmax = np.quantile(field_val, q=0.95)
    m = Basemap(width=1700000,height=3000000, projection='aea',lat_1=80.,lat_2=40, lon_0=-45, lat_0=71, resolution='l')
    m.drawcoastlines()
    m.fillcontinents()
    m.scatter(
        data[lon_field].to_numpy(),
        data[lat_field].to_numpy(),
        c=field_val,
        vmin=vmin,
        vmax=vmax,
        s=5,
        latlon=True,
        cmap=plt.cm.plasma,
        zorder=10
    )
    cbar = m.colorbar(location='right')
    cbar.ax.set_ylabel(field_label,rotation=270,labelpad=15)
    return m
