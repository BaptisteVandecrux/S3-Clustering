#!/bin/sh

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

echo $DIR
# Clément Miège. 2018. Spatial extent of Greenland firn aquifer detected by airborne radars, 2010-2014. Arctic Data Center. doi:10.18739/A2985M.
# https://arcticdata.io/catalog/view/doi:10.18739/A2985M
if [ ! -d "$DIR/../../data" ]
then
  mkdir "$DIR/../../data"
fi

if [ ! -d "$DIR/../../data/firn_aquifer" ]
then
  mkdir "$DIR/../../data/firn_aquifer"
fi

wget --no-check-certificate https://arcticdata.io/metacat/d1/mn/v2/packages/application%2Fbagit-097/resource_map_doi%3A10.18739%2FA2985M -O "$DIR/../../data/firn_aquifer/data.zip"

cd "$DIR/../../data/firn_aquifer"
echo "$DIR/"
unzip data.zip

# Greenland Ice Slabs Data by Michael MacFerrin
# https://figshare.com/articles/Greenland_Ice_Slabs_Data/8309777
if [ ! -d "$DIR/../../data/ice_slabs" ]
then
  mkdir "$DIR/../../data/ice_slabs"
fi

wget --no-check-certificate https://ndownloader.figshare.com/articles/8309777/versions/1 -O "$DIR/../../data/ice_slabs/data.zip"

cd "$DIR/../../data/ice_slabs"
unzip data.zip
