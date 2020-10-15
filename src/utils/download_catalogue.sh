# This script downloads the full catalogue for Sentinel-3 A&B

wget --no-check-certificate -r -np -R "index.html*" -nH -P "data" "https://scihub.copernicus.eu/catalogueview/S3A/"
wget --no-check-certificate -r -np -R "index.html*" -nH -P "data" "https://scihub.copernicus.eu/catalogueview/S3B/"
