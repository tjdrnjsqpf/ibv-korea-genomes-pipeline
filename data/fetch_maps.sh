#!/bin/bash
# Fetch the two boundary files that Figure 1b needs. They are third-party data under their own
# terms, so they are not redistributed here.
#
#   skorea-municipalities-2013-geo.json  si/gun boundaries, Statistics Korea 2013, via southkorea-maps
#   skorea-provinces-geo.json            province boundaries, GADM (academic use; redistribution is
#                                        not permitted, so download it yourself)
set -euo pipefail
cd "$(dirname "$0")"
curl -fsSL -o skorea-municipalities-2013-geo.json \
  https://raw.githubusercontent.com/southkorea/southkorea-maps/master/kostat/2013/json/skorea_municipalities_geo_simple.json
echo "skorea-municipalities-2013-geo.json fetched"
echo
echo "Province boundaries: download South Korea level-1 GeoJSON from https://gadm.org/download_country.html"
echo "and save it here as skorea-provinces-geo.json (features need the NAME_1 property)."
