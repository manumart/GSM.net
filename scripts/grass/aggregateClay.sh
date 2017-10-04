#!/bin/bash - 
#===============================================================================
#
#          FILE: aggregateClay.sh
# 
#         USAGE: ./aggregateClay.sh 
# 
#   DESCRIPTION:  == GNU General Public License
#   COPYRIGHT:   ## This file is part of Infosol Digital Soil Mapping (IDSM) Software Copyright(c) 2007, INRA
#       OPTIONS: ---
#  REQUIREMENTS: ---
#          BUGS: ---
#         NOTES: ---
#        AUTHOR: Manuel Martin, manuel.martin@inra.fr
#  ORGANIZATION: INRA
#       CREATED: 04/10/2017 11:03
#      REVISION:  ---
#===============================================================================

set -o nounset                              # Treat unset variables as an error


g.mapset location=GSM90L93 mapset=500m

for wd in 0 5 15
do
for var in lowerpred predicted upperpred 
do
	r.import --overwrite \
	  input=/home/manuel/Y/communs_infosol/Projets/GlobalSoilMap.net/GSM500m/CLAY/clay_l$wd'_'$var'_mean.tif' output=clay_l$wd'_'$var
done
done

#averaged C content in the top 3 layers for 23 cm
r.mapcalc --overwrite  expression='Clay23cm=(clay_l0_predicted*5/23)+\
                                   (clay_l5_predicted*10/23)+(clay_l15_predicted*8/23)'

g.mapset location=GSM90L93 mapset=1km

## aggregates to 1km
r.resamp.stats --o input=Clay23cm@500m output=Clay23cm method=average

