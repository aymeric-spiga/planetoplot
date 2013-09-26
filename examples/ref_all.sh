#! /bin/bash

cd ppclass_reference

./aire.py
./function.py
./hodograph.py
./les_psfc.py
./minimal.py
./simple.py
./zonalmean.py
./anomaly.py
./gw.py
./tide.py
./windspeed.py
./gw_v2.py
./meso_profile.py
./scatter.py
./vector.py
./zonalcontour.py

tar czvf /home/aymeric/Dropbox/Public/demo_data.tar.gz *.ppobj

mv *.ppobj ../../demo_data/
\rm *.png

