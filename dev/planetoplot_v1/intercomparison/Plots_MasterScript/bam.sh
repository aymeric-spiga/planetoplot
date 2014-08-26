#!/bin/bash

########## Constants
mramsfiles=/home/local_home/colaitis/intercomp_test/converted_files
lmdfiles=/planeto/acolmd/SwRI_files/lmd/HOLDEN_nests_mrams
lat=-26. # landing site
lon=-34.   # landing site
#lon=-7.919998  #central longitude to test synchro
#lon=0.
filemrams='_2024-06-17_07:30:00'
filemramsa='_2024-06-16_07:30:00'
filemramsb='_2024-06-18_07:30:00'
filelmd='_2024-06-17_06:00:00'
filelmda='_2024-06-16_06:00:00'
filelmdb='_2024-06-18_06:00:00'
########## ---------
if [ $1 == boundaries ]
then
echo 'Boundaries'
# ------------------
if [ $2 == temp ]
then
for i in 1 2 3
do
/san0/acolmd/SVN/trunk//UTIL/PYTHON/pp.py -f $mramsfiles/mramsout_d0${i}$filemrams --div 20 --axtime lt --time 12 --redope edge_x1 -v tk --title "West-most T boundary MRAMS Nest $i at 12:00" -O mrams_lat1tbound_d0$i -S png -d --res 200 -i 4 -l 0,40,30 -c spectral -m 100 -M 250
/san0/acolmd/SVN/trunk//UTIL/PYTHON/pp.py -f $mramsfiles/mramsout_d0${i}$filemrams --div 20 --axtime lt --time 12 --redope edge_x2 -v tk --title "East-most T boundary MRAMS Nest $i at 12:00" -O mrams_lat2tbound_d0$i -S png -d --res 200 -i 4 -l 0,40,30 -N -c spectral -m 100 -M 250
/san0/acolmd/SVN/trunk//UTIL/PYTHON/pp.py -f $mramsfiles/mramsout_d0${i}$filemrams --div 20 --axtime lt --time 12 --redope edge_y1 -v tk --title "South-most T boundary MRAMS Nest $i at 12:00" -O mrams_lon1tbound_d0$i -S png -d --res 200 -i 4 -l 0,40,30 -N -c spectral -m 100 -M 250
/san0/acolmd/SVN/trunk//UTIL/PYTHON/pp.py -f $mramsfiles/mramsout_d0${i}$filemrams --div 20 --axtime lt --time 12 --redope edge_y2 -v tk --title "North-most T boundary MRAMS Nest $i at 12:00" -O mrams_lon2tbound_d0$i -S png -d --res 200 -i 4 -l 0,40,30 -N -c spectral -m 100 -M 250
/san0/acolmd/SVN/trunk//UTIL/PYTHON/pp.py -f $lmdfiles/wrfout_d0${i}$filelmd --div 20 --axtime lt --time 12 --redope edge_x1 -v tk --title "West-most T boundary LMD_MMM Nest $i at 12:00" -O lmd_lat1tbound_d0$i -S png -d --res 200 -i 4 -l 0,40,30 -c spectral -m 100 -M 250
/san0/acolmd/SVN/trunk//UTIL/PYTHON/pp.py -f $lmdfiles/wrfout_d0${i}$filelmd --div 20 --axtime lt --time 12 --redope edge_x2 -v tk --title "East-most T boundary LMD_MMM Nest $i at 12:00" -O lmd_lat2tbound_d0$i -S png -d --res 200 -i 4 -l 0,40,30 -N -c spectral -m 100 -M 250
/san0/acolmd/SVN/trunk//UTIL/PYTHON/pp.py -f $lmdfiles/wrfout_d0${i}$filelmd --div 20 --axtime lt --time 12 --redope edge_y1 -v tk --title "South-most T boundary LMD_MMM Nest $i at 12:00" -O lmd_lon1tbound_d0$i -S png -d --res 200 -i 4 -l 0,40,30 -N -c spectral -m 100 -M 250
/san0/acolmd/SVN/trunk//UTIL/PYTHON/pp.py -f $lmdfiles/wrfout_d0${i}$filelmd --div 20 --axtime lt --time 12 --redope edge_y2 -v tk --title "North-most T boundary LMD_MMM Nest $i at 12:00" -O lmd_lon2tbound_d0$i -S png -d --res 200 -i 4 -l 0,40,30 -N -c spectral -m 100 -M 250
montage mrams_lat1tbound_d0${i}_200.png mrams_lat2tbound_d0${i}_200.png mrams_lon1tbound_d0${i}_200.png mrams_lon2tbound_d0${i}_200.png lmd_lat1tbound_d0${i}_200.png lmd_lat2tbound_d0${i}_200.png lmd_lon1tbound_d0${i}_200.png lmd_lon2tbound_d0${i}_200.png -mode concatenate -tile 4x2 bounding_d0${i}_temp.png
rm -rf mrams_lat1tbound_d0${i}_200.png mrams_lat2tbound_d0${i}_200.png mrams_lon1tbound_d0${i}_200.png mrams_lon2tbound_d0${i}_200.png lmd_lat1tbound_d0${i}_200.png lmd_lat2tbound_d0${i}_200.png lmd_lon1tbound_d0${i}_200.png lmd_lon2tbound_d0${i}_200.png mrams_lat1tbound_d0${i}.sh mrams_lat2tbound_d0${i}.sh mrams_lon1tbound_d0${i}.sh mrams_lon2tbound_d0${i}.sh lmd_lat1tbound_d0${i}.sh lmd_lat2tbound_d0${i}.sh lmd_lon1tbound_d0${i}.sh lmd_lon2tbound_d0${i}.sh
done
fi
# ------------------
if [ $2 == uv ]
then
for i in 1 2 3
do
/san0/acolmd/SVN/trunk//UTIL/PYTHON/pp.py -f $mramsfiles/mramsout_d0${i}$filemrams --div 20 --axtime lt --time 12 --redope edge_x1 -v uvmet --title "West-most UV boundary MRAMS Nest $i at 12:00" -O mrams_lat1uvbound_d0$i -S png -d --res 200 -i 4 -l 0,40,30 -c spectral -m 0 -M 50
/san0/acolmd/SVN/trunk//UTIL/PYTHON/pp.py -f $mramsfiles/mramsout_d0${i}$filemrams --div 20 --axtime lt --time 12 --redope edge_x2 -v uvmet --title "East-most UV boundary MRAMS Nest $i at 12:00" -O mrams_lat2uvbound_d0$i -S png -d --res 200 -i 4 -l 0,40,30 -N -c spectral -m 0 -M 50
/san0/acolmd/SVN/trunk//UTIL/PYTHON/pp.py -f $mramsfiles/mramsout_d0${i}$filemrams --div 20 --axtime lt --time 12 --redope edge_y1 -v uvmet --title "South-most UV boundary MRAMS Nest $i at 12:00" -O mrams_lon1uvbound_d0$i -S png -d --res 200 -i 4 -l 0,40,30 -N -c spectral -m 0 -M 50
/san0/acolmd/SVN/trunk//UTIL/PYTHON/pp.py -f $mramsfiles/mramsout_d0${i}$filemrams --div 20 --axtime lt --time 12 --redope edge_y2 -v uvmet --title "North-most UV boundary MRAMS Nest $i at 12:00" -O mrams_lon2uvbound_d0$i -S png -d --res 200 -i 4 -l 0,40,30 -N -c spectral -m 0 -M 50
/san0/acolmd/SVN/trunk//UTIL/PYTHON/pp.py -f $lmdfiles/wrfout_d0${i}$filelmd --div 20 --axtime lt --time 12 --redope edge_x1 -v uvmet --title "West-most UV boundary LMD_MMM Nest $i at 12:00" -O lmd_lat1uvbound_d0$i -S png -d --res 200 -i 4 -l 0,40,30 -c spectral -m 0 -M 50
/san0/acolmd/SVN/trunk//UTIL/PYTHON/pp.py -f $lmdfiles/wrfout_d0${i}$filelmd --div 20 --axtime lt --time 12 --redope edge_x2 -v uvmet --title "East-most UV boundary LMD_MMM Nest $i at 12:00" -O lmd_lat2uvbound_d0$i -S png -d --res 200 -i 4 -l 0,40,30 -N -c spectral -m 0 -M 50
/san0/acolmd/SVN/trunk//UTIL/PYTHON/pp.py -f $lmdfiles/wrfout_d0${i}$filelmd --div 20 --axtime lt --time 12 --redope edge_y1 -v uvmet --title "South-most UV boundary LMD_MMM Nest $i at 12:00" -O lmd_lon1uvbound_d0$i -S png -d --res 200 -i 4 -l 0,40,30 -N -c spectral -m 0 -M 50
/san0/acolmd/SVN/trunk//UTIL/PYTHON/pp.py -f $lmdfiles/wrfout_d0${i}$filelmd --div 20 --axtime lt --time 12 --redope edge_y2 -v uvmet --title "North-most UV boundary LMD_MMM Nest $i at 12:00" -O lmd_lon2uvbound_d0$i -S png -d --res 200 -i 4 -l 0,40,30 -N -c spectral -m 0 -M 50
montage mrams_lat1uvbound_d0${i}_200.png mrams_lat2uvbound_d0${i}_200.png mrams_lon1uvbound_d0${i}_200.png mrams_lon2uvbound_d0${i}_200.png lmd_lat1uvbound_d0${i}_200.png lmd_lat2uvbound_d0${i}_200.png lmd_lon1uvbound_d0${i}_200.png lmd_lon2uvbound_d0${i}_200.png -mode concatenate -tile 4x2 bounding_d0${i}_uv.png
rm -rf mrams_lat1uvbound_d0${i}_200.png mrams_lat2uvbound_d0${i}_200.png mrams_lon1uvbound_d0${i}_200.png mrams_lon2uvbound_d0${i}_200.png lmd_lat1uvbound_d0${i}_200.png lmd_lat2uvbound_d0${i}_200.png lmd_lon1uvbound_d0${i}_200.png lmd_lon2uvbound_d0${i}_200.png mrams_lat1uvbound_d0${i}.sh mrams_lat2uvbound_d0${i}.sh mrams_lon1uvbound_d0${i}.sh mrams_lon2uvbound_d0${i}.sh lmd_lat1uvbound_d0${i}.sh lmd_lat2uvbound_d0${i}.sh lmd_lon1uvbound_d0${i}.sh lmd_lon2uvbound_d0${i}.sh
done
fi
# ------------------
fi
########## ---------
if [ $1 == sync ]
then
echo 'Synchronization'
for i in [1,2,3]
do
/san0/acolmd/SVN/trunk//UTIL/PYTHON/pp.py -f $mramsfiles/mramsout_d0${i}$filemrams --div 20 --axtime lt --time 12 -v tsurf --title "Tsurf MRAMS Nest $i at 12:00" -O mrams_tsurf_d0$i -S png -d --res 200 ---proj moll --mark $lon,$lat
/san0/acolmd/SVN/trunk//UTIL/PYTHON/pp.py -f $lmdfiles/wrfout_d0${i}$filelmd --div 20 --axtime lt --time 12 -v tsurf --title "Tsurf LMD_MMM Nest $i at 12:00" -O lmd_tsurf_d0$i -S png -d --res 200 --proj moll --mark $lon,$lat
done
montage mrams_tsurf_d01_200.png mrams_tsurf_d02_200.png mrams_tsurf_d03_200.png lmd_tsurf_d01_200.png lmd_tsurf_d02_200.png lmd_tsurf_d03_200.png tsurf_sync.png
rm -f mrams_tsurf_d01_200.png mrams_tsurf_d02_200.png mrams_tsurf_d03_200.png lmd_tsurf_d01_200.png lmd_tsurf_d02_200.png lmd_tsurf_d03_200.png mrams_tsurf_d01_200.sh mrams_tsurf_d02_200.sh mrams_tsurf_d03_200.sh lmd_tsurf_d01_200.sh lmd_tsurf_d02_200.sh lmd_tsurf_d03_200.sh

/san0/acolmd/SVN/trunk//UTIL/PYTHON/pp.py -f $mramsfiles/mramsout_d03$filemrams,$lmdfiles/wrfout_d03$filelmd  --lon $lon --lat $lat --axtime lt -v tsurf --ylabel "Surface Temperature" --xlabel "Local time (h)" --title "NEST 3 SYNCHRONIZATION" -O sync_tsurfnest3 -S png -d --res 200 --labels "MRAMS","LMD_MMM" 
rm -f sync_tsurfnest3.sh

fi
########## ---------

if [ $1 == nests ]
then
echo 'Safe Reference Landing Site'

if [ $2 == power ]
then
echo 'Power spectrum'
/san0/acolmd/SVN/trunk//UTIL/PYTHON/pp.py -f $mramsfiles/mramsout_d01$filemrams,$mramsfiles/mramsout_d02$filemrams,$mramsfiles/mramsout_d03$filemrams --div 20  --lon $lon --lat $lat --axtime lt --time 8,20 0 -i 4 -l 0.1,10,999 -v uv --analysis fft --logy --logx --xlabel "Wavelength (m)" --ylabel "Spectrum amplitude" --title "MRAMS Nest $i Daytime Power spectrum of horizontal winds at landing site (0-10km AGL)" -O mrams_power_nests -S png -d --res 200 --xmin 10 --xmax 10000 --ymin 1 --ymax 10000 --labels "MRAMS Nest 1","MRAMS Nest 2","MRAMS Nest 3"
/san0/acolmd/SVN/trunk//UTIL/PYTHON/pp.py -f $lmdfiles/wrfout_d01$filelmd,$lmdfiles/wrfout_d02$filelmd,$lmdfiles/wrfout_d03$filelmd --div 20  --lon $lon --lat $lat --axtime lt --time 8,20 0 -i 4 -l 0.1,10,999 -v uv --analysis fft --logy --logx --xlabel "Wavelength (m)" --ylabel "Spectrum amplitude" --title "LMD_MMM Nest $i Daytime Power spectrum of horizontal winds at landing site (0-10km AGL)" -O lmd_power_nests -S png -d --res 200 --xmin 10 --xmax 10000 --ymin 1 --ymax 10000 --labels "LMD_MMM Nest 1","LMD_MMM Nest 2","LMD_MMM Nest 3"
rm -f mrams_power_nests.sh lmd_power_nests.sh
fi

if [ $2 == histo ]
then
echo 'Wind histograms'
/san0/acolmd/SVN/trunk//UTIL/PYTHON/pp.py -f $mramsfiles/mramsout_d01$filemrams,$mramsfiles/mramsout_d02$filemrams,$mramsfiles/mramsout_d03$filemrams -v uv --analysis histo --lon $lon --lat $lat --div 20 --labels "MRAMS Nest 1","MRAMS Nest 2","MRAMS Nest 3" --title "MRAMS horizontal wind histogram at landing site for all time and altitudes" -d -S png -O mrams_histo2d_nests
/san0/acolmd/SVN/trunk//UTIL/PYTHON/pp.py -f $lmdfiles/wrfout_d01$filelmd,$lmdfiles/wrfout_d02$filelmd,$lmdfiles/wrfout_d03$filelmd -v uv --analysis histo --lon $lon --lat $lat --div 20 --labels "LMD_MMM Nest 1","LMD_MMM Nest 2","LMD_MMM Nest 3" --title "LMD_MMM horizontal wind histogram at landing site for all time and altitudes" -d -S png -O lmd_histo2d_nests
rm -f mrams_histo2d_nests.sh lmd_histo2d_nests.sh
fi
fi

if [ $1 == srlsn3 ]
then
echo 'Safe Reference Landing Site inner nest'

if [ $2 == histo ]
then

echo 'WARNING'
echo 'WARNING'
echo 'WARNING'
echo 'THOSE HIST ARE FOR ALL TIMES, CHECK THAT SYNCHRONIZATION IS OK OR SPECIFY COMMON LOCAL TIMES'
echo 'THOSE HIST ARE FOR ALL TIMES, CHECK THAT SYNCHRONIZATION IS OK OR SPECIFY COMMON LOCAL TIMES'
echo 'THOSE HIST ARE FOR ALL TIMES, CHECK THAT SYNCHRONIZATION IS OK OR SPECIFY COMMON LOCAL TIMES'

/san0/acolmd/SVN/trunk//UTIL/PYTHON/pp.py -f $mramsfiles/mramsout_d03$filemrams,$lmdfiles/wrfout_d03$filelmd -v uv --analysis histo --lon $lon --lat $lat --div 20 --labels "MRAMS","LMD_MMM" --title "Inner-Nest Horizontal wind histogram at landing site for all time and altitudes" -d -S png -O comp_histo2d_srlsn3 --ylabel "Density of probability" --xlabel "Horizontal wind speed (m/s)"
/san0/acolmd/SVN/trunk//UTIL/PYTHON/pp.py -f $mramsfiles/mramsout_d03$filemrams,$lmdfiles/wrfout_d03$filelmd -v uv --analysis histo --lon $lon --lat $lat --div 20 --labels "MRAMS","LMD_MMM" --title "Hor. wind histogram / 10m-600m AGL range: POWERED DESCENT" -d -S png -O comp_histo1d_10m_600m_AGL_srlsn3 --vert 0,10 -i 4 -l 0.01,0.6,20 --ylabel "Density of probability" --xlabel "Horizontal wind speed (m/s)"
/san0/acolmd/SVN/trunk//UTIL/PYTHON/pp.py -f $mramsfiles/mramsout_d03$filemrams,$lmdfiles/wrfout_d03$filelmd -v uv --analysis histo --lon $lon --lat $lat --div 20 --labels "MRAMS","LMD_MMM" --title "Hor. wind histogram / 1km-4km AML range: MAIN PARACHUTE DESCENT" -d -S png -O comp_histo1d_1km_4km_AML_srlsn3 --vert 0,20 -i 3 -l 1,4,20 --ylabel "Density of probability" --xlabel "Horizontal wind speed (m/s)"
/san0/acolmd/SVN/trunk//UTIL/PYTHON/pp.py -f $mramsfiles/mramsout_d03$filemrams,$lmdfiles/wrfout_d03$filelmd -v uv --analysis histo --lon $lon --lat $lat --div 20 --labels "MRAMS","LMD_MMM" --title "Hor. wind histogram / 4km-10km AML range: DROGUE AND PARACHUTE DEPLOYMENT" -d -S png -O comp_histo1d_4km_10km_AML_srlsn3 --vert 0,20 -i 3 -l 4,10,20 --ylabel "Density of probability" --xlabel "Horizontal wind speed (m/s)"
/san0/acolmd/SVN/trunk//UTIL/PYTHON/pp.py -f $mramsfiles/mramsout_d03$filemrams,$lmdfiles/wrfout_d03$filelmd -v uv --analysis histo --lon $lon --lat $lat --div 20 --labels "MRAMS","LMD_MMM" --title "Hor. wind histogram / 10km-40km AML range: DECELERATION PHASE" -d -S png -O comp_histo1d_10km_40km_AML_srlsn3 --vert 0,20 -i 3 -l 10,40,20 --ylabel "Density of probability" --xlabel "Horizontal wind speed (m/s)"
montage comp_histo1d_10m_600m_AGL_srlsn3_200.png comp_histo1d_1km_4km_AML_srlsn3_200.png comp_histo1d_4km_10km_AML_srlsn3_200.png comp_histo1d_10km_40km_AML_srlsn3_200.png -mode concatenate -tile 2x2 comp_histo1d_decomposed_srlsn3.png
rm -f comp_histo2d_srlsn3.sh comp_histo1d_10m_600m_AGL_srlsn3.sh comp_histo1d_1km_4km_AML_srlsn3.sh comp_histo1d_4km_10km_AML_srlsn3.sh comp_histo1d_10km_40km_AML_srlsn3.sh comp_histo1d_10m_600m_AGL_srlsn3_200.png comp_histo1d_1km_4km_AML_srlsn3_200.png comp_histo1d_4km_10km_AML_srlsn3_200.png comp_histo1d_10km_40km_AML_srlsn3_200.png

#/san0/acolmd/SVN/trunk//UTIL/PYTHON/pp.py -f $mramsfiles/mramsout_d03$filemrams,$lmdfiles/wrfout_d03$filelmd -v tsurf --lon $lon --lat $lat --div 20 --labels "MRAMS","LMD_MMM" --title "Inner-Nest Surface temperature histogram" -d -S png -O comp_histo1d_tsurf_srlsn3
#rm -f comp_histo1d_tsurf_srlsn3.sh
fi

if [ $2 == power ]
then
/san0/acolmd/SVN/trunk//UTIL/PYTHON/pp.py -f $mramsfiles/mramsout_d03$filemrams,$lmdfiles/wrfout_d03$filelmd  --lon $lon --lat $lat --axtime lt --time 8,20 0 -i 4 -l 0.1,10,999 -v uv --analysis fft --logy --logx --xlabel "Wavelength (m)" --ylabel "Spectrum amplitude" --title "Inner-Nest Daytime Power spectrum of horizontal winds at landing site (0-10km AGL)" -O comp_power_srlsn3 -S png -d --res 200 --xmin 10 --xmax 10000 --ymin 1 --ymax 10000 --labels "MRAMS","LMD_MMM" 
rm -f comp_power_srlsn3.sh
fi

if [ $2 == series ]
then

/san0/acolmd/SVN/trunk//UTIL/PYTHON/pp.py -f $mramsfiles/mramsout_d03$filemrams,$lmdfiles/wrfout_d03$filelmd -v UV --vert 0 --lon $lon --lat $lat --div 20 --labels "MRAMS","LMD_MMM" --title "Horizontal winds at 10m AGL (final drop and touchdown)" -d -S png -O comp_uv10m_serie_srlsn3 -i 4 -l 0.01 --xlabel "Local time (h)" --axtime lt --ylabel "Horizontal wind speed (m/s)"
rm -f comp_uv10m_serie_srlsn3.sh

/san0/acolmd/SVN/trunk//UTIL/PYTHON/pp.py -f $mramsfiles/mramsout_d03$filemrams,$lmdfiles/wrfout_d03$filelmd -v tk --vert 0 --lon $lon --lat $lat --div 20 --labels "MRAMS","LMD_MMM" --title "Temperature at 10m AGL (final drop and touchdown)" -d -S png -O comp_temp10m_serie_srlsn3 -i 4 -l 0.01 --xlabel "Local time (h)" --axtime lt --ylabel "Temperature (K)"
rm -f comp_temp10m_serie_srlsn3.sh

/san0/acolmd/SVN/trunk//UTIL/PYTHON/pp.py -f $mramsfiles/mramsout_d03$filemrams,$lmdfiles/wrfout_d03$filelmd -v UV --lon $lon --lat $lat --div 20 --title "MRAMS","LMD_MMM" -d -S png -O comp_uv_serie2d_srlsn3 -i 4 -l 0.01,40,40 --xlabel "Local time (h)" --axtime lt --ylabel "Altitude above surface (km)" --xmin 5 --xmax 24 -m 0 -M 75
rm -f comp_uv_serie2d_srlsn3.sh

/san0/acolmd/SVN/trunk//UTIL/PYTHON/pp.py -f $mramsfiles/mramsout_d03$filemrams,$lmdfiles/wrfout_d03$filelmd -v tk --lon $lon --lat $lat --div 20 --title "MRAMS","LMD_MMM" -d -S png -O comp_temp_serie2d_srlsn3 -i 4 -l 0.01,40,40 --xlabel "Local time (h)" --axtime lt --ylabel "Altitude above surface (km)" -m 105 -M 250 -c spectral --xmin 5 --xmax 24
rm -f comp_temp_serie2d_srlsn3.sh

/san0/acolmd/SVN/trunk//UTIL/PYTHON/pp.py -f $mramsfiles/mramsout_d03$filemrams,$lmdfiles/wrfout_d03$filelmd -v hodograph --vert 0 --lon $lon --lat $lat --div 20 --labels "MRAMS","LMD_MMM" --title "Hodograph of horizontal winds at 10m AGL (final drop and touchdown)" -d -S png -O comp_hodograph_srlsn3 -i 4 -l 0.01 --xlabel "U component (m/s)" --axtime lt --ylabel "V component (m/s)" --xmin -15 --xmax 15 --ymin -15 --ymax 15
rm -f comp_hodograph_srlsn3.sh

/san0/acolmd/SVN/trunk//UTIL/PYTHON/pp.py -f $mramsfiles/mramsout_d03$filemrams,$lmdfiles/wrfout_d03$filelmd -v hodograph_2 --vert 0 --lon $lon --lat $lat --div 20 --labels "MRAMS","LMD_MMM" --title "Horiztonal winds rotation at 10m AGL (final drop and touchdown)" -d -S png -O comp_hodograph_2_srlsn3 -i 4 -l 0.01 --xlabel "Local Time (h)" --axtime lt --ylabel "Arbitrary wind rotation (degrees)"
rm -f comp_hodograph_2_srlsn3.sh

fi

if [ $2 == profiles ]
then

/san0/acolmd/SVN/trunk//UTIL/PYTHON/pp.py -f $mramsfiles/mramsout_d03$filemrams,$lmdfiles/wrfout_d03$filelmd -v UV --lon $lon --lat $lat --div 20 --labels "MRAMS","LMD_MMM" --title "Horizontal winds at LT 13:00" -d -S png -O comp_uv_profile_lt13_srlsn3 -i 4 -l 0.01,40,40 --xlabel "Horizontal wind speed (m/s)" --axtime lt --ylabel "Altitude above surface (km)" --axtime lt --time 13
rm -f comp_uv_profile_lt13_srlsn3.sh

/san0/acolmd/SVN/trunk//UTIL/PYTHON/pp.py -f $mramsfiles/mramsout_d03$filemrams,$lmdfiles/wrfout_d03$filelmd -v tk --lon $lon --lat $lat --div 20 --labels "MRAMS","LMD_MMM" --title "Temperature serie at LT 13:00" -d -S png -O comp_temp_profile_lt13_srlsn3 -i 4 -l 0.01,40,40 --xlabel "Temperature (K)" --axtime lt --ylabel "Altitude above surface (km)" --axtime lt --time 13
rm -f comp_temp_profile_lt13_srlsn3.sh

fi

if [ $2 == surface ]
then

/san0/acolmd/SVN/trunk//UTIL/PYTHON/pp.py -f $mramsfiles/mramsout_d03$filemrams,$lmdfiles/wrfout_d03$filelmd -v PSFC --axtime lt --title "MRAMS Inner Nest surface pressure","LMD_MMM Inner Nest surface pressure" -d -S png -O comp_pressure_srlsn3 --blat -20 --div 20 --res 200 --lon $lon --lat $lat --xlabel "Local Time (h)" --ylabel "Surface Pressure (Pa)"
rm -f comp_pressure_srlsn3.sh
/san0/acolmd/SVN/trunk//UTIL/PYTHON/pp.py -f $mramsfiles/mramsout_d03$filemrams,$lmdfiles/wrfout_d03$filelmd -v CDU --axtime lt --title "MRAMS Inner Nest surface stress","LMD_MMM Inner Nest surface stress" -d -S png -O comp_surfstress_srlsn3 --blat -20 --div 20 --res 200 --lon $lon --lat $lat --xlabel "Local Time (h)" --ylabel "Surface stress (m/s)"
rm -f comp_surfstress_srlsn3.sh
/san0/acolmd/SVN/trunk//UTIL/PYTHON/pp.py -f $mramsfiles/mramsout_d03$filemrams,$lmdfiles/wrfout_d03$filelmd -v HFX --axtime lt --title "MRAMS Inner Nest sensib. heat flux ","LMD_MMM Inner Nest sensib. heat flux" -d -S png -O comp_sensibflux_srlsn3 --blat -20 --div 20 --res 200 --lon $lon --lat $lat --xlabel "Local Time (h)" --ylabel "Sensible heat flux (K.m/s)"
rm -f comp_sensibflux_srlsn3.sh

fi

fi

########## ---------

if [ $1 == maps ]
then

if [ $2 == winds ]
then

/san0/acolmd/SVN/trunk//UTIL/PYTHON/pp.py -f $mramsfiles/mramsout_d03$filemrams -v UV --axtime lt --time 14 -i 4 -l 0.01 --title "MRAMS horizontal winds at 10m AGL" -d -S png -W --vert 0 -O mrams_winds_10mAGL_maps --proj spstere --blat -20 --div 20 --res 200 -s 5 --mark $lon,$lat
/san0/acolmd/SVN/trunk//UTIL/PYTHON/pp.py -f $mramsfiles/mramsout_d03$filemrams  -v UV --axtime lt --time 14 -i 4 -l 0.6 --title "MRAMS horizontal winds at 600m AGL" -d -S png -W --vert 0 -O mrams_winds_600mAGL_maps --proj spstere --blat -20 --div 20 --res 200 -s 5 --mark $lon,$lat
/san0/acolmd/SVN/trunk//UTIL/PYTHON/pp.py -f $mramsfiles/mramsout_d03$filemrams  -v UV --axtime lt --time 14 -i 3 -l 1 --title "MRAMS horizontal winds at 1km AML" -d -S png -W --vert 0 -O mrams_winds_1kmAML_maps --proj spstere --blat -20 --div 20 --res 200 -s 5 --mark $lon,$lat
/san0/acolmd/SVN/trunk//UTIL/PYTHON/pp.py -f $mramsfiles/mramsout_d03$filemrams  -v UV --axtime lt --time 14 -i 3 -l 4 --title "MRAMS horizontal winds at 4km AML" -d -S png -W --vert 0 -O mrams_winds_4kmAML_maps --proj spstere --blat -20 --div 20 --res 200 -s 5 --mark $lon,$lat
/san0/acolmd/SVN/trunk//UTIL/PYTHON/pp.py -f $mramsfiles/mramsout_d03$filemrams -v UV --axtime lt --time 14 -i 3 -l 10 --title "MRAMS horizontal winds at 10km AML" -d -S png -W --vert 0 -O mrams_winds_10kmAML_maps --proj spstere --blat -20 --div 20 --res 200 -s 5 --mark $lon,$lat
/san0/acolmd/SVN/trunk//UTIL/PYTHON/pp.py -f $mramsfiles/mramsout_d03$filemrams  -v UV --axtime lt --time 14 -i 3 -l 20 --title "MRAMS horizontal winds at 20km AML" -d -S png -W --vert 0 -O mrams_winds_20kmAML_maps --proj spstere --blat -20 --div 20 --res 200 -s 5 --mark $lon,$lat

/san0/acolmd/SVN/trunk//UTIL/PYTHON/pp.py -f $lmdfiles/wrfout_d03$filelmd -v UV --axtime lt --time 14 -i 4 -l 0.01 --title "LMD_MMM horizontal winds at 10m AGL" -d -S png -W --vert 0 -O lmd_winds_10mAGL_maps --proj spstere --blat -20 --div 20 --res 200 -s 5 --mark $lon,$lat
/san0/acolmd/SVN/trunk//UTIL/PYTHON/pp.py -f $lmdfiles/wrfout_d03$filelmd  -v UV --axtime lt --time 14 -i 4 -l 0.6 --title "LMD_MMM horizontal winds at 600m AGL" -d -S png -W --vert 0 -O lmd_winds_600mAGL_maps --proj spstere --blat -20 --div 20 --res 200 -s 5 --mark $lon,$lat
/san0/acolmd/SVN/trunk//UTIL/PYTHON/pp.py -f $lmdfiles/wrfout_d03$filelmd  -v UV --axtime lt --time 14 -i 3 -l 1 --title "LMD_MMM horizontal winds at 1km AML" -d -S png -W --vert 0 -O lmd_winds_1kmAML_maps --proj spstere --blat -20 --div 20 --res 200 -s 5 --mark $lon,$lat
/san0/acolmd/SVN/trunk//UTIL/PYTHON/pp.py -f $lmdfiles/wrfout_d03$filelmd  -v UV --axtime lt --time 14 -i 3 -l 4 --title "LMD_MMM horizontal winds at 4km AML" -d -S png -W --vert 0 -O lmd_winds_4kmAML_maps --proj spstere --blat -20 --div 20 --res 200 -s 5 --mark $lon,$lat
/san0/acolmd/SVN/trunk//UTIL/PYTHON/pp.py -f $lmdfiles/wrfout_d03$filelmd -v UV --axtime lt --time 14 -i 3 -l 10 --title "LMD_MMM horizontal winds at 10km AML" -d -S png -W --vert 0 -O lmd_winds_10kmAML_maps --proj spstere --blat -20 --div 20 --res 200 -s 5 --mark $lon,$lat
/san0/acolmd/SVN/trunk//UTIL/PYTHON/pp.py -f $lmdfiles/wrfout_d03$filelmd  -v UV --axtime lt --time 14 -i 3 -l 20 --title "LMD_MMM horizontal winds at 20km AML" -d -S png -W --vert 0 -O lmd_winds_20kmAML_maps --proj spstere --blat -20 --div 20 --res 200 -s 5 --mark $lon,$lat

montage mrams_winds_10mAGL_maps_200.png lmd_winds_10mAGL_maps_200.png mrams_winds_600mAGL_maps_200.png lmd_winds_600mAGL_maps_200.png mrams_winds_1kmAML_maps_200.png lmd_winds_1kmAML_maps_200.png -mode concatenate -tile 2x3 comp_winds_maps_10m_1km.png
montage mrams_winds_4kmAML_maps_200.png lmd_winds_4kmAML_maps_200.png mrams_winds_10kmAML_maps_200.png lmd_winds_10kmAML_maps_200.png mrams_winds_20kmAML_maps_200.png lmd_winds_20kmAML_maps_200.png -mode concatenate -tile 2x3 comp_winds_maps_4km_20km.png
rm -f mrams_winds_10mAGL_maps.sh lmd_winds_10mAGL_maps_200.png mrams_winds_600mAGL_maps_200.png lmd_winds_600mAGL_maps_200.png mrams_winds_1kmAML_maps_200.png lmd_winds_1kmAML_maps_200.png mrams_winds_4kmAML_maps_200.png lmd_winds_4kmAML_maps_200.png mrams_winds_10kmAML_maps_200.png lmd_winds_10kmAML_maps_200.png mrams_winds_20kmAML_maps_200.png lmd_winds_20kmAML_maps_200.png mrams_winds_10mAGL_maps.sh lmd_winds_10mAGL_maps.sh mrams_winds_600mAGL_maps.sh lmd_winds_600mAGL_maps.sh mrams_winds_1kmAML_maps.sh lmd_winds_1kmAML_maps.sh mrams_winds_4kmAML_maps.sh lmd_winds_4kmAML_maps.sh mrams_winds_10kmAML_maps.sh lmd_winds_10kmAML_maps.sh mrams_winds_20kmAML_maps.sh lmd_winds_20kmAML_maps.sh 

fi

if [ $2 == surface ]
then
#/san0/acolmd/SVN/trunk//UTIL/PYTHON/pp.py -f $mramsfiles/mramsout_d03$filemrams,$lmdfiles/wrfout_d03$filelmd -v PSFC --axtime lt --time 14 --title "MRAMS Inner Nest surface pressure","LMD_MMM Inner Nest surface pressure" -d -S png -O comp_pressure_maps --proj spstere --blat -20 --div 20 --res 200 --mark $lon,$lat
#rm -f comp_pressure_maps.sh
#/san0/acolmd/SVN/trunk//UTIL/PYTHON/pp.py -f $mramsfiles/mramsout_d03$filemrams,$lmdfiles/wrfout_d03$filelmd -v CDU --axtime lt --time 14 --title "MRAMS Inner Nest surface stress","LMD_MMM Inner Nest surface stress" -d -S png -O comp_surfstress_maps --proj spstere --blat -20 --div 20 --res 200 --mark $lon,$lat
#rm -f comp_surfstress_maps.sh
#/san0/acolmd/SVN/trunk//UTIL/PYTHON/pp.py -f $mramsfiles/mramsout_d03$filemrams,$lmdfiles/wrfout_d03$filelmd -v HFX --axtime lt --time 14 --title "MRAMS Inner Nest sensib. heat flux ","LMD_MMM Inner Nest sensib. heat flux" -d -S png -O comp_sensibflux_maps --proj spstere --blat -20 --div 20 --res 200 --mark $lon,$lat
#rm -f comp_sensibflux_maps.sh

/san0/acolmd/SVN/trunk//UTIL/PYTHON/pp.py -f $lmdfiles/wrfout_d03$filelmd -v HFX --axtime lt --time 14 --title "MRAMS Inner Nest sensib. heat flux ","LMD_MMM Inner Nest sensib. heat flux" --proj spstere --blat -20 --div 20 --mark $lon,$lat

fi

fi

########## ---------

if [ $1 == variability ]
then

if [ $2 == histo ]
then
/san0/acolmd/SVN/trunk//UTIL/PYTHON/pp.py -f $mramsfiles/mramsout_d03$filemramsa,$mramsfiles/mramsout_d03$filemrams,$mramsfiles/mramsout_d03$filemramsb -v uv --analysis histo --lon $lon --lat $lat --div 20 --labels "MRAMS day-1","MRAMS day 0","MRAMS day+1" --title "MRAMS Inner-Nest Horizontal wind histogram at landing site for all time and altitudes" -d -S png -O mrams_histo2d_variability --ylabel "Density of probability" --xlabel "Horizontal wind speed (m/s)" --ymin 0 --ymax 0.05
/san0/acolmd/SVN/trunk//UTIL/PYTHON/pp.py -f $mramsfiles/mramsout_d03$filemramsa,$mramsfiles/mramsout_d03$filemrams,$mramsfiles/mramsout_d03$filemramsb -v uv --analysis histo --lon $lon --lat $lat --div 20 --labels "MRAMS day-1","MRAMS day 0","MRAMS day+1" --title "Hor. wind histogram / 10m-600m AGL range: POWERED DESCENT" -d -S png -O mrams_histo1d_10m_600m_AGL_variability --vert 0,10 -i 4 -l 0.01,0.6,20 --ylabel "Density of probability" --xlabel "Horizontal wind speed (m/s)"
/san0/acolmd/SVN/trunk//UTIL/PYTHON/pp.py -f $mramsfiles/mramsout_d03$filemramsa,$mramsfiles/mramsout_d03$filemrams,$mramsfiles/mramsout_d03$filemramsb -v uv --analysis histo --lon $lon --lat $lat --div 20 --labels "MRAMS day-1","MRAMS day 0","MRAMS day+1" --title "Hor. wind histogram / 1km-4km AML range: MAIN PARACHUTE DESCENT" -d -S png -O mrams_histo1d_1km_4km_AML_variability --vert 0,20 -i 3 -l 1,4,20 --ylabel "Density of probability" --xlabel "Horizontal wind speed (m/s)"
/san0/acolmd/SVN/trunk//UTIL/PYTHON/pp.py -f $mramsfiles/mramsout_d03$filemramsa,$mramsfiles/mramsout_d03$filemrams,$mramsfiles/mramsout_d03$filemramsb -v uv --analysis histo --lon $lon --lat $lat --div 20 --labels "MRAMS day-1","MRAMS day 0","MRAMS day+1" --title "Hor. wind histogram / 4km-10km AML range: DROGUE AND PARACHUTE DEPLOYMENT" -d -S png -O mrams_histo1d_4km_10km_AML_variability --vert 0,20 -i 3 -l 4,10,20 --ylabel "Density of probability" --xlabel "Horizontal wind speed (m/s)"
/san0/acolmd/SVN/trunk//UTIL/PYTHON/pp.py -f $mramsfiles/mramsout_d03$filemramsa,$mramsfiles/mramsout_d03$filemrams,$mramsfiles/mramsout_d03$filemramsb -v uv --analysis histo --lon $lon --lat $lat --div 20 --labels "MRAMS day-1","MRAMS day 0","MRAMS day+1" --title "Hor. wind histogram / 10km-40km AML range: DECELERATION PHASE" -d -S png -O mrams_histo1d_10km_40km_AML_variability --vert 0,20 -i 3 -l 10,40,20 --ylabel "Density of probability" --xlabel "Horizontal wind speed (m/s)"
montage mrams_histo1d_10m_600m_AGL_variability_200.png mrams_histo1d_1km_4km_AML_variability_200.png mrams_histo1d_4km_10km_AML_variability_200.png mrams_histo1d_10km_40km_AML_variability_200.png -mode concatenate -tile 2x2 mrams_histo1d_decomposed_variability.png
rm -f mrams_histo1d_10m_600m_AGL_variability_200.png mrams_histo1d_1km_4km_AML_variability_200.png mrams_histo1d_4km_10km_AML_variability_200.png mrams_histo1d_10km_40km_AML_variability_200.png mrams_histo1d_10m_600m_AGL_variability.sh mrams_histo1d_1km_4km_AML_variability.sh mrams_histo1d_4km_10km_AML_variability.sh mrams_histo1d_10km_40km_AML_variability.sh

/san0/acolmd/SVN/trunk//UTIL/PYTHON/pp.py -f $lmdfiles/wrfout_d03$filelmda,$lmdfiles/wrfout_d03$filelmd,$lmdfiles/wrfout_d03$filelmdb -v uv --analysis histo --lon $lon --lat $lat --div 20 --labels "LMD_MMM day-1","LMD_MMM day 0","LMD_MMM day+1" --title "LMD_MMM Inner-Nest Horizontal wind histogram at landing site for all time and altitudes" -d -S png -O lmd_histo2d_variability --ylabel "Density of probability" --xlabel "Horizontal wind speed (m/s)" --ymin 0 --ymax 0.05
/san0/acolmd/SVN/trunk//UTIL/PYTHON/pp.py -f $lmdfiles/wrfout_d03$filelmda,$lmdfiles/wrfout_d03$filelmd,$lmdfiles/wrfout_d03$filelmdb -v uv --analysis histo --lon $lon --lat $lat --div 20 --labels "LMD_MMM day-1","LMD_MMM day 0","LMD_MMM day+1" --title "Hor. wind histogram / 10m-600m AGL range: POWERED DESCENT" -d -S png -O lmd_histo1d_10m_600m_AGL_variability --vert 0,10 -i 4 -l 0.01,0.6,20 --ylabel "Density of probability" --xlabel "Horizontal wind speed (m/s)"
/san0/acolmd/SVN/trunk//UTIL/PYTHON/pp.py -f $lmdfiles/wrfout_d03$filelmda,$lmdfiles/wrfout_d03$filelmd,$lmdfiles/wrfout_d03$filelmdb -v uv --analysis histo --lon $lon --lat $lat --div 20 --labels "LMD_MMM day-1","LMD_MMM day 0","LMD_MMM day+1" --title "Hor. wind histogram / 1km-4km AML range: MAIN PARACHUTE DESCENT" -d -S png -O lmd_histo1d_1km_4km_AML_variability --vert 0,20 -i 3 -l 1,4,20 --ylabel "Density of probability" --xlabel "Horizontal wind speed (m/s)"
/san0/acolmd/SVN/trunk//UTIL/PYTHON/pp.py -f $lmdfiles/wrfout_d03$filelmda,$lmdfiles/wrfout_d03$filelmd,$lmdfiles/wrfout_d03$filelmdb -v uv --analysis histo --lon $lon --lat $lat --div 20 --labels "LMD_MMM day-1","LMD_MMM day 0","LMD_MMM day+1" --title "Hor. wind histogram / 4km-10km AML range: DROGUE AND PARACHUTE DEPLOYMENT" -d -S png -O lmd_histo1d_4km_10km_AML_variability --vert 0,20 -i 3 -l 4,10,20 --ylabel "Density of probability" --xlabel "Horizontal wind speed (m/s)"
/san0/acolmd/SVN/trunk//UTIL/PYTHON/pp.py -f $lmdfiles/wrfout_d03$filelmda,$lmdfiles/wrfout_d03$filelmd,$lmdfiles/wrfout_d03$filelmdb -v uv --analysis histo --lon $lon --lat $lat --div 20 --labels "LMD_MMM day-1","LMD_MMM day 0","LMD_MMM day+1" --title "Hor. wind histogram / 10km-40km AML range: DECELERATION PHASE" -d -S png -O lmd_histo1d_10km_40km_AML_variability --vert 0,20 -i 3 -l 10,40,20 --ylabel "Density of probability" --xlabel "Horizontal wind speed (m/s)"
montage lmd_histo1d_10m_600m_AGL_variability_200.png lmd_histo1d_1km_4km_AML_variability_200.png lmd_histo1d_4km_10km_AML_variability_200.png lmd_histo1d_10km_40km_AML_variability_200.png -mode concatenate -tile 2x2 lmd_histo1d_decomposed_variability.png
rm -f lmd_histo1d_10m_600m_AGL_variability_200.png lmd_histo1d_1km_4km_AML_variability_200.png lmd_histo1d_4km_10km_AML_variability_200.png lmd_histo1d_10km_40km_AML_variability_200.png lmd_histo1d_10m_600m_AGL_variability.sh lmd_histo1d_1km_4km_AML_variability.sh lmd_histo1d_4km_10km_AML_variability.sh lmd_histo1d_10km_40km_AML_variability.sh

montage mrams_histo2d_variability_200.png lmd_histo2d_variability_200.png -mode concatenate -tile 2x1 comp_histo2d_variability.png
rm -f mrams_histo2d_variability_200.png lmd_histo2d_variability_200.png mrams_histo2d_variability.sh lmd_histo2d_variability.sh
fi

if [ $2 == density ]
then
/san0/acolmd/SVN/trunk//UTIL/PYTHON/pp.py -f $mramsfiles/mramsout_d03$filemramsa,$mramsfiles/mramsout_d03$filemramsb -v uv --lon $lon --lat $lat --div 60 --labels "dummy","dummy","MRAMS day-1/0","dummy","dummy","MRAMS day+1/0" --title "MRAMS day-to-day horizontal winds variability at Landing Site" -d -S png -O mrams_histodensity_uv_variability --ope - --fref $mramsfiles/mramsout_d03$filemrams  --analysis histodensity --xmin -15 --xmax 15 --ymin 0 --ymax 0.35
/san0/acolmd/SVN/trunk//UTIL/PYTHON/pp.py -f $lmdfiles/wrfout_d03$filelmda,$lmdfiles/wrfout_d03$filelmdb -v uv --lon $lon --lat $lat --div 20 --labels "dummy","dummy","LMD_MMM day-1/0","dummy","dummy","LMD_MMM day+1/0" --title "LMD_MMM day-to-day horizontal winds variability at Landing Site" -d -S png -O lmd_histodensity_uv_variability --ope - --fref $lmdfiles/wrfout_d03$filelmd  --analysis histodensity  --xmin -15 --xmax 15 --ymin 0 --ymax 0.35
montage mrams_histodensity_uv_variability_200.png lmd_histodensity_uv_variability_200.png -mode concatenate -tile 2x1 comp_histodensity_uv_variability.png
rm -f mrams_histodensity_uv_variability_200.png lmd_histodensity_uv_variability_200.png mrams_histodensity_uv_variability.sh lmd_histodensity_uv_variability.sh

/san0/acolmd/SVN/trunk//UTIL/PYTHON/pp.py -f $mramsfiles/mramsout_d03$filemramsa,$mramsfiles/mramsout_d03$filemramsb -v tk --lon $lon --lat $lat --div 200 --labels "dummy","dummy","MRAMS day-1/0","dummy","dummy","MRAMS day+1/0" --title "MRAMS day-to-day temperature variability at Landing Site" -d -S png -O mrams_histodensity_tk_variability --ope - --fref $mramsfiles/mramsout_d03$filemrams --analysis histodensity --xmin -5 --xmax 5 --ymin 0 --ymax 0.8
/san0/acolmd/SVN/trunk//UTIL/PYTHON/pp.py -f $lmdfiles/wrfout_d03$filelmda,$lmdfiles/wrfout_d03$filelmdb -v tk --lon $lon --lat $lat --div 40 --labels "dummy","dummy","LMD_MMM day-1/0","dummy","dummy","LMD_MMM day+1/0" --title "LMD_MMM day-to-day temperature variability at Landing Site" -d -S png -O lmd_histodensity_tk_variability --ope - --fref $lmdfiles/wrfout_d03$filelmd --analysis histodensity  --xmin -5 --xmax 5 --ymin 0 --ymax 0.8
montage mrams_histodensity_tk_variability_200.png lmd_histodensity_tk_variability_200.png -mode concatenate -tile 2x1 comp_histodensity_tk_variability.png
rm -f mrams_histodensity_tk_variability_200.png lmd_histodensity_tk_variability_200.png mrams_histodensity_tk_variability.sh lmd_histodensity_tk_variability.sh 
fi

fi
