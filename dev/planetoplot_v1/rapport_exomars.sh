base=/planeto/acolmd/RUN/MESOSCALE/MERIDIANI_EXOMARS
file=wrfout_d01_2024-09-07_23:02:00
for i in tau05 tau1 tau2 tau5
do
##
##_TSURFLT01_Th_Vs_Conv
##
# pp.py --mark -6.13,-1.88 --mope -2 --Mope 2 -f $base/M_E_Th_$i/$file --fref $base/M_E_Conv_$i/$file -w HGT --operation - -v TSURF --time 1 --axtime lt --div 40 -m 180 -M 260 -c Spectral_r -d -S png --res 300 -O ${i}_TSURFLT01_Th_Vs_Conv --title "<TH+RiSL+MY4> Surface Temperature LT01 (K)" --titleref "<Convadj> Surface Temperature LT01 (K)"
##
##_TSURFLT14_Th_Vs_Conv
##
# pp.py --mark -6.13,-1.88 --mope -2 --Mope 2 -f $base/M_E_Th_$i/$file --fref $base/M_E_Conv_$i/$file -w HGT --operation - -v TSURF --time 14 --axtime lt --div 40 -m 260 -M 310 -c Spectral_r -d -S png --res 500 -O ${i}_TSURFLT14_Th_Vs_Conv --title "<TH+RiSL+MY4> Surface Temperature LT14 (K)" --titleref "<Convadj> Surface Temperature LT 14 (K)"
##
##_TSURFLT14
##
# pp.py --mark -6.13,-1.88 -f $base/M_E_Th_$i/$file -w HGT -v TSURF --time 14 --axtime lt --div 40 -m 260 -M 310 -c Spectral_r -d -S png --res 300 -O ${i}_TSURFLT14 --title "<TH+RiSL+MY4> Surface Temperature LT14 (K)"
##
##_TSURFLT16_Th_Vs_Conv
##
# pp.py --mark -6.13,-1.88 --mope -2 --Mope 2 -f $base/M_E_Th_$i/$file --fref $base/M_E_Conv_$i/$file -w HGT --operation - -v TSURF --time 16 --axtime lt --div 40 -m 240 -M 300 -c Spectral_r -d -S png --res 300 -O ${i}_TSURFLT16_Th_Vs_Conv --title "<TH+RiSL+MY4> Surface Temperature LT16 (K)" --titleref "<Convadj> Surface Temperature LT 16 (K)"
##
##_UV5LT15_Th_Vs_Conv
##
# pp.py --mark -6.13,-1.88 --mope -2 --Mope 2 -f $base/M_E_Th_$i/$file --fref $base/M_E_Conv_$i/$file -w HGT --operation - -v uvmet --vert 10 -i 4 -l 0.01 --time 15 --axtime lt --div 40 -m 0 -M 35 -c Spectral_r -d -S png --res 300 -O ${i}_UVLT14_Th_Vs_Conv --title "<TH+RiSL+MY4> Horizontal winds at 10m LT15 (K)" --titleref "<Convadj> First layer horizontal winds LT 15 (K)"
##
##_CHECK_BOUNDING_Th_Vs_Conv
##
# pp.py --axtime lt -f $base/M_E_Th_$i/$file --fref $base/M_E_Conv_$i/$file --operation - -v tk --div 30 -c Spectral_r -m 190 -M 255 --lon -24.5 --lat -20 -i 4 -l 0.0005,30,60 --mope -10 --Mope 10 --ymin 0 --ymax 30-d -S png --res 300 -O ${i}_CHECK_BOUNDING_Th_Vs_Conv --title "<TH+RiSL+MY4> Temperature Serie (K)" --titleref "<Convadj> Temperature Serie (K)"
##
##_TSERIEHIGH_Th_Vs_Conv
##
# pp.py --axtime lt -f $base/M_E_Th_$i/$file --fref $base/M_E_Conv_$i/$file --operation - -v tk --div 30 -c Spectral_r -m 190 -M 255 --lon -6.13 --lat -1.88 -i 4 -l 0.0005,30,60 --mope -10 --Mope 10 --ymin 0 --ymax 30-d -S png --res 300 -O ${i}_TSERIEHIGH_Th_Vs_Conv --title "<TH+RiSL+MY4> Temperature Serie (K)" --titleref "<Convadj> Temperature Serie (K)"
##
##_TSERIE_Th_Vs_Conv
##
# pp.py --axtime lt -f $base/M_E_Th_$i/$file --fref $base/M_E_Conv_$i/$file --operation - -v tk --div 30 -c Spectral_r -m 205 -M 255 --lon -6.13 --lat -1.88 -i 4 -l 0.0005,8,60 --mope -3 --Mope 3 --ymin 0 --ymax 8-d -S png --res 300 -O ${i}_TSERIE_Th_Vs_Conv --title "<TH+RiSL+MY4> Temperature Serie (K)" --titleref "<Convadj> Temperature Serie (K)"
##
##_TSERIE_Th
##
# pp.py --ylabel "Altitude above local surface (km)" --xlabel "Temperature (K)" --axtime lt -f $base/M_E_Th_$i//$file -v tk --div 30 -c Spectral_r -m 205 -M 255 --lon -6.13 --lat -1.88 -i 4 -l 0.005,10,60 --ymin 0 --ymax 10-d -S png --res 300 -O ${i}_TSERIE_Th --title "<TH+RiSL+MY4> Temperature Serie (K)"
##
##_UVSERIE_Th_Vs_Conv
##
# pp.py --axtime lt -f $base/M_E_Th_$i/$file --fref $base/M_E_Conv_$i/$file --operation - -v uvmet --div 30 -c Spectral_r -m 0 -M 35 --lon -6.13 --lat -1.88 -i 4 -l 0.0005,30,60 --mope -5 --Mope 5 --ymin 0 --ymax 30-d -S png --res 300 -O ${i}_UVSERIE_Th_Vs_Conv --title "<TH+RiSL+MY4> Horizontal Wind Serie (m/s)" --titleref "<Convadj> Horizontal Wind Serie (K)"
##
##_UVSERIE_Th
##
 #pp.py --axtime lt --ylabel "Altitude above Mola reference (km)" --xlabel "Horizontal Wind Velocity (m/s)" -f $base/M_E_Th_$i//$file -v uvmet --div 20 -c Spectral_r -m 0 -M 35 --lon -6.13 --lat -1.88 -i 3 -l -2,30,60 --ymin -2--ymax 22-d -S png --res 300 -O ${i}_UVSERIE_Th --title "<TH+RiSL+MY4> Horizontal Wind Serie"
##
##_UVSERIEZOOM_Th
##
# pp.py --axtime lt --ylabel "Altitude above local surface (km)" --xlabel "Horizontal Wind Velocity (m/s)" -f $base/M_E_Th_$i//$file -v uvmet --div 30 -c Spectral_r -m 0 -M 26 --lon -6.13 --lat -1.88 -i 4 -l 0.005,10,60 --ymin 0 --ymax 10-d -S png --res 300 -O ${i}_UVSERIEZOOM_Th --title "<TH+RiSL+MY4> Horizontal Wind Serie"
##
##_UVSERIEZOOMx2_Th
##
# pp.py --axtime lt --ylabel "Altitude above local surface (km)" --xlabel "Horizontal Wind Velocity (m/s)" -f $base/M_E_Th_$i//$file -v uvmet --div 30 -c Spectral_r -m 0 -M 16 --lon -6.13 --lat -1.88 -i 4 -l 0.0035,0.01,20 --ymin 0 --ymax 10 -d -S png --res 300 -O ${i}_UVSERIEZOOMx2_Th --title "<TH+RiSL+MY4> Horizontal Wind Serie"
##
##_WSERIE_Th
##
#pp.py --axtime lt --ylabel "Altitude above Mola reference (km)" --xlabel "Vertical Wind Velocity (m/s)" -f $base/M_E_Th_$i//$file -v W --div 31 -c Spectral_r -m -0.6 -M 0.6 --lon -6.13 --lat -1.88 -i 3 -l -2,30,60 --ymin -2--ymax 22-d -S png --res 300 -O ${i}_WSERIE_Th --title "<TH+RiSL+MY4> Vertical Wind Serie"
##########################
### NEAR SURFACE STUFF ###
##########################
##
##_DTMAP1mLT15_Th_Vs_Conv
##
#pp.py -w HGT -i 4 --mark -6.13,-1.88 -l 0.005,1,5 --mope -5 --Mope 5 -m -15 -M 70 --axtime lt -f $base/M_E_Th_$i/$file --fref $base/M_E_Conv_$i/$file --operation - -v DELTAT --div 30 -c Spectral_r --time 15 -d -S png --res 500 -O ${i}_DTMAP1mLT15_Th_Vs_Conv --title "<TH+RiSL+MY4> Temperature Delta at LT15 (m/s)" --titleref "<Convadj> Temperature Delta at LT15 (m/s)"
##
##_UVMAP1mLT15_Th_Vs_Conv
##
# pp.py -w HGT --mark -6.13,-1.88 --vert 0 --axtime lt -f $base/M_E_Th_$i/$file --fref $base/M_E_Conv_$i/$file --operation - -v U_OUT1 --div 30 -c Spectral_r --time 15 -m 0 -M 35 --mope -5 --Mope 5 -d -S png --res 500 -O ${i}_UVMAP1mLT15_Th_Vs_Conv --title "<TH+RiSL+MY4> Horizontal Winds at 1m and LT15 (m/s)" --titleref "<Convadj> Horizontal Winds at 1m and LT15 (m/s)"
##
##_TMAP1mLT15_Th_Vs_Conv
##
# pp.py --mark -6.13,-1.88 -w HGT --mope -15 --Mope 15 --axtime lt -f $base/M_E_Th_$i/$file --fref $base/M_E_Conv_$i/$file --operation - -v T_OUT1 --div 30 -c Spectral_r --time 15 -d -S png --res 500 -O ${i}_TMAP1mLT15_Th_Vs_Conv --title "<TH+RiSL+MY4> Temperature at 1m and LT15 (m/s)" --titleref "<Convadj> Temperature at 1m and LT15 (m/s)"
###
##_UVMAPLT15_Th_Vs_Conv
##
# pp.py -i 4 -l 0.004 -w HGT --mark -6.13,-1.88 --axtime lt -f $base/M_E_Th_$i/$file --fref $base/M_E_Conv_$i/$file --operation - -v uvmet --div 30 -c Spectral_r --time 15 -m 0 -M 35 --vert 4 --mope -5 --Mope 5 -d -S png --res 500 -O ${i}_UVMAPLT15_Th_Vs_Conv --title "<TH+RiSL+MY4> Horizontal Winds at 5m and LT15 (m/s)" --titleref "<Convadj> Horizontal Winds at 5m and LT15 (m/s)"
##
##_UVMAPLT14_500_Th & LT16
##
# pp.py -s 5 -i 4 -l 0.5 --axtime lt -W --mark -6.13,-1.88 -f $base/M_E_Th_$i/$file -v uvmet --div 30 -c Spectral_r --time 14 -m 0 -M 35 --vert 500 --mope -5 --Mope 5 -d -S png --res 300 -O ${i}_UVMAPLT14_500_Th --title "<TH+RiSL+MY4> Horizontal Winds at 500m AGL and LT14 (m/s)"
# pp.py -s 5 -N -i 4 -l 0.5 --axtime lt -W --mark -6.13,-1.88 -f $base/M_E_Th_$i/$file -v uvmet --div 30 -c Spectral_r --time 16 -m 0 -M 35 --vert 500 --mope -5 --Mope 5 -d -S png --res 300 -O ${i}_UVMAPLT16_500_Th --title "<TH+RiSL+MY4> Horizontal Winds at 500m AGL and LT16 (m/s)"
##
##_UVMAPLT14_1000_Th & LT16
##
#pp.py -s 5 -i 3 -l 1. --axtime lt -W --mark -6.13,-1.88 -f $base/M_E_Th_$i/$file -v uvmet --div 30 -c Spectral_r --time 14 -m 0 -M 35 --vert 1--mope -5 --Mope 5 -d -S png --res 300 -O ${i}_UVMAPLT14_1000_Th --title "<TH+RiSL+MY4> Horizontal Winds at 1km AMR and LT14 (m/s)"
#pp.py -s 5 -N -i 3 -l 1. --axtime lt -W --mark -6.13,-1.88 -f $base/M_E_Th_$i/$file -v uvmet --div 30 -c Spectral_r --time 16 -m 0 -M 35 --vert 1--mope -5 --Mope 5 -d -S png --res 300 -O ${i}_UVMAPLT16_1000_Th --title "<TH+RiSL+MY4> Horizontal Winds at 1km AMR and LT16 (m/s)"
##
##_UVMAPLT14_4000_Th & LT16
##
#pp.py -s 5 -i 3 -l 4. --axtime lt -W --mark -6.13,-1.88 -f $base/M_E_Th_$i/$file -v uvmet --div 30 -c Spectral_r --time 14 -m 0 -M 35 --vert 4--mope -5 --Mope 5 -d -S png --res 300 -O ${i}_UVMAPLT14_4000_Th --title "<TH+RiSL+MY4> Horizontal Winds at 4km AMR and LT14 (m/s)"
#pp.py -s 5 -N -i 3 -l 4. --axtime lt -W --mark -6.13,-1.88 -f $base/M_E_Th_$i/$file -v uvmet --div 30 -c Spectral_r --time 16 -m 0 -M 35 --vert 4--mope -5 --Mope 5 -d -S png --res 300 -O ${i}_UVMAPLT16_4000_Th --title "<TH+RiSL+MY4> Horizontal Winds at 4km AMR and LT16 (m/s)"
##
##_UVMAPLT14_6000_Th & LT16
##
#pp.py -s 5 -i 3 -l 6. --axtime lt -W --mark -6.13,-1.88 -f $base/M_E_Th_$i/$file -v uvmet --div 30 -c Spectral_r --time 14 -m 0 -M 35 --vert 6--mope -5 --Mope 5 -d -S png --res 300 -O ${i}_UVMAPLT14_6000_Th --title "<TH+RiSL+MY4> Horizontal Winds at 6km AMR and LT14 (m/s)"
#pp.py -s 5 -N -i 3 -l 6. --axtime lt --mark -6.13,-1.88 -W -f $base/M_E_Th_$i/$file -v uvmet --div 30 -c Spectral_r --time 16 -m 0 -M 35 --vert 6--mope -5 --Mope 5 -d -S png --res 300 -O ${i}_UVMAPLT16_6000_Th --title "<TH+RiSL+MY4> Horizontal Winds at 6km AMR and LT16 (m/s)"
##
##_UVMAPLT14_10000_Th & LT16
##
#pp.py -s 5 -i 3 -l 10. --axtime lt -W --mark -6.13,-1.88 -f $base/M_E_Th_$i/$file -v uvmet --div 30 -c Spectral_r --time 14 -m 0 -M 35 --vert 10--mope -5 --Mope 5 -d -S png --res 300 -O ${i}_UVMAPLT14_10000_Th --title "<TH+RiSL+MY4> Horizontal Winds at 10km AMR and LT14 (m/s)"
#pp.py -s 5 -N -i 3 -l 10; --axtime lt -W --mark -6.13,-1.88 -f $base/M_E_Th_$i/$file -v uvmet --div 30 -c Spectral_r --time 16 -m 0 -M 35 --vert 10--mope -5 --Mope 5 -d -S png --res 300 -O ${i}_UVMAPLT16_10000_Th --title "<TH+RiSL+MY4> Horizontal Winds at 10km AMR and LT16 (m/s)"
##
##_UVMAPLT14_20000_Th & LT16
##
# pp.py  -s 5 -i 3 -l 20. --axtime lt -W --mark -6.13,-1.88 -f $base/M_E_Th_$i/$file -v uvmet --div 30 -c Spectral_r --time 14 -m 0 -M 35 --vert 20--mope -5 --Mope 5 -d -S png --res 300 -O ${i}_UVMAPLT14_20000_Th --title "<TH+RiSL+MY4> Horizontal Winds at 20km AMR and LT14 (m/s)"
# pp.py -s 5  -i 3 -l 20. --axtime lt -W --mark -6.13,-1.88 -f $base/M_E_Th_$i/$file -v uvmet --div 30 -c Spectral_r --time 16 -m 0 -M 35 --vert 20--mope -5 --Mope 5 -d -S png --res 300 -O ${i}_UVMAPLT16_20000_Th --title "<TH+RiSL+MY4> Horizontal Winds at 20km AMR and LT16 (m/s)"
##
##_UVMAPLT14_1_Th & LT16
##
#pp.py -W --mark -6.13,-1.88 -s 7 --axtime lt -f $base/M_E_Th_$i/$file -v U_OUT1 --div 30 -c Spectral_r --time 14 -m 0 -M 25 --vert 0 -d -S png --res 300 -O ${i}_UVMAPLT14_1_Th --title "<TH+RiSL+MY4> Horizontal Winds interpolated at 1m AGL and LT14"
#pp.py -W --mark -6.13,-1.88 -s 7 --axtime lt -f $base/M_E_Th_$i/$file -v U_OUT1 --div 30 -c Spectral_r --time 16 -m 0 -M 25 --vert 0 -d -S png --res 300 -O ${i}_UVMAPLT16_1_Th --title "<TH+RiSL+MY4> Horizontal Winds interpolated at 1m AGL and LT16"
##
##_TMAPLT14_1_Th & LT16
##
## note: the T_OUT1_MOD variable is a special mod for pp.py that takes max(T_1,T_OUT1).
## This mod is not updated in the svn repository because pbl_parameters.F has been modified so that it does it.
##
#pp.py -m 230 -M 275 -i 4 --mark -6.13,-1.88 -l 0.004 --axtime lt -f $base/M_E_Th_$i/$file -v T_OUT1_MOD --div 30 -c Spectral_r --time 14 -d -S png --res 300 -O ${i}_TMAPLT14_1_Th --title "<TH+RiSL+MY4> Temperature interpolated at 1m AGL and LT14"
#pp.py -m 230 -M 275  -i 4 --mark -6.13,-1.88 -l 0.004 --axtime lt -f $base/M_E_Th_$i/$file -v T_OUT1_MOD --div 30 -c Spectral_r --time 16 -d -S png --res 300 -O ${i}_TMAPLT16_1_Th --title "<TH+RiSL+MY4> Temperature interpolated at 1m AGL and LT16"
#####################################
######### Vertical Winds ############
#####################################
#pp.py -s 5 -i 3 -l 20. --axtime lt --mark -6.13,-1.88 -f $base/M_E_Th_$i/$file -v W --div 61 -c Spectral_r --time 1 -m -0.6 -M 0.6 --vert 20--mope -5 --Mope 5 -d -S png --res 300 -O ${i}_WMAPLT14_20000_Th --title "<TH+RiSL+MY4> Vertical Winds at 20km AMR and LT14 (m/s)"
#pp.py -s 5 -N -i 3 -l 20. --axtime lt -f $base/M_E_Th_$i/$file -v W --div 61 -c Paired --time 10 -m -0.6 -M 0.6 --vert 20-d -S png --res 300 -O ${i}_WMAPLT10_20000_Th --title "<TH+RiSL+MY4> Vertical Winds at 20km AMR and LT10 (m/s)"
#pp.py -s 5 -i 3 -l 10 --axtime lt --mark -6.13,-1.88 -f $base/M_E_Th_$i/$file -v W --div 30 -c Spectral_r --time 13 -m -0.6 -M 1 --vert 10--mope -5 --Mope 5 -d -S png --res 300 -O ${i}_WMAPLT13_10000_Th --title "<TH+RiSL+MY4> Vertical Winds at 10km AMR and LT13 (m/s)"
#pp.py -s 5 -i 3 -l 6 --axtime lt --mark -6.13,-1.88 -f $base/M_E_Th_$i/$file -v W --div 30 -c Spectral_r --time 13 -m -0.6 -M 1 --vert 6--mope -5 --Mope 5 -d -S png --res 300 -O ${i}_WMAPLT13_6000_Th --title "<TH+RiSL+MY4> Vertical Winds at 6km AMR and LT13 (m/s)"
#pp.py -s 5 -i 3 -l 20 --axtime lt --mark -6.13,-1.88 -f $base/M_E_Th_$i/$file -v W --fref $base/M_E_Conv_$i//$file --operation - --div 30 -c Spectral_r --time 13 --titleref "<Convadj> Vertical Winds at 6km AMR and LT13 (m/s)" -m -0.6 -M 1 --vert 20--mope -0.5 --Mope 0.5 -d -S png --res 300 -O ${i}_WMAPLT13_20000_Th_vs_Conv --title "<TH+RiSL+MY4> Vertical Winds at 20km AMR and LT13 (m/s)"
#pp.py -s 5 -i 3 -l 10 --axtime lt --mark -6.13,-1.88 -f $base/M_E_Th_$i//$file -v W --div 30 -c Spectral_r --time 13 -m -0.6 -M 1 --vert 10--fref $base/M_E_Conv_$i//$file --titleref "<Convadj> Vertical Winds at 6km AMR and LT13 (m/s)" --operation - --mope -0.5 --Mope 0.5 -d -S png --res 300 -O ${i}_WMAPLT13_10000_Th_vs_Conv --title "<TH+RiSL+MY4> Vertical Winds at 10km AMR and LT13 (m/s)"
#pp.py -s 5 -i 3 -l 6 --axtime lt --mark -6.13,-1.88 -f $base/M_E_Th_$i//$file -v W --div 30 -c Spectral_r --time 13 -m -0.6 -M 1 --vert 6--mope -0.5 --Mope 0.5 -d -S png --res 300 --titleref "<Convadj> Vertical Winds at 6km AMR and LT13 (m/s)" -O ${i}_WMAPLT13_6000_Th_vs_Conv --fref $base/M_E_Conv_$i//$file --operation - --title "<TH+RiSL+MY4> Vertical Winds at 6km AMR and LT13 (m/s)"
#####################################
######## Convection #################
#####################################
##
##_WSTARMAPLT13_Th
##
# pp.py -w HGT --axtime lt -f $base/M_E_Th_$i/$file -v WSTAR --mark -6.13,-1.88 --div 30 -c jet --time 13 -m 0 -M 6 --res 300 -O ${i}_WSTARMAPLT13_Th --title "<TH+RiSL+MY4> Free convection velocity at LT13 (m/s)" -S png -d 
##
##_ZMAXMAPLT13_Th
##
# pp.py -w HGT --axtime lt -f $base/M_E_Th_$i/$file -v ZMAX_TH --mark -6.13,-1.88 --div 30 -c jet --time 13 -m 0 -M 9--res 300 -O ${i}_ZMAXMAPLT13_Th --title "<TH+RiSL+MY4> Boundary Layer Depth (km)" -S png -d
echo done
done
##
##_ALL_TSURFSERIE_Th_Vs_Conv
##
# pp.py -S png -d --lstyle "-b","--b","-r","--r","-y","--y","-g","--g" --labels "<TH+RiSL+MY4> tau=0.5","<Convadj> tau=0.5","<TH+RiSL+MY4> tau=1","<Convadj> tau=1","<TH+RiSL+MY4> tau=2","<Convadj> tau=2","<TH+RiSL+MY4> tau=5","<Convadj> tau=5" --axtime lt -f $base/M_E_Th_tau05/$file,$base/M_E_Conv_tau05/$file,$base/M_E_Th_tau1/$file,$base/M_E_Conv_tau1/$file,$base/M_E_Th_tau2/$file,$base/M_E_Conv_tau2/$file,$base/M_E_Th_tau5/$file,$base/M_E_Conv_tau5/$file -v TSURF --lon -6.13 --lat -1.88 --res 200 -O ALL_TSURFSERIE_Th_Vs_Conv --title "Surface Temperature Series (K)" --ymax 310 --ylabel "Temperature (K)" --xlabel "Local time (h)"
##
##ALL_ZMAXSERIE_Th
##
# pp.py -S png -d --lstyle "-b","-r","-y","-g" --labels "<TH+RiSL+MY4> tau=0.5","<TH+RiSL+MY4> tau=1","<TH+RiSL+MY4> tau=2","<TH+RiSL+MY4> tau=5" --axtime lt -f $base/M_E_Th_tau05/$file,$base/M_E_Th_tau1/$file,$base/M_E_Th_tau2/$file,$base/M_E_Th_tau5/$file -v ZMAX_TH --lon -6.13 --lat -1.88 --res 200 -O ALL_ZMAXSERIE_Th --title "Boundary Layer Height (km)" --ylabel "Boundary Layer Height (km)" --xlabel "Local time (h)"
##
##ALL_WSTARSERIE_Th
##
# pp.py -S png -d --lstyle "-b","-r","-y","-g" --labels "<TH+RiSL+MY4> tau=0.5","<TH+RiSL+MY4> tau=1","<TH+RiSL+MY4> tau=2","<TH+RiSL+MY4> tau=5" --axtime lt -f $base/M_E_Th_tau05/$file,$base/M_E_Th_tau1/$file,$base/M_E_Th_tau2/$file,$base/M_E_Th_tau5/$file -v WSTAR --lon -6.13 --lat -1.88 --res 200 -O ALL_WSTARSERIE_Th --title "Typical Mean Updraft Vertical Velocity (m/s)" --ylabel "Free Convection Velocity (m/s)" --xlabel "Local time (h)"
##
##ALL_TPROF_Th_Vs_Conv_LT15_zoom
##
# pp.py -S png -d --lstyle "-b","--b","-r","--r","-y","--y","-g","--g" --labels "<TH+RiSL+MY4> tau=0.5","<Convadj> tau=0.5","<TH+RiSL+MY4> tau=1","<Convadj> tau=1","<TH+RiSL+MY4> tau=2","<Convadj> tau=2","<TH+RiSL+MY4> tau=5","<Convadj> tau=5" --axtime lt -f $base/M_E_Th_tau05//$file,$base/M_E_Conv_tau05//$file,$base/M_E_Th_tau1//$file,$base/M_E_Conv_tau1//$file,$base/M_E_Th_tau2//$file,$base/M_E_Conv_tau2//$file,$base/M_E_Th_tau5//$file,$base/M_E_Conv_tau5//$file -v tk -i 4 -l 0.005,10,120 --lon -6.13 --lat -1.88 --time 15 --xmin 220 --xmax 270 --ymin 0 --ymax 10--res 200 -O ALL_TPROFZOOM_Th_Vs_Conv_LT15 --title "Temperature Profiles LT 15 (K)" --xlabel "Temperature (K)" --ylabel "Altitude Above Local Surface (km)"
##
##ALL_TPROF_Th_Vs_Conv_LT15
##
# pp.py -S png -d --lstyle "-b","--b","-r","--r","-y","--y","-g","--g" --labels "<TH+RiSL+MY4> tau=0.5","<Convadj> tau=0.5","<TH+RiSL+MY4> tau=1","<Convadj> tau=1","<TH+RiSL+MY4> tau=2","<Convadj> tau=2","<TH+RiSL+MY4> tau=5","<Convadj> tau=5" --axtime lt -f $base/M_E_Th_tau05/$file,$base/M_E_Conv_tau05/$file,$base/M_E_Th_tau1/$file,$base/M_E_Conv_tau1/$file,$base/M_E_Th_tau2/$file,$base/M_E_Conv_tau2/$file,$base/M_E_Th_tau5/$file,$base/M_E_Conv_tau5/$file -v tk -i 3 -l -3,30,60 --lon -6.13 --lat -1.88 --time 15 --xmin 190 --xmax 260 --ymin -2--ymax 30--res 200 -O ALL_TPROF_Th_Vs_Conv_LT15 --title "Temperature Profiles LT 15 (K)" --xlabel "Temperature (K)" --ylabel "Altitude Above Mola Reference (km)"
##
##ALL_TPROF_Th_Vs_Conv_LT06
##
# pp.py  -S png -d --lstyle "-b","--b","-r","--r","-y","--y","-g","--g" --labels "<TH+RiSL+MY4> tau=0.5","<Convadj> tau=0.5","<TH+RiSL+MY4> tau=1","<Convadj> tau=1","<TH+RiSL+MY4> tau=2","<Convadj> tau=2","<TH+RiSL+MY4> tau=5","<Convadj> tau=5" --axtime lt -f $base/M_E_Th_tau05/$file,$base/M_E_Conv_tau05/$file,$base/M_E_Th_tau1/$file,$base/M_E_Conv_tau1/$file,$base/M_E_Th_tau2/$file,$base/M_E_Conv_tau2/$file,$base/M_E_Th_tau5/$file,$base/M_E_Conv_tau5/$file -v tk -i 3 -l -3,30,60 --lon -6.13 --lat -1.88 --time 6 --xmin 190 --xmax 260 --ymin -2--ymax 30--res 200 -O ALL_TPROF_Th_Vs_Conv_LT06 --title "Temperature Profiles LT 06 (K)" --xlabel "Temperature (K)" --ylabel "Altitude Above Mola Reference (km)"
##
##ALL_TPROF_Th_Vs_Conv_LT00
##
# pp.py  -S png -d --lstyle "-b","--b","-r","--r","-y","--y","-g","--g" --labels "<TH+RiSL+MY4> tau=0.5","<Convadj> tau=0.5","<TH+RiSL+MY4> tau=1","<Convadj> tau=1","<TH+RiSL+MY4> tau=2","<Convadj> tau=2","<TH+RiSL+MY4> tau=5","<Convadj> tau=5" --axtime lt -f $base/M_E_Th_tau05/$file,$base/M_E_Conv_tau05/$file,$base/M_E_Th_tau1/$file,$base/M_E_Conv_tau1/$file,$base/M_E_Th_tau2/$file,$base/M_E_Conv_tau2/$file,$base/M_E_Th_tau5/$file,$base/M_E_Conv_tau5/$file -v tk -i 3 -l -3,30,60 --lon -6.13 --lat -1.88 --time 0 --xmin 190 --xmax 260 --ymin -2--ymax 30--res 200 -O ALL_TPROF_Th_Vs_Conv_LT00 --title "Temperature Profiles LT 00 (K)" --xlabel "Temperature (K)" --ylabel "Altitude Above Mola Reference (km)"
##
##ALL_TPROF_Th_Vs_Conv_LT20
##
# pp.py  -S png -d --lstyle "-b","--b","-r","--r","-y","--y","-g","--g" --labels "<TH+RiSL+MY4> tau=0.5","<Convadj> tau=0.5","<TH+RiSL+MY4> tau=1","<Convadj> tau=1","<TH+RiSL+MY4> tau=2","<Convadj> tau=2","<TH+RiSL+MY4> tau=5","<Convadj> tau=5" --axtime lt -f $base/M_E_Th_tau05/$file,$base/M_E_Conv_tau05/$file,$base/M_E_Th_tau1/$file,$base/M_E_Conv_tau1/$file,$base/M_E_Th_tau2/$file,$base/M_E_Conv_tau2/$file,$base/M_E_Th_tau5/$file,$base/M_E_Conv_tau5/$file -v tk -i 3 -l -3,30,60 --lon -6.13 --lat -1.88 --time 20 --xmin 190 --xmax 260 --ymin -2--ymax 30--res 200 -O ALL_TPROF_Th_Vs_Conv_LT20 --title "Temperature Profiles LT 20 (K)" --xlabel "Temperature (K)" --ylabel "Altitude Above Mola Reference (km)"
##
##ALL_UVPROFZOOM_Th_Vs_Conv_LT15
##
# pp.py -S png -d --lstyle "-b","--b","-r","--r","-y","--y","-g","--g" --labels "<TH+RiSL+MY4> tau=0.5","<Convadj> tau=0.5","<TH+RiSL+MY4> tau=1","<Convadj> tau=1","<TH+RiSL+MY4> tau=2","<Convadj> tau=2","<TH+RiSL+MY4> tau=5","<Convadj> tau=5" --axtime lt -f $base/M_E_Th_tau05//$file,$base/M_E_Conv_tau05//$file,$base/M_E_Th_tau1//$file,$base/M_E_Conv_tau1//$file,$base/M_E_Th_tau2//$file,$base/M_E_Conv_tau2//$file,$base/M_E_Th_tau5//$file,$base/M_E_Conv_tau5//$file -v uvmet -i 4 -l 0.005,7,120 --lon -6.13 --lat -1.88 --time 15 --xmin 0. --xmax 30. --ymin -0 --ymax 7--res 200 -O ALL_UVPROFZOOM_Th_Vs_Conv_LT15 --title "Horizontal Wind Profiles LT 15 (m/s)" --xlabel "Horizontal Winds (m/s)" --ylabel "Altitude Above Ground Level (km)"
##
##ALL_UVPROFZOOM_Th_Vs_Conv_LT15_t051
## 
# pp.py -S png -d --lstyle "-b","--b","-r","--r" --labels "<TH+RiSL+MY4> tau=0.5","<Convadj> tau=0.5","<TH+RiSL+MY4> tau=1","<Convadj> tau=1"  --axtime lt -f $base/M_E_Th_tau05//$file,$base/M_E_Conv_tau05//$file,$base/M_E_Th_tau1//$file,$base/M_E_Conv_tau1//$file -v uvmet -i 4 -l 0.005,7,120 --lon -6.13 --lat -1.88 --time 15 --xmin 0. --xmax 30. --ymin -0 --ymax 7--res 200 -O ALL_UVPROFZOOM_Th_Vs_Conv_LT15_t051 --title "Horizontal Wind Profiles LT 15 (m/s)" --xlabel "Horizontal Winds (m/s)" --ylabel "Altitude Above Ground Level (km)"
##
##ALL_UVPROFZOOM_Th_Vs_Conv_LT15_t25
## 
# pp.py -S png -d --lstyle "-y","--y","-g","--g" --labels "<TH+RiSL+MY4> tau=2","<Convadj> tau=2","<TH+RiSL+MY4> tau=5","<Convadj> tau=5"  --axtime lt -f $base/M_E_Th_tau2//$file,$base/M_E_Conv_tau2//$file,$base/M_E_Th_tau5//$file,$base/M_E_Conv_tau5//$file -v uvmet -i 4 -l 0.005,7,120 --lon -6.13 --lat -1.88 --time 15 --xmin 0. --xmax 30. --ymin -0 --ymax 7--res 200 -O ALL_UVPROFZOOM_Th_Vs_Conv_LT15_t25 --title "Horizontal Wind Profiles LT 15 (m/s)" --xlabel "Horizontal Winds (m/s)" --ylabel "Altitude Above Ground Level (km)"
##
##ALL_UVPROFZOOM_Th_Vs_Conv_LT13_t051
## 
# pp.py -N -S png -d --lstyle "-b","--b","-r","--r" --labels "<TH+RiSL+MY4> tau=0.5","<Convadj> tau=0.5","<TH+RiSL+MY4> tau=1","<Convadj> tau=1" --axtime lt -f $base/M_E_Th_tau05//$file,$base/M_E_Conv_tau05//$file,$base/M_E_Th_tau1//$file,$base/M_E_Conv_tau1//$file -v uvmet -i 4 -l 0.005,7,120 --lon -6.13 --lat -1.88 --time 13 --xmin 0. --xmax 30. --ymin -0 --ymax 7 --res 200 -O ALL_UVPROFZOOM_Th_Vs_Conv_LT13_t051 --title "Horizontal Wind Profiles LT 13 (m/s)" --xlabel "Horizontal Winds (m/s)" --ylabel "Altitude Above Ground Level (km)"
##
##ALL_UVPROFZOOM_Th_Vs_Conv_LT13_t25
## 
# pp.py -S png -d --lstyle "-y","--y","-g","--g" --labels "<TH+RiSL+MY4> tau=2","<Convadj> tau=2","<TH+RiSL+MY4> tau=5","<Convadj> tau=5" --axtime lt -f $base/M_E_Th_tau2//$file,$base/M_E_Conv_tau2//$file,$base/M_E_Th_tau5//$file,$base/M_E_Conv_tau5//$file -v uvmet -i 4 -l 0.005,7,120 --lon -6.13 --lat -1.88 --time 13 --xmin 0. --xmax 30. --ymin -0 --ymax 7--res 200 -O ALL_UVPROFZOOM_Th_Vs_Conv_LT13_t25 --title "Horizontal Wind Profiles LT 13 (m/s)" --xlabel "Horizontal Winds (m/s)" --ylabel "Altitude Above Ground Level (km)"
##
##ALL_TSERIEFIRST_Th_Vs_Conv
##
#pp.py -S png -d --lstyle "-b","--b","-r","--r","-y","--y","-g","--g" --labels "<TH+RiSL+MY4> tau=0.5","<Convadj> tau=0.5","<TH+RiSL+MY4> tau=1","<Convadj> tau=1","<TH+RiSL+MY4> tau=2","<Convadj> tau=2","<TH+RiSL+MY4> tau=5","<Convadj> tau=5" --axtime lt -f $base/M_E_Th_tau05/$file,$base/M_E_Conv_tau05/$file,$base/M_E_Th_tau1/$file,$base/M_E_Conv_tau1/$file,$base/M_E_Th_tau2/$file,$base/M_E_Conv_tau2/$file,$base/M_E_Th_tau5/$file,$base/M_E_Conv_tau5/$file -v T_OUT1_MOD -i 4 -l 0.0035 --lon -6.13 --lat -1.88 --vert 3.5 --xmin -2 --xmax 22 --res 200 -O ALL_TSERIEFIRST_Th_Vs_Conv --title "Temperature Series at 1m AGL" --ylabel "Temperature (K)" --xlabel "Local Time (h)"
##
##ALL_TSERIEFIRST_Th_Vs_Conv_t051
##
# pp.py --ymin 190 --ymax 270 -S png -d --lstyle "-b","--b","-r","--r" --labels "<TH+RiSL+MY4> tau=0.5","<Convadj> tau=0.5","<TH+RiSL+MY4> tau=1","<Convadj> tau=1" --axtime lt -f $base/M_E_Th_tau05/$file,$base/M_E_Conv_tau05/$file,$base/M_E_Th_tau1/$file,$base/M_E_Conv_tau1/$file -v T_OUT1_MOD -i 4 -l 0.0035 --lon -6.13 --lat -1.88 --vert 3.5 --xmin -2 --xmax 22 --res 200 -O ALL_TSERIEFIRST_Th_Vs_Conv_t051 --title "Temperature Series at 1m AGL" --ylabel "Temperature (K)" --xlabel "Local Time (h)"
##
##ALL_TSERIEFIRST_Th_Vs_Conv_t25
##
# pp.py --ymin 190 --ymax 270 -S png -d --lstyle "-y","--y","-g","--g" --labels "<TH+RiSL+MY4> tau=2","<Convadj> tau=2","<TH+RiSL+MY4> tau=5","<Convadj> tau=5" --axtime lt -f $base/M_E_Th_tau2/$file,$base/M_E_Conv_tau2/$file,$base/M_E_Th_tau5/$file,$base/M_E_Conv_tau5/$file -v T_OUT1_MOD -i 4 -l 0.0035 --lon -6.13 --lat -1.88 --vert 3.5 --xmin -2 --xmax 22 --res 200 -O ALL_TSERIEFIRST_Th_Vs_Conv_t25 --title "Temperature Series at 1m AGL" --ylabel "Temperature (K)" --xlabel "Local Time (h)"
##
##ALL_UVSERIE_Th_Vs_Conv
##
#pp.py -S png -d --lstyle "-b","--b","-r","--r","-y","--y","-g","--g" --labels "<TH+RiSL+MY4> tau=0.5","<Convadj> tau=0.5","<TH+RiSL+MY4> tau=1","<Convadj> tau=1","<TH+RiSL+MY4> tau=2","<Convadj> tau=2","<TH+RiSL+MY4> tau=5","<Convadj> tau=5" --axtime lt -f $base/M_E_Th_tau05/$file,$base/M_E_Conv_tau05/$file,$base/M_E_Th_tau1/$file,$base/M_E_Conv_tau1/$file,$base/M_E_Th_tau2/$file,$base/M_E_Conv_tau2/$file,$base/M_E_Th_tau5/$file,$base/M_E_Conv_tau5/$file -v uvmet -i 4 -l 0.005 --lon -6.13 --lat -1.88 --vert 5 --ymin 0 --ymax 15 --xmin -2 --xmax 22 --res 200 -O ALL_UVSERIE_Th_Vs_Conv --title "Horizontal Wind Series at 5m (m/s)" --ylabel "Horizontal Wind Velocity (m/s)" --ylabel "Local Time (h)"
##
##ALL_UVSERIE_Th
##
#pp.py -S png -d --lstyle "-b","-r","-y","-g" --labels "<TH+RiSL+MY4> tau=0.5","<TH+RiSL+MY4> tau=1","<TH+RiSL+MY4> tau=2","<TH+RiSL+MY4> tau=5" --axtime lt -f $base/M_E_Th_tau05/$file,$base/M_E_Th_tau1/$file,$base/M_E_Th_tau2/$file,$base/M_E_Th_tau5/$file -v U_OUT1 --lon -6.13 --lat -1.88 --vert 0 --ymin 0 --ymax 15 --xmin -2 --xmax 22 --res 200 -O ALL_UVSERIE_Th --title "Horizontal Winds interpolated at 1m (m/s)" --ylabel "Horizontal Wind Velocity (m/s)" --ylabel "Local Time (h)"
##
##ALL_UVSERIE_10m_Th
##
#pp.py -S png -d --lstyle "-b","-r","-y","-g" --labels "<TH+RiSL+MY4> tau=0.5","<TH+RiSL+MY4> tau=1","<TH+RiSL+MY4> tau=2","<TH+RiSL+MY4> tau=5" --axtime lt -f $base/M_E_Th_tau05/$file,$base/M_E_Th_tau1/$file,$base/M_E_Th_tau2/$file,$base/M_E_Th_tau5/$file -v uvmet -i 4 -l 0.01 --lon -6.13 --lat -1.88 --vert 10 --ymin 0 --ymax 20 --xmin -2 --xmax 22 --res 200 -O ALL_UVSERIE_10m_Th --title "Horizontal Winds at 10m AGL (m/s)" --ylabel "Horizontal Wind Velocity (m/s)" --ylabel "Local Time (h)"
##
##ALL_TSERIEint_Th
##
# pp.py  -S png -d --lstyle "-b","-r","-y","-g" --labels "<TH+RiSL+MY4> tau=0.5","<TH+RiSL+MY4> tau=1","<TH+RiSL+MY4> tau=2","<TH+RiSL+MY4> tau=5" --axtime lt -f $base/M_E_Th_tau05/$file,$base/M_E_Th_tau1/$file,$base/M_E_Th_tau2/$file,$base/M_E_Th_tau5/$file -v T_OUT1_MOD -i 4 -l 0.004 --lon -6.13 --lat -1.88 --vert 0 --xmin -2 --xmax 22 --res 200 -O ALL_TSERIEint_Th --title "Interpolated temperature at 1m (m/s)" --ylabel "Temperature (K)" --ylabel "Local Time (h)"
##
##ALL_UVSERIE_Th_Vs_Conv_t051
##
#pp.py -S png -d --lstyle "-b","--b","-r","--r" --labels "<TH+RiSL+MY4> tau=0.5","<Convadj> tau=0.5","<TH+RiSL+MY4> tau=1","<Convadj> tau=1" --axtime lt -f $base/M_E_Th_tau05/$file,$base/M_E_Conv_tau05/$file,$base/M_E_Th_tau1/$file,$base/M_E_Conv_tau1/$file -v U_OUT1 --lon -6.13 --lat -1.88 --vert 0 --ymin 0 --ymax 15 --xmin -2 --xmax 22 --res 200 -O ALL_UVSERIE_Th_Vs_Conv_t051 --title "Horizontal Winds interpolated at 1m (m/s)" --ylabel "Horizontal Wind Velocity (m/s)" --ylabel "Local Time (h)"
##
##ALL_UVSERIE_Th_Vs_Conv_t25
##
#pp.py -S png -d --lstyle "-y","--y","-g","--g" --labels "<TH+RiSL+MY4> tau=2","<Convadj> tau=2","<TH+RiSL+MY4> tau=5","<Convadj> tau=5" --axtime lt -f $base/M_E_Th_tau2/$file,$base/M_E_Conv_tau2/$file,$base/M_E_Th_tau5/$file,$base/M_E_Conv_tau5/$file -v uvmet -i 4 -l 0.005,10,20 --lon -6.13 --lat -1.88 --vert 0 --ymin 0 --ymax 15 --xmin -2 --xmax 22 --res 200 -O ALL_UVSERIE_Th_Vs_Conv_t25 --title "Horizontal Wind Series at 5m (m/s)" --ylabel "Horizontal Wind Velocity (m/s)" --ylabel "Local Time (h)"
##
##ALL_UVPROF_Th_Vs_Conv_LT00
##
# pp.py -S png -d --lstyle "-b","--b","-r","--r","-y","--y","-g","--g" --labels "<TH+RiSL+MY4> tau=0.5","<Convadj> tau=0.5","<TH+RiSL+MY4> tau=1","<Convadj> tau=1","<TH+RiSL+MY4> tau=2","<Convadj> tau=2","<TH+RiSL+MY4> tau=5","<Convadj> tau=5" --axtime lt -f $base/M_E_Th_tau05/$file,$base/M_E_Conv_tau05/$file,$base/M_E_Th_tau1/$file,$base/M_E_Conv_tau1/$file,$base/M_E_Th_tau2/$file,$base/M_E_Conv_tau2/$file,$base/M_E_Th_tau5/$file,$base/M_E_Conv_tau5/$file -v uvmet -i 3 -l -3,30,60 --lon -6.13 --lat -1.88 --time 0 --xmin 0. --xmax 50. --ymin -2--ymax 30--res 200 -O ALL_UVPROF_Th_Vs_Conv_LT00 --title "Horizontal Wind Profiles LT 15 (m/s)" --xlabel "Horizontal Winds (m/s)" --ylabel "Altitude Above Mola Reference (km)"
##
##ALL_UVPROF_Th_Vs_Conv_LT06
##
# pp.py  -S png -d --lstyle "-b","--b","-r","--r","-y","--y","-g","--g" --labels "<TH+RiSL+MY4> tau=0.5","<Convadj> tau=0.5","<TH+RiSL+MY4> tau=1","<Convadj> tau=1","<TH+RiSL+MY4> tau=2","<Convadj> tau=2","<TH+RiSL+MY4> tau=5","<Convadj> tau=5" --axtime lt -f $base/M_E_Th_tau05/$file,$base/M_E_Conv_tau05/$file,$base/M_E_Th_tau1/$file,$base/M_E_Conv_tau1/$file,$base/M_E_Th_tau2/$file,$base/M_E_Conv_tau2/$file,$base/M_E_Th_tau5/$file,$base/M_E_Conv_tau5/$file -v uvmet -i 3 -l -3,30,60 --lon -6.13 --lat -1.88 --time 6 --xmin 0. --xmax 50. --ymin -2--ymax 30--res 200 -O ALL_UVPROF_Th_Vs_Conv_LT06 --title "Horizontal Wind Profiles LT 15 (m/s)" --xlabel "Horizontal Winds (m/s)" --ylabel "Altitude Above Mola Reference (km)"
##
##ALL_UVPROF_Th_Vs_Conv_LT20
##
# pp.py  -S png -d --lstyle "-b","--b","-r","--r","-y","--y","-g","--g" --labels "<TH+RiSL+MY4> tau=0.5","<Convadj> tau=0.5","<TH+RiSL+MY4> tau=1","<Convadj> tau=1","<TH+RiSL+MY4> tau=2","<Convadj> tau=2","<TH+RiSL+MY4> tau=5","<Convadj> tau=5" --axtime lt -f $base/M_E_Th_tau05/$file,$base/M_E_Conv_tau05/$file,$base/M_E_Th_tau1/$file,$base/M_E_Conv_tau1/$file,$base/M_E_Th_tau2/$file,$base/M_E_Conv_tau2/$file,$base/M_E_Th_tau5/$file,$base/M_E_Conv_tau5/$file -v uvmet -i 3 -l -3,30,60 --lon -6.13 --lat -1.88 --time 20 --xmin 0. --xmax 50. --ymin -2--ymax 30--res 200 -O ALL_UVPROF_Th_Vs_Conv_LT20 --title "Horizontal Wind Profiles LT 15 (m/s)" --xlabel "Horizontal Winds (m/s)" --ylabel "Altitude Above Mola Reference (km)"
