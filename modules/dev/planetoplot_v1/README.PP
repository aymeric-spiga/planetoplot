**************************************
**************************************
**************************************
    PLANETOPLOT TUTORIAL EXAMPLES
**************************************
         Authors : AC + AS
**************************************
  DON'T FORGET YOUR BEST FRIEND IS
     pp.py -h [or] pp.py --help
**************************************
**************************************
**************************************

The program pp.py works on virtually any NETCDF file. 
If not, please contact us, it is usually an easy fix to make it work.

NB: To use the movie capability, you will have to install mencoder.
NB: To use advance interpolation methods for GCM and mesoscale files, you will have to compile fortran routines with f2py (or compile those independently).

************************************
SIMPLE EXAMPLES on a SAMPLE GCM FILE
************************************
Goal: The simplest, most minimal example. Mapping topography.
pp.py -f diagfired.nc

Goal: I would like finer contours.
pp.py -f diagfired.nc --div 30

Goal: I would like wind vectors.
pp.py -f diagfired.nc -W

Goal: I would like more vectors [i.e. lower the stride].
pp.py -f diagfired.nc -W -s 1

Goal: I want to map a given field (surface temperature).
pp.py -f diagfired.nc -v tsurf

Goal: I want to map two fields next to one another (topography and tauice).
pp.py -f diagfired.nc -v phisinit,tauice

Goal: I want to map two fields, tauice shaded, topography contoured, same plot.
pp.py -f diagfired.nc -v tauice -w phisinit

Goal: I want to map a field but projected on the sphere.
pp.py -f diagfired.nc -v tauice -p ortho

Goal: I want to redefine the minimum and maximum values shown.
pp.py -f diagfired.nc -v tauice -m 0.2 -M 0.9

Goal: I want to insert holes wherever values are lower than 0.2 and higher than 0.9
pp.py -f diagfired.nc -v tauice -m 0.2 -M 0.9 -H

Goal: I want to fill holes with an background image of Mars [you have to be connected to Internet]
pp.py -f diagfired.nc -v tauice -m 0.2 -M 0.9 -H -b vishires

Goal: I want the same map, but projected on the sphere
pp.py -f diagfired.nc -v tauice -m 0.2 -M 0.9 -H -b vishires -p ortho

Goal: I want the same map, but projected on the sphere and centered on longitude 100°
pp.py -f diagfired.nc -v tauice -m 0.2 -M 0.9 -H -b vishires -p ortho --blon -100

Goal: I want the same map, but projected with north polar stereographic view
pp.py -f diagfired.nc -v tauice -m 0.2 -M 0.9 -H -b vishires -p npstere

Goal: I want the same map, but with a transparent field to see background image 
pp.py -f diagfired.nc -v tauice -m 0.2 -M 0.9 -H -b vishires -p npstere --trans 0.6

Goal: I want to save this in PNG format
pp.py -f diagfired.nc -v tauice -m 0.2 -M 0.9 -H -b vishires -p ortho -S png

Goal: I want to animate this along time axis
pp.py -f diagfired.nc -v tauice -m 0.2 -M 0.9 -H -b vishires -p ortho -S avi

Goal: I want to animate this along time axis with fps=12
pp.py -f diagfired.nc -v tauice -m 0.2 -M 0.9 -H -b vishires -p ortho --rate 12

Goal: I want to plot results from two simulation files next to one another
pp.py -f diagfired.nc,diagfired.nc -v tsurf 

Goal: I want to plot results for two different times in the file next to one another
pp.py -f diagfired.nc -v tsurf --time 0.5 --time 0.9

Goal: I want to plot averaged results in the file from one time to another time
pp.py -f diagfired.nc -v tsurf --time 0.5,0.9

Goal: I want to plot the minimum/maximum value over times stored in the file
pp.py -f diagfired.nc -v tsurf --redope mint
pp.py -f diagfired.nc -v tsurf --redope maxt

Goal: I want to plot a section of temperature at longitude 0
pp.py -f diagfired.nc --var temp --lon 0 --time 0.5

Goal: I want to plot a section of zonally-averaged temperature
pp.py -f diagfired.nc -v temp --lon -180,180 --time 0.5

Goal: I want to plot a section of zonally-averaged temperature with contours of zonally-averaged zonal wind
pp.py -f diagfired.nc -v temp -w u --lon -180,180 --time 0.5

Goal: I want to plot a section of zonally-averaged temperature and zonally-averaged zonal wind
pp.py -f diagfired.nc -v temp,u --lon -180,180 --time 0.5

Goal: I want to plot a globally-averaged 1D temperature profile
pp.py -f diagfired.nc -v temp --time 0.5 --lat -90,90 --lon -180,180

Goal: I want to overplot few globally-averaged 1D temperature profiles at different times
pp.py -f diagfired.nc -v temp --time 0.5 --time 0.9 --lat -90,90 --lon -180,180

Goal: I want to overplot lat=0,lon=0 1D temperature profiles at different times
pp.py -f diagfired.nc -v temp --time 0.5 --time 0.9 --lat 0 --lon 0 --ymax 20 --xmin 180.

******************** [specific mesoscale]

Goal: I want to plot results for two different LOCAL times in the file next to one another
pp.py -f wrfout_d01_2024-05-03_01:00:00 -v TSURF --time 4 --time 7 --axtime lt

Goal: I want to plot a section of temperature at longitude -120 between AMR altitudes 0 km and 30 km with 50 levels
pp.py -f wrfout_d01_2024-05-03_01:00:00 -v tk --time 4 --axtime lt --lon -120 -i 3 -l 0,30,50

Goal: I want to plot a section of temperature at longitude -120 between ALS altitudes 0 km and 30 km with 50 levels
pp.py -f wrfout_d01_2024-05-03_01:00:00 -v tk --time 4 --axtime lt --lon -120 -i 4 -l 0,30,50

***********************************************************************************
EXAMPLE : The classic mountain GW plot
***********************************************************************************
pp.py -f wrfout_d01_9999-09-09_09:00:00 -v W,tpot --lat 60 --time 15 -i 4 -l 30,130,100 --div 50
***********************************************************************************

***********************************************************************************
EXAMPLE : Plotting surface temperature maps from TES
***********************************************************************************
pp.py -f TES/TES.MappedClimatology.nadir.MY25.nc -v Tsurf_day -p npstere --blat 60 --trans 0.5 -b vishires --time 50
***********************************************************************************
pp.py -f /d5/emlmd/TES/TES.MappedClimatology.nadir.MY25.nc,/d5/emlmd/TES/TES.MappedClimatology.nadir.MY26.nc   -p npstere -v Tsurf_day --time 90,105 --time 105,120 --time 120,135 --time 135,150 --blat 70 --div 30 -m 190 -M 270 -b vishires --trans 0.7 
***********************************************************************************
pp.py -f /d5/emlmd/TES/TES.MappedClimatology.nadir.MY25.nc,/d5/emlmd/TES/TES.MappedClimatology.nadir.MY26.nc -p npstere -v Tsurf_day --time 110,130 --blat 70 --div 30 -m 190 -M 270 -b vishires --trans 0.7 -c onebar --title "Surface temperature (K)"
***********************************************************************************
dust from TES
pp.py -f TES.MappedClimatology.nadir.MY25.nc -v tau_dust -p cyl -c Oranges -m 0.0 -M 2.0 --div 20 --trans 0.8 -b molabw -S avi -T

***********************************************************************************
COMMENTED EXAMPLE : The globe with surface temperature and winds
***********************************************************************************
pp.py -f diagfired.nc -v tsurf -w phisinit -m 120 -M 320 --div 20 -W -s 1 --vert 0 -p ortho --blat 20 --blon -80 -S html -t $W
***********************************************************************************
See results here: http://www.lmd.jussieu.fr/~aslmd/EXAMPLES/LMD_GCM_movie_tsurf_UV/anim.html
***********************************************************************************
pp.py -f diagfired.nc
	OK. You probably get that one.
-v tsurf -w phisinit
	Shade surface temperature. Contour topography.
-m 120 -M 320 --div 20
	Surface temperature is shown with bounds 120K to 320K. Use 20 levels for shading.
-W -s 1
	Include wind vectors. Prescribe a stride of 1: vectors are shown at every grid point.
--vert 0
	Show fields in the first (lowermost) level.
-p ortho --blat 20 --blon -80
	Use orthographic projection ('whole sphere' view). Center view on lon -80E and lat 20N.
-S html
	Make nice webpage with animation and controls.
-t /u/aslmd/WWW/EXAMPLES
	Move resulting plot files to the given folder.

***********************************************************************************
COMMENTED EXAMPLE : The dust storm section movie
***********************************************************************************
pp.py -f wrfout_d01_2024-05-30_12:00:00,wrfout_d01_2024-05-30_18:00:00,wrfout_d01_2024-05-31_00\:00\:00 --operation cat -v QDUST --lat -3. -i 3 -l -1,37,100 --div 30 -c Oranges_r -m 0. -M 5.e-5 -t $W --rate 12 --xmin=5 --xmax=115
***********************************************************************************
See results here: http://www.lmd.jussieu.fr/~aslmd/EXAMPLES/LMD_MMM_d1_10km_movie_QDUST_-1000m-AMR_lat_-3_Ls134.8/anim.html
***********************************************************************************
pp.py 
	OK. You probably get that one.
-f wrfout_d01_2024-05-30_12:00:00,wrfout_d01_2024-05-30_18:00:00,wrfout_d01_2024-05-31_00\:00\:00 --operation cat 
	Mesoscale outputs are splitted in several files. In that case, this was 1 file per 6 simulated hours. 
	The above options allow to concatenate files along time axis for 1D time series or animated movies 
	[if --operation cat is omitted, this is a multiplot call, with one subplot per files in -f].
-v QDUST 
	Choose to plot dust mass mixing ratio.
-i 3 -l -1,37,100 
	Set a call to vertical interpolator [compiled with f2py, thereby being embedded as a Python routine] for each of the files in the -f instance. 
	-i sets the kind of interpolation, 3 means Above MOLA Reference Altitude. -l sets the range for altitude levels: from -1 km to 37 km with 100 levels.
--div 30 
	The number of contours used for shaded plots. Higher value means smoother appearance.
-c Oranges_r 
	Choose a colorbar adapted to display a dust storm.
-m 0. -M 5.e-5 
	Choose bounds for the plotted field. This one is adapted to show dust mass mixing ratio.
-t $W
	Put resulting figure or movie in another destination folder. 
	Personally I have an environnement variable W which is somewhere in my system where the file automatically appears on the web, 
		hence is easy to see from a remote place.
--rate 12 --lat -3.
	Define prescribed axis. A section in latitude -3°N. A time animation with 12 frame per seconds. 
	So the displayed field will be an altitude/longitude section. Alternative: "-S avi" instead of "--rate 12" creates a default 8 fps movie. 
	Alternative II: "-S html" instead of --rate 12 creates a nice webpage.
--xmin=5 --xmax=115
	Define limits for the displayed section. Here we just want to get rid of transition rows where atmospheric fields are 
		relaxed towards prescribed GCM fields.


***********************************************************************************
Simple 2D plot: Zonal mean.
**********************************************************************************

Goal:

Plot the zonal mean temperature from a netcdf field representing one month.

Command:

pp.py -f POLAR_NIGHT_RUN/diagfi16.nc --var temp --lon 180,-180 --time 0,65

Note:

The --time, --lat, --lon and --vert command takes in input values corresponding to the unit stored in the netcdf file, and not indices ! For example, if the "Time" unit is in sol (which is common for a gcm output), --time 2 means sol 2.0 and not index 2 along the time direction. Consequently, one can ask --time 2.5 for temperatures at 12:00 on sol 2.

Means are easy to perform by specifying a range. Here: --time 0,65 and --lat -180,180.

***********************************************************************************
Vertical interpolation of the field.
***********************************************************************************

Goal: 

Calls to zrecast and api are built-in the python functions. One can call them using -i with the appropriate argument (see meso.py -h or gcm.py -h). Here is an example that re-interpolates data using zrecast before plotting it in a 2D contour.

Command:

pp.py -f POLAR_NIGHT_RUN/diagfi16.nc --var temp --lon -180,180 --time 0,65 -i 4

Note:

All interpolation modes available in zrecast and api (pressure, AGL, distance from planet center, etc...) are theoretically possible, but may not be coded yet in the routine. See gcm.py -h or meso.py -h.

For this example, the default behavior of zrecast for -i 4 is to interpolate in (m) from the local surface, between 0 and 150 km. The command will generate a reinterpolated netcdf file "POLAR_NIGHT_RUN/diagfi16_S.nc" with only the requested field, which is not deleted afterward.


***********************************************************************************
2D plot of the difference between two files.
***********************************************************************************

Goal:

Comparing two .nc files with similar dimension axis can be done in a single command, by specifying which files to compare and the comparison operator (i.e. is it a difference, an addition, etc...). When comparing data along a vertical axis, it can be wise to also ask for an interpolation of the fields to make sure the comparison is correct.

Command:

pp.py -f POLAR_NIGHT_RUN/stats16.nc --var temp --lon -180,180 --time 1 -i 4 --fref POLAR_NIGHT_REF/stats16.nc --operation - --mope -2 --Mope 2 --title "Polar temperatures with new parametrizations" --titleref "Reference run"

Note:

The command will output 3 plots: the field from file 1, the field from file 2, and the comparison between the two. One can specify specific names for the title of these plots by using --title and --titleref for the titles of file 1 and file 2, and can specify different plotting range for the normal field (-m -M) and the compared field (--mope --Mope).

One can combine this command with projections and means, so that for example, to compare co2 depletion at the south pole:

pp.py -f POLAR_NIGHT_RUN/start16.nc --var co2 --vert 0,150 --proj spstere --time 1 -i 4 --fref POLAR_NIGHT_REF/start16.nc --operation - --mope -0.5 --Mope 0.5 --title "Polar co2 with new parametrizations" --titleref "Reference run"

***********************************************************************************
2D plot of data with missing values, along a pressure axis (decreasing with height).
***********************************************************************************

Goal:

By default, python will force the y-axis of a 2D plot to be ordered by increasing values. Here is how to force it otherwise.

Command:

pp.py -f TES.MappedClimatology.nadir.MY25.nc --var T_nadir_day --lat -76. --time 90 --ymin 500 --ymax 1 -m 128 -M 148 -H

Note:

The axis reversal is done by specifying ymin and ymax in the right order. One can also simply use --inverty and not specify (ymin,ymax). Missing values (out of range values) are replaced by holes by the option "-H".
