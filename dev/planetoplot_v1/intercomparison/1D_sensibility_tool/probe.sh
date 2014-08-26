#!/bin/bash
taus="0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0"
TAU=($taus)

if [ $1 == run ]
then

for i in 0 1 2 3 4 5 6 7 8 9
do
rm -rf r$i
mkdir r$i
ln -sf ../run.def r$i/.
ln -sf ../z2sig.def r$i/.
ln -sf ../traceur.def r$i/.
ln -sf ../testphys1d.e r$i/.
ln -sf ../profile r$i/profile
cp callphys.def r$i/.
sed -e s/'tauvis=999'/'tauvis='${TAU[ $i ]}/g r$i/callphys.def > r$i/callphys.def.tmp ; \mv r$i/callphys.def.tmp r$i/callphys.def
done

for i in 0 1 2 3 4 5 6 7 8 9
do
cd r$i/.
echo 'running case '$i
./testphys1d.e > a.out
cd ..
done

fi

if [ $1 == analyse -o $1 == run ]
then
  if [ $2 == eps ]
  then
    format='eps'
  else
    format='png'
  fi
pp.py -f r0/diagfi.nc,r1/diagfi.nc,r2/diagfi.nc,r3/diagfi.nc,r4/diagfi.nc,r5/diagfi.nc,r6/diagfi.nc,r7/diagfi.nc,r8/diagfi.nc,r9/diagfi.nc -v fluxsurf_lw --labels "tau=0.1","tau=0.2","tau=0.3","tau=0.4","tau=0.5","tau=0.6","tau=0.7","tau=0.8","tau=0.9","tau=1.0" --lstyle '-k','--k','-r','--r','-y','--y','-g','--g','-b','--b' -S $format -d -O fluxsurf_lw_taus --title "Incident LW flux at surface"
pp.py -f r0/diagfi.nc,r1/diagfi.nc,r2/diagfi.nc,r3/diagfi.nc,r4/diagfi.nc,r5/diagfi.nc,r6/diagfi.nc,r7/diagfi.nc,r8/diagfi.nc,r9/diagfi.nc -v fluxsurf_sw --labels "tau=0.1","tau=0.2","tau=0.3","tau=0.4","tau=0.5","tau=0.6","tau=0.7","tau=0.8","tau=0.9","tau=1.0" --lstyle '-k','--k','-r','--r','-y','--y','-g','--g','-b','--b' -S $format -d -O fluxsurf_sw_taus --title "Incident SW flux at surface"
pp.py -f r0/diagfi.nc,r1/diagfi.nc,r2/diagfi.nc,r3/diagfi.nc,r4/diagfi.nc,r5/diagfi.nc,r6/diagfi.nc,r7/diagfi.nc,r8/diagfi.nc,r9/diagfi.nc -v tsurf --labels "tau=0.1","tau=0.2","tau=0.3","tau=0.4","tau=0.5","tau=0.6","tau=0.7","tau=0.8","tau=0.9","tau=1.0" --lstyle '-k','--k','-r','--r','-y','--y','-g','--g','-b','--b' -S $format -d -O tsurf_taus --title "Surface temperature"
pp.py -f r0/diagfi.nc,r1/diagfi.nc,r2/diagfi.nc,r3/diagfi.nc,r4/diagfi.nc,r5/diagfi.nc,r6/diagfi.nc,r7/diagfi.nc,r8/diagfi.nc,r9/diagfi.nc -v "rh-1" --labels "tau=0.1","tau=0.2","tau=0.3","tau=0.4","tau=0.5","tau=0.6","tau=0.7","tau=0.8","tau=0.9","tau=1.0" --lstyle '-k','--k','-r','--r','-y','--y','-g','--g','-b','--b' -S $format -d -O invrh_taus --title "Heat aerodyn conductance"

if [ $format == 'png' ]
then
  montage fluxsurf_lw_taus_200.png fluxsurf_sw_taus_200.png tsurf_taus_200.png invrh_taus_200.png -mode concatenate -tile 2x2 tau_sensitivity.png
  \mv fluxsurf_lw_taus_200.png fluxsurf_sw_taus_200.png tsurf_taus_200.png invrh_taus_200.png figures/.
  \rm fluxsurf_lw_taus.sh fluxsurf_sw_taus.sh tsurf_taus.sh invrh_taus.sh
  cp tau_sensitivity.png figures/.
fi
\rm fluxsurf_lw_taus.sh fluxsurf_sw_taus.sh tsurf_taus.sh invrh_taus.sh
fi

if [ $1 == clean ]
then

rm -rf figures
mkdir figures
rm -f *.png

fi
