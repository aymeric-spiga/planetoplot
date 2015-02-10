#/bin/bash
touch cdo-current.tar.gz
rm -rf cdo-current.tar.gz
wget https://code.zmaw.de/attachments/download/3604/cdo-current.tar.gz --no-check-certificate

touch cdo-1.5.6.1
rm -rf cdo-1.5.6.1
tar -xvzf cdo-current.tar.gz
cd cdo-1.5.6.1 

./configure --with-netcdf=/donnees/emlmd/netcdf64-4.0.1_pgi --prefix=/u/acolmd/san0/cdo/cdo-1.5.6.1 CPPFLAGS="-DNDEBUG -DpgiFortran" CC=pgcc CFLAGS="-O2 -Msignextend -Mipa" CXX=pgCC CXXFLAGS="-O2 -Msignextend -Mipa" FC=pgf90 FFLAGS="-O2 -Mipa"

make

make install

