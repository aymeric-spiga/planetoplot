#! /bin/bash

./pp.py  /home/aymeric/Big_Data/DATAPLOT/diagfired.nc
./pp.py  /home/aymeric/Big_Data/DATAPLOT/diagfired.nc -v u -v v -t 0.5 -z 10.
./pp.py  /home/aymeric/Big_Data/DATAPLOT/diagfired.nc -v u -t 0.5 -z 10. -y 2
./pp.py  /home/aymeric/Big_Data/DATAPLOT/diagfired.nc -v u -t 0.5 -z 10. -y 2 -m 1000 -L --
./pp.py  /home/aymeric/Big_Data/DATAPLOT/diagfired.nc -v u -t 0.5 -z 10. -T Yeah -C spectral
./pp.py  /home/aymeric/Big_Data/DATAPLOT/diagfired.nc -v u -t 0.5 -z 10. -T Yeah -P ortho
./pp.py  /home/aymeric/Big_Data/DATAPLOT/diagfired.nc -v u -t 0.5 -z 10. -T Yeah -P moll
./pp.py  /home/aymeric/Big_Data/DATAPLOT/diagfired.nc -v u -t 0.5 -z 10. -T Yeah -P npstere
./pp.py  /home/aymeric/Big_Data/DATAPLOT/diagfired.nc -v u -t 0.5 -z 10. -T Yeah -P spstere
./pp.py  /home/aymeric/Big_Data/DATAPLOT/diagfired.nc -v u -t 0.5 -z 10. -T Yeah -P cyl -A Tharsis
./pp.py  /home/aymeric/Big_Data/DATAPLOT/diagfired.nc -v u -t 0.5 -z 10. -T Yeah -P lcc -A Tharsis
./pp.py  /home/aymeric/Big_Data/DATAPLOT/diagfired.nc -v u -t 0.5 -z 10. -T Yeah -P laea -A Tharsis
./pp.py  /home/aymeric/Big_Data/DATAPLOT/diagfired.nc -v ps -t 0.5 -z 10. -C jet -A Tropics -P merc
./pp.py  /home/aymeric/Big_Data/DATAPLOT/diagfired.nc -v u -c phisinit -t 0.5 -z 10. -T Yeah -P laea -A Tharsis
./pp.py  /home/aymeric/Big_Data/DATAPLOT/diagfired.nc -v ps -t 0.5 -z 10. -C jet -P merc -A Whole_No_High
./pp.py  /home/aymeric/Big_Data/DATAPLOT/diagfired.nc -v ps -t 0.5 -z 10. -C jet -P robin

./pp.py  /home/aymeric/Big_Data/GALE/wrfout_d03_2024-06-09_00\:00\:00 
./pp.py  /home/aymeric/Big_Data/GALE/wrfout_d03_2024-06-09_00\:00\:00 -v HGT -B vishires
./pp.py  /home/aymeric/Big_Data/GALE/wrfout_d03_2024-06-09_00\:00\:00 -v HGT -t 1. -B vishires -H 0.5
./pp.py  /home/aymeric/Big_Data/GALE/wrfout_d03_2024-06-09_00\:00\:00 -v HGT -t 1. -P ortho -B vishires

./pp.py  /home/aymeric/Big_Data/DATAPLOT/diagfired.nc -v ps -t 0.5 -z 10. -C jet -P robin
./pp.py  /home/aymeric/Big_Data/DATAPLOT/diagfired.nc -v u -t 0.5 -x 10 -C jet -P robin
./pp.py  /home/aymeric/Big_Data/DATAPLOT/diagfired.nc -v ps -t 0.5 -z 10. -C jet -P robin -H 0.5 -B vishires

./pp.py  /home/aymeric/Big_Data/DATAPLOT/diagfired.nc -v icetot -P robin
./pp.py  /home/aymeric/Big_Data/DATAPLOT/diagfired.nc -v phisinit -P ortho -I -120.
./pp.py  /home/aymeric/Big_Data/DATAPLOT/diagfired.nc -v phisinit -P ortho -I -120. -v phisinit -I 120
./pp.py  /home/aymeric/Big_Data/DATAPLOT/diagfired.nc -v phisinit -P ortho -I -120. -v phisinit
./pp.py  /home/aymeric/Big_Data/DATAPLOT/diagfired.nc -v icetot -t 1. -P robin -I -120. -v icetot -I 120.

./pp.py  /home/aymeric/Big_Data/DATAPLOT/diagfired.nc -v icetot -t 1. -P moll
./pp.py  /home/aymeric/Big_Data/DATAPLOT/diagfired.nc -v icetot -t 1. -P moll -i u -j v
./pp.py  /home/aymeric/Big_Data/DATAPLOT/diagfired.nc -v icetot -t 1. -P moll -i u -j v -z 10.
./pp.py  /home/aymeric/Big_Data/DATAPLOT/diagfired.nc -v icetot -t 1. -P moll -i u -j v -z 10000.
./pp.py  /home/aymeric/Big_Data/DATAPLOT/diagfired.nc -v icetot -t 1. -P cyl -i u -j v -z 10000. --verbose
./pp.py  /home/aymeric/Big_Data/DATAPLOT/diagfired.nc -v icetot -t 1. -P cyl -i u -j v -z 10000. -c phisinit
./pp.py  /home/aymeric/Big_Data/DATAPLOT/diagfired.nc -v icetot -t 1. -P cyl -z 10000. -c temp

./pp.py  wrfout_d03_2024-06-09_00\:00\:00_z -v tk -c HGT -t 1 -z 0 -P lcc -i Um -j Vm --verbose
