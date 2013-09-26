#! /usr/bin/env python
from ppclass import pp
tes = pp()
tes.file = [\
'/home/aymeric/Big_Data/tes/TES.MappedClimatology.nadir.MY24.nc',\
'/home/aymeric/Big_Data/tes/TES.MappedClimatology.nadir.MY25.nc',\
'/home/aymeric/Big_Data/tes/TES.MappedClimatology.nadir.MY26.nc']
tes.var='T_nadir_day'
tes.x = "-180,180"
tes.y = 60.
tes.z = 50.
tes.superpose = True
tes.getdefineplot()
tes.p[0].legend = "MY24"
tes.p[1].legend = "MY25"
tes.p[2].legend = "MY26"
tes.makeplot()
