#! /usr/bin/env python
from ppclass import pp

les = pp(file="/home/aymeric/Big_Data/LES_dd/psfc_f18.nc",var="PSFC",y=["50,250"],x=["50","100"])
les.verbose = True
les.getplot()

les.x = None
les.y = ["200","250"]
les.t = "10"
les.getplot()

les.x = None
les.y = None
les.t = ["10"]
#les.stridex = 3
#les.stridey = 3
les.filename = "les_psfc"
les.getplot()
