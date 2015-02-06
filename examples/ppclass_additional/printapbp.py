#! /usr/bin/env python
from ppclass import pp


ff = "ref64_p20b_diagfi.nc"
ff = "ref80_p20b_diagfi.nc"

apbp = open("apbp.txt", "w")
prof = open("temp_profile.txt","w")

ap = pp(file=ff,x=0,y=0,var="ap").getf()
bp = pp(file=ff,x=0,y=0,var="bp").getf()
temp = pp(file=ff,x=0,y=0,t=1e10,var="temp").getf()

for nn in range(len(ap)):
   apbp.write("%12.5e%12.5e\n"%(ap[nn],bp[nn]))

for nn in range(len(temp)):
   prof.write("%12.5f \n"%(temp[nn]))

apbp.close()
prof.close()

