#! /usr/bin/env python
from ppclass import pp

########
## ppclass allows for getting fields very easily from a netcdf file
## ... this is contained in the "f" attributes of the pp object
## ... and the "l" attributes of the pp object contains useful information
## ... two methods getf() and getfl() are helpful shortcuts
###########################################

## 1 --> an unique and straightforward request
##   --> very easy ! see also minimal_field.py
##############################################
icetot = pp(file="/home/aymeric/Big_Data/DATAPLOT/diagfired.nc",var="icetot",t=0.5).getf()
print("icetot", icetot[10,10])

## 2 --> simple multiple request, no labelling
##############################################
test = pp(file="/home/aymeric/Big_Data/DATAPLOT/diagfired.nc",t=0.5)
test.var = ["mtot","icetot"]
allf = test.getf() # or allf = test.get().f
##
mtot = allf[0]
icetot = allf[1]
print("mtot", mtot[10,10])
print("icetot", icetot[10,10])

## 3 --> complex multiple requests and labelling
################################################
test = pp(file="/home/aymeric/Big_Data/DATAPLOT/diagfired.nc")
test.var = ["mtot","icetot"]
test.t = [0.4,0.5]
allf,lab = test.getfl()
##
icetot04 = allf[lab.index("_v=icetot_t=0.4_")]
mtot04 = allf[lab.index("_v=mtot_t=0.4_")]
icetot05 = allf[lab.index("_v=icetot_t=0.5_")]
mtot05 = allf[lab.index("_v=mtot_t=0.5_")]
print("mtot04", mtot04[10,10])
print("icetot04", icetot04[10,10])
print("mtot05", mtot05[10,10])
print("icetot05", icetot05[10,10])

## 4 --> an example of complete labelling
## .... a rather unlikely example ....
## .... but shows label ordering ....
#########################################
test = pp()
test.file = ["/home/aymeric/Big_Data/DATAPLOT/diagfired.nc","/home/aymeric/Big_Data/DATAPLOT/diagfired.nc"]
test.var = ["u","v"]
test.x = [10.,20.]
test.y = [10.,20.]
test.z = [10.,20.]
test.t = [0.4,0.5]
print("... please wait. this one is a bit stupid...")
allf,lab = test.getfl()
## note label ordering: file -- var -- x -- y -- z -- t
l1 = "_f=#2_v=u_x=10.0_y=10.0_z=20.0_t=0.4_"
l2 = "_f=#2_v=v_x=10.0_y=10.0_z=20.0_t=0.4_"
u_example = allf[lab.index(l1)]
v_example = allf[lab.index(l2)]
print(l1, u_example)
print(l2, v_example)
