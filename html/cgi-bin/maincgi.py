#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import cgi, cgitb 

import sys
sys.path.insert(0, "../planetoplot/modules")
import ppplot
import ppclass




########################################

import numpy as np
xx = np.arange(25)
yy = 3.*xx

fig = ppplot.figuref(x=8,y=6)

pl = ppplot.plot1d()
pl.fig = fig # have to send to figure
pl.f = yy
pl.x = xx
pl.make()


######################################## more sophisticated example
## RETRIEVE DATA
#from ppclass import pp
#fifi = "/home/aspiga/soft/mcd_python/minimal_server/cgi-bin/wrfout_d01_2024-10-04_06z00z00_zabg"
#ff,xx,yy,zz,tt = pp(file=fifi,var="HGT",z=0,t=0).getfd()
#xx = pp(file=fifi,var="XLONG",z=0,t=0).getf()
#yy = pp(file=fifi,var="XLAT",z=0,t=0).getf()
#uu = pp(file=fifi,var="Um",z=0,t=0).getf()
#vv = pp(file=fifi,var="Vm",z=0,t=0).getf()
#
## PLOT
#pl = ppplot.plot2d()
#pl.fig = fig # have to send to figure
#pl.f = ff
#pl.x = xx
#pl.y = yy
#pl.vx = uu
#pl.vy = vv
#pl.legend = "yorgl"
#pl.marker = None
#pl.nyticks = 20
#pl.ylabel = "YAARGL"
#pl.proj = "laea"
#pl.make()
########################################

# create figure
ppplot.sendagg(fig,filename='webapp.png', dpi=150)

# for debugging in web browser
cgitb.enable()

## Create instance of FieldStorage 
#form = cgi.FieldStorage() 

##### NOW WRITE THE HTML PAGE TO USER
print "Content-type:text/html;charset=utf-8\n"
print     #Apache needs a space after content-type
header="""<html><head><title>Mars Climate Database: The Web Interface</title></head><body>"""
print header
print "THIS IS A TEST!"
print "<img src='../webapp.png'><br />"
bottom = "</body></html>"
print bottom

