#! /usr/bin/env python
from ppclass import pp
from    netCDF4               import    Dataset
from	numpy		      import	*
import  numpy                 as        np
import  matplotlib.pyplot     as        mpl

## Author: Tanguy Bertrand

## Exemple d'histogramme
#
# INPUT : 
#
# file
# nb_dataset: number of dataset that is number of superposed histograms
# Time range
# altitude range
# longitude range
# latitude range
# step: gap on xaxis between two ticks. Ex: step = 1 leads to one bar for each unit of the selected variable  
# variable : u, v, w, icetot, tsurf ...

############################
filename="diagfi1.nc"
nb_dataset=3
tint=[["0.25,1.5"],["0.25,5.5"],["0.25,10.5"]] #Time must be as written in the input file
zint=[["0.1,10"],["0.1,10"],["0.1,10"]] #alt in km
xarea="-180,180"
yarea="-90,90"
step=5  
var="u" #variable
############################

x=np.zeros(nb_dataset,dtype='object') #object of all datasets

bornemax=0. # initialisation boundary values
bornemin=0.

for i in range(nb_dataset):

	myvar = pp(file=filename,var=var,t=tint[i],z=zint[i],x=xarea,y=yarea,compute="nothing").getf()  # get data to be changed according to selected variable
	data=np.ravel(myvar)
	        
	x[i]=np.array(data)
	
        #upper lower bounds:    
        maxval=np.amax(myvar)       
        romax=round(maxval,0)
        if abs(romax/maxval) < 1: romax=romax+1	

        minval=np.amin(myvar)       
        romin=round(minval,0)
        if abs(romin/minval) < 1: romin=romin-1	

	bornemax=np.amax([romax,bornemax])
        bornemin=np.amin([romin,bornemin])
        #
        
        print('minval,maxval=',minval,romax)
        print('romin,romax=',romin,romax)
	print('bornemin,bornemax=',bornemin,bornemax)
        
# bins definition:
bins=np.arange(bornemin,bornemax+1,step)



vect=[x[0],x[1],x[2]] #vecteur of datasets to be changed according to the number of desired datasets

## PLOT
# legend:
zelab=['time: '+str(tint[0])+' - alt: '+str(zint[0])+' km','time: '+str(tint[1])+' - alt: '+str(zint[1])+' km','time: '+str(tint[2])+' - alt: '+str(zint[2])+' km']
# Histogram : normed, probability
mpl.hist(vect,bins,normed=1,histtype='bar',align='mid',rwidth=0.8,label=zelab)
mpl.legend(prop={'size':20})
mpl.title('U wind distribution',fontsize=20)
mpl.xlabel('Horizontal Wind (m.s-1)',fontsize=20)
mpl.ylabel('Probability',fontsize=20)
mpl.xticks(bins,fontsize=16)
mpl.yticks(fontsize=16)
mpl.grid(True)

mpl.show()


mpl.figure(5)
 
