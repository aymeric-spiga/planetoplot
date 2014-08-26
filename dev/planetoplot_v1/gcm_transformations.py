### Wrapper for all gcm utilities.
### 
### A. Colaitis -- LMD -- 08/11/2011 <-- Added ZRECAST for altitude mode, and MCS and TES modes for Pressure
### A. Colaitis -- LMD -- 08/11/2011 <-- Added HRECAST for MCS and TES predefined grids
###                                      Added wrapper to streamfunction.e
###                                      Added wrapper to localtime.e

def call_zrecast (  interp_mode   = '4', \
                    input_name      = None, \
                    fields  = 'all', \
                    limites = None, \
                    predefined = None, \
                    verbose = True):

    import numpy as np
    from myplot import separatenames
    from os import system
    pressure_axis_tes=[1658.152,1291.37,1005.72,783.2555,610.,475.0685,369.9837,288.1436,224.4065,174.7679,136.1094,106.0021,82.55452,64.29353,50.07185,38.99599,30.37011,23.65227,18.4204,14.34583,11.17254]
    pressure_axis_mcs=[1878.9, 1658.2, 1463.3, 1291.4, 1139.6, 1005.7, 887.54, 783.26,691.22, 610., 538.32, 475.07, 419.25, 369.98, 326.51, 288.14, 254.29,224.41, 198.04, 174.77, 154.23, 136.11, 120.12, 106., 93.547, 82.555,72.854, 64.294, 56.739, 50.072, 44.188, 38.996, 34.414, 30.37, 26.802,23.652, 20.873, 18.42, 16.256, 14.346, 12.66, 11.173, 9.8597, 8.7012,7.6788, 6.7765, 5.9802, 5.2775, 4.6574, 4.1101, 3.6272, 3.201, 2.8249,2.4929, 2.2, 1.9415, 1.7134, 1.512, 1.3344, 1.1776, 1.0392, 0.9171,0.80934, 0.71424, 0.63031, 0.55625, 0.49089, 0.43321, 0.3823, 0.33738,0.29774, 0.26275, 0.23188, 0.20463, 0.18059, 0.15937, 0.14064, 0.12412,0.10953, 0.096661, 0.085303, 0.07528, 0.066434, 0.058628, 0.051739,0.04566, 0.040294, 0.03556, 0.031381, 0.027694, 0.02444, 0.021568,0.019034, 0.016797, 0.014824, 0.013082, 0.011545, 0.010188, 0.0089909,0.0079345, 0.0070021, 0.0061794, 0.0054533, 0.0048125, 0.004247]
    system("rm -f zrecast.auto.def")
    system("touch zrecast.auto.def")
    indicefile=0
    outputfilename=""

    for zfile in input_name:
        f = open('zrecast.auto.def', 'w')
        f.write(zfile+"\n")
        for zvar in fields:
             f.write(zvar+"\n")
        f.write("\n")
        if interp_mode == 4:
             append="_S"
             f.write("3"+"\n")
             f.write("yes"+"\n")
             f.write("0 150000"+"\n")
        elif interp_mode == 2:
             append="_P"
             if predefined in ["TES","tes"]:
                print "Using TES pressure levels"
                f.write("1"+"\n")
                f.write("no"+"\n")
                f.write(str(len(pressure_axis_tes))+"\n")
                for zp in pressure_axis_tes:
                    f.write(str(zp)+"\n")
                f.write("\n")
             elif predefined in ["MCS","mcs"]:
                print "Using MCS pressure levels"
                f.write("1"+"\n")
                f.write("no"+"\n")
                f.write(str(len(pressure_axis_mcs))+"\n")
                for zp in pressure_axis_mcs:
                    f.write(str(zp)+"\n")
                f.write("\n")
             else:
                f.write("1"+"\n")
                f.write("yes"+"\n")
                if limites[0] != -9999.:  f.write(str(limites[0])+" "+str(limites[-1])+"\n")
                else:                     f.write("370 0.1"+"\n")
                #f.write("1000000 100"+"\n")
                f.write("20"+"\n")
        else:
             print "zrecast interp option unsupported for now. Exiting."
             exit()
        f.close()
        if verbose:
           system("zrecast.e < zrecast.auto.def")
        else:
           system("zrecast.e < zrecast.auto.def > /dev/null")

        if indicefile is 0:
              outputfilename=input_name[indicefile][0:len(input_name[0])-3]+append+".nc"
        else:
              outputfilename=outputfilename+","+input_name[indicefile][0:len(input_name[0])-3]+append+".nc"
        indicefile = indicefile +1

        for i in range(len(outputfilename)):
            output_name = separatenames(outputfilename)

    return output_name

def call_hrecast (  input_name = None, \
                    fields  = 'all', \
                    predefined = None,\
                    verbose = True):

    import numpy as np
    from myplot import separatenames
    from os import system
    # Beware that latitude axis in mcs and tes files is reversed
    # When comparing to MCS or TES data, one has to reverse the latitude axis and data along latitude axis.
    # This is for exemple automatically done in mcs.py.
    # This could also be done in a future tes.py.
    latitudes_reversed_tes_mcs=[88.5, 85.5, 82.5, 79.5, 76.5, 73.5, 70.5, 67.5, 64.5, 61.5, 58.5, 55.5, 52.5, 49.5, 46.5, 43.5, 40.5, 37.5, 34.5, 31.5, 28.5, 25.5, 22.5, 19.5, 16.5, 13.5, 10.5, 7.5, 4.5, 1.5, -1.5, -4.5, -7.5, -10.5, -13.5, -16.5, -19.5, -22.5, -25.5, -28.5, -31.5, -34.5, -37.5, -40.5, -43.5, -46.5, -49.5, -52.5, -55.5, -58.5, -61.5, -64.5, -67.5, -70.5, -73.5, -76.5, -79.5, -82.5, -85.5, -88.5]
    longitudes_tes_mcs=[-176.25, -168.75, -161.25, -153.75, -146.25, -138.75, -131.25,-123.75, -116.25, -108.75, -101.25, -93.75, -86.25, -78.75, -71.25, -63.75, -56.25, -48.75, -41.25, -33.75, -26.25, -18.75, -11.25, -3.75, 3.75, 11.25, 18.75, 26.25, 33.75, 41.25, 48.75, 56.25, 63.75, 71.25, 78.75, 86.25, 93.75, 101.25, 108.75, 116.25, 123.75, 131.25, 138.75, 146.25, 153.75, 161.25, 168.75, 176.25]
    system("rm -f hrecast.auto.def")
    system("touch hrecast.auto.def")
    indicefile=0
    outputfilename=""
    append="_h"
    f = open('hrecast.auto.def', 'w')
    for zfile in input_name:
        f.write(zfile+"\n")
        for zvar in fields:
             f.write(zvar+"\n")
        f.write("\n")
    if predefined in ["MCS","mcs","TES","tes"]:
        f.write(str(len(longitudes_tes_mcs))+"\n")
        for lon in longitudes_tes_mcs:
            f.write(str(lon)+"\n")
        f.write(str(len(latitudes_reversed_tes_mcs))+"\n")
        for lat in latitudes_reversed_tes_mcs:
            f.write(str(lat)+"\n")
        f.close()
        if verbose:
           system("hrecast.e < hrecast.auto.def")
        else:
           system("hrecast.e < hrecast.auto.def > /dev/null")
           
        if indicefile is 0:
              outputfilename=input_name[indicefile][0:len(input_name[0])-3]+append+".nc"
        else:
              outputfilename=outputfilename+","+input_name[indicefile][0:len(input_name[0])-3]+append+".nc"
        indicefile = indicefile +1

        for i in range(len(outputfilename)):
            output_name = separatenames(outputfilename)

    return output_name

def call_localtime (input_name = None, \
                    fields  = 'all', \
                    times = None, \
                    verbose = True):

    import numpy as np
    from myplot import separatenames
    from os import system
    system("rm -f localtime.auto.def")
    system("touch localtime.auto.def")
    outputfilename=""
    append="_LT"
    output_name = []
    #print input_name
    for zfile in input_name:
        print 'LOCALTIME FOR FILE: '+zfile
        f = open('localtime.auto.def', 'w')
        #print zfile
        f.write(zfile+"\n")
        for zvar in fields:
           f.write(zvar+"\n")
        f.write("\n")
        f.write(str(len(times))+"\n")
        for t in times:
           f.write(str(t)+"\n")
        f.close()
        if verbose:
          system("localtime.e < localtime.auto.def")
        else:
          system("localtime.e < localtime.auto.def >  /dev/null ")
        #print input_name
        zeoutput_name = zfile.replace('.nc','_LT'+str(times[0])+'.nc')
    
        system("mv -f "+zfile.replace('.nc','_LT.nc')+" "+zeoutput_name)
        
        output_name = np.append(output_name,zeoutput_name)

    return output_name
    

def call_streamfunction (  input_name = None):
    import numpy as np
    from myplot import separatenames
    from os import system
    system("rm -f streamfunction.auto.def")
    system("touch streamfunction.auto.def")
    indicefile=0
    outputfilename=""
    append="_stream"
    f = open('streamfunction.auto.def', 'w')
    for zfile in input_name:
        f.write(zfile+"\n")
    f.close()
    system("streamfunction.e < streamfunction.auto.def")

    if indicefile is 0:
          outputfilename=input_name[indicefile][0:len(input_name[0])-3]+append+".nc"
    else:
          outputfilename=outputfilename+","+input_name[indicefile][0:len(input_name[0])-3]+append+".nc"
    indicefile = indicefile +1

    for i in range(len(outputfilename)):
        output_name = separatenames(outputfilename)

    return output_name


def call_concat (  input_name = None, \
                    fields  = 'all', \
                    begin = 0, \
                    verbose = True):

    import numpy as np
    from os import system
    system("rm -f concatnc.auto.def")
    system("touch concatnc.auto.def")
        
    f = open('concatnc.auto.def', 'w')
    for zfile in input_name:
        f.write(zfile+"\n")
    
    f.write("\n")
    f.write(str(begin)+"\n")
    f.write("sol"+"\n")
    
    for zvar in fields:
        f.write(zvar+"\n")
    f.write("\n")

    f.close()
    if verbose:
       system("concatnc.e < concatnc.auto.def")
    else:
       system("concatnc.e < concatnc.auto.def > /dev/null")
    
    
    #output_name = 'concat' + input_name[0].replace('diagfi','').replace('.nc','') \
    #              + '_' + input_name[len(input_name)-1].replace('diagfi','').replace('.nc','') + '.nc'

    #system('mv -f concat.nc '+output_name)

    output_name = 'concat.nc'

    return [output_name]



def call_lslin ( input_name = None, \
                 begin = 0,\
                 timestep  = None,\
                 isconcat = True, \
                 verbose = True):

    import numpy as np
    from os import system
    system("rm -f lslin.auto.def")
    system("touch lslin.auto.def")
        
    f = open('lslin.auto.def', 'w')
    
    f.write(input_name[0]+"\n")
    if isconcat:
       f.write("y"+"\n")
    else:
       f.write("n"+"\n")

    f.write(str(begin)+"\n")
    if timestep is None:
       f.write("y"+"\n")
    else:
       f.write("n"+"\n")
       f.write(str(timestep)+"\n")
       
       
    f.close()
    if verbose:
       system("lslin.e < lslin.auto.def")
    else:
       system("lslin.e < lslin.auto.def > /dev/null")
    
    output_name = input_name[0].replace('.nc','_Ls.nc')

    return [output_name]




