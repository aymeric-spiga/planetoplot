#!/u/acolmd/san0/Python/64bits/epd-7.1-2-rh5-x86_64/bin/python
if __name__ == "__main__":
    import subprocess
    import sys
    from os import system,unlink
    from glob import glob
    from optparse import OptionParser 
    from myplot import errormess
    from ic import check_file, merge_atm_sfc
    #############################
    ### Get options and variables
    parser = OptionParser()
    parser.add_option('-g','--generate',  action='store_true',dest='generate',          default=False,help='Generate netcdf files')
    parser.add_option('-c','--clean',     action='store_true',dest='clean',            default=False,help='Clean mrams_ files')
    parser.add_option('--deepclean',      action='store_true',dest='deepclean',            default=False,help='Clean everything')
    (opt,args) = parser.parse_args()
    if opt.deepclean:
       for zfile in ['converted_files','figures']:
           subprocess.call(['touch',zfile],shell=False)
           subprocess.call(['rm','-rf',zfile],shell=False)
           subprocess.call(['mkdir',zfile],shell=False)
    if opt.clean:
       subprocess.call(['touch','converted_files/mramsout_d_dummy'],shell=False)
       for f in glob ('converted_files/mramsout_d*'):
           unlink(f)
    if opt.generate:
       opath='./converted_files/'
       ipath='./vis/'
#       for grid in ['g1','g2','g3','g4']: ## INDEX 4 OF TEST CASE HOLDEN DOES NOT WORK (CDO ERROR)
       for grid in ['g1','g2','g3']:
           outfileatm='holden_ls150-atm-S-'+grid+'.nc'
           outfilesfc='holden_ls150-sfc-S-'+grid+'.nc'
           outfilecat='mramsout_d'
           infileatm='holden_ls150-atm-S-'+grid+'.ctl'
           infilesfc='holden_ls150-sfc-S-'+grid+'.ctl'
           print '- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -'
           err=check_file(opath+outfileatm)
           print '=> Begining conversion for File '+infileatm+' -> in -> '+outfileatm
           if err == 1: subprocess.call(['cdo','-f','nc','import_binary',ipath+infileatm,opath+outfileatm],shell=False)
           print '=> File '+infileatm+' -> converted in -> '+outfileatm
           err=check_file(opath+outfilesfc)
           print '=> Begining conversion for File '+infilesfc+' -> in -> '+outfilesfc
           if err == 1: subprocess.call(['cdo','-f','nc','import_binary',ipath+infilesfc,opath+outfilesfc],shell=False)
           print '=> File '+infilesfc+' -> converted in -> '+outfilesfc
           print '- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -'
           print '=> Begining concatenation for Files '+outfilesfc+' and '+outfileatm
           merge_atm_sfc(opath+outfilesfc,opath+outfileatm,opath+outfilecat)
           print '=> Files concatenated !'
           print '- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -'

#   subprocess.call(cmdstring,shell=False)

# Check if a file is present and readable
def check_file(filename):
    from os import path, access, R_OK  # W_OK for write permission. 
    PATH=filename 
    if path.exists(PATH) and path.isfile(PATH) and access(PATH, R_OK): 
       return 0  # file exists and is readable
    else: 
       return 1  # file is not present and/or readable

# Find lat lon center of the grid
def find_center(lon,lat,xlon,xlat):
    import numpy as np
    i_lon=0 ; i_lat = 0 ; i_lon2 = 0 ; i_lat2=0
    if len(lon)%2!=0:i_lon=np.floor(len(lon)/2) # no +1 because arrays start from 0
    if len(lat)%2!=0:i_lat=np.floor(len(lat)/2) # no +1 because arrays start from 0
    if i_lon == 0:
         i_lon=len(lon)/2 -1
         i_lon2=len(lon)/2
    if i_lat == 0:
         i_lat=len(lat)/2 -1
         i_lat2=len(lat)/2
    if i_lon2 == 0 and i_lat2 == 0:
         cen_lon = xlon[i_lat,i_lon]
         cen_lat = xlat[i_lat,i_lon]
    if i_lon2 != 0 and i_lat2 != 0:
         cen_lon = (xlon[i_lat,i_lon] + xlon[i_lat2,i_lon2] + xlon[i_lat,i_lon2] + xlon[i_lat2,i_lon])/4.
         cen_lat = (xlat[i_lat,i_lon] + xlat[i_lat2,i_lon2] + xlat[i_lat,i_lon2] + xlat[i_lat2,i_lon])/4.
    if i_lon2 != 0 and i_lat2 == 0:
         cen_lon = (xlon[i_lat,i_lon] + xlon[i_lat,i_lon2])/2.
         cen_lon = (xlat[i_lat,i_lon] + xlat[i_lat,i_lon2])/2.
    if i_lon2 == 0 and i_lat2 != 0:
         cen_lon = (xlon[i_lat,i_lon] + xlon[i_lat2,i_lon])/2.
         cen_lon = (xlat[i_lat,i_lon] + xlat[i_lat2,i_lon])/2.
    return cen_lon,cen_lat

# Merge atm and sfc netcdf. cdo cat does not work for that. (neither does ncecat)
def merge_atm_sfc(fsfc,fatm,fmerged):
    from os import system
    from netCDF4 import Dataset
    import numpy as np
    import subprocess

    # Calendar
    path_calendar='/san0/acolmd/SVN/trunk/MESOSCALE/LMD_MM_MARS/SIMU/calendar'
    gcm_sol=[]
    gcm_ls=[]
    mmm_date=[]
    calendar = open(path_calendar, 'r')
    calendar.readline() #get rid of the first line
    calendar_l=calendar.readlines()
    for line in calendar_l:
        s=str.split(line)
        gcm_sol=np.append(gcm_sol,np.float(s[1]))
        gcm_ls=np.append(gcm_ls,np.float(s[2]))
        mmm_date=np.append(mmm_date,s[3])
    calendar.close()
    calendar_l=0.

    # Constants
    grav=3.72

    # Grid properties
    ls_simu_start=147.  #ls of start simu
#    lt_start="07:30:00" #local time of start simu at cen_lon # read in .ctl
#   NOTE: I ASSUME THE LT_START IN .CTL FOR MRAMS IS LOCAL TIME OF START AT LONGITUDE 0
    dx=[240000.,80000.,26667,8889]
    dy=dx
   
    # Date management
    idx=np.abs(gcm_ls - ls_simu_start).argmin()
    date_start=mmm_date[idx]
    # USER: fill the date you want to analyse:
    days=[16,17,18]

    # Grid id
    if 'g1' in fatm:
        gid=1
    elif 'g2' in fatm:
        gid=2
    elif 'g3' in fatm:
        gid=3
    else:
        gid=4

    # Load datasets   
    atm=Dataset(fatm) ; sfc=Dataset(fsfc)

    # reconstruct dimensions
    lon = atm.variables['lon'] ; nlon = len(lon)
    lat = atm.variables['lat'] ; nlat = len(lat)
    lev = atm.variables['lev'] ; nlev = len(lev)
    time = atm.variables['time'] ; ntime = len(time)
    units_time=time.units[12:len(time.units)]

    # date from mrams file:
    year_mrams, month_mrams, day_mrams, hour_mrams, minute_mrams, second_mrams = date_to_time(units_time)
    lt_start = hour_mrams+':'+minute_mrams+':'+second_mrams
    # or date from calendar with initial Ls:
    year, month, day, hour, minute, second = date_to_time(date_start,lt_start)
    dt=(time[1]-time[0])*60.
    daylength=3600.*24./(dt*60.)

    i=0
    imax=np.floor(len(time)/daylength)
    while (i < imax+1):
        i=i+1
        dayslab=[int((i-1)*daylength),int(daylength*i)]
        if i==imax+1: dayslab=[int((i-1)*daylength),len(time)]
        fmerged_slab=fmerged+'0'+str(gid)+'_'+time_to_date(mmm_date[idx+i-1],lt_start)
        err=check_file(fmerged_slab)
        year, month, day, hour, minute, second = date_to_time(mmm_date[idx+i-1],lt_start)
        if int(day) in days:
           print '  -->> hyperslab : ',dayslab,' => mramsout_d0'+str(gid)+'_'+time_to_date(mmm_date[idx+i-1],lt_start)
        else: 
           err=0
           print '  -->> hyperslab : ',dayslab,' => mramsout_d0'+str(gid)+'_'+time_to_date(mmm_date[idx+i-1],lt_start)+' --> Skipped'
        if err == 1:
            # prepare merged file
            subprocess.call(['touch',fmerged_slab],shell=False)
            subprocess.call(['rm','-rf',fmerged_slab],shell=False)
            merged_file = Dataset(fmerged_slab, 'w', format='NETCDF3_64BIT')
            merged_file.description = 'MRAMS merged atm and sfc converted (netcdf) file.'
            ## Dimensions
            merged_file.createDimension('Time', None)
            merged_file.createDimension('bottom_top',nlev)
            merged_file.createDimension('south_north',nlat)
            merged_file.createDimension('west_east',nlon)
            merged_file.createDimension('bottom_top_stag',nlev+1) # very important for PHTOT in API
#            merged_file.createDimension('south_north_stag',nlat+1)
#            merged_file.createDimension('west_east_stag',nlon+1)
            ## Variables for dimensions
            times = merged_file.createVariable('Time', 'f', ('Time',))
            times[:] = time[dayslab[0]:dayslab[1]]
            altitudes = merged_file.createVariable('bottom_top', 'f', ('bottom_top',))
            altitudes[:] = lev[:]
            latitudes = merged_file.createVariable('south_north', 'f', ('south_north',))
            latitudes[:] = lat[:]
            longitudes = merged_file.createVariable('west_east', 'f', ('west_east',))
            longitudes[:] = lon[:]
            ## Set attributes
            setattr(merged_file, 'TITLE', 'OUTPUT FROM MRAMS MODEL')
            setattr(merged_file, 'START_DATE', time_to_date(mmm_date[idx],lt_start))
            setattr(merged_file, 'SIMULATION_START_DATE', time_to_date(mmm_date[idx],lt_start))
            setattr(merged_file, 'WEST-EAST_GRID_DIMENSION', nlon+1)
            setattr(merged_file, 'SOUTH-NORTH_GRID_DIMENSION', nlat+1)
            setattr(merged_file, 'BOTTOM-TOP_GRID_DIMENSION', nlev+1)
            setattr(merged_file, 'DX', np.float(dx[gid-1])) 
            setattr(merged_file, 'DY', np.float(dy[gid-1]))
            setattr(merged_file, 'DT', np.float(dt))
            setattr(merged_file, 'GRID_TYPE', 'C') #
            setattr(merged_file, 'GRID_ID', gid) 
            setattr(merged_file, 'PARENT_ID', gid-1) 
            setattr(merged_file, 'MAP_PROJ', 2) #
            setattr(merged_file, 'JULYR', year) 
            setattr(merged_file, 'JULDAY', gcm_sol[idx+i-1])
            setattr(merged_file, 'WEST-EAST_PATCH_START_UNSTAG' , 1)
            setattr(merged_file, 'WEST-EAST_PATCH_END_UNSTAG' , nlon)
            setattr(merged_file, 'WEST-EAST_PATCH_START_STAG' , 1)
            setattr(merged_file, 'WEST-EAST_PATCH_END_STAG' , nlon) ###
            setattr(merged_file, 'SOUTH-NORTH_PATCH_START_UNSTAG' , 1)
            setattr(merged_file, 'SOUTH-NORTH_PATCH_END_UNSTAG' , nlat)
            setattr(merged_file, 'SOUTH-NORTH_PATCH_START_STAG' , 1)
            setattr(merged_file, 'SOUTH-NORTH_PATCH_END_STAG' , nlat) ###
            setattr(merged_file, 'BOTTOM-TOP_PATCH_START_UNSTAG' , 1)
            setattr(merged_file, 'BOTTOM-TOP_PATCH_END_UNSTAG' , nlev)
            setattr(merged_file, 'BOTTOM-TOP_PATCH_START_STAG' , 1)
            setattr(merged_file, 'BOTTOM-TOP_PATCH_END_STAG' , nlev+1)
            ## Set grid variables
            xlon = merged_file.createVariable('XLONG', 'f', ('Time','south_north','west_east',))
            xlat = merged_file.createVariable('XLAT', 'f', ('Time','south_north','west_east',))
            print '    glon -> XLON' ; xlon[:,:,:] = sfc.variables['glon'][dayslab[0]:dayslab[1],:,:]
            print '    glat -> XLAT' ; xlat[:,:,:] = sfc.variables['glat'][dayslab[0]:dayslab[1],:,:]
            cen_lon, cen_lat = find_center(lon,lat,xlon[0,:,:],xlat[0,:,:])
            setattr(merged_file, 'CEN_LAT', cen_lat)
            setattr(merged_file, 'CEN_LON', cen_lon)
            ## Renaming
            print '    topo -> HGT' ; hgt = merged_file.createVariable('HGT', 'f', ('Time','south_north','west_east',))
            hgt[:,:,:] =  sfc.variables['topo'][dayslab[0]:dayslab[1],:,:]
            print '    tempk -> tk' ; tk = merged_file.createVariable('tk', 'f', ('Time','bottom_top','south_north','west_east',))
            tk[:,:,:,:] = atm.variables['tempk'][dayslab[0]:dayslab[1],:,:,:]
            print '    press -> PTOT' ; ptot = merged_file.createVariable('PTOT', 'f', ('Time','bottom_top','south_north','west_east',))
            ptot[:,:,:,:] = atm.variables['press'][dayslab[0]:dayslab[1],:,:,:]
            u = merged_file.createVariable('U', 'f', ('Time','bottom_top','south_north','west_east',))
            v = merged_file.createVariable('V', 'f', ('Time','bottom_top','south_north','west_east',))
            w = merged_file.createVariable('W', 'f', ('Time','bottom_top','south_north','west_east',))
            print '    u_avg -> U' ; u[:,:,:,:] = atm.variables['u_avg'][dayslab[0]:dayslab[1],:,:,:]
            print '    v_avg -> V' ; v[:,:,:,:] = atm.variables['v_avg'][dayslab[0]:dayslab[1],:,:,:]
            print '    w_avg -> W' ; w[:,:,:,:] = atm.variables['w_avg'][dayslab[0]:dayslab[1],:,:,:]
            tke = merged_file.createVariable('TKE', 'f', ('Time','bottom_top','south_north','west_east',))
            print '    sgs_tke -> TKE' ; tke[:,:,:,:] = atm.variables['sgs_tke'][dayslab[0]:dayslab[1],:,:,:]
            tau_dust = merged_file.createVariable('TAU_DUST', 'f', ('Time','south_north','west_east',))
            print '    dust_sfc_od_vis -> TAU_DUST' ; tau_dust[:,:,:] = atm.variables['dust_sfc_od_vis'][dayslab[0]:dayslab[1],:,:]
            z_lyrmid_agl = merged_file.createVariable('z_lyrmid_agl', 'f', ('Time','bottom_top','south_north','west_east',))
            print '    z_lyrmid_agl -> z_lyrmid_agl' ; z_lyrmid_agl[:,:,:,:] = atm.variables['z_lyrmid_agl'][dayslab[0]:dayslab[1],:,:,:]
            # Compute PHTOT for API
            print '    Computing PHTOT'
            phtot = merged_file.createVariable('PHTOT', 'f', ('Time','bottom_top_stag','south_north','west_east',))
            phtot[:,0:nlev,:,:]=z_lyrmid_agl[:,:,:,:]*grav ; phtot[:,nlev,:,:]=0.
            for l in np.arange(nlev):
                phtot[:,l,:,:]=phtot[:,l,:,:] + grav*hgt[:,:,:]

            merged_file.close()

def date_to_time(date,ltstart=None):
    year=date[0:4]
    month=date[5:7]
    day=date[8:10]
    hour=date[11:13]
    minute=date[14:16]
    second=date[17:19]
    if ltstart is not None:
       hour=ltstart[0:2]
       minute=ltstart[3:5]
       second=ltstart[6:8]
    return year,month,day,hour,minute,second

def time_to_date(date,ltstart=None):
    if ltstart is not None: return date[0:11]+ltstart
    else: return date
