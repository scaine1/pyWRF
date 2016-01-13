#module: pyWRF
"""	This is a python  package for those iterested in working with 
	the Weather Research and Forecasting model output. 

	Users are required to add PYWRFPATH to their .profile
	For example:

	export PYWRFPATH=/path/to/pyWRF

	you should also add PYWRFPATH to your PYTHONPATH so you can import 
	pyWRF from any directory you want. i.e. 

	export PYTHONPATH=$PYTHONPATH:$PYWRFPATH

	The pyWRF module comes with a library directory which contains 
	information required for plotting,microphysics calculations etc.
	DO NOT REMOVE THIS DIRECTORY.

	Example python scripts can be found in the example_scripts directory.
	"""
import sys
import os
PYWRFPATH=os.environ['PYWRFPATH']
sys.path.append(PYWRFPATH+'/library/')
import Nio
import wrf_user_unstagger
import perturbation_variables
import numpy as n
import pylab as pl
import scipy as s
#from mpl_toolkits.basemap import Basemap, shiftgrid
from mpl_toolkits import basemap
import coltbls
import dircache

class calc_vars():
	""" This is where routines to calcalate variables should be stored. 
	If you add something please document the innput and output variables
        """
	def __init__(self):
	    pass


	def compute_lat_lon_nmm(self,var):
	    """convert lat or lon in degrees (from radians)
	    Inputs:
	    -------
	    GLAT or GLON from NMM core

	    Returns:
	    -------
	    XLAT or XLONG in degrees

	    """

	    if self.wrf_core == 'NMM':
		print 'converting your input array from radians to degrees'
		outvar=var*57.2957795
	    else:
		print 'you are not using the nmm core and should not need this, returning input variable'
		outvar=var

	    return outvar



	def compute_tk(self,pressure=None,theta=None):
	    """Compute temperature in kelvin
	    Inputs:
	    -------
	    Pressure  (pa I think)
	    Potential Temperature  (Kelvin)
	    If called without inputs, function will try and get the variables for you.
	    
	    Returns:
	    -------
	    Temperature in Kelvin
	    """
	    if (pressure == None):	    
		try:
		   pressure=self.variable_dict['PRES']
		except:
		   pressure=self.get_var('PRES')
	
	    if (theta == None):	    
		try:
		   theta=self.variable_dict['THETA']
		except:
		   theta=self.get_var('THETA')

	    pressure=pressure*100.
	    print 'calculating temperature in kelvin\n'

	    p1000mb=100000.
	    r_d=287. 
	    cp = (7/2.)*r_d
	
	    pi=(pressure/p1000mb)**(r_d/cp)
	    tk=pi*theta
	    return tk
	
	def compute_rh(self,qvapor=None,pressure=None,temperature=None):
	    """Compute the relative humidity
	    Inputs:
	    ------
	    Qvapor (kg/kg)
	    pressure (pa)
	    Temperature (K) ---> This should be actual temperature not potential temperature I think
	
	    ***if called without input args function will try to get the variables for you

	    Returns:
	    --------
	    Relative Humidity
	
	    """

	    if (qvapor == None):
		try:
		   qvapor=self.variable_dict['QVAPOR']
		except:
		   qvapor=self.get_var('QVAPOR')

	    if (pressure == None):	    
		try:
		   pressure=self.variable_dict['PRES']
		except:
		   pressure=self.get_var('PRES')
	
	    if (temperature == None):	    
		try:
		   temperature=self.variable_dict['TEMP']
		except:
		   temperature=self.get_var('TEMP')


	    pressure=pressure*100.
	    print 'calculating relative humidity\n'
	    svp1=0.6112
	    svp2=17.67
	    svp3=29.65
	    svpt0=273.15
	
	    r_d=287.
	    r_v=461.6
	    ep_2=r_d/r_v
	    ep_3=0.622
	
	    es	=10.*svp1*n.exp(svp2*(temperature-svpt0)/(temperature-svp3))
	    qvs	=ep_3*es/(0.01 * pressure- (1-ep_3)*es)
	
	    qvap_on_qvs=qvapor/qvs
	
	    inner_part=n.where(qvap_on_qvs >= 1.0, 1.0,qvap_on_qvs)
	    rh=100*n.where(inner_part <=0,0,inner_part)
	    return rh
	
	def compute_vtmk(self,qvapor=None,tmk=None):
	    """calculate the virtual temperature
	    Inputs:
	    -------
	    None, will obtain variables for you.

	    Optional Inputs:
	    -------	    	    
	    QVAPOR (kg/kg):- use qvapor=
	    Temperature (kevlin i think): use tmk=
	    
	    Returns:
	    -------
	    Virtual emperature in Kelvin

	    Exammple:
	    -------
	    wf=pyWRF.wrf_file('wrfoutxxxx')
	    wf.compute_vtmk() 

	    or 

	    qvp=wf.get_var('QVAPOR')
	    tmk=wf.get_var('TEMP')
	    wf.compute_vtmk(qvpaor=qvp, tmk=tmk)

	
	    Notes:
	    ------

	    doing this quick, could be errors,check later simon

	
	    """
	    if (qvapor == None):
		try:
		   qvapor=self.variable_dict['QVAPOR']
		except:
		   qvapor=self.get_var('QVAPOR')

	    if (tmk == None):	    
		try:
		   tmk=self.variable_dict['TEMP']
		except:
		   tmk=self.get_var('TEMP')

	    print 'calculating virtual temperature in kelvin\n'
    
	    qvapor_temp=qvapor*0.001
	    eps=0.622

	    vtmk=tmk * (eps+qvapor_temp)/(eps*(1.+qvapor_temp))

	
	    return vtmk

	
	def compute_td(self,qvapor=None,pressure=None):
	    """calculate the dew point temperature
	    Inputs:
	    -------
	    QVAPOR (kg/kg)
	    Pressure  (mb I think)
	    
	    Returns:
	    -------
	    Dew point temperature in Kelvin
	
	
	    Notes:
	    ------
	    Need to modify this so that it does not matter if pressure is in pa or hpa
	
	    """
	    if (qvapor == None):
		try:
		   qvapor=self.variable_dict['QVAPOR']
		except:
		   qvapor=self.get_var('QVAPOR')

	    if (pressure == None):	    
		try:
		   pressure=self.variable_dict['PRES']
		except:
		   pressure=self.get_var('PRES')

	    print 'calculating dew point temperature in kelvin\n'
#	    pressure=pressure/100.				#I think pressure needs to be in mb
	    qv 		= n.where(qvapor <0.,0.,qvapor)
	    vapor_pres  = (qv*pressure)/(0.622+qv)	
	    vapor_pres  = n.where(vapor_pres < 0.001, 0.001,vapor_pres)	#avoid problems near zero
	    td		= (243.5*n.log(vapor_pres)-440.8)/(19.48-n.log(vapor_pres))
	
	    return td
	
	
	
	def compute_mslp(self):
	    """calculate mslp using the hypsometric equation in the following form:
	    MSLP =Psfc*exp^g*dz/(R*T)

	    Inputs:
	    -------
	    Z (distance in meters from the surface to mean sea-level)
	    g is gravity
	    R is the dry gas constant
	    T is the mean layer temperature from the surface to sea-level

	    
	    Returns:
	    -------
	    Mean Sea Level Pressure  (in mb)
	
	
	    Notes:
	    ------
	    """

	    if (self.wrf_core == 'ARW') or (self.wrf_core == 'NMM'):
#		try:
#		    psfc=self.variable_dict['PSFC']/100.
#		except:
#		    psfc=self.get_var('PSFC')/100.
		try:
		    pres=self.variable_dict['PRES']
		except:
		    pres=self.get_var('PRES')


		try:
		    Z=self.variable_dict['Z']
		except:
		    Z=self.get_var('Z')

		try:
		    T=self.variable_dict['TEMP']
		except:
		    T=self.get_var('TEMP')

	


	    rgas=287.04
	    grav=9.81

	    exponent=(grav*Z[0,0,:,:])/(rgas*T[0,0,:,:])
#	    mslp = psfc[0,:,:]*n.exp(exponent)
	    mslp = pres[0,0,:,:]*n.exp(exponent)


	    return mslp
	
	
	def compute_sph(self,qvapor=None):
	    """calculate the specific humidity
	    Inputs:
	    -------
	    QVAPOR (kg/kg)

	    Returns:
	    -------
	    Specific humidity in kg/kg

	    """

	    print 'calculating specific humidity in kg/kg'
	    if (qvapor == None):
		try:
		   qvapor=self.variable_dict['QVAPOR']
		except:
		   qvapor=self.get_var('QVAPOR')

	    sph=qvapor/(1+qvapor)

	    return sph



	def compute_sph(self,qvapor=None):
	    """calculate the specific humidity
	    Inputs:
	    -------
	    QVAPOR (kg/kg)

	    Returns:
	    -------
	    Specific humidity in kg/kg

	    """

	    print 'calculating specific humidity in kg/kg'
	    if (qvapor == None):
		try:
		   qvapor=self.variable_dict['QVAPOR']
		except:
		   qvapor=self.get_var('QVAPOR')

	    sph=qvapor/(1+qvapor)

	    return sph

	
	def get_ij_lat_long(self,lat_array,long_array,user_lat,user_lon):
	    """This function is designed to return the closest grid point to your chosen lat, long.
	    Inputs:
	    ------
		lat_array: (i.e. XLAT)
		lon_array: (i.e. XLONG)
		user_lat:  The lat you want to find
		user_lon: The long you want to find.
	    Returns:
	    --------
		i,j coordinate of the WRF grid point closest to your chosen lat,long.
		where i corresponds to ns dimension and j corresponds to the we dimenion (I think)
		
	
	    Notes:
	    ------
		lat_array and lon_array currently needs to be 3 dimensional.
		the first dimension is time, which wrf creates for no good 
		reason unles you have a moving nest.
		will try and make this script smarter so that it doesnt die
		if xlat and xlon is 2d.
		"""


            if (user_lat < n.min(lat_array)) | (user_lat > n.max(lat_array)) |  (user_lon < n.min(long_array)) | (user_lon > n.max(long_array)):
                print 'point outside array bounds, skipping'
                return n.nan,n.nan


	    del_lat=del_lon=2.0   #2 degrees! this should easily be enough for all possible cases
	    i=j=[0,0]

	    loop=1
	    div_factor_lat=2.0
	    div_factor_lon=2.0
	    loop_count=0.
	    while loop == 1:
		while (len(i) > 1) | ( len(j) >1):
		    upper_lat=user_lat+del_lat
		    lower_lat=user_lat-del_lat
		    upper_lon=user_lon+del_lon
		    lower_lon=user_lon-del_lon


		    if (len(n.shape(lat_array)) == 3):
			lat_ij=n.where((lat_array[0,:,:] < upper_lat) & (lat_array[0,:,:] >= lower_lat),1,0)
			lon_ij=n.where((long_array[0,:,:] > lower_lon) & (long_array[0,:,:] <= upper_lon),1,0)

		    elif (len(n.shape(lat_array)) == 2):
			lat_ij=n.where((lat_array[:,:] < upper_lat) & (lat_array[:,:] >= lower_lat),1,0)
			lon_ij=n.where((long_array[:,:] > lower_lon) & (long_array[:,:] <= upper_lon),1,0)
		    else:
			pass


		    i,j = n.where((lat_ij == 1) & (lon_ij==1))

		    del_lat_old=del_lat
		    del_lon_old=del_lon
		    del_lat=del_lat/div_factor_lat
		    del_lon=del_lon/div_factor_lon

#		    print 'len i j ', len(i),len(j)
#
#		    print 'del lat',del_lat, 'del lon',del_lon
#
#
#		    print 'user lat', user_lat, 'user lon', user_lon
#
#		    print 'lower lat', lower_lat, 'upper lat', upper_lat
#
#		    print 'lower lon', lower_lon, 'upper lon', upper_lon
#
#		    print n.sum(lat_ij), n.sum(lon_ij),loop_count
		    if (loop_count >= 2000):
			print 'loop greater than 2000'
			try:
			    #print i[0],j[0]
			    #print lat_array[0,j[0],i[0]]
			    #print lon_array[0,j[0],i[0]]
                            print 'error is', user_lat - i[0], user_long - j[0]
			    return i[0], j[0]
			except:
			    print 'passing'
			    pass


#
#		    print 'div factor lat', div_factor_lat
#		    print 'div factor lon', div_factor_lon
#		    print ''

		    loop_count+=1
#		    if loop_count == 70:
#			raw_input('pause')


		if (n.shape(i)[0] == 0):
		    del_lat=del_lat_old*div_factor_lat
		    del_lat=del_lat*(99/100.)

#		    div_factor_lat=div_factor_lat/1.001
		    i=[0,0]
		if (n.shape(j)[0] == 0): 
		    del_lon=del_lon_old*div_factor_lon
		    del_lon=del_lon*(99/100.)
#		    div_factor_lon=div_factor_lon/1.001
		    j=[0,0]


#		print loop_count
		if (n.shape(i)[0] == 1) & (n.shape(j)[0] == 1): 
		    loop=0.

	    print 'finding the closest grid point to coordinate (',user_lat,',' ,user_lon,')'
	    print i,j
	    if (len(n.shape(lat_array)) == 3):
	        print 'result of the above coordinate points is     (',lat_array[0,i[0],j[0]],',' ,long_array[0,i[0],j[0]],')'
	    elif (len(n.shape(lat_array)) == 2):
	        print 'result of the above coordinate points is     (',lat_array[i[0],j[0]],',' ,long_array[i[0],j[0]],')'
            print ''
#	    print del_lat_old
#	    print del_lon_old
	    return i[0],j[0]
	    #return point_i,point_j


	def interp_to_height(self,height_levels, height_on_model_levels,input_variable):
	    """interpolate your data from model levels to height above ground
	    Inputs:
	    -------
	    height_levels (in meters, 1-D vector)
	    height_on_model_levels  (in meters, basically geopotental height / 9.81) i.e. get_var('Z')
	    input_variable  (which ever variable you wish to interpolate)

	    Returns:
	    -------
	    Your variable interpolated to constant height above ground
	
	
	    Notes:
	    ------
	    Assuming the following dimensions  
	    4 dimensions = (time, bottom_top,sound_north, west_east)
	
	    """


	    if (len(n.shape(input_variable)) == 4):
		ti_dim=n.shape(input_variable)[0]
		hi_dim=n.shape(input_variable)[1]
		sn_dim=n.shape(input_variable)[2]
		we_dim=n.shape(input_variable)[3]
		output_variable=n.zeros((ti_dim,len(height_levels),sn_dim,we_dim),dtype=n.float32)

		for the_time in range(ti_dim):
		    for i in range(sn_dim):
			for j in range(we_dim):
			    output_variable[the_time,:,i,j] = s.interp(height_levels[:],height_on_model_levels[the_time,:,i,j],input_variable[the_time,:,i,j])


	    if (len(n.shape(input_variable)) == 1):
		hi_dim=n.shape(input_variable)[0]
		output_variable=n.zeros((len(height_levels)),dtype=n.float32)

		output_variable[:] = s.interp(height_levels[:],height_on_model_levels[:],input_variable[:])

	    return output_variable


	def interp_to_pressure(self,pres_levels, pres_on_model_levels,input_variable):
	    """interpolate your data from model levels to a given set of pressure levels
	    Note that this is just a test and might be buggy, check the code and use at own risk
	    simon made this when tired.
	    Inputs:
	    -------
	    pres_levels (in meters, 1-D vector)
	    prest_on_model_levels  
	    input_variable  (which ever variable you wish to interpolate)

	    Returns:
	    -------
	    Your variable interpolated to constant given pressure
	
	
	    Notes:
	    ------
	    Assuming the following dimensions  
	    4 dimensions = (time, bottom_top,sound_north, west_east)
	
	    """


	    if (len(n.shape(input_variable)) == 4):
		ti_dim=n.shape(input_variable)[0]
		hi_dim=n.shape(input_variable)[1]
		sn_dim=n.shape(input_variable)[2]
		we_dim=n.shape(input_variable)[3]
		output_variable=n.zeros((ti_dim,len(pres_levels),sn_dim,we_dim),dtype=n.float32)
		input_variable_rev=n.zeros((ti_dim,hi_dim,sn_dim,we_dim),dtype=n.float32)
		pres_on_model_levels_rev=n.zeros((ti_dim,hi_dim,sn_dim,we_dim),dtype=n.float32)


		for ti in range(ti_dim):
		    for k in range(hi_dim):
			pres_on_model_levels_rev[ti,hi_dim -1 - k,:,:] = pres_on_model_levels[ti,k,:,:]
			input_variable_rev[ti,hi_dim -1 - k,:,:] = input_variable[ti,k,:,:]
	    
		for the_time in range(ti_dim):
		    for i in range(sn_dim):
			for j in range(we_dim):
			    output_variable[the_time,:,i,j] = s.interp(pres_levels[:],pres_on_model_levels_rev[the_time,:,i,j],input_variable_rev[the_time,:,i,j])


		


	    if (len(n.shape(input_variable)) == 1):
		hi_dim=n.shape(input_variable)[0]
		output_variable=n.zeros((len(pres_levels)),dtype=n.float32)
		output_variable[:] = s.interp(pres_levels[:],pres_on_model_levels[:],input_variable[:])

	    return output_variable
	

	def calculate_dbz_lin(self):
	    """Calculate simulated radar reflectivity based on Lin et al microphysics
	    Inputs:
	    -------
		None: will get the variables it needs itself
	    Returns:
	    -------
		reflectivity fields
	
	    Notes:
	    ------
		The actual convesion routine is based on a fortran routine stolen from grads
		If this does not work you may need to recomiple the libary object for your 
		system. Use the program f2py
		example: f2py -c -m dbzcalc_lin_py dbzcalc_lin_py.f

		Further note: you need to check in the dbzcalc_lin_py.f that the intercept variables
		are correct for your version of the lin code. I.e. Check the actual microphysics
		routine in wrf.
	
	    """

	    try:
	        import dbzcalc_lin_py
	    except:
		print 'cannot import dbzcalc_lin_py, please run f2py on dbzcalc_lin_py.f'
		print 'reflectivity functions for Lin scheme will not be available untill you do this'


	    try:
		qvp=self.variable_dict['QVAPOR']
	    except:
		qvp=self.get_var('QVAPOR')
	    try:
		qra=self.variable_dict['QRAIN']
	    except:
		qra= self.get_var('QRAIN')
	    try:
		qsn=self.variable_dict['QSNOW']
	    except:
		qsn=self.get_var('QSNOW')
	    try:
		qgr=self.variable_dict['QRGAUP']
	    except:
		qgr=self.get_var('QGRAUP')
	    try:
		tk=self.variable_dict['TEMP']
	    except:
		tk=self.get_var('TEMP')
	    try:
		prs=self.variable_dict['PRES']
	    except:
		prs=self.get_var('PRES')



	    in0r=0
	    in0s=0
	    in0g=0
	    iliqskin=0
	    if (len(n.shape(qvp)) == 4):
		ntimes=n.shape(qvp)[0]
		mkzh=n.shape(qvp)[1]
		mjx=n.shape(qvp)[2]  #NORTH-SOUTH
		miy=n.shape(qvp)[3]  #WEST-EAST



	    #need to reorient z x y to x y z
	    
	    qvp_r=n.zeros((miy,mjx,mkzh),dtype=n.float32)
	    qra_r=n.zeros((miy,mjx,mkzh),dtype=n.float32)
	    qsn_r=n.zeros((miy,mjx,mkzh),dtype=n.float32)
	    qgr_r=n.zeros((miy,mjx,mkzh),dtype=n.float32)
	    tk_r =n.zeros((miy,mjx,mkzh),dtype=n.float32)
	    prs_r=n.zeros((miy,mjx,mkzh),dtype=n.float32)
	    qvp_r=n.zeros((miy,mjx,mkzh),dtype=n.float32)
	    dbz=n.zeros((ntimes,miy,mjx,mkzh),dtype=n.float32)

	    
	

	    
	    for ti in range(ntimes):
		for k in range(mkzh):
		    qvp_r[:,:,k]=n.transpose(qvp[ti,k,:,:])
		    qra_r[:,:,k]=n.transpose(qra[ti,k,:,:])
		    qsn_r[:,:,k]=n.transpose(qsn[ti,k,:,:])
		    qgr_r[:,:,k]=n.transpose(qgr[ti,k,:,:])
		    tk_r[:,:,k]= n.transpose(tk[ti,k,:,:])
		    prs_r[:,:,k]=n.transpose(prs[ti,k,:,:])
		    #qvp_r[:,:,k]=n.transpose(qvp[ti,k,:,:])


		dbz[ti,:,:,:]=dbzcalc_lin_py.dbzcalc(qvp_r[:,:,:],qra_r[:,:,:],qsn_r[:,:,:],qgr_r[:,:,:],tk_r[:,:,:],prs_r[:,:,:],dbz[ti,:,:,:],in0r,in0s,in0g,iliqskin,miy=miy,mjx=mjx,mkzh=mkzh)


	    #transpose back


	    dbz_lin=n.zeros((ntimes,mkzh,mjx,miy),dtype=n.float32)

	    for ti in range(ntimes):
		for k in range(mkzh):
		    dbz_lin[ti,k,:,:]=n.transpose(dbz[ti,:,:,k])

	    self.variable_dict.update({'dbz_lin':dbz_lin})

	    return dbz_lin

	def calculate_dbz_thompson(self):
	    """Calculate simulated radar reflectivity based on Thompson et al microphysics
	    Inputs:
	    -------
		None: will get the variables it needs itself
	    Returns:
	    -------
		reflectivity fields
	
	    Notes:
	    ------
	        simon needs to write something here
	
	    """

	    try:
	        #import dbzcalc_thompson_py_ext
	        import dbzcalc_thompson_py
	    except:
		print 'cannot import dbzcalc_thompson_py_ext please run f2py on fortran file'
		print 'reflectivity functions for Thompson scheme will not be available untill you do this'


	    try:
		qvp=self.variable_dict['QVAPOR']
	    except:
		qvp=self.get_var('QVAPOR')
	    try:
		qra=self.variable_dict['QRAIN']
	    except:
		qra= self.get_var('QRAIN')
	    try:
		qsn=self.variable_dict['QSNOW']
	    except:
		qsn=self.get_var('QSNOW')
	    try:
		qgr=self.variable_dict['QRGAUP']
	    except:
		qgr=self.get_var('QGRAUP')
	    try:
		tk=self.variable_dict['TEMP']
	    except:
		tk=self.get_var('TEMP')
	    try:
		prs=self.variable_dict['PRES']
	    except:
		prs=self.get_var('PRES')

	    def compute_rho_dry(rho,tmk, pressure):
	        rgas=287.04
	        rho = pressure*100.0/(rgas*tmk)

	        return rho



	    def compute_rho(rho,tmk, pressure, qvp,nx,ny,nz):
	        rgas=287.04
		rho = pressure*100.0/(rgas*(tmk*(0.622+qvp)/(0.622*(1.+qvp))))

		return rho





	    in0r=0
	    in0s=0
	    in0g=0
	    iliqskin=0
	    if (len(n.shape(qvp)) == 4):
		ntimes=n.shape(qvp)[0]
		mkzh=n.shape(qvp)[1]
		mjx=n.shape(qvp)[2]  #NORTH-SOUTH
		miy=n.shape(qvp)[3]  #WEST-EAST



	    #need to reorient z x y to x y z
	    
	    qvp_r=n.zeros((miy,mjx,mkzh),dtype=n.float32)
	    qnr_r=n.zeros((miy,mjx,mkzh),dtype=n.float32)
	    qra_r=n.zeros((miy,mjx,mkzh),dtype=n.float32)
	    qsn_r=n.zeros((miy,mjx,mkzh),dtype=n.float32)
	    qgr_r=n.zeros((miy,mjx,mkzh),dtype=n.float32)
	    tk_r =n.zeros((miy,mjx,mkzh),dtype=n.float32)
	    prs_r=n.zeros((miy,mjx,mkzh),dtype=n.float32)
	    rho_r=n.zeros((miy,mjx,mkzh),dtype=n.float32)
	    dbz=n.zeros((ntimes,miy,mjx,mkzh),dtype=n.float32)

	    rho=n.zeros((ntimes,miy,mjx,mkzh),dtype=n.float32)

	    rho=compute_rho(rho,tk,prs,qvp,miy,mjx,mkzh)
	    #rho=compute_rho_dry(rho,tk,prs)
	    
	

	    
	    for ti in range(ntimes):
		for k in range(mkzh):
		    rho_r[:,:,k]=n.transpose(rho[ti,k,:,:])
		    qra_r[:,:,k]=n.transpose(qra[ti,k,:,:])
		    qsn_r[:,:,k]=n.transpose(qsn[ti,k,:,:])
		    qgr_r[:,:,k]=n.transpose(qgr[ti,k,:,:])
		    tk_r[:,:,k]= n.transpose(tk[ti,k,:,:])
		    prs_r[:,:,k]=n.transpose(prs[ti,k,:,:])
		    #qvp_r[:,:,k]=n.transpose(qvp[ti,k,:,:])




		print miy,mjx,mkzh
	        print n.shape(qra_r)
		print n.shape(qnr_r)
	        print n.shape(qsn_r)
	        print n.shape(qgr_r)
		print n.shape(tk_r)
	        print n.shape(rho_r)
		print n.shape(dbz[ti,:,:,:])

#		dbz[ti,:,:,:]=dbzcalc_lin_py      .dbzcalc(qvp_r[:,:,:],qra_r[:,:,:],qsn_r[:,:,:],qgr_r[:,:,:],tk_r[:,:,:],prs_r[:,:,:],dbz[ti,:,:,:],in0r,in0s,in0g,iliqskin,miy=miy,mjx=mjx,mkzh=mkzh)
		#dbz[ti,:,:,:]=dbzcalc_thompson_py_ext.dbzcalc2(qra_r[:,:,:],qnr_r[:,:,:],qsn_r[:,:,:],qgr_r[:,:,:],tk_r[:,:,:],rho_r[:,:,:],dbz[ti,:,:,:],miy=miy,mjx=mjx,mkzh=mkzh)   #dbz will be top-down
		#poo=dbzcalc_thompson_py_ext.dbzcalc2(qra_r[:,:,:],qnr_r[:,:,:],qsn_r[:,:,:],qgr_r[:,:,:],tk_r[:,:,:],rho_r[:,:,:],dbz[ti,:,:,:],miy=miy,mjx=mjx,mkzh=mkzh)   #dbz will be top-down
		poo=dbzcalc_thompson_py.dbzcalc2(qra_r[:,:,:],qnr_r[:,:,:],qsn_r[:,:,:],qgr_r[:,:,:],tk_r[:,:,:],rho_r[:,:,:],dbz[ti,:,:,:],miy=miy,mjx=mjx,mkzh=mkzh)   #dbz will be top-down
		print n.shape(poo)
	    #transpose back
	    
#	    dbz_test=poo[6]
	    dbz_test=poo#[6]
	    print n.shape(dbz_test)

	    dbz_thomp=n.zeros((ntimes,mkzh,mjx,miy),dtype=n.float32)

	    for ti in range(ntimes):
		for k in range(mkzh):
		    dbz_thomp[ti,k,:,:]=n.transpose(dbz_test[:,:,k])

	    self.variable_dict.update({'dbz_thomp':dbz_thomp})

	    return dbz_thomp




	def calculate_dbz_ferrier(self):
	    """Calculate simulated radar reflectivity based on Ferrier microphysics
	    Inputs:
	    -------
		None: will get the variables it needs itself
	    Returns:
	    -------
		reflectivity fields
	
	    Notes:
	    ------
		The actual convesion routine is based on a fortran routine stolen from the CALMICT.f90 routine (Brad Ferrier)
		If this does not work you may need to recomiple the libary object for your system. 
		Use the program f2py,
		example: f2py -c -m dbzcalc_ferrier_py dbzcalc_ferrier_py.f90

	    """
	    
	    try:
		import dbzcalc_ferrier_py
	    except:
		print 'cannot import dbzcalc_ferrier_py, please run f2py on dbzcalc_ferrier_py.f'
		print 'reflectivity functions for ferrier scheme will not be available untill you do this'



	    try:
		P1D=self.variable_dict['PRES']*100 # Pressure (PA)
	    except:
		P1D=self.get_var('PRES')*100
	    try:
		T1D=self.variable_dict['TEMP'] # Temperature (K)
	    except:
		T1D=self.get_var('TEMP')
	    try:
		Q1D=self.variable_dict['SPH'] # Specific Humidity (kg/kg)
	    except:
		Q1D=self.get_var('SPH')
	    try:
		C1D=self.variable_dict['CWM'] # Total Condensate (CWM, kg/kg)
	    except:
		C1D=self.get_var('CWM')
	    try:
		QW1=self.variable['QCLOUD'] # QCLOUD
	    except:
		QW1=self.get_var('QCLOUD')
	    try:
		QR1=self.variable['QRAIN'] # QRAIN
	    except:
		QR1=self.get_var('QRAIN')
	    try:
		QI1=self.variable['QSNOW'] # QSNOW
	    except:
		QI1=self.get_var('QSNOW')

	    if (self.wrf_core == 'ARW'):
		try:
		    FI1D=self.variable_dict['F_ICE_PHY'] # F_ice (fraction of condensate in form of ice)
		except:
		    FI1D=self.get_var('F_ICE_PHY')
		try:
		    FR1D=self.variable_dict['F_RAIN_PHY'] # F_rain (fraction of liquid water in form of rain
		except:
		    FR1D=self.get_var('F_RAIN_PHY')
		try:
		    FS1D=self.variable_dict['F_RIMEF_PHY'] # F_RimeF ("Rime Factor", ratio of total ice growth to deposition growth)
		except:
		    FS1D=self.get_var('F_RIMEF_PHY')

	    elif (self.wrf_core == 'NMM'):
		try:
		    FI1D=self.variable_dict['F_ICE'] # F_ice (fraction of condensate in form of ice)
		except:
		    FI1D=self.get_var('F_ICE')
		try:
		    FR1D=self.variable_dict['F_RAIN'] # F_rain (fraction of liquid water in form of rain
		except:
		    FR1D=self.get_var('F_RAIN')
		try:
		    FS1D=self.variable_dict['F_RIMEF'] # F_RimeF ("Rime Factor", ratio of total ice growth to deposition growth)
		except:
		    FS1D=self.get_var('F_RIMEF')
	    else:
		print 'microhpysics variables not found, check for F_ICE, FRAIN,F_RIMEF'
		return None



	    # Read in MASSR and MASSI arrays from text files
	    # The text files are from ETAMPNEW_DATA 
	    # These are used in the MICROINIT.f subroutine of CALMICT.f
	    fid_massr=open(PYWRFPATH+'/library/eta_tables_massr.txt','r')
	    fid_massi=open(PYWRFPATH+'/library/eta_tables_massi.txt','r')

	    MASSR=n.array(fid_massr.readline().split(','),dtype=n.float32)
	    MASSI=n.array(fid_massi.readline().split(','),dtype=n.float32)


	    # Set up grid (vertical levels = mkzh, NSdim = y = miy, WEdim = x = mjx, ) 
	    if (len(n.shape(FS1D)) == 4):
		ntimes=n.shape(FS1D)[0]
		mkzh=n.shape(FS1D)[1]
		mjx=n.shape(FS1D)[2]   # NS-DIMENSION
		miy=n.shape(FS1D)[3]   # WE-DIMENSION

		
	    QS1		=n.zeros((mjx,miy),dtype=n.float32) 	# "Snow" (precipitation ice) mixing ratio (kg/kg)
	    Tdbz	=n.zeros((ntimes,mkzh,mjx,miy),dtype=n.float32) 	# output array (currently DBZ1)	
	    Rdbz	=n.zeros((ntimes,mkzh,mjx,miy),dtype=n.float32) 	# output array (currently DBZR)
	    Idbz	=n.zeros((ntimes,mkzh,mjx,miy),dtype=n.float32) 	# output array (currently DBZI)
	    DBZ1	=n.zeros((mjx,miy),dtype=n.float32) 	# Equivalent radar reflectivity factor in dBZ; ie., 10*Log10(Z)
	    DBZR1	=n.zeros((mjx,miy),dtype=n.float32) 	# Equivalent radar reflectivity factor from rain in dBZ
	    DBZI1	=n.zeros((mjx,miy),dtype=n.float32) 	# Equivalent radar reflectivity factor from ice (all forms) in dBZ
	    DBZC1	=n.zeros((mjx,miy),dtype=n.float32) 	# Equivalent radar reflectivity factor from parameterized convection in dBZ
	    NLICE1	=n.zeros((mjx,miy),dtype=n.float32)	# Time-averaged number concentration of large ice
	    CUREFL	=n.zeros((mjx,miy),dtype=n.float32)	# Radar reflectivity contribution from convection (mm**6/m**3)

	    im = mjx
	    jm = miy

	    for ti in range(ntimes):
		for k in range(mkzh):
		     #Tdbz[ti,k,:,:], Rdbz[ti,k,:,:], Idbz[ti,k,:,:] = dbzcalc_ferrier_py.calmict(P1D[ti,k,:,:],T1D[ti,k,:,:],Q1D[ti,k,:,:],C1D[ti,k,:,:],FI1D[ti,k,:,:],FR1D[ti,k,:,:],FS1D[ti,k,:,:],CUREFL,QW1[ti,k,:,:],QI1[ti,k,:,:],QR1[ti,k,:,:],QS1,DBZ1,DBZR1,DBZI1,DBZC1,NLICE1,miy,mjx,MASSR,MASSI)

		     Tdbz[ti,k,:,:] = dbzcalc_ferrier_py.calmict(P1D[ti,k,:,:],T1D[ti,k,:,:],Q1D[ti,k,:,:],C1D[ti,k,:,:],FI1D[ti,k,:,:],FR1D[ti,k,:,:],FS1D[ti,k,:,:],CUREFL,QW1[ti,k,:,:],QI1[ti,k,:,:],QR1[ti,k,:,:],QS1,DBZ1,DBZR1,DBZI1,DBZC1,NLICE1,miy,mjx,MASSR,MASSI,im,jm)

	    self.variable_dict.update({'dbz_ferrier':Tdbz})
	    
	    return Tdbz #, Rdbz, Idbz


	def compute_height(self):
	    """Calculate geopotential height (3d) in meters from surace geopotient 
	    this is based on the hypsometer equation which is
	     Z_2 = Z_1 + (R_T_v/g)*[log(p_1) - log(p_2)] 
	     and only technicaly valid for hydrostatic balance
	     use this only for the NMM core case they do not give you geopotential height.
	     

	    Inputs:
	    -------
		None: will get the variables it needs itself
	    Returns:
	    -------
		geopotential height in meters (I think)
	
	    Notes:
	    ------
		Only technicaly valid for hydrostatic balance, I may add further.
		Simon change this when you add non-hydro support.
		(even if non-hydro just use it for now cause it will be close enough)
	    """
	    
	    try:
		#print 'getting the unstaggered version of PINT'
		prs=self.niofile.variables['PINT'].get_value()
	    except:
		print 'cannt fint your variable PINT'
	    try:
		ter=self.variable['FIS']/9.81 # Terrain height = surface geopotential height/9.81
	    except:
		ter=self.get_var('FIS')/9.81

	    try:
		vtmk=self.compute_vtmk()
	    except:
		print 'could not calculate virtual temperature'
		return  None


	    #using the hypsometer equation which is
	    # Z_2 - Z_1 = (R*T_v/g)*log(p_1/p_2)
	    # Z_2 = Z_1 + (R_T_v/g)*[log(p_1) - log(p_2)] 
	    
	    #fix up geopotential height at broken points
	    ter=n.where(ter < 0,0,ter)
	    


	    # Set up grid (vertical levels = mkzh, NSdim = y = miy, WEdim = x = mjx, ) 
	    if (len(n.shape(prs)) == 4):
		ntimes=n.shape(prs)[0]
		mkzh=n.shape(prs)[1]-1
		mjx=n.shape(prs)[2]   # NS-DIMENSION
		miy=n.shape(prs)[3]   # WE-DIMENSION


	    log_p=n.zeros((ntimes,mkzh+1,mjx,miy), dtype="float")
	    base_level=n.zeros((ntimes,mjx,miy), dtype="float")
	    ght=n.zeros((ntimes,mkzh+1,mjx,miy), dtype="float")
	

	    for ti in range(ntimes):
	        for k in range(mkzh+1):
		   log_p[ti,k,:,:]=n.log(prs[ti,k,:,:])


	    rgas=287.04
	    grav=9.81


	    for ti in range(ntimes):
		ght[ti,0,:,:] = ter[ti,:,:]
	        for k in range(1, mkzh+1):
		    tv_average=(vtmk[ti,k-1,:,:]) # not doing average for some reason
		    ght[ti,k,:,:]  = ght[ti,k-1,:,:] + ((rgas*tv_average)/grav)*(log_p[ti,k-1,:,:] - log_p[ti,k,:,:])


	    #now we have ght on staggared levels, we must destagger

	    ght_destag=n.zeros((ntimes,mkzh,mjx,miy), dtype="float")

	    for ti in range(ntimes):
		for k in range(mkzh):
		    ght_destag[ti,k,:,:] = 0.5*( ght[ti,k,:,:] + ght[ti,k+1,:,:] )

		    
	    self.variable_dict.update({'GHT':ght_destag})
	    
	    return ght_destag 


	def compute_extrema(self,mat,mode='wrap',window=10):
	    from scipy.ndimage.filters import maximum_filter, minimum_filter
	    """find the indicies of local extrema (min and max
	    in the local array. """
	    mn=minimum_filter(mat,size=window,mode=mode)
	    mx=maximum_filter(mat,size=window,mode=mode)
	    return n.nonzero(mat == mn), n.nonzero(mat==mx)



	def write_netcdf_file(self,input_variable,var_name,directory='.',filename='new_output', var_dim=('T','BT','SN','WE'), var_type='f'):
	    """This is a function to write out a selected field to a netcdf file
	    Inputs:
	    -------
	    input_data	(this is an N-dimensional array you want to write)
	    var_name	(the name you want to call your variable)

	    Optional Inputs:
	    ----------------
	    directory (the directory where your file will be produced, default is '.'
	    filename  (the name of the netcdf file you want to write, defile is new_output.nc')
	    var_type = 'type' (default = float ='f'), for string use 's'
	    var_dim  = a type of any lenght, but use the index 'T' for Time, 'BT' for height or bottom_top
			'SN' for south-north or y-dimension and 'WE' for west-east or x-dimenion'

		      for example var_dim = ('T','BT','SN','WE') for 4-d array (default)
					    ('T','SN','WE') for 2-d array at multiple times
					    ('BT','SN','WE') for 3-d array, single time

	    Returns:
	    -------
	
	
	    Notes:
	    ------
	    Assuming the following dimensions as default when no var_dim is called 
	    4 dimensions = (time, bottom_top,south_north, west_east)

	    will add stuff for different dimensions as required

	    This is simple and does not produce as much meta-data as it should, will add later
	
	    """

	   

	    dir_exist=os.path.isdir(directory)
	    if (dir_exist == False):
		print 'chosen directory does not exist, making it now'
		os.system('mkdir -p '+directory)

	    list_files = dircache.listdir(directory)

	    create_option='c'

	    for file in list_files:
		if file == filename+'.nc':
		    print 'your file already exists'
		    for try_count in range(20):
			create_option = raw_input('add variable '+var_name +'  to file, will overwrite variable if it exists?\n')
			if (create_option == 'y') | (create_option =='n'):
			    break
			else:
			    print 'type y or n'

	    if (create_option == 'n'):
		print 'quitting, doing nothing then'
		return
		    
	    
	    master_var_dim=('T','BT','SN','WE')
	    master_dim_name=('Time', 'bottom_top','south_north','west_east')

	    dimension_dict={}
	    dimension_names=[]
	    dimension_numbers=[]
	    for dim_n in range(len(var_dim)):
		
		if( len(var_dim) != len(n.shape(input_variable) )):
		    print 'Error, your array and your dimension shape differ'
		    return
		
		dimension_dict.update({var_dim[dim_n]:n.shape(input_variable)[dim_n]} )
		dimension_names.append(master_dim_name[master_var_dim.index(var_dim[dim_n])])
		dimension_numbers.append(dimension_dict[var_dim[dim_n]])
#	    return tuple(dimension_names),tuple(dimension_numbers)

		
	    var_names=tuple(dimension_names)  #var_name is the wrf way of describing dimension names i.e. ('Time','bottom_top','south_north','west_east')
	    var_nums=tuple(dimension_numbers)



        

	    if (create_option=='c'):
		new_file=Nio.open_file(directory+'/'+filename,mode='w',format='nc')

    
		try:
		    current_dims=new_file.dimensions
		except:
		    current_dims={}


		if (len(var_names) == 1) & (var_names[0] == 'Time'):
		    if (current_dims.has_key('DateStrLen') == False):
			 # get DateStrLen for Times
			self.get_dims()
			DateStrLen=self.dims['DateStrLen']
			new_file.create_dimension('DateStrLen',DateStrLen)
			var_names=('Time','DateStrLen')
			var_type='S1'

		for dim_n in range(len(var_names)):
		    if ( current_dims.has_key(var_names[dim_n]) == False):
			new_file.create_dimension(var_names[dim_n],var_nums[dim_n])

#
	        #the new variable will be called "var_name" in the netcdf file and has dimensions var_names
		new_var=new_file.create_variable(var_name,var_type,var_names)
		new_var.assign_value(input_variable)
		new_file.close()

	    if (create_option=='y'):
		print directory
		if filename[-3:] != '.nc':
		    filename=filename+'.nc'
		print filename
		new_file=Nio.open_file(directory+'/'+filename,mode='rw')
		
		try:
		    old_var=new_file.variables[var_name]	# look for the variable in the existing file
		except:

		    if (len(var_names) == 1) & (var_names[0] == 'Time'):
			try:
			# get DateStrLen for Times
			    self.get_dims()
			    DateStrLen=self.dims['DateStrLen']
			    new_file.create_dimension('DateStrLen',DateStrLen)
			except:
			    pass
			var_names=('Time','DateStrLen')
			var_type='S1'
			
	
		    
		    old_var=new_file.create_variable(var_name,var_type,var_names)	# if the variable does not exist, create it
	
		old_var.assign_value(input_variable)
		new_file.close()



class wrf_plots():
    """ This is where generic plotting functions should go"
    """
    def __init__(self):
	self.bmap()
	pass

    def bmap(self,lat=None,lon=None,proj='lcc'):
	"""create map projection and transform lat and long values to x and y
	Inputs:
	 -------
	     None. Function will get the variable it needs from your variable dict or file'
	Optional Inputs:
	     lat : This is useful if you want specify a smaller domain (i.e. not plot the whole wrf domain)
	     lon : This is useful if you want specify a smaller domain (i.e. not plot the whole wrf domain)



	 Returns:
	  -------
	     This function also returns nothing. However the important variables become associated
	     with your file object.
	     I.e. map projection becomes        wrf_file.map_proj
	          transformed lat array becomes wrf_file.map_x
	          transformed lan array becomes wrf_file.map_y


	     m or wrf_file.map_proj is the basemap projection itself, once createed you can use this to plot your
	     data nicely.

	Example:
	  -------
	   import pyWRF
	   wrf_file=pyWRF.wrf_file('input_file/wrfout_d01_2007-12-18_00:00:00.nc')
	   wrf_file.bmap()  this creates your projection and transforms your arrays

	   Can do the following for a once off
	        wrf_file.map_proj.contour(map_x,map_y, variable)

	   I prefer the following.
	        x=wrf_file.map_x
	        y=wrf_file.map_y
	        m=wrf_file.map_proj

		m.contour(x,y,variable)
		i.e. rain = wrf_file.get_var('RAINNC')
		m.contourf(x,y,rain[0,:,:])     <--- must be 2d, first dimensions is time
		m.drawcoastlines()
		imrort pylab as pl (should already be done if you are using a script)
		pl.show()


		Can add meridians and paralllels and states and lots

	"""
	#try:
	#    lon=self.variable_dict['XLONG']
	#    lat=self.variable_dict['XLAT']
	#except KeyError:

	if (lat == None) | (lon == None):
	    if self.wrf_core == 'ARW':
		try:
		    lon=self.get_var('XLONG')
		    lat=self.get_var('XLAT')
		except ValueError:
		    print 'Cant find your lat/long data, checking if you are using met_em files'
		    lon=self.get_var('XLONG_M')
		    lat=self.get_var('XLAT_M')
		    print 'yes you are using met_em files'
	    elif self.wrf_core == 'NMM':
 		try:
		    glon=self.get_var('GLON')
		    glat=self.get_var('GLAT')
		    lon=glon*57.2957795
		    lat=glat*57.2957795
		    print "found GLON AND GLAT"
		except ValueError:
		    print "GLON AND GLAT NOT FOUND"

        cen_lat	=self.niofile.CEN_LAT[0]
	cen_lon	=self.niofile.CEN_LON[0]
	map_proj=self.niofile.MAP_PROJ[0]

	if (map_proj == 1):
	    print 'your input data has map projection of Lambert Conformal'
	    #proj = 'lcc'									#testing, default call with nothing is llc, else proj is user defined
	    print 'using plotting projection of ' +proj

	elif (map_proj == 2):
	    print 'your input data has amap projection of Polar stereographic'
	    print 'assuming north pole (not tested this projection type, fix if you use it'
	    proj = 'npstere'										#testing, default call with nothing is llc, else proj is user defined
	    #print 'using projection ' +proj

	elif (map_proj == 3):
	    print 'your input data has a map projection of Mercator'
	    #proj = 'merc'										#testing, default call with nothing is llc, else proj is user defined
	    print 'using plotting projection of ' +proj


	elif (map_proj == 203):
	    print 'your input data has a  map projection of rotated lat/long'
	    print 'using plotting  projection of '+proj
	    #proj = 'lcc'											#testing, default call with nothing is llc, else proj is user defined

	else:
	    print 'map projection not known, quitting'
	    return


	m=basemap.Basemap(projection=proj,llcrnrlon=lon[0,0,0],llcrnrlat=lat[0,0,0],urcrnrlon=lon[0,-1,-1],urcrnrlat=lat[0,-1,-1],lat_0=cen_lat,lon_0=cen_lon,resolution='h')
	x , y = m(lon[0,:,:],lat[0,:,:])
	self.map_proj=m
	self.map_x=x
	self.map_y=y
	#return m,x,y

    def map_lines(self,map_proj=None,nlines=10.,coast=True):
	"""This is intended to be a semi-smart function that will determine how
	many lat/lon lines should be draw on your map (and plots coast lines)
	Inputs:
	 -------
	     None. Function uses basemap, so will call bmap() if it is not realdy called'
	Optional inputs:
	     Modified so that you can input basemap projection m
	     i.e.e self.bmap() --> self.map_proj
	     this is included because calculating this takes a lot of work, thus if
	     you are looping through files, you want to calculate this only once at
	     the start, then input your value to the function later to avoid doubling up
	     on work

	     coast=True or False. default True,
	     determines whether or not you want the coast lines drawn

	 Returns:
	  -------
	     This function also returns nothing. However lat/lon lines and coast lines will be plotted
	     on your current figure

	Example:
	"""
	#has basemap been called before?
	if (map_proj == None):
	    try:
		m=self.map_proj
	    except:
		self.bmap()
		m=self.map_proj
	else:
	    m=map_proj




	#need to determine how what the grid spacing should be for your map
	if self.wrf_core == 'ARW':
	    lat_var='XLAT'
	    lon_var='XLONG'
	elif self.wrf_core == 'NMM':
	    lat_var='GLAT'
	    lon_var='GLON'
	    

	try:
	    lat_min=n.min(self.variable_dict[lat_var])
	except:
	    self.get_var(lat_var)
	    lat_min=n.min(self.variable_dict[lat_var])

	try:
	    lon_min=n.min(self.variable_dict[lon_var])
	except:
	    self.get_var(lon_var)
	    lon_min=n.min(self.variable_dict[lon_var])

	lat_max=n.max(self.variable_dict[lat_var])
	lon_max=n.max(self.variable_dict[lon_var])

	lat_range= lat_max-lat_min
	lon_range= lon_max-lon_min

	if (self.wrf_core == 'NMM'):
	    lat_range=lat_range*57.2957795
	    lon_range=lon_range*57.2957795


    
	#going to assume we should have no more than 10 lines in each direction
	#and that we want lines either 1 deg, 5deg, 10deg, 20deg etc
	del_lat=lat_range/nlines
	del_lon=lon_range/nlines


	spacing_array=n.array([0.5,1.0,5.0,10.0,20.0],dtype='float')
	
	lat_diff= n.sqrt((spacing_array-del_lat)**2)
	lon_diff= n.sqrt((spacing_array-del_lon)**2)

	lat_pos = n.where(lat_diff == n.min(lat_diff))
	lon_pos = n.where(lon_diff == n.min(lon_diff))

	delat=spacing_array[lat_pos]
	delon=spacing_array[lon_pos]



	if coast:
	    m.drawcoastlines()
	circles = n.arange(0.,90.+delat,delat).tolist()+\
	    n.arange(-delat,-90.-delat,-delat).tolist()
	m.drawparallels(circles,labels=[1,0,0,0],fontsize=10.)
	meridians = n.arange(10.,360.,delon)
	m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10.)


    def plot_skewt(self,i,j,imagename='skewt.png', title='',timestep=0):
	"""plot a skewt log p diagram
	Inputs:
	    ------
		i  : (grid point in the i direction) (see notes)
		j  : (grid point in the j direction) (see notes)

	Optional input:
	     --------
	     imagename = 'your file name.png', defaults to skewt.png'
	     title = 'your title' defaults to none
	     timestep = integer. If there are multiple time periods per file then
	     choose your timestep, defaults to 0'
	    Returns:
	    --------
		skewT diagram
	
	    Notes:
	    ------	
	    if you want to produce a skewt at a specific lat/long location,
	    use  get_ij_lat_long function to find the closest i ,j point to 
	    your lat/long

	"""
	fig_skewt=pl.figure(figsize=(6,8))
	import skewt_old as skewt
	try:
	    pres=self.variable_dict['PRES']
	except:
	    pres=self.get_var('PRES')
	try:
	    height=self.variable_dict['Z']
	except:
	    height=self.get_var('Z')
	try:
	    temp_c=self.variable_dict['TEMP']- 273.15
	except:
	    temp_c=self.get_var('TEMP')-273.15
	try:
	    td=self.variable_dict['TD']
	except:
	    td= self.get_var('TD')
	try:
	    u=self.variable_dict['U']
	except:
	    u=self.get_var('U')
	try:
	    v=self.variable_dict['V']
	except:
	    v=self.get_var('V')

	skewt.draw_skewt(pres[timestep,:,i,j],height[timestep,:,i,j],temp_c[timestep,:,i,j],td[timestep,:,i,j],u[timestep,:,i,j],v[timestep,:,i,j],imagename=self.plot_directory+'/'+imagename,title=title)

	fig_skewt.clf()

    def plot_vapor(self,imagename='vapor.png', title='',timestep=0):
	"""plot  qvapor
	Optional input:
	     --------
	     imagename = 'your file name.png', defaults to vapor.png'
	     title = 'your title' defaults to none
	     timestep = integer. If there are multiple time periods per file then
	     choose your timestep, defaults to 0'

	    Returns:
	    --------
		contour plot of qvapor
	
	"""

	pl.figure(2)
	pl.clf()
	try:
	    vapor_full=self.variable_dict['QVAPOR']
	except:
	    vapor_full=self.get_var('QVAPOR')

	qvapor=vapor_full[timestep,0,:,:].copy()


	custom_map=pl.get_cmap('jet',lut=22)
	try:
	    m=self.map_proj
	    x=self.map_x
	    y=self.map_y
	except:
	    self.bmap()
	    m=self.map_proj
	    x=self.map_x
	    y=self.map_y

	
	custom_map=pl.get_cmap('jet',lut=10)
	contour_levs=n.arange(0.012,0.025, 0.001)
	vapor_plot=m.contourf(x,y,qvapor[:,:],contour_levs,cmap=custom_map, title=title)#,cmap=custom_map)
	vapor_plot.set_clim((0.012,0.025))

	self.map_lines()

	pl.colorbar(orientation='horizontal')
	pl.title(title)

	pl.savefig(self.plot_directory+'/'+imagename)


    def plot_winds(self,U,V,strip=10,imagename='winds.png',clearfig='yes'):
	"""plot  2-d wind barbs 
	Inputs:
	---------
	    U (X-wind component)
	    V (Y-wind component)

	Optional input:
	     --------
	     strip = integer (this means only some of the barbs will be plotted)
	     i.e. if strip = 10 then every 10th barb will be shown (default)
	     imagename = 'your file name.png', defaults to winds.png'
	     title = 'your title' defaults to none
	     clearfig = 'yes'/'no', the default is yes

	    Returns:
	    --------
		2-d wind barb plot

	    Notes:
	    -------
	        U and V must be two dimensional, i.e. you have already chosen
		which vertical level you wish to plot. 
		This can be model, pressure, or height levels (depends on your input data).
		i.e. if you choose level 10 in  U and V from the wrf file it will be on model level 10
		if you interploate to height levels then choose it is blah blah
		example:
		wrf_file.plot_winds(U[0,10,:,:],V[0,10,:,:])"""



	sn_points=n.shape(U)[0]
	we_points=n.shape(V)[1]

	u_stripped=n.zeros((sn_points,we_points),dtype=n.float32)
	v_stripped=n.zeros((sn_points,we_points),dtype=n.float32)
	u_stripped[:,:]='nan'
	v_stripped[:,:]='nan'


	for i in range(0,sn_points,strip):
	    for j in range(0,we_points,strip):
		u_stripped[i,j]=U[i,j]
		v_stripped[i,j]=V[i,j]



	try:
	    m=self.map_proj
	except:
	    self.bmap()
	    m=self.map_proj


	try:
	    lat=self.variable_dict['XLAT']
	    lon=self.variable_dict['XLONG']
	except KeyError:
	    lat=self.get_var('XLAT')
	    lon=self.get_var('XLONG')


	urot,vrot,x,y = m.rotate_vector(u_stripped,v_stripped,lon[0,:,:],lat[0,:,:], returnxy=True)

	print x
	#pl.figure(11,figsize=(16,8))
	pl.figure(11)
	if (clearfig == 'yes'):
	   pl.clf()
	
	Q=m.quiver(x,y,urot,vrot, units='width', scale=500.0, width=0.001,headwidth=3, color='r')
	QK=pl.quiverkey(Q,0.91,1.05,5,'5 m/s', labelpos='W', fontproperties={'weight':'bold'})


	self.map_lines()

	pl.savefig(self.plot_directory+'/'+imagename)


	

    def plot_cappi(self,level,MPscheme,dbz_array=None,dx=0.,dy=0.,min_dbz=-10.,radar_range=[-150,150,-150,150],imagename=None,title=None):
	"""This function will plot a cappi at a user defined level.
	Please note that the level given will actually be an array coordinate.
	i.e. if you choose level 10 and your  ata is on pressure levels then this function  will plot data for pressure level 10.
	To do a real CAPPI (constant altitude) you should first interpolate to height (i.e. wrf_file.interp_to_height)

	This functions assumes standard (time, height, lenght, width) 

	Inputs:
	    MPscheme	- microphysics scheme (lin, ferrier)
	    Level 	- The level number you wish to plot

	Optional Inputs:
	    dx, dy  - The grid spacing of your data (will assume same as the wrf spacing unless specified otherwise'

	    radar_rage = [ xmin,xmax,ymin,ymax] in km, please use a square domain, will assume 150 if not specified


	"""

	if (dbz_array == None):

	    if (MPscheme == 'lin'):
		try:
		    dbz_data=self.variable_dict['dbz_lin']
		    print 'getting lin data from variable dict'
		except:
		    dbz_data=self.calculate_dbz_lin()
	    if (MPscheme == 'ferrier'):
		try:
		    dbz_data=self.variable_dict['dbz_ferrier']
		    print 'getting ferrier data from variable dict'
		except:
		    dbz_data=self.calculate_dbz_ferrier()

	else:
	    dbz_data=dbz_array



	if (len(n.shape(dbz_data)) == 4):
	    #assuming (time,height,..,...)
	    ntimes		=n.shape(dbz_data)[0]
	    nlevels		=n.shape(dbz_data)[1]
	    xgridpoints		=n.shape(dbz_data)[2]
	    ygridpoints		=n.shape(dbz_data)[3]



	if (dx==0.) | (dy==0.):
	    print 'you did not specify a grid spacing so assuming it is the same as your wrf data'
	    delx=self.niofile.DX[0]
	    dely=self.niofile.DY[0]


	dbz_plane=pl.zeros((xgridpoints,ygridpoints),dtype=n.float32)




	xmina= radar_range[0]
	ymina= radar_range[2]

        xar=pl.zeros(xgridpoints,dtype=n.float32)
        yar=pl.zeros(ygridpoints,dtype=n.float32)
	for ii in range(xgridpoints):
	    xar[ii]=xmina+(ii)*delx
        for jj in range(ygridpoints):
	    yar[jj]=ymina+(jj)*dely

  

        times=self.get_var('Times')

	print 'xgridpoint', xgridpoints, len(xar)
	print 'ygridpoints', ygridpoints, len(yar)


        for ti in range(ntimes):


	    dbz_plane[:,:]=dbz_data[ti,level,:,:]

	    dbz_mask=n.ma.masked_less(dbz_plane, min_dbz)

	    print 'dbz_plane'
	    print n.shape(dbz_plane)
	    print 'dbz_mask'
	    print n.shape(dbz_mask)

	    pl.figure()
	    pl.clf()


	    cma=pl.get_cmap('jet',lut=10)
	    if (dbz_mask.count() == 0):
		dbz_mask[0,0]=0.0
	        #cappi=pl.pcolor(xar-delx/2,yar-dely/2,dbz_mask[:,:],cmap=cma)
	        cappi=pl.pcolor(dbz_mask[:,:],cmap=cma)
		cappi.set_clim((10,60))
    
	    if (dbz_mask.count() >0):
		dbz_mask[0,0]=0.0
		#cappi=pl.pcolor(xar-delx/2,yar-dely/2,dbz_mask[:,:],cmap=cma)
		cappi=pl.pcolor(dbz_mask[:,:],cmap=cma)
		cappi.set_clim((10,60))
		pl.colorbar()

	    #pl.axis(radar_range,'equal')

	    if title == None:
		the_title=times[ti]
	    else:
		the_title=title

	    pl.title(the_title,fontsize=25)

	    if (imagename == None):
		outimagename=times[ti]+'_dbz.png'
	    else:
		outimagename=n.str(times[ti])+'_'+imagename+'_dbz.png'


	    pl.savefig(self.plot_directory+'/'+outimagename)
	    del(cappi)


	
    

    def plot_precip(self,imagename='precip.png', title='',timestep=0):
	"""plot  accumulated precipitation
	Optional input:
	     --------
	     imagename = 'your file name.png', defaults to vapor.png'
	     title = 'your title' defaults to none
	     timestep = integer. If there are multiple time periods per file then
	     choose your timestep, defaults to 0'

	    Returns:
	    --------
		contour plot of accumulated precipitaiton
	"""

	pl.figure(6)
	pl.clf()
	try:
	    rainnc_full=self.variable_dict['RAINNC']
	except ValueError:
	    print 'you do not have any non-convective rain data'
	    rainnc_full = None
	except KeyError:
	    rainnc_full=self.get_var('RAINNC')


	try:
	    rainc_full=self.variable_dict['RAINC']
	except KeyError:
	    rainc_full=self.get_var('RAINC')
	except ValueError:
	    print 'you do not have any convective rain data'
	    rainc_full = None


	times=self.get_var('Times')



	if (timestep == 0):
	    start_time =0
	    end_time = len(times)
	else:
	    start_time=timestep
	    end_time=timestep+1

	try:
	    m=self.map_proj
	    x=self.map_x
	    y=self.map_y
	except:
	    self.bmap()
	    m=self.map_proj
	    x=self.map_x
	    y=self.map_y

	


	for the_time in range(start_time,end_time):
	    time=times[the_time]
	    pl.clf()
	    if (rainnc_full == None) & (rainc_full == None):
		print 'you have no precipitation data at all, quitting'
		return None

	    if (rainnc_full == None) & (rainc_full != None):
		rain_total = rainc_full[the_time,:,:]

	    if (rainc_full == None) & (rainnc_full != None):
		rain_total = rainnc_full[the_time,:,:]
    
	    
	    if (rainc_full != None) & (rainnc_full != None):
		rain_total=rainnc_full[the_time,:,:] + rainc_full[the_time,:,:]  


	
	    custom_map=pl.get_cmap('jet',lut=22)

	    PCP_LEVELS = [0.01,0.03,0.05,0.10,0.15,0.20,0.25,0.30,0.40,0.50,0.60,0.70,0.80,0.90,1.00,1.25,1.50,1.75,2.00,2.50]

	    max_precip = n.max(rain_total)
	    scaling_factor = n.ceil(max_precip/PCP_LEVELS[-1])


	    if scaling_factor > 0.0:
		for cont in range(len(PCP_LEVELS)):
		    PCP_LEVELS[cont]=PCP_LEVELS[cont]*scaling_factor

	    m.contourf(x,y,rain_total[:,:],PCP_LEVELS,cmap = coltbls.precip1()) 


	    self.map_lines()

	    pl.colorbar(orientation='horizontal')
	    pl.title(title)

	    pl.savefig(self.plot_directory+'/'+imagename)


    def plot_mslp(self,imagename='mslp.png', title='',timestep=0):
	"""plot  mean sea level pressure
	Optional input:
	     --------
	     imagename = 'your file name.png', defaults to mslp.png'
	     title = 'your title' defaults to none
	     timestep = integer. If there are multiple time periods per file then
	     choose your timestep, defaults to 0'

	    Returns:
	    --------
		contour plot of mean sea level pressures with Highs and Lows
		labeled 
	"""

	pl.figure(7)
	pl.clf()

	times=self.get_var('Times')

	if (timestep == 0):
	    start_time =0
	    end_time = len(times)
	else:
	    start_time=timestep
	    end_time=timestep+1

	try:
	    m=self.map_proj
	    x=self.map_x
	    y=self.map_y
	except:
	    self.bmap()
	    m=self.map_proj
	    x=self.map_x
	    y=self.map_y

	mslp=self.compute_mslp()


	local_min,local_max=self.compute_extrema(mslp,mode='constant',window=50)

	#    out_dir='/home/simon/model_simulations/plots/plots_arw/pres/'
	#    if (os.path.isdir(out_dir) == False):
	#	os.system("mkdir -p "+out_dir)
	
	#    the_time= wrf_out[11:]
	clevs=n.arange(900,1100,5)
        cs=  m.contour(x,y,mslp[:,:],clevs,colors='k', linewidths=1.0)
	pl.clabel(cs,inline=1,fontsize=8,fmt='%4.0f')
        self.map_lines(map_proj=m,coast=False)
	m.fillcontinents(color='orange',alpha=0.3)
        m.drawcoastlines(color='orange',linewidth=0.5)

	xlows  	=x[local_min]
        xhighs 	=x[local_max]
	ylows  	=y[local_min]
        yhighs	=y[local_max]
        lowvals	=mslp[local_min]
        highvals	=mslp[local_max]
        #plot lows as blue L's with min pressure value underneith
        xyplotted=[]
        #dont plot if there is already a L or H within dmin meters
        yoffset=0.022*(m.ymax-m.ymin)
	dmin=yoffset
        xlen=(m.xmax-m.xmin)*0.022
        ylen=(m.ymax-m.ymin)*0.022
	dmin=max(xlen,ylen)

        #need to order low values and high values so that the lowest and the highest get plotted

        switch=1
	while switch == 1:
	    switch =0
	    for lval in range(len(lowvals)-1):
		if lowvals[lval] > lowvals[lval+1]:
		    new_lval=lowvals[lval+1]
		    new_xval=xlows[lval+1]
		    new_yval=ylows[lval+1]

		    lowvals[lval+1]=lowvals[lval]
		    xlows[lval+1]=xlows[lval]
		    ylows[lval+1]=ylows[lval]

		    lowvals[lval]=new_lval
		    xlows[lval]=new_xval
		    ylows[lval]=new_yval
		    switch = 1



        switch=1
	while switch == 1:
	    switch =0
	    for hval in range(len(highvals)-1):
		if highvals[hval] < highvals[hval+1]:
		    new_hval=highvals[hval+1]
		    new_xval=xhighs[hval+1]
		    new_yval=yhighs[hval+1]

		    highvals[hval+1]=highvals[hval]
		    xhighs[hval+1]=xhighs[hval]
		    yhighs[hval+1]=yhighs[hval]

		    highvals[hval]=new_hval
		    xhighs[hval]=new_xval
		    yhighs[hval]=new_yval
		    switch = 1


        for xx,yy,pp in zip(xlows,ylows,lowvals):
	    if xx < m.xmax and xx > m.xmin and yy < m.ymax and yy > m.ymin:
		dist=[n.sqrt((xx-x0)**2+(yy-y0)**2) for x0,y0, in xyplotted]
		if not dist or min(dist) > dmin:
		    pl.text(xx,yy,'L',fontsize=14,fontweight='bold',ha='center',va='center',color='b')
		    pl.text(xx,yy-yoffset,repr(int(pp)), fontsize=9,ha='center',va='top',color='b',bbox= dict(boxstyle='square',ec='None',fc=(1,1,1,0.5)))
		    xyplotted.append((xx,yy))

        #plot highs as red H's with max pressure value underneight
    #    xyplotted=[]
	for xx,yy,pp in zip(xhighs,yhighs,highvals):
	    if xx<m.xmax-dmin and xx>m.xmin+dmin and yy < m.ymax-dmin and yy >m.ymin+dmin:
		dist=[n.sqrt((xx-x0)**2+(yy-y0)**2) for x0,y0 in xyplotted]
		if not dist or min(dist) > dmin:
		    pl.text(xx,yy,'H',fontsize=14,fontweight='bold',ha='center',va='center',color='r')
		    pl.text(xx,yy-yoffset,repr(int(pp)), fontsize=9,ha='center',va='top',color='r',bbox= dict(boxstyle='square',ec='None',fc=(1,1,1,0.5)))
		    xyplotted.append((xx,yy))


	pl.title(title)

	pl.savefig(self.plot_directory+'/'+imagename)





		    
	
class wrf_file(calc_vars, wrf_plots):
    """Is used to read in wrf files and extract variables 
    Variables should automatically have perturbations and base-states combined, returning the full fields 
    Variables that require it should automatically be destaggared on to a normal grid. 
    Some variables that are not standard output from the WRF model will be calculated behind the seens.
    Check to see if WRF NMM or ARW core is being used, this will then change how some of the later modules work.
    Adding grib support. I currently assume that if you are using grb files, they have been generated with UPP
    and therefore already destagged and you have calculated variables that you want.
    """
    def __init__(self,filename=None):
	self.filename=filename
	self.variable_dict={}
	self.plot_directory='.'
	import sys as sys
	sys.path.append('wrf_plots/')

	if (filename.count('wrfout')==1) & (filename[-3:] != '.nc'):
	    self.filename=str(self.filename+'.nc')
	    self.file_type='nc'
	if (filename.count('dbz_d0')==1) & (filename[-3:] != '.nc'):
	    self.filename=str(self.filename+'.nc')
	    self.file_type='nc'
	if (filename.count('WRFPRS')==1) & (filename[-4:] != '.grb'):
	    self.filename=str(self.filename+'.grb')
	    self.file_type='grb'
	    self.wrf_core="NULL"#should be indepeneant of NMM or WRF as it has been processed
	    print 'Found UPP processed grib file, assuming all data has been processed'
	if (filename.count('wrfprs')==1) & (filename[-4:] != '.grb'):
	    self.filename=str(self.filename+'.grb')
	    self.file_type='grb'
	    self.wrf_core="NULL"  #should be indepeneant of NMM or WRF as it has been processed
	    print 'Found UPP processed grib file, assuming all data has been processed'
	if (filename.count('wrfwnd')==1) & (filename[-4:] != '.grb'):
	    self.filename=str(self.filename+'.grb')
	    self.file_type='grb'
	    self.wrf_core="NULL"  #should be indepeneant of NMM or WRF as it has been processed
	    print 'Found UPP processed grib file, assuming all data has been processed'
	if (filename.count('WRFWND')==1) & (filename[-4:] != '.grb'):
	    self.filename=str(self.filename+'.grb')
	    self.file_type='grb'
	    self.wrf_core="NULL"  #should be indepeneant of NMM or WRF as it has been processed
	    print 'Found UPP processed grib file, assuming all data has been processed'





	if (self.filename.find('wrfout') == 0) | (self.filename.find('wrfprs') == 0) | (self.filename.find('WRFPRS') == 0) | (self.filename.find('dbz_d') == 0):
	    self.wrf_directory = './'
	else:
	    if self.filename.count('wrfout') == 1:
		self.wrf_directory = str(self.filename[0:self.filename.find('wrfout')])
	    elif self.filename.count('dbz_d0') == 1:
		self.wrf_directory = str(self.filename[0:self.filename.find('dbz_d0')])
	try:
	    wrf_file=Nio.open_file(self.filename,'r')
	    self.niofile=wrf_file

	    if (self.file_type=='nc'):
		if (self.niofile.variables.keys().count('DX_NMM') >=1):
		    print 'Found variable DX_NMM, assuming you are running the NMM core'
		    self.wrf_core='NMM'
		else:
		    print 'DX_NMM not found, assuming you are running the ARW core'
		    self.wrf_core='ARW'

	except:
	    print '\n \n \n \n SOMETHING WENT WRONG.\n This is likely because your file path is incorrect \n   please try again'
	    self.niofile=None

	pass

    def close(self):
	"""
	This function is used to close any wrf file that is open
	"""
	self.niofile.close()


    def list_var(self):
	"""this should list the variables available for you to use in get_var
	Note that some of these variables will be derived"""

	if self.niofile:
	    #print self.wrf_core
	    if self.wrf_core=='ARW':
	        pert_vars=perturbation_variables.pert_variable_dict_ARW
		calc_vars=perturbation_variables.calc_variable_dict_ARW
	    elif self.wrf_core=='NMM': 
	        pert_vars=perturbation_variables.pert_variable_dict_NMM
		calc_vars=perturbation_variables.calc_variable_dict_NMM
	    elif self.wrf_core=='NULL':
	        pert_vars={}
		calc_vars={}

		pass
	    else:
		print "wrf core type not found, you should be using either ARW or NMM"
		exit()


	    print "you can extract the following variables\n"
	    var_list=self.niofile.variables.keys()#+calc_vars.keys()+pert_vars.keys()
	    var_list.sort()

    
	    var_it=iter(var_list)
	    line=''

	    max_len_varname=14

	    for variable in range(len(var_list)):
		var_done=0
		var=var_it.next()

		if (len(var_list[variable]) > max_len_varname):
		    max_len_varname = len(var_list[variable])+2


	        if (variable % 4  == 3):
		    line = '\n'
		else:
		    line = ''

	
		for vvar in range(len(pert_vars.keys())):
		    if (var in pert_vars.values()[vvar][0]) | (var in str(pert_vars.values()[vvar][1])):
		        #var_color='\033[3;9m'+var.ljust(14)+'\033[m'
		        var_half='\033[0;41m'+var+'\033[m'
		        var_color=var_half.ljust(24)
		        print var_color  , line,
		        var_done=1
		        break

		if var_done == 0:
		    var_color=var.ljust(max_len_varname)
		    print var_color  , line,


	    print '\n'
	    print 'The following variables are created from perterbation and base state variables \n'

	    count=0
	    for cvar in pert_vars:
	        if (count % 4  == 3):
		    line = '\n'
		else:
		    line = ''

		var_color=''
		var_color='\033[0;31m'+cvar.ljust(14)+'\033[m'
		count+=1
		print var_color  , line,



	    print '\n'
	    print 'The following variables will be calculated for you\n'

	    count=0
	    for cvar in calc_vars:
	        if (count % 4  == 3):
		    line = '\n'
		else:
		    line = ''

		var_color=''
		var_color='\033[0;32m'+cvar.ljust(14)+'\033[m'
		count+=1
		print var_color  , line,




#    def _full_variable(self,var_tuple,fullvar):
#	'this routine combines base and perterbation variables'
#	print fullvar, '--->',var_tuple[0],'+', var_tuple[1] ,'\n'
#
#	if (type(var_tuple[0]) == type('')):
#	    invar1=self.niofile.variables[var_tuple[0]]
#	    stag1=invar1.stagger
#	    array1=invar1.get_value()
#	elif (type(var_tuple[0]) == type(1.)):
#	    array1=var_tuple[0] #does not have to be an array, in fact probably isnt
#	    stag1=''
#	else:
#	    print 'error'
#
#
#
#	if (type(var_tuple[1]) == type('')):
#	    invar2=self.niofile.variables[var_tuple[1]]
#	    stag2=invar2.stagger	
#	    array2=invar2.get_value()
#	elif (type(var_tuple[1]) == type(1.)):
#	    array2=var_tuple[1] #does not have to be an array, in fact probably isnt
#	    stag2=''
#	else:
#	    print 'error'
#
#	if (len(stag1) != len(stag2)):
#	    print 'warning, both variables are staggared the same'
#	else:
#	    stag=stag1
#
#	full_array=array1+array2
#
#	return full_array, stag
	


    def _full_variable(self,var_tuple,fullvar):
	'this routine combines base and perterbation variables'
	print fullvar, '--->',var_tuple[0],'+', var_tuple[1] ,'\n'

	if (type(var_tuple[0]) == type('')):
	    invar1=self.get_var(var_tuple[0])
	    stag1='-'
#	    array1=invar1.get_value()
	elif (type(var_tuple[0]) == type(1.)):
	    invar1=var_tuple[0] #does not have to be an array, in fact probably isnt
	    stag1='-'
	else:
	    print 'error'



	if (type(var_tuple[1]) == type('')):
	    invar2=self.get_var(var_tuple[1])
	    stag2='-'
	elif (type(var_tuple[1]) == type(1.)):
	    invar2=var_tuple[1] #does not have to be an array, in fact probably isnt
	    stag2='-'
	else:
	    print 'error'

	if (len(stag1) != len(stag2)):
	    print 'warning, both variables are staggared the same'
	else:
	    stag=stag1

	full_array=invar1+invar2

	return full_array, stag
	

    def get_var(self,var,multi=False):


	"""get WRF variables from the file, destagger if necessary
	can use list_var to see which variables your file contains.
	
	Arguements: 
	-----------
	Var:  Variable Name
	-------
	Returns:
	numpy.ndarray

	See Also
	--------
	list_var()
	will show you which variables can be obtained

	Example:
	--------
	get_var('U') 
	returns a numpy array containing the U wind data

	"""
	#import calculate_variables
	#import perturbation_variables


        if self.wrf_core=='ARW':
	    pert_vars=perturbation_variables.pert_variable_dict_ARW
	    calc_vars=perturbation_variables.calc_variable_dict_ARW
	elif self.wrf_core=='NMM': 
            pert_vars=perturbation_variables.pert_variable_dict_NMM
	    calc_vars=perturbation_variables.calc_variable_dict_NMM
	elif self.wrf_core == 'NULL':
            pert_vars={}
	    calc_vars={}

	else:
	    print "wrf core type not found, you should be using either ARW or NMM"
	    exit()


	if (var in pert_vars.keys()):
	    data_array, stag = self._full_variable(pert_vars[var],var)
	elif (var in calc_vars.keys()):
	    stag='-'
	else:
	    
	    try:
		input_var  =  self.niofile.variables[var]

	    except KeyError:
		print 'your variable ' + var+ ' does not exist'
		#raise ValueError('your WRF variable does not exist')
		input_var=None
		return input_var
 	    except ValueError:
		print 'Error ValueError'

	    vartype    =  input_var.typecode()
	    numDims    =  input_var.rank
	    dimSizes   =  input_var.shape
	    dimNames   =  input_var.dimensions
	    try:
		varAttsandVals=input_var.attributes
	    except:
		pass

	    try:
		stag	   = input_var.stagger
	    except:
		stag = '-'

	    if stag == "":
		stag = '-'
    
	    if self.wrf_core == 'NMM':
	#	print 'ASSUMING ALL NMM HORIZONTAL VARIABLES ARE NOTSTAGGERED, NOT SURE IF THIS IS CORRECT'
		if stag !='Z':
		    stag="-"
		


	    #try:
	    if (var != 'Times'):
		    #cpu_check=int(self.filename[-7:-3])
		    #print 'cpu_check', cpu_check
		    if (self.filename[-8] == '_') & (multi == False):
			print 'I think you are using multiple files, i.e. 1 per cpu, going to automatically combine for you'
			data_array , stag= self.get_var_multi_cpu(self.wrf_directory,self.filename,var,stag=stag)
			if (stag != '-'):
			    print 'Your variable', var, ' is now being de-staggered\n'
			    if self.wrf_core == 'ARW':
				data_array=    wrf_user_unstagger.wrf_user_unstagger_ARW(data_array,stag)
			    elif self.wrf_core == 'NMM':
				data_array=    wrf_user_unstagger.wrf_user_unstagger_NMM(data_array,stag)
			    else:
				print 'Error: WRF core not known, should be either ARW or NMM'

			return data_array
		    elif(self.filename[-8] == '_') & (multi == True):
			pass
		    else:
			pass
#


	    data_array=input_var.get_value()

        if (multi == True):
	    return data_array



	if (stag != '-'):
	    print 'Your variable', var, ' is now being de-staggered\n'
	    if self.wrf_core == 'ARW':
		data_array=    wrf_user_unstagger.wrf_user_unstagger_ARW(data_array,stag)
	    elif self.wrf_core == 'NMM':
		data_array=    wrf_user_unstagger.wrf_user_unstagger_NMM(data_array,stag)
	    else:
		print 'Error: WRF core not known, should be either ARW or NMM'


	
	#variables should alredy be destaggared by this point because calc works on variable already recieved.
	if (var in calc_vars.keys()):
	    for dependant_variable in calc_vars[var]:
		if (self.variable_dict.has_key(dependant_variable) == False):
		    print 'also need to extract/calculate', dependant_variable, ' if you want to calculate',var, ' doing this now\n'
		    data_array = self.get_var(dependant_variable,multi=False)
		else:
		    data_array = self.variable_dict[dependant_variable]

	    if (var == 'TEMP') & (self.wrf_core == 'ARW'):
		data_array= self.compute_tk(self.variable_dict[calc_vars[var][0]],self.variable_dict[calc_vars[var][1]])
	    elif(var == 'TEMP') & (self.wrf_core == 'NMM'):
		data_array=self.variable_dict['T']

	    if (var == 'RH'):
		data_array= self.compute_rh(self.variable_dict[calc_vars[var][0]],self.variable_dict[calc_vars[var][1]],self.variable_dict[calc_vars[var][2]])
	    if (var == 'TD'):
		data_array= self.compute_td(self.variable_dict[calc_vars[var][0]],self.variable_dict[calc_vars[var][1]])
	    if (var == 'PRES'):

		data_array= data_array/100.
	    if (var == 'Z'):
		if self.wrf_core == 'ARW':
		    print 'dividing geopotential by 9.81 to give height in m'
		    data_array= data_array/9.81
		elif self.wrf_core == 'NMM':
		    data_array= self.compute_height()
		else:
		    print 'cant find your wrf core (should be ARW or NMM)'

	    if (var == 'SPH'):
		data_array= self.compute_sph(self.variable_dict[calc_vars[var][0]])
	  
	    if (var == 'VTMK'):
		data_array= self.compute_vtmk()

	    if (self.wrf_core == 'NMM'):
		if (var == 'XLONG'):
		    print 'you asked for XLONG but are using the NMM core, I will give you GLON as XLONG so that code still works'
		    data_array=self.get_var('GLON')
		    data_array=data_array*57.2957795

		if (var == 'XLAT'):
		    print 'you asked for XLAT but are using the NMM core, I will give you GLAT as XLAT so that code still works'
		    data_array=self.get_var('GLAT')
		    data_array=data_array*57.2957795
	

	if (var == 'Times'):
	    try:
		len(time_list)
	    except:
		time_list=[]
	    for ti in range(len(data_array)):
		time_list.append(''.join(data_array[ti]))

	    del(data_array)
	    data_array=time_list


    
	if (multi == False):
	    self.variable_dict.update({var:data_array})



	return data_array




    def get_var_multi_cpu(self,multi_wrf_directory, multi_file_filename,var_desired,stag='-',all=False):
	import pyWRF
	"""create full array from multiple files. This is requed when io_form_history = 102.
	with this setting, a seperate file is created for each cpu. 
	please note that this function uses the previously defined get_var function,
	so your variables should be automatically destaggared and perterbation fields combined.
	This function is a bit lame as it requires wrf file to have already been loaded
	
	Arguements: 
	-----------
	Var:  Variable Name
	wrf_file_directory: the directory where all your wrf files are stored, assuming files
	of the form. wrfout_d0x_YYYY-MM-DD_HH:MM:SS_CPUID where CPUID = xxxx (i.e. four units long)
	-------
	Returns:
	numpy.ndarray

	See Also
	--------
	list_var()
	will show you which variables can be obtained

	Example:
	--------
	get_var('U') 
	returns a numpy array containing the U wind data

	"""

	multi_wrf_directory=str(multi_wrf_directory)
	multi_file_filename=str(multi_file_filename)


	multi_file_filename=multi_file_filename[self.filename.find('wrfout'):]
    
	list_files_temp=dircache.listdir(str(multi_wrf_directory))
	list_files=[]
	for i in range(len(list_files_temp)):
	    list_files.append(str(list_files_temp[i]))
	##############################################################
	############determine if we have one file per CPU#############
	##############################################################
	unique={}


	for file in list_files:
	    prefix=str(file[0:-5])
	    if (prefix != str(multi_file_filename[:-8])):
		continue

	    if unique.has_key(prefix) == False:
		cpu_number=file[-4:]
		unique.update({prefix:cpu_number})
	    if unique.has_key(prefix) == True:
		cpu_number=file[-4:]
		if (unique[prefix] < cpu_number):
		    unique.update({prefix:cpu_number})
	
	#print unique

	unique.keys().sort()
	for wrf_file_split in unique.keys():
	    first_time = True
	    print 'getting variable '+var_desired+' from your '+ str(int(unique[wrf_file_split])+1) +' seperate files'
	    for cpu_id in range(int(unique[wrf_file_split])+1):
		str_cpu=str(cpu_id)
		str_cpu=str_cpu.rjust(4,'0')
	

		wrf_file_m=pyWRF.wrf_file(multi_wrf_directory+wrf_file_split+'_'+str_cpu)
		wrf_file_m.get_tiles()
		if (first_time == True):
		    SN_TOTAL=wrf_file_m.tile['SN-TOTAL']-1
		    WE_TOTAL=wrf_file_m.tile['WE-TOTAL']-1
		    BT_TOTAL=wrf_file_m.tile['BT-TOTAL']-1
		    NTIMES=wrf_file_m.tile['NTIMES']

		    if (stag == 'X'):
			WE_TOTAL=WE_TOTAL+1
		    if (stag == 'Y'):
			SN_TOTAL=SN_TOTAL+1
		    if (stag == 'Z'):
			BT_TOTAL=BT_TOTAL+1
	    
		    variable_dims=int(len(n.shape(wrf_file_m.niofile.variables[var_desired].get_value())))
		    #variable_dims=len(n.shape(wrf_file_m.get_var(var_desired,multi=True)))

		    if (variable_dims == 4):
			total_array=n.zeros((NTIMES,BT_TOTAL,SN_TOTAL,WE_TOTAL),dtype=n.float32)
		    elif (variable_dims == 3):
			total_array=n.zeros((NTIMES,SN_TOTAL,WE_TOTAL),dtype=n.float32)
		    else:
			print ' need to add code to use this type of variable'
			print n.shape(variable_dims)

		    first_time = False

	        if (stag == 'X'):
		    SN_START=wrf_file_m.tile['SN-START']-1
		    SN_END=wrf_file_m.tile['SN-END']-1
		    WE_START=wrf_file_m.tile['WE-START_STAG']-1
		    WE_END=wrf_file_m.tile['WE-END_STAG']-1
		elif(stag == 'Y'):
		    SN_START=wrf_file_m.tile['SN-START_STAG']-1
		    SN_END=wrf_file_m.tile['SN-END_STAG']-1
		    WE_START=wrf_file_m.tile['WE-START']-1
		    WE_END=wrf_file_m.tile['WE-END']-1
	        else:
		    SN_START=wrf_file_m.tile['SN-START']-1
		    SN_END=wrf_file_m.tile['SN-END']-1
		    WE_START=wrf_file_m.tile['WE-START']-1
		    WE_END=wrf_file_m.tile['WE-END']-1

		if (variable_dims == 4):
		    total_array[:,:,SN_START:SN_END+1,WE_START:WE_END+1] = wrf_file_m.get_var(var_desired,multi=True)
		elif (variable_dims == 3):
		    total_array[:,SN_START:SN_END+1,WE_START:WE_END+1] = wrf_file_m.get_var(var_desired,multi=True)
		wrf_file_m.close()


	self.variable_dict.update({var_desired:total_array})


	return total_array, stag



    def get_dims(self):
	"""get WRF dimensions from the file
	
	Arguements: 
	-----------
	None:
	-------
	Returns:
	None: However, dimension dictionary now assosicated with your wrf_file object
	dictionary is named  yourfile.dims

	"""
	dims=self.niofile.dimensions
	self.dims=dims





    def get_att(self):
	"""get WRF global attributes
	
	Arguements: 
	-----------
	None:
	-------
	Returns:
	None: However, dimension dictionary now assosicated with your wrf_file object
	dictionary is named yourfile.atts

	"""
	atts=self.niofile.attributes
	self.atts=atts


    def get_tiles(self):
	"""get info needed to combine tiled data
	i.e. when WRF has written 1 output file per CPU.
	
	Arguements: 
	-----------
	None:
	-------
	Returns:
	None: However, dimension dictionary now assosicated with your wrf_file object
	dictionary is named yourfile.tile

	"""
	atts=self.niofile.attributes
	tile={}

	tile.update({'WE-START':atts['WEST-EAST_PATCH_START_UNSTAG'][0]})
	tile.update({'SN-START':atts['SOUTH-NORTH_PATCH_START_UNSTAG'][0]})
	tile.update({'WE-END':atts['WEST-EAST_PATCH_END_UNSTAG'][0]})
	tile.update({'SN-END':atts['SOUTH-NORTH_PATCH_END_UNSTAG'][0]})
	tile.update({'WE-TOTAL':atts['WEST-EAST_GRID_DIMENSION'][0]})
	tile.update({'SN-TOTAL':atts['SOUTH-NORTH_GRID_DIMENSION'][0]})
	tile.update({'BT-TOTAL':atts['BOTTOM-TOP_GRID_DIMENSION'][0]})

	tile.update({'WE-START_STAG':atts['WEST-EAST_PATCH_START_STAG'][0]})
	tile.update({'SN-START_STAG':atts['SOUTH-NORTH_PATCH_START_STAG'][0]})
	tile.update({'WE-END_STAG':atts['WEST-EAST_PATCH_END_STAG'][0]})
	tile.update({'SN-END_STAG':atts['SOUTH-NORTH_PATCH_END_STAG'][0]})


	self.get_dims()
	tile.update({'NTIMES':self.dims['Time']})

	self.tile=tile




