from astropy.io.votable import parse
from astropy.io.votable.tree import VOTableFile, Resource, Table, Field
from astroquery.gaia import Gaia
import pickle
import numpy as np

"""
Function to import data, either from a pre existing table or from a .vot file
After importing data from vot file a table file is saved for later import
This is used to save importing time, but is not very memory efficient
Table name is set for later reference 
s_ is a PATH to the file 
"""
def importData( s_ ):
	
	s = str( s_ )
	try:
		with open( s +'_table', 'rb' ) as input:
			table = pickle.load( input )
			print( "Data was loaded from an existing table." )
	except:
		votable = parse(s)
		table = votable.get_first_table()
		with open( s + '_table', 'wb' ) as output:
			pickle.dump(table,output,pickle.HIGHEST_PROTOCOL)
			
	table.name = s 
	return table
"""
Reads desired files from table object
Groups data in a tuple (dataTuple) with 5 smaller arrays: DataK, labelsK, dataErrorK, dataL , labelsL
dataK contains data about position and velocity : ra, dec, pmra, pmdec, parallax, radial_velocity(optional) which is used for clustering 
dataerrorK contains error for data in dataK
dataL contains data all other data: bp_rp, phot_g_mean_val, teff_val, lum_val, phot_bp_rp_excess_factor, a_g_val, e_bp_min_rp_val, radius_val and further desired data can be added
labels store the names of the data used in the same order they appear inside array
include_rad_vel is used if reading already processed files which may not contain radial_velocity at all
NOTE: Throughout code dataL is only referenced from front index, so extra data can be added to the end without danger of breaking old code, same is not true for dataK in which radial_velocity is assumed to be at the end
NOTE: DataTuple is used as file transfer through most of code 
NOTE: Error imported is absolute error but is later transformed to relative error, it is up to user to know which error is currently in dataTuple
"""
def getData( table , include_rad_vel = True ):
	
	# dataK	
	ra = table.array['ra'].tolist()
	dec =  table.array['dec'].tolist()
	pmra =  table.array['pmra'].tolist()
	pmdec = table.array['pmdec'].tolist()
	parallax =  table.array['parallax'].tolist()
	

	# dataL
	bp_rp = np.array( table.array['bp_rp'], dtype = 'float' )
	phot_g_mean_mag = np.array( table.array['phot_g_mean_mag'], dtype = 'float' )
	teff_val = np.array( table.array['teff_val'], dtype = 'float' )
	lum_val = np.array( table.array['lum_val'] , dtype = 'float' )
	phot_bp_rp_excess_factor = np.array ( table.array['phot_bp_rp_excess_factor'] , dtype = 'float' )
	a_g_val = np.array ( table.array['a_g_val'] , dtype = 'float' )
	e_bp_min_rp_val = np.array ( table.array['e_bp_min_rp_val'] , dtype = 'float' )
	radius_val = np.array ( table.array['radius_val'] , dtype = 'float' )
	
	# Get error for  dataK,  data_errorK 
	ra_error = np.array ( table.array['ra_error'] , dtype = 'float' )
	dec_error = np.array ( table.array['dec_error'] , dtype = 'float' )
	pmra_error = np.array ( table.array['pmra_error'] , dtype = 'float' )
	pmdec_error = np.array ( table.array['pmdec_error'] , dtype = 'float' )
	parallax_error = np.array ( table.array['parallax_error'] , dtype = 'float' )
		
	# Generate data arrays from individual lists of each parameter 
	if(include_rad_vel):
		radial_velocity = table.array['radial_velocity'].tolist()
		radial_velocity_error = table.array['radial_velocity_error'].tolist()
	
		dataK_raw = np.array( [ra, dec, pmra, pmdec, parallax, radial_velocity] ).T
		data_errorK_raw = np.array( [ra_error, dec_error, pmra_error, pmdec_error, parallax_error, radial_velocity_error] ).T
		labelsK = ["ra", "dec", "pmra", "pmdec", "parallax", "radial_velocity" ]

	else:
		dataK_raw = np.array( [ra, dec, pmra, pmdec, parallax] ).T
		labelsK = ["ra", "dec", "pmra", "pmdec", "parallax",]
		data_errorK_raw = np.array( [ra_error, dec_error, pmra_error, pmdec_error, parallax_error] ).T

	dataL_raw = np.array( [bp_rp, phot_g_mean_mag, teff_val, lum_val, phot_bp_rp_excess_factor, radius_val, a_g_val, e_bp_min_rp_val ] ).T
	labelsL = ["bp_rp", "phot_g_mean_mag", "teff_val", "lum_val", "phot_bp_rp_excess_factor", "radius_val", "a_g_val", "e_bp_min_rp_val"]	

	return ( dataK_raw, labelsK, data_errorK_raw, dataL_raw, labelsL )

"""
export data as a VOTable file (.xml)
dataTuple is the tuple of data which is being saved
name is used to specify PATH and file name
"""
def exportData( dataTuple ,name):
	
	dataK, labelK, data_errorK, dataL, labelL = dataTuple
	
	# concatenate all arrays from dataTuple and add labels for errors
	data_full = np.concatenate( (np.concatenate( ( dataK, data_errorK ) , axis = 1 ), dataL), axis = 1 )
	error_label = []
	for i in labelK:
		error_label += [i + "_error"]

	label_full = labelK + error_label + labelL
		
	votable = VOTableFile()
	resource = Resource()
	votable.resources.append(resource)
	table = Table(votable)
	resource.tables.append(table)
	
	fields = []
	for i in range( data_full.shape[1] ):
		fields += [Field( votable, name = label_full[i], datatype='float' )]
	
	table.fields.extend( fields )
	table.create_arrays( data_full.shape[0] )
	
	for i in range( data_full.shape[0] ):    
		table.array[i] = tuple( data_full[i,:] )   
	
	votable.to_xml(name +".xml")     

"""
Downloads data from Gaia database
Any query can be set, current one includes radius_val as well
Data is searched at (ra,dec) with radius of radius and saved under ./FINAL/Data/name_
ra, dec and radius are in degrees
"""
def downloadData( ra , dec, radius, name_ ):
	Gaia.launch_job_async( "SELECT gaia_source.source_id,gaia_source.ra,gaia_source.ra_error,gaia_source.dec,gaia_source.dec_error,gaia_source.parallax,gaia_source.parallax_error,gaia_source.parallax_over_error,gaia_source.pmra,gaia_source.pmra_error,gaia_source.pmdec,gaia_source.pmdec_error,gaia_source.phot_g_n_obs,gaia_source.phot_g_mean_flux,gaia_source.phot_g_mean_flux_error,gaia_source.phot_g_mean_flux_over_error,gaia_source.phot_g_mean_mag,gaia_source.phot_bp_n_obs,gaia_source.phot_bp_mean_flux,gaia_source.phot_bp_mean_flux_error,gaia_source.phot_bp_mean_flux_over_error,gaia_source.phot_bp_mean_mag,gaia_source.phot_rp_n_obs,gaia_source.phot_rp_mean_flux,gaia_source.phot_rp_mean_flux_error,gaia_source.phot_rp_mean_flux_over_error,gaia_source.phot_rp_mean_mag,gaia_source.phot_bp_rp_excess_factor,gaia_source.bp_rp,gaia_source.radial_velocity,gaia_source.radial_velocity_error,gaia_source.rv_template_logg,gaia_source.rv_template_fe_h,gaia_source.phot_variable_flag,gaia_source.teff_val,gaia_source.teff_percentile_lower,gaia_source.teff_percentile_upper,gaia_source.a_g_val,gaia_source.e_bp_min_rp_val,gaia_source.flame_flags,gaia_source.radius_val,gaia_source.radius_percentile_lower,gaia_source.radius_percentile_upper,gaia_source.lum_val,gaia_source.lum_percentile_lower,gaia_source.datalink_url,gaia_source.epoch_photometry_url FROM gaiadr2.gaia_source WHERE CONTAINS(POINT('ICRS',gaiadr2.gaia_source.ra,gaiadr2.gaia_source.dec),CIRCLE('ICRS',%f,%f,%f))=1" %(ra,dec,radius), dump_to_file=True, output_file="FINAL/Data/%s" %name_ )
	print( "Finished downloading %s" %name_ )

"""
Reads a table of objects where columns are separated by tabs
Lines which contain %filter at a given position are filtered out and their ra, dec, size and name are saved
Sized is multiplied by a factor (in this case 4 )
All numbers are transformed to degrees 
NOTE: This function is highly dependent on file which is being read, so array indexes may need to be changed for different files 
"""
def readFile( filename, filters = "Oc" ):
	file = open( filename,"r" )
	lines = file.readlines()
	output = []
	factor = 4
	
	for i in lines:
		fields = i.split( "\t" )
		if ( fields[2].strip() == filters ):
			name = fields[0]
			objecttype = fields[2]
			size = float( fields[4] ) * factor / 60
			hours, mins = fields[6].split( " " ) 
			ra = ( float( hours[:-1] ) + ( np.sign( float( hours[:-1] ) )* float( mins[:-1] )/60 ) )* 15
			deg, minutes = fields[7].split()
			dec = float( deg[:-1] ) + ( np.sign( float(deg[:-1]) ) ) * float( minutes[:-1] ) / 60
			output += [ [name,size,ra,dec] ]

	return output
	







