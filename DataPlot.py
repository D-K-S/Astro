import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as clr
from mpl_toolkits import mplot3d
import math as m
from astropy import units as u
from astropy.coordinates import SkyCoord
from adjustText import adjust_text # manually installed library 

"""
Plots right ascension vs declination graph with points coloured according to their cluster 
if color is specified point color is chosen from colorlist, otherwise all points are coloured black
Colorlist colors are chosen from biggest to smallest clusters, noise is always grey and clusters of size < 50 are not plotted at all 
Legend shows cluster labels from biggest being 1 to smallest, these do not correspond to labels given by algortihm, these labels are clusterlabels[i] 
if name is specified the graph will be specified and not shown 
NOTE: colorlist has finite lenght so an error may occur, but that is higly unlikely. If the algorithm finds more than 120 clusters there are bigger issues than list length 
"""
def plotRaDec(dataTuple,color = (0,0,0), name = ""):
	
	dataK, labelsK, data_errorK, dataL, labelsL = dataTuple
	
	colorlist = [(0,0,1,1),(0,1,0,1),(1,0,0,1),(0.5,0.5,0,1),(0.5,0,0.5,1),(0,0.5,0.5,1)] + [(0.5,0,0,1),(0,0.5,0,1),(0,0,0.5,1),(0.5,0.5,0,1),(0.5,0,0.5,1),(0,0.5,0.5,1)] + [(1,0.5,0,1),(0.5,1,0,1),(1,0,0.5,1),(0.5,0,1,1),(0,1,0.5,1),(0,0.5,1,1)] + [(0,0,0,0.5)] * 100

	ra = dataK[:,0]
	dec = dataK[:,1]

	plt.figure(figsize = (20,12))
	plt.title("Declination plotted against right ascension", fontsize = 24, weight = 'bold' )
	plt.xlabel("Right ascension/° ", fontsize = 24, weight = 'bold' )
	plt.ylabel("Declination/° ", fontsize = 24, weight = 'bold' )
	
	colorlabels = list( set(color ) )
	colorlabels.sort(key = color.tolist().count,reverse = True)
	if( color != ( 0 , 0 , 0 ) ):
		
		for i in range( len( colorlabels) ):
			if ( color.tolist().count( colorlabels[i] )  < 50):
				break
			if (colorlabels[i] != -1):
				# -1 is added to colorlist because noise is expected to be first in the colorlabels list so not to waste colour 
				plt.scatter( ra[color == colorlabels[i]],dec[color == colorlabels[i]], c = [colorlist[i - 1]] * color.tolist().count( colorlabels[i] ) , s = 10, label = 'Cluster %i' %i, lw = 0)
			else:
				plt.scatter(ra[color == -1], dec[color==-1], c = [colorlist[-1]] * color.tolist().count( colorlabels[i] ) , lw = 0, s = 7)
			
		plt.legend(ncol = 1, frameon = True, fontsize = 12, handlelength = 2, borderpad = 1.8, handletextpad = 1, scatterpoints = 3,loc ='lower left')
	
	else:
		plt.scatter(ra,dec,c = color, s = 5)
	
		
	if (name != ""):
		plt.savefig(name + "_decvsra.png",dpi = 600)
		plt.close()
	else:
		plt.show()

"""
plots three HR diagram:
#1: showing only cluster stars, colorcoded with error
#2: showing all stars, colorcoded with error
#3: showing only cluster stars, colorcoded with teff_val and size coded with radius_val (this means that surface of the point is proportional to radius of star), radius to radius was tried but points tend to be too big
scaling factor determines how big a default point will be
if name is specified plots are saved separately without showing 
"""
def plotHR(dataTuple_cluster, dataTuple_out, name = "" ):
	
	dataK_cluster, labelsK, data_errorK_cluster, dataL_cluster, labelsL = dataTuple_cluster
	#dataK_out, labelsK, data_errorK_out, dataL_out, labelsL = dataTuple_out
	
	phot_mag_in = dataL_cluster[:,1] 
	starcolor_in = dataL_cluster[:,0] 
	
	#phot_mag_out = dataL_out[:,1] 
	#starcolor_out = dataL_out[:,0] 
	"""
	Generate a colormap based on relative errors in proper motion and parallax
	Root is used to magnify the effect of error
	Points with lower errors are darker (higher alpha)
	Points with relative error higher than 1 will have alpha 0 and not be visible
	"""
	"""
	color_cluster = []
	for i in range( data_errorK_cluster.shape[0] ):
		relpmra  = data_errorK_cluster[i,2]
		relpmdec = data_errorK_cluster[i,3]
		relpar   = data_errorK_cluster[i,4]
		if ( relpmra > 1 or relpmdec > 1 or relpar > 1 ):
			alpha = 0
			color_cluster += [( 0, 0, 0, 0 )]
		else:
			alpha = ( 1 - relpmra**0.5 ) * ( 1 - relpmdec**0.5 ) * ( 1 - relpar**0.5 )
			color_cluster += [( relpmra**0.5, relpmdec**0.5, relpar**0.5, alpha )]	
	"""
	
	"""
	color_out = [] 
	for i in range( data_errorK_out.shape[0] ):
		relpmra  = data_errorK_out[i,2]
		relpmdec = data_errorK_out[i,3]
		relpar   = data_errorK_out[i,4]
		if ( relpmra > 1 or relpmdec > 1 or relpar > 1 ):
			alpha = 0
			color_out += [( 0, 0, 0, 0 )]
		else:
			alpha = ( 1 - relpmra**0.5 ) * ( 1 - relpmdec**0.5 ) * ( 1 - relpar**0.5 )
			color_out += [( relpmra**0.5, relpmdec**0.5, relpar**0.5, alpha )]	
	"""
	#1 
	plt.figure( figsize = ( 20 , 12 ) )
	plt.title("Apparent magnitude plotted against colour, cluster stars only") 
	plt.xlabel("Colour(B-R) [mag]")
	plt.ylabel("Apparent magnitude [mag]")
	plt.gca().invert_yaxis()
	plt.scatter( starcolor_in, phot_mag_in , s = 5 , c = [(0,0,0)] * len(phot_mag_in.tolist())  )
		
	if ( name != "" ):
		plt.savefig( name + "_HR_in.png" , dpi = 600 )
		plt.close()
	"""
	#2 
	plt.figure( figsize = ( 20 , 12 ) )
	plt.title("Apparent magnitude plotted against colour, all stars ") 
	plt.xlabel("Colour(B-R) [mag]")
	plt.ylabel("Apparent magnitude [mag]")
	plt.gca().invert_yaxis() 
	plt.scatter( starcolor_out , phot_mag_out , s = 5 , c = color_out , alpha = 0.05 )
	plt.scatter( starcolor_in, phot_mag_in, s = 5 , c = color_cluster, alpha = 1 )
	
	if (name != ""):
		plt.savefig( name+"_HR_all.png", dpi = 600 )
		plt.close()
	"""

	"""
	#3 
	#NORM used is LOGNORM() to get better colors from cmap
	plt.figure( figsize = ( 20 , 12 ) )
	temperature = np.array( dataL_cluster[:,2], dtype = 'float' )
	scaling_factor = 5
	size = list( map( lambda x:scaling_factor*x, np.array( dataL_cluster[:,5], dtype = 'float' ) ) )	
	plt.title("Apparent magnitude plotted against colour, cluster stars only", fontsize = 24, weight = 'bold') 
	plt.xlabel("Colour(B-R) [mag]", fontsize = 24, weight = 'bold' )
	plt.ylabel("Apparent magnitude [mag]", fontsize = 24, weight = 'bold' )
	plt.gca().invert_yaxis()
	plt.scatter( starcolor_in, phot_mag_in, s = size, c = temperature, cmap = 'RdYlBu', vmax = 15000 , vmin = 3000, norm = clr.LogNorm() )
	cb = plt.colorbar()
	cb.set_label('Temperature', size = 15, weight = 'bold' )
	
	l1, = plt.plot( [], [] ,'or', markersize = 0.5*scaling_factor )
	l2, = plt.plot( [], [], 'or', markersize = scaling_factor )
	l3, = plt.plot( [], [], 'or', markersize = scaling_factor*3 )
	
	legendlabels = ['0.5 Solar radii', '1 solar radius', '3 Solar radii']
	
	leg = plt.legend( [l1, l2, l3], legendlabels, ncol = 1, frameon = True, fontsize = 20, handlelength = 2, borderpad = 1.8, handletextpad = 1, scatterpoints = 1,loc ='lower left')
	
	if ( name != "" ):
		plt.savefig( name + "_HR_all_colour.png", dpi = 600 )
		plt.close()	
	else:
		plt.show()

	"""
	
"""
plots a histogram of relative errors
by default parallax error is show, as that is biggest error
if name is specified plot is saved without showin
"""
def plotHist( dataTuple, index = 4, name = "" ):
	dataK, labelsK, data_errorK, dataL, labelsL = dataTuple
	
	# 16 bins because data will already be from 0 to 0.8 relative error
	parallaxerror = data_errorK[:,index]
	plt.hist( parallaxerror,bins=16,rwidth = 0.8 )
	plt.title("%s relative error histogram" %labelsK[index])
	plt.xlabel("relative parallax error")
	plt.ylabel("count")
	
	if ( name != "" ):
		plt.savefig( name + "_hist.png" )
		plt.close()
	else:
		plt.show()

"""
plots a 3d plot in Cartesian coordinates
it transforms ra, dec and parallax from spherical to cartesian coordinates
"""
def plot3Dcoord(ra,dec,parallax, color = (0,0,1)):
	ax = plt.axes(projection='3d')
	ax.set_xlabel('x')
	ax.set_ylabel('y')
	ax.set_zlabel('z')
	ax.set_xlim(-2,2)
	ax.set_ylim(-2,2)
	ax.set_zlim(-2,2)

	z = list( map( lambda r,fi:(1/r)*m.sin(fi*m.pi/180), parallax, dec) )
	y = list( map( lambda r,fi,th:(1/r)*m.cos(fi*m.pi/180)*m.sin(th*m.pi/180),parallax,dec,ra ) )
	x = list( map( lambda r,fi,th:(1/r)*m.cos(fi*m.pi/180)*m.cos(th*m.pi/180),parallax,dec,ra ) )

	ax.scatter3D(x,y,z, s = 3, c = color,cmap = 'Spectral')
	plt.show()

"""
Takes in a list of clusters and plots their galactic coordinates in mollweide coordinate system
faillist constains list of all failed clusters
if the cluster was NOT found it will be painted red, otherwise it is painted blue 
NOTE: It uses an external library for adjust_text 
"""
def plotGalactic( clusters ):
	
	faillist = ['M16','M26','M29','M39']
	
	hold = np.array( clusters )
	coord = np.array( hold[:,2:4], dtype = 'float' )
	pointlabels = hold[:,0]
	ra = coord[:,0]
	dec = coord[:,1]
	
	coord = SkyCoord( ra = ra*u.degree, dec = dec*u.degree )
	l = np.radians ( coord.galactic.l.degree ) 
	b = np.radians( coord.galactic.b.degree )
	# Convert to range from -pi to pi 
	l = (l > np.pi) * (-2*np.pi) + l
	b = (b > np.pi) * (-2*np.pi) + b 
	
	flagarray = []
	for i in hold[:,0]:
		if (i in faillist):
			flagarray += [True]
		else:
			flagarray += [False]
	
	flagarray = np.array(flagarray)
	
	fig = plt.figure( figsize = ( 20 , 12 ) )
	plt.subplot( 111, projection='mollweide' )
	# Good looking shade of blue 
	plt.scatter( l[flagarray], b[flagarray], c = [(1,0,0)] * flagarray.tolist().count(True), label = 'Not deteceted')
	plt.scatter( l[flagarray == False], b[flagarray == False], c = [(0,0,1)] * flagarray.tolist().count(False), label = 'Detected')
	plt.xlabel("Galactic longitude: l/° ", fontsize = 24, weight = 'bold' )
	plt.ylabel("Galactic latitude: b/° ", fontsize = 24, weight = 'bold' )
	plt.title("Galactic position of Messier open clusters", fontsize = 24, weight = 'bold' )
	plt.grid()
	 
	texts = []
	for i in range( len(l) ):
		texts.append( plt.text( l[i] , b[i] , pointlabels[i] ,fontsize = 12) )
	
	adjust_text(texts, arrowprops = dict(arrowstyle = "-", color = 'b', lw = 0.6) )	
	
	leg = plt.legend( ncol = 1, frameon = True, fontsize = 24, handlelength = 2, borderpad = 1.8, handletextpad = 1, scatterpoints = 1,loc ='lower left')
	
	
	plt.show()

















