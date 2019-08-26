from sklearn.cluster import DBSCAN
from sklearn.preprocessing import StandardScaler
from sklearn.neighbors import NearestNeighbors
import pickle
import numpy as np
import matplotlib.pyplot as plt 
import statistics

# Directory to which files are saved and loaded
savedir = ""

"""
sets the savedirectory for clustering algorithm
"""
def setDirectory( s_ ):
	global savedir  
	s = str(s_)
	print( "Save directory was set to %s" %s )
	savedir = s 

"""
runs the DBScan clustering algorithm on data from dataTuple
epsilon needs to be set and min_samples = 2* %number_of_dimensions
data is scaled before running the algorithm
Algorithm is called on epsilon, epsilon-interval and epsilon+interval, this was added due to rough estimation of epsilon with knn distance graph 
If the algorithm was already run with that epsilon it is loaded from a file, epsilon is saved to three digit precision
Function returns a list of labels found by algorithm for inputed data
It also prints epsilon value used, number of clusters found (not including noise) and number of stars for all clusters which have >= 50 points
50 is chosen as boundary in order to eliminate some noise from data, typical size of cluster is 100+ 
NOTE: Make sure to change directories when changing data, otherwise data might become corrupted due to load/save method
"""
def runDBScan( dataTuple, epsilon,interval = 0.0):

	dataK, labelsK, data_errorK, dataL, labelsL = dataTuple

	try:
		with open( savedir +'.DBScan_%.3f_est' %epsilon, 'rb' ) as input:
			db,y_pred = pickle.load( input )
		
	except:
			data_scaled = StandardScaler().fit_transform( dataK )
			db = DBSCAN( eps = epsilon, min_samples = 2*dataK.shape[1] )
			db.fit( data_scaled )

			y_pred = db.labels_
			with open( savedir +'.DBScan_%.3f_est' %epsilon , 'wb' ) as output:
				pickle.dump( ( db,y_pred ) , output, pickle.HIGHEST_PROTOCOL )

	print( "Finished DBScan with epsilon = %.3f." %epsilon )

	cluster_number = len( set( y_pred ) )
	
	print( "Cluster number: %i " %( cluster_number - 1 ) )  	 
	
	normal_cluster_number = 0

	for i in set( db.labels_ ):
		if i != -1:
			occurance = db.labels_.tolist().count( i )
			if occurance >= 50:
				print( "Star number with label %i is %i" %( i, list( db.labels_ ).count( i ) ) )
				normal_cluster_number += 1
	
	print("There are %i clusters with 50 or more stars" %normal_cluster_number)
	
	if ( abs( interval ) > 0.001 ):
		print( "Running DBScan for upper and lower bound" )
		runDBScan(dataTuple,epsilon+interval,interval = 0.0)
		runDBScan(dataTuple,epsilon-interval,interval = 0.0)
		
	return y_pred

"""
plots a K nearest neighbors graph which is used to estimate epsilon for DBScan
K is 2*%number_of_dimensions by default 
if name is specified graph is saved, as name_knn.png, without showing it, otherwise it will be shown and not saved
#NOTE: It is recommended to look at knn graph and save it manually, as it needs to be zoomed to an area which cannot be determined a priori
"""
def knnHelp( dataTuple, name = "" ):
	
	data, labelsK, data_errorK, dataL, labelsL = dataTuple
	
	data_scaled = StandardScaler().fit_transform(data)
	nbrs = NearestNeighbors( n_neighbors = data_scaled.shape[1] * 2 ).fit( data_scaled )
	distances, indicies = nbrs.kneighbors( data_scaled )
	distanceLast = sorted( distances[:,-1] , reverse = True)
	
	plt.figure( figsize = ( 20, 12 ) )
	plt.plot( np.arange( 0, len( distanceLast ), 1 ) , distanceLast )
	
	if ( name != "" ):
		plt.savefig( name + "_knn.png", dpi = 600 )
		plt.close()
	else:
		plt.show()

"""
it determines the average position of the data provided
Position is calculated as a regular mean of data points
Function prints out the name of object observed ( if specified ) and its right ascension in hh/mm and declination in deg/mm
NOTE: This is mostly used to simplify mapping  between labels and cluster colours in the plot 
"""
def averagePosition( dataTuple, name = "" ):
	
	data, labelsK, data_errorK, dataL, labelsL = dataTuple
	ra = statistics.mean( data[:,0] )
	dec = statistics.mean( data[:,1] )
	
	decd = int(dec)
	decm = ( dec - decd ) * 60
	rah = int (ra / 15)
	ram = ( ra / 15 - rah ) * 60
	
	print ( name + ( " @ Ra: %ih %im, Dec: %.i° %i'" %( rah, ram, decd, decm ) ) )






