import numpy as np

"""s
takes in a dataTuple and returns it filtered out
it filters out all points in dataK with missing data 
radial_velocity is assumed to be missing most of data, therefore if more than threshold (default 80%) of data is missing it is completely removed
"""
def filterMissing( dataTuple,threshold = 0.8 ):

	dataK, labelsK, data_errorK, dataL, labelsL = dataTuple
	
	count = 0 
	for i in dataK[:,-1]:
		if not isinstance( i, float ):
			count +=1
		
	length = dataK.shape[0]	
	
	if ( 1 - count / length < threshold ):
		dataK = dataK[:,0:-1]
		data_errorK = data_errorK[:,0:-1] 
		labelsK = labelsK[0:-1]
		print( "radial_velocity is NOT included in the data, as only %.2f of stars are not missing data" %( 1- count / length ) )
	else:
		print( "radial_velocity IS included in the data" )

	flagarray = []
	for i in range( dataK.shape[0] ):
		for j in dataK[i]:
			if not isinstance( j, float ):
				flagarray += [False]
				break
		else:
			flagarray += [True]

	
	dataK = dataK[flagarray,:]

	data_errorK = data_errorK[flagarray,:]
	
	dataL = dataL[flagarray,:]

	print( "There are %i stars left" %dataK.shape[0] )
	
	return ( dataK, labelsK, data_errorK, dataL, labelsL )
	
"""
Computes relative error from absolute and replaces it in dataTuple 
It prints out that replacement has happened, as a warning to the user 
"""
def computeError( dataTuple ):
	dataK, labelsK, data_errorK, dataL, labelsL = dataTuple
	
	rel_error = []
	for i in range( dataK.shape[0] ):
		hold = []
		for j in range( dataK.shape[1] ):
			hold += [abs( data_errorK[i,j] / dataK[i,j] )]
		rel_error +=[hold]

	rel_error = np.array(rel_error)
	
	print("Absolute error array is replaced with relative error array")

	return (dataK, labelsK, rel_error, dataL, labelsL)

"""
Filters out data with relative error higher than threshold at data_errorK[:,index]
It returns filtered out dataTuple 
NOTE: This can only be applied if the dataTuple contains relative rather than apsolute error
"""
def filterRelativeError( dataTuple, threshold, index ):
	
	dataK, labelsK, data_errorK, dataL, labelsL = dataTuple

	flagarray = []
	for i in range( dataK.shape[0] ):
		if ( data_errorK[i][index] > threshold ):
				flagarray += [False]
		else:
			flagarray += [True]
	
	dataK = dataK[flagarray,:]
	data_errorK = data_errorK[flagarray,:]	
	dataL = dataL[flagarray,:]
	
	print("%i stars are filtered out and %i are left, due to %s error higher than %.2f" %( flagarray.count( False ), flagarray.count( True ),labelsK[index], threshold ) )
	
	return ( dataK, labelsK, data_errorK, dataL, labelsL )

def filterParallaxError(dataTuple,threshold=0.8):
	return filterRelativeError(dataTuple,threshold,4)

def filterPmraError(dataTuple,threshold=0.8):
	return filterRelativeError(dataTuple,threshold,2)
	
def filterPmdecError(dataTuple,threshold=0.8):
	return filterRelativeError(dataTuple,threshold,3)
	
"""
Filters out stars with magnitude higher than 20, as they are too dark and will have high error 
"""
def filterHighMagnitude( dataTuple ):
	dataK, labelsK, data_errorK, dataL, labelsL = dataTuple
	
	flagarray = []
	for i in range( dataL.shape[0] ):
		if ( dataL[i,1] > 20 ):
			flagarray += [False]
		else:
			flagarray += [True]
	
	dataK = dataK[flagarray,:]
	dataL = dataL[flagarray,:]
	data_errorK = data_errorK[flagarray,:]	
	
	print("%i stars are filtered out (too dark) and %i are left" %( flagarray.count( False ), flagarray.count( True ) ) )
	
	return ( dataK,labelsK, data_errorK, dataL, labelsL )

"""
Separates data which is inside a cluster with label = mark from other data 
inside separates weather data inside (True) or outside (False) should be taken (default True)
"""
def filterCluster(dataTuple, filt ,inside = True,mark = 0): 
	
	dataK, labelsK, data_errorK, dataL, labelsL = dataTuple

	# Split the data on stars inside and stars outside the cluster
	dataL_cluster = dataL[filt == mark]
	dataL_out = dataL[filt != mark]	
	dataK_cluster = dataK[filt == mark]
	dataK_out = dataK[filt != mark]
	data_errorK_cluster = data_errorK[filt == mark]
	data_errorK_out = data_errorK[filt != mark]
	
	if ( inside ):
		return ( dataK_cluster, labelsK, data_errorK_cluster, dataL_cluster, labelsL )
	else:
		return ( dataK_out, labelsK, data_errorK_out, dataL_out, labelsL )

"""
applies all filters in correct order with default values 
this is used to speed up process, since this is applied to the data at the very beginning 
"""
def filterAll (datatuple):
	return   filterHighMagnitude( filterPmraError( filterPmdecError( filterParallaxError( computeError( filterMissing( datatuple ) ) ) ) ) ) 

