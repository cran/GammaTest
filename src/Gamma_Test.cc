#include <math.h>		// math routines
#include "ANN.h"        // ANN headers
#include <R.h>			// R headers

//-----------------------------------------------------------------------------
//			Main Gamma Test Program for R
//-----------------------------------------------------------------------------
extern "C"
{
    void Gamma_Test_Main(double *data, int *mask, int *sumMask, int *k, int *dim, 
	                 int *m_pts, double *eps, double *delta, double *gamma)
    {
	int				d;			// Actual Dimension
	int				M;			// Number of Data points
	double			error_bound;// enough said!
	int				numNN;		// Max. num of NN
	ANNpointArray	data_pts;	// Data points
	ANNidxArray		nn_idx;		// Near neighbor indices
	ANNdistArray	dists;		// Near neighbor distances
	ANNkd_tree		*the_tree;	// Search structure

	d			= *dim;
	M			= *m_pts;
	numNN	    = *k;
	error_bound = *eps;

	double output_pts[M];					// Allocate query point
	data_pts   = annAllocPts(M, *sumMask);	// Allocate data points
	nn_idx 	   = new ANNidx[numNN];			// Allocate near neigh indices
	dists 	   = new ANNdist[numNN];		// Allocate near neighbor dists

	float sumDists[numNN];		       // Sum the delta values
	float sumOutDists[numNN];	       // Sum the Gamma values
	float outDist[numNN];    	       // Calculate each gamma value

	//---------------- Put data into ANN format -------------------------
	int incOutputData 	= (d-1)*M;
	int d_ptr[d-1];

	for(int i = 0; i < d-1; i++)
	{
	    d_ptr[i] = 0;
	    d_ptr[i] = i*M;
	}
	
	for(int i = 0; i < M; i++)
	{
	    int myJ = 0;
	    for(int j = 0; j < d-1; j++)
	    {
		int temp = d_ptr[j];

		if(mask[j] == 1)
		{
		    data_pts[i][myJ] = 00.00;
		    data_pts[i][myJ] = data[temp];
		    myJ++;
		}

		d_ptr[j] = 0;
		d_ptr[j] = temp + 1;
	    }
	    
	    output_pts[i] = 0000.000;
	    output_pts[i] = data[incOutputData];
	    incOutputData++; 
	}

	//-------------------- Build kd-tree -----------------------------

	the_tree = new ANNkd_tree(	
		    data_pts,	// The data points
		    M,		// Number of points
		    *sumMask);	// Dimension of space

	for (int j = 0; j < numNN; j++) 	// Initialize the neccessary arrays
	{
	    dists[j] 		= 0.0000;
	    sumDists[j] 	= 0.0000;
	    outDist[j] 		= 0.0000;
	    sumOutDists[j] 	= 0.0000;
	}

	for(int i = 0; i < M; i++)	// read query points
	{
	    the_tree->annkSearch(	// search
		data_pts[i],			// query point
		numNN,					// number of near neighbors
		nn_idx,					// nearest neighbors (returned)
		dists,					// distance (returned)
		error_bound);			// error bound

	    for (int j = 0; j < numNN; j++)
	    {
	    	dists[j] = sqrt(dists[j])*sqrt(dists[j]); // unsquare distance
	    	outDist[j] = (output_pts[nn_idx[j]] - output_pts[i])*(output_pts[nn_idx[j]] - output_pts[i]);

	    	sumDists[j] += dists[j];
	    	sumOutDists[j] += outDist[j];
	    }
	}

	for(int i = 0; i < numNN; i++)
	{
	    delta[i] = 0.00000000;
	    gamma[i] = 0.00000000;
	    delta[i] = double(sumDists[i] / M);
	    gamma[i] = double(sumOutDists[i] / (2*M));
	}

	// Do a little bit of memory management......
	annDeallocPts(data_pts);
	delete [] nn_idx;
	delete [] dists;
	delete the_tree;
    }
}
