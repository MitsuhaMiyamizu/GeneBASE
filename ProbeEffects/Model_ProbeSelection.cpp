/*
 *  Model_ProbeSelection.cpp
 *  ProbeEffects
 *
 *  Created by Karen Kapur on 2/26/08.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */


#include "Model_ProbeSelection.h"
#include "Model_AffyLiWong.h"
#include "cluster.h"
#include "stdafx.h"
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>



void Model_ProbeSelection::UpdateExpressionEstimates( gsl_matrix* i_pInten, gsl_vector* o_pUpdatedEstimates, set<int>& i_SelectedProbeIndexes ) {
	int numsamples = i_pInten->size1;
	int numprobes = i_SelectedProbeIndexes.size();
	gsl_matrix* pSelectedProbeInten = gsl_matrix_alloc( numsamples, numprobes );
	int colctr = 0;
	set<int>::iterator itBegin = i_SelectedProbeIndexes.begin();
	set<int>::iterator itEnd = i_SelectedProbeIndexes.end();
	for (; itBegin != itEnd; ++itBegin ) {
		int idx = *itBegin;
		gsl_vector_view colview = gsl_matrix_column (i_pInten, idx);
		gsl_matrix_set_col (pSelectedProbeInten, colctr, &colview.vector);
		++colctr;
	}
	Model_AffyLiWong* pModel = new Model_AffyLiWong( pSelectedProbeInten );
	gsl_vector_memcpy( o_pUpdatedEstimates, pModel->GetTheta() );
	delete pModel;
	gsl_matrix_free( pSelectedProbeInten );
}


void Model_ProbeSelection::GetInitialEstimates( gsl_matrix* i_pInten, gsl_vector* o_pPrevEstimates, set<int>& o_SelectedProbeIndexes ) {
	
	int nrows = i_pInten->size1;
	int ncolumns = i_pInten->size2;
	double** pData = &(i_pInten->data); // will this work?  
										// mask a matrix[nrows][ncolumns] with values set to 1.  (or at least not 0). will not be used
										// can set mask to null.  
										// weight double[nrows] = 1 weights the samples, will not be used.  
	int transpose = 1; // cluster probes
	char dist = 'c';
	char method = 'a';
	
	double** pDistmatrix = new double*[ ncolumns ]; // calculate distmatrix ourselves, to avoid errors if sd == 0 for any probe
	for ( int i = 0; i < ncolumns; ++i ) {
		pDistmatrix[i] = new double[ ncolumns ]; // do we need to set the i,i'th element to 1???
		for ( int j = 0; j < i; j++ ) {
			// calculate the correlation of the ith and jth probe
			// get the sds of both probes.  
			double sdi = gsl_stats_sd (&i_pInten->data[i], ncolumns, nrows);
			double sdj = gsl_stats_sd (&i_pInten->data[j], ncolumns, nrows);
			if ( sdi == 0 || sdj == 0 ) {
				pDistmatrix[i][j] = 1;
			} else {
				double covval = gsl_stats_covariance (&i_pInten->data[i], ncolumns, &i_pInten->data[j], ncolumns, nrows);
				double corval = covval / ( sdi * sdj);
				pDistmatrix[i][j] = 1 - corval;
			}
		}
	}
	// Do hierarchical clustering
	// Call treecluster with mask and weights set to NULL.  
	// Ok for average linkage clustering with distance matrix known.  
	Node* pTree = treecluster( nrows, ncolumns, pData, NULL, NULL, transpose, dist, method, pDistmatrix );
	if ( pTree == NULL ) {
		// problem with treecluster
		Parameters::LOGFILE << "Error:  treecluster failed." << endl;
		gsl_vector_set_all( o_pPrevEstimates, -1000 );
		// don't add any indexes to o_SelectedProbeIndexes
	} else {
		int* clusterid = new int[ ncolumns ];
		gsl_vector* pSubclusterSize = gsl_vector_calloc( ncolumns ); 
		cuttreebydistance (ncolumns, pTree, this->m_cut_height, clusterid); 
		// cut the tree at the height m_cut_height and find the largest subcluster
		for ( int i = 0; i < ncolumns; ++i ) {
			double val = gsl_vector_get( pSubclusterSize, clusterid[i] );
			gsl_vector_set( pSubclusterSize, clusterid[i], (val + 1));
		}
		int maxindex = gsl_stats_max_index (pSubclusterSize->data, 1, ncolumns);
		int numSelectedProbes = gsl_vector_get( pSubclusterSize, maxindex );
		gsl_matrix* pSelectedProbeInten = gsl_matrix_alloc( nrows, numSelectedProbes );
		int colctr = 0;
		for ( int j = 0; j < ncolumns; ++j ) {
			if ( clusterid[j] == maxindex ) {
				gsl_matrix_set_col (pSelectedProbeInten, colctr, &(gsl_matrix_column (i_pInten, j).vector) );
				o_SelectedProbeIndexes.insert( j );
				++colctr;
			}
		}
		// Get initial li-wong gene expr estimates
		// store in pPrevEstimates
		if ( o_SelectedProbeIndexes.size() > 0 ) {
			Model_AffyLiWong* pModel = new Model_AffyLiWong( pSelectedProbeInten );
			gsl_vector_memcpy (o_pPrevEstimates, pModel->GetTheta());
			delete pModel;
		} else {
			gsl_vector_set_all( o_pPrevEstimates, 0);
		}
		gsl_matrix_free( pSelectedProbeInten );
		gsl_vector_free( pSubclusterSize );
		delete [] clusterid;
	}
	for ( int i = 0; i < ncolumns; ++i ) { // free distmatrix allocated memory
		delete [] pDistmatrix[i];
	}
	delete [] pDistmatrix;
	free( pTree );
}

Model_ProbeSelection::Model_ProbeSelection(gsl_matrix* i_pInten ) : m_corr_cutoff(0.70), m_cut_height(0.50), m_pExpression(NULL)
{
	//***
	// Step 1:  Get Li-Wong estimates
	//***
	// i_pInten is an I by J matrix of probe intensities, 
	// where I = num samples, J = num probes
	int numsamples = i_pInten->size1;
	int numprobes = i_pInten->size2;
	this->m_pExpression = gsl_vector_alloc( numsamples );
	gsl_vector* pLiWongExpression = gsl_vector_alloc( numsamples );
	Model_AffyLiWong* pModel = new Model_AffyLiWong( i_pInten );
	gsl_vector_memcpy( pLiWongExpression, pModel->GetTheta() );
	delete pModel;
	set<int> selectedProbeIndexes;
	set<int> finalProbeIndexes;
	this->GetInitialEstimates( i_pInten, m_pExpression, selectedProbeIndexes );
	int numSelectedProbes = selectedProbeIndexes.size();
	if ( numSelectedProbes < 6) {
		// copy Li-Wong estimates to m_pExpression
		gsl_vector_memcpy( m_pExpression, pLiWongExpression );
	} else { // start iterative method
		int numIter = 0;
		while ( numIter <= 50 && numSelectedProbes >= 6 ) {			
			double expressionsd = gsl_stats_sd (m_pExpression->data, 1, numsamples);
			// calculate the pearson correlation between probe intensity and estimates
			for ( int j = 0; j < numprobes; ++j ) {
				double probesd = gsl_stats_sd (&i_pInten->data[j], numprobes, numsamples);
				if ( expressionsd != 0 && probesd != 0 ) {
					double covval = gsl_stats_covariance ( &i_pInten->data[j], numprobes, m_pExpression->data, 1, numsamples);
					double corval = covval / ( probesd * expressionsd );
					if ( corval > 0.7 ) {
						finalProbeIndexes.insert( j ); // if correlation is greater than 0.7, add the probe to final_probes
					}			
				}
			}
			// if final_probes == selected_probes, break (use stl equal algorithm)
			// else set selected_probes = final_probes, 
			bool setequality = ( selectedProbeIndexes.size() == finalProbeIndexes.size() );
			if ( setequality == true ) {
				setequality = equal( selectedProbeIndexes.begin(), selectedProbeIndexes.end(), finalProbeIndexes.begin() );
			}
			if ( setequality == true ) {
				break;
			} else {
				selectedProbeIndexes = finalProbeIndexes; // will this work? 
				numSelectedProbes = selectedProbeIndexes.size();
				if ( numSelectedProbes >= 6 ) {
					UpdateExpressionEstimates( i_pInten, m_pExpression, selectedProbeIndexes ); 
					// should we add .5 to each theta estimate? 
				} else {
					// replace m_pExpression with pLiWongExpression
					gsl_vector_memcpy( m_pExpression, pLiWongExpression );
					break;
				}
			}
			++numIter;
		}
	}
	gsl_vector_free( pLiWongExpression );
}

Model_ProbeSelection::~Model_ProbeSelection(void) {
	if ( m_pExpression != NULL ) {
		gsl_vector_free( m_pExpression );
	}
}

gsl_vector* Model_ProbeSelection::GetExpression()
{
	return this->m_pExpression;
}



